import sys
import os
from collections import Counter, defaultdict
from contaminants_cython import test_for_contaminant
from general_sequence_tools import dna_rev_comp


def get_fastqc_contaminant_list(
        contaminant_file,
        max_contaminant_length=None,
        include_rev_comp=True
        ):
    """
    Reads given contaminant file and returns a list of tuples with contaminant name and sequence.
    """
    contaminant_list = []
    for line in open(contaminant_file):
        if line[0] == '#':
            continue
        var = line.strip().split('\t')
        name = var[0]
        seq = var[-1]
        if not name:
            continue  # Blank lines
        if 'A' not in seq and 'C' not in seq and 'G' not in seq and 'T' not in seq:
            sys.exit('Non-DNA contaminant, %s: %s' % (name, seq))

        contaminant_list.append((name, seq[:max_contaminant_length]))
        if include_rev_comp:
            contaminant_list.append(
                (name + ' RevComp', dna_rev_comp(seq[:max_contaminant_length])))

    if not contaminant_list:
        sys.exit('No contaminants found in ' + contaminant_file)

    return contaminant_list


def build_contaminant_kmers_dict(contaminant_list, k):
    """
    Builds dict with kmers as keys and a list of all contaminant names containing the given kmer as
    values for all contaminants in given contaminant list and for given value of k.
    """
    # Build the contaminant_kmers dict as a defaultdict of sets for uniqueness
    contaminant_kmers = defaultdict(set)
    for name, seq in contaminant_list:
        for i in xrange(len(seq)-k+1):
            contaminant_kmers[seq[i:i+k]].add(name)
    # Clean the dict by making it a regular dict and converting sets to lists
    output_dict = {}
    for kmer, contaminant_name_set in contaminant_kmers.items():
        output_dict[kmer] = list(contaminant_name_set)
    return output_dict


def get_fastqc_contaminant_kmers_dict(contaminant_file, k):
    """
    Wrapper around build_contaminant_kmers_dict given contaminant filename instead of contaminant
    list.
    """
    contaminant_list = get_fastqc_contaminant_list(contaminant_file)
    return build_contaminant_kmers_dict(contaminant_list, k)


def get_fastqc_contaminant_dict(contaminant_file):
    """
    Builds dict of seq given name for contaminants in given file name.
    """
    contaminant_list = get_fastqc_contaminant_list(contaminant_file)
    contaminant_dict = {}
    for name, seq in contaminant_list:
        contaminant_dict[name] = seq
    return contaminant_dict


def output_contaminant_removal_statistics(
        total_reads,
        trimmed_reads,
        deleted_reads,
        contaminants_found,
        log_file_handle=sys.stdout,
        contaminant_label='Contaminant',
        ):
    """
    Write a summary of the contaminants removed to log_file_handle (default stdout).
    """
    contaminant_label = contaminant_label[0].upper() + contaminant_label[1:]
    log_file_handle.write('%s removal statistics:\n' % contaminant_label)
    log_file_handle.write('Total reads: %d\n' % total_reads)
    if trimmed_reads:
        log_file_handle.write('Trimmed reads: %d, %f%%\n' %
                              (trimmed_reads, 100*float(trimmed_reads)/total_reads))
    if deleted_reads:
        log_file_handle.write('Deleted reads: %d, %f%%\n' %
                              (deleted_reads, 100*float(deleted_reads)/total_reads))
    if contaminants_found:
        log_file_handle.write('%ss found:\n' % contaminant_label)
        for name, quantity in reversed(sorted(contaminants_found.items(), key=lambda tup: tup[1])):
            log_file_handle.write('  %12d  %s\n' % (quantity, name))


def remove_paired_read_contaminants(
        fname1,
        fname2,
        k=11,                         # k-mer to use exact matches to contaminants
        alignment_score_cutoff=14,    # Minimum Smith-Waterman score required to delete a read
        log_file_handle=sys.stdout
        ):
    outname1 = os.path.splitext(fname1)[0] + '_decontaminated_%d.fastq' % alignment_score_cutoff
    outname2 = os.path.splitext(fname2)[0] + '_decontaminated_%d.fastq' % alignment_score_cutoff

    contaminant_file = '/home/hawkjo/python_src/sequence_tools/contaminant_list.txt'
    contaminant_list = get_fastqc_contaminant_list(contaminant_file)
    contaminant_dict = get_fastqc_contaminant_dict(contaminant_file)
    contaminant_kmers_dict = build_contaminant_kmers_dict(contaminant_list, k)

    f1 = open(fname1)
    f2 = open(fname2)
    out1 = open(outname1, 'w')
    out2 = open(outname2, 'w')

    total_reads = 0
    trimmed_reads = 0
    deleted_reads = 0
    contaminants_found = Counter()
    while True:
        defline1 = f1.readline().strip()
        defline2 = f2.readline().strip()
        if not defline1 or not defline2:
            if not defline1 and not defline2:
                break  # End of file
            else:
                sys.exit('Unequal number of lines in paired files')

        total_reads += 1

        seqline1 = f1.readline().strip()
        plusline1 = f1.readline().strip()
        qualline1 = f1.readline().strip()

        seqline2 = f2.readline().strip()
        plusline2 = f2.readline().strip()
        qualline2 = f2.readline().strip()

        contaminant1_name = test_for_contaminant(seqline1,
                                                 alignment_score_cutoff,
                                                 k,
                                                 contaminant_kmers_dict,
                                                 contaminant_dict)
        if contaminant1_name:
            deleted_reads += 1
            contaminants_found[contaminant1_name] += 1
            continue

        contaminant2_name = test_for_contaminant(seqline2,
                                                 alignment_score_cutoff,
                                                 k,
                                                 contaminant_kmers_dict,
                                                 contaminant_dict)
        if contaminant2_name:
            deleted_reads += 1
            contaminants_found[contaminant2_name] += 1
            continue

        # Non-deleted read
        out1.write('\n'.join([defline1, seqline1, plusline1, qualline1]) + '\n')
        out2.write('\n'.join([defline2, seqline2, plusline2, qualline2]) + '\n')

    f1.close()
    f2.close()
    out1.close()
    out2.close()

    output_contaminant_removal_statistics(
        total_reads,
        trimmed_reads,
        deleted_reads,
        contaminants_found,
        log_file_handle,
        )

    return outname1, outname2


def remove_single_read_contaminants(
        fname,
        k=11,                         # k-mer to use exact matches to contaminants
        alignment_score_cutoff=14,    # Minimum Smith-Waterman score required to delete a read
        log_file_handle=sys.stdout
        ):
    outname = os.path.splitext(fname)[0] + '_decontaminated_%d.fastq' % alignment_score_cutoff

    contaminant_file = '/home/hawkjo/python_src/sequence_tools/contaminant_list.txt'
    contaminant_list = get_fastqc_contaminant_list(contaminant_file)
    contaminant_dict = get_fastqc_contaminant_dict(contaminant_file)
    contaminant_kmers_dict = build_contaminant_kmers_dict(contaminant_list, k)

    f = open(fname)
    out = open(outname, 'w')

    total_reads = 0
    trimmed_reads = 0
    deleted_reads = 0
    contaminants_found = Counter()
    while True:
        defline = f.readline().strip()
        if not defline:
            break  # End of file
        total_reads += 1

        seqline = f.readline().strip()
        plusline = f.readline().strip()
        qualline = f.readline().strip()

        contaminant_name = test_for_contaminant(seqline,
                                                alignment_score_cutoff,
                                                k,
                                                contaminant_kmers_dict,
                                                contaminant_dict)
        if contaminant_name:
            deleted_reads += 1
            contaminants_found[contaminant_name] += 1
        else:
            # Non-deleted read
            out.write('\n'.join([defline, seqline, plusline, qualline]) + '\n')

    f.close()
    out.close()

    output_contaminant_removal_statistics(
        total_reads,
        trimmed_reads,
        deleted_reads,
        contaminants_found,
        log_file_handle,
        )

    return outname
