import sys
import os
from collections import Counter
from adapters_cython import find_adapter_positions
from fastqtools import get_GSAF_barcode
from contaminants import get_fastqc_contaminant_list, output_contaminant_removal_statistics

"""
Let's lay out a schematic of the full read. The pieces are:

    R[12]:  Read 1/2
    SP[12]: Sequence Primer 1/2
    P[57]:  P5/7 adapters
    Index:  Demultiplexing index
    dA/dT:  Poly-dA/T tails, dA only if readthrough
    Ns:     Unread bases
    *:      Reverse Complement
    |>:     Wall attachment during amplification

Then the read structure is:

        dT      P5      SP1     R1      Ns      R2*     SP2*    Index   P7*     dA      
    |> ------- ------- ------- ------- ------- ------- ------- ------- ------- ----
          ---- ------- ------- ------- ------- ------- ------- ------- ------- ------- <|
        dA      P5*     SP1*    R1*     Ns*     R2      SP2     Index*  P7      dT
"""

tru_seq_SP1_rc = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGA'
tru_seq_SP2_rc = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'

# Previous generation of adapters
paired_end_SP1_rc = tru_seq_SP1_rc
paired_end_SP2_rc = 'AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG'

P5_rc = 'TCTCGGTGGTCGCCGTATCATT'
P7_rc = 'ATCTCGTATGCCGTCTTCTGCTTG'

A_tail = 'A' * 10


def build_adapters(index_sequence='', max_length=None, primer_type='tru_seq'):
    if primer_type == 'tru_seq':
        SP1_rc = tru_seq_SP1_rc
        SP2_rc = tru_seq_SP2_rc
    elif primer_type == 'PE':
        SP1_rc = paired_end_SP1_rc
        SP2_rc = paired_end_SP2_rc
    else:
        raise ValueError('Invalid primer type: {0}'.format(primer_type))
    adapter_in_R1 = SP2_rc + index_sequence + P7_rc + A_tail
    adapter_in_R2 = SP1_rc + P5_rc + A_tail
    truncated_slice = slice(None, max_length)
    return adapter_in_R1[truncated_slice], adapter_in_R2[truncated_slice]


def single_end_adapter_position(read_seq,
                                adapter_seq,
                                min_comparison_length,
                                max_distance):
    read_positions = find_adapter_positions(read_seq,
                                            adapter_seq,
                                            min_comparison_length,
                                            max_distance)
    read_positions = set(read_positions)
    if read_positions:
        return min(read_positions)
    else:
        return None


def paired_end_adapter_position(R1_seq,
                                R2_seq,
                                adapter_in_R1,
                                adapter_in_R2,
                                min_comparison_length,
                                max_distance):
    R1_positions = find_adapter_positions(R1_seq,
                                          adapter_in_R1,
                                          min_comparison_length,
                                          max_distance)
    R2_positions = find_adapter_positions(R2_seq,
                                          adapter_in_R2,
                                          min_comparison_length,
                                          max_distance)

    R1_positions = set(R1_positions)
    R2_positions = set(R2_positions)
    common_positions = R1_positions & R2_positions
    if common_positions:
        return min(common_positions)
    else:
        return None


def get_fastqc_paired_read_adapters(fname1, fname2, max_comparison_length=16):
    adapter_file = '/home/hawkjo/python_src/sequence_tools/adapter_list.txt'
    adapter_list = get_fastqc_contaminant_list(adapter_file,
                                               max_contaminant_length=max_comparison_length,
                                               include_rev_comp=False)
    adapter_tuples = [(name, adapter, adapter) for name, adapter in adapter_list]
    return adapter_tuples


def get_GSAF_paired_read_adapter_tuples(fname1, fname2, max_comparison_length=16):
    barcode = get_GSAF_barcode(fname1, fname2)
    adapter_in_R1, adapter_in_R2 = build_adapters(index_sequence=barcode,
                                                  max_length=max_comparison_length)
    adapter_tuples = [('GSAF (TruSeq)', adapter_in_R1, adapter_in_R2)]
    return adapter_tuples


def trim_paired_read_adapters(fname1, fname2,  # Files to trim
                              adapter_tuples_func=get_fastqc_paired_read_adapters,
                              min_read_len=25,           # Minimum length of output read
                              min_comparison_length=9,   # Minimum length to compare
                              max_comparison_length=16,  # Maximum length to compare
                              max_mismatches=2,          # Max hamming distance
                              log_file_handle=sys.stdout):

    # Adapter_tuples_func generates a list of tuples:
    #       (adapter_name, adapter_in_R1, adapter_in_R2)
    adapter_tuples = adapter_tuples_func(fname1, fname2, max_comparison_length)

    if 'GSAF' in adapter_tuples_func.__name__:
        adapter_type = 'GSAF'
    elif 'fastqc' in adapter_tuples_func.__name__:
        adapter_type = 'fastqc'
    else:
        adapter_type = adapter_tuples_func.__name__

    outname1 = os.path.splitext(fname1)[0] + '_deadaptered_%s_%d-%d_%d.fastq' \
        % (adapter_type, min_comparison_length, max_comparison_length, max_mismatches)
    outname2 = os.path.splitext(fname2)[0] + '_deadaptered_%s_%d-%d_%d.fastq' \
        % (adapter_type, min_comparison_length, max_comparison_length, max_mismatches)

    f1 = open(fname1)
    f2 = open(fname2)
    o1 = open(outname1, 'w')
    o2 = open(outname2, 'w')

    total_reads = 0
    trimmed_reads = 0
    deleted_reads = 0
    adapters_found = Counter()
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

        for adapter_name, adapter_in_R1, adapter_in_R2 in adapter_tuples:
            adapter_position = paired_end_adapter_position(seqline1,
                                                           seqline2,
                                                           adapter_in_R1,
                                                           adapter_in_R2,
                                                           min_comparison_length,
                                                           max_mismatches)

            if adapter_position is not None:
                adapters_found[adapter_name] += 1
                if adapter_position < min_read_len:
                    deleted_reads += 1
                    break
                else:
                    trimmed_reads += 1
                    seqline1 = seqline1[:adapter_position]
                    qualline1 = qualline1[:adapter_position]
                    seqline2 = seqline2[:adapter_position]
                    qualline2 = qualline2[:adapter_position]

        else:
            # Did not delete the read
            o1.write('\n'.join([defline1, seqline1, plusline1, qualline1]) + '\n')
            o2.write('\n'.join([defline2, seqline2, plusline2, qualline2]) + '\n')

    f1.close()
    f2.close()
    o1.close()
    o2.close()

    output_contaminant_removal_statistics(
        total_reads,
        trimmed_reads,
        deleted_reads,
        adapters_found,
        log_file_handle,
        contaminant_label='Adapter')

    return (outname1, outname2)


def trim_single_read_adapters(
        fname,                       # File to trim
        min_read_len=25,           # Minimum length of output read
        min_comparison_length=12,  # Minimum length to compare
        max_comparison_length=16,  # Maximum length to compare
        max_mismatches=1,          # Max hamming distance
        log_file_handle=sys.stdout):
    # Default min/max comparison len values of 12/16 have 2e-6 and 1e-8 probability of randomly
    # appearing within distance 1 of a given sequence.
    adapter_file = '/home/hawkjo/python_src/sequence_tools/adapter_list.txt'
    adapter_list = get_fastqc_contaminant_list(adapter_file,
                                               max_contaminant_length=max_comparison_length,
                                               include_rev_comp=False)
    outname = os.path.splitext(fname)[0] + '_deadaptered_fastqc_%d-%d_%d.fastq' \
        % (min_comparison_length, max_comparison_length, max_mismatches)
    f = open(fname)
    out = open(outname, 'w')

    total_reads = 0
    trimmed_reads = 0
    deleted_reads = 0
    adapters_found = Counter()
    while True:
        defline = f.readline().strip()
        if not defline:
            break  # End of file
        total_reads += 1

        seqline = f.readline().strip()
        plusline = f.readline().strip()
        qualline = f.readline().strip()

        for adapter_name, adapter_seq in adapter_list:
            adapter_position = single_end_adapter_position(seqline,
                                                           adapter_seq,
                                                           min_comparison_length,
                                                           max_mismatches)

            if adapter_position is not None:
                adapters_found[adapter_name] += 1
                if adapter_position < min_read_len:
                    deleted_reads += 1
                    break
                else:
                    trimmed_reads += 1
                    seqline = seqline[:adapter_position]
                    qualline = qualline[:adapter_position]
        else:
            # Non-deleted read
            out.write('\n'.join([defline, seqline, plusline, qualline]) + '\n')

    f.close()
    out.close()

    output_contaminant_removal_statistics(
        total_reads,
        trimmed_reads,
        deleted_reads,
        adapters_found,
        log_file_handle,
        contaminant_label='Adapter')

    return outname
