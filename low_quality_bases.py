import sys
import os
import numpy as np
from itertools import izip
from fastq_tools import determine_phred_offset
from contaminants import output_contaminant_removal_statistics


def bookend_qual_subread(qualline, phred_offset, qual_cutoff, min_len=25):
    """
    Accepts part of read starting from the first position in which read (or both reads) are higher
    than 30 and ending on the last position where both are higher than 30.
    """
    if type(qualline) is str:
        high_qual_inds = [i for i, c in enumerate(qualline.strip())
                          if ord(c) >= phred_offset + qual_cutoff]
    elif type(qualline) is list:
        if len(qualline) != 2:
            sys.exit('Unexpected qualline input to bookend_qual_subreads: ' + str(qualline))
        high_qual_inds = [i for i, (c1, c2) in enumerate(izip(qualline[0].strip(),
                                                              qualline[1].strip()))
                          if min(ord(c1), ord(c2)) >= phred_offset + qual_cutoff]
    else:
        sys.exit('Unexpected qualline input to bookend_qual_subreads: ' + str(qualline))

    if not high_qual_inds:
        return None
    start = high_qual_inds[0]
    end = high_qual_inds[-1] + 1
    if end - start < min_len:
        return None

    return (start, end)


def longest_subread_above_thresh(qualline, phred_offset, qual_cutoff, min_len=25):
    """
    Finds the longest subread where all bases are of quality above the cutoff, and returns the
    start and end indices. If shorter than min_len, returns None.
    """
    if type(qualline) is str:
        low_qual_inds = [i for i, c in enumerate(qualline.strip())
                         if ord(c) < phred_offset + qual_cutoff]
    elif isinstance(qualline, list):
        if len(qualline) != 2:
            sys.exit('Unexpected qualline input to bookend_qual_subreads: ' + str(qualline))
        low_qual_inds = [i for i, (c1, c2) in enumerate(izip(qualline[0].strip(),
                                                             qualline[1].strip()))
                         if min(ord(c1), ord(c2)) < phred_offset + qual_cutoff]

    # Add index before first and after last as 'bad bases' to consider as non-included end-points
    low_qual_inds.insert(0, -1)
    low_qual_inds.append(len(qualline))

    max_len = 0
    for ind_i, ind_ip1 in izip(low_qual_inds[:-1], low_qual_inds[1:]):
        if ind_ip1 - ind_i > max_len:
            start = ind_i + 1  # Half-open, not including bad bases
            end = ind_ip1
            max_len = end - start
    if max_len < min_len:
        return None
    return (start, end)


def max_score_subread(qualline, phred_offset, score_zero, min_len=25):
    """
    Finds the max-score subread given a fastq score line, phred offset, and phred value to set as
    zero. Optional argument sets the min acceptable length (default 25).

    If qualline is a tuple of two quallines, the scores are added together before processing.
    """
    # We create an upper diagonal score matrix which stores all possible contiguous subread scores.
    # This is easily accomplished by adding the score of the ith base to the submatrix in the
    # upper-right submatrix from the (i,i)th element of the score matrix.
    #
    # We then zero out the scores of all subreads along the first [min_len-1] diagonals,
    # corresponding to the reads shorter than the min length, thus enforcing the minimum length,
    # since only positive scores are accepted.

    if type(qualline) is str:
        base_scores = [ord(c) - phred_offset - score_zero for c in qualline]
    elif type(qualline) is list:
        if len(qualline) != 2:
            sys.exit('Unexpected qualline input to max_score_subreads: ' + str(qualline))
        base_scores = [ord(c1) + ord(c2) - 2*(phred_offset + score_zero)
                       for c1, c2 in izip(qualline[0], qualline[1])]
    else:
        sys.exit('Unexpected qualline input to max_score_subreads: ' + str(qualline))

    score_mat = np.zeros((len(qualline), len(qualline)))
    for i in xrange(len(base_scores)):  # Calculate scores
        score_mat[0:i+1, i:len(base_scores)] += base_scores[i]
    for i in xrange(len(base_scores)):  # Enforce min length
        score_mat[i, i:i+min_len-1] = 0
    max_score = score_mat.max()
    if max_score <= 0:
        return None

    max_bases = 0
    start = 0
    end = 0
    max_inds = np.where(score_mat == max_score)
    for a, b in izip(max_inds[0], max_inds[1]):
        tmp_len = b-a+1
        if tmp_len > max_bases:
            max_bases = tmp_len
            start = a
            end = b+1
    return (start, end)


def trim_single_read_low_quality_bases(
        fname,
        subread_func=bookend_qual_subread,      # Function which finds where to trim
        cutoff=10,                              # Quality cutoff, used differently in different fns
        min_len=25,                             # Min length of output read
        log_file_handle=sys.stdout,
        ):
    """
    trim_single_read_low_quality_bases trims off the low quality bases according to the
    subread_func and cutoff given.
    """
    phred_offset = determine_phred_offset(fname)
    outname = os.path.splitext(fname)[0] + '_%s_%s.fastq' % \
        (subread_func.__name__, cutoff)

    total_reads = 0
    trimmed_reads = 0
    deleted_reads = 0
    with open(outname, 'w') as out:
        with open(fname) as infile:
            while True:
                defline = infile.readline().strip()
                if not defline:
                    break
                total_reads += 1

                seqline = infile.readline().strip()
                plusline = infile.readline().strip()
                qualline = infile.readline().strip()

                subread_ends = subread_func(qualline, phred_offset, cutoff, min_len)
                if not subread_ends:
                    deleted_reads += 1
                    continue
                start, end = subread_ends
                if start != 0 or end != len(seqline):
                    trimmed_reads += 1

                out.write(defline + '\n')
                out.write(seqline[start:end] + '\n')
                out.write(plusline + '\n')
                out.write(qualline[start:end] + '\n')

    output_contaminant_removal_statistics(
        total_reads,
        trimmed_reads,
        deleted_reads,
        contaminants_found=None,
        log_file_handle=log_file_handle,
        contaminant_label='Low Quality Base',
        )


def trim_paired_read_low_quality_bases(
        fname1,
        fname2,
        subread_func=bookend_qual_subread,      # Function which finds where to trim
        cutoff=10,                              # Quality cutoff, used differently in different fns
        min_len=25,                             # Min length of output read
        log_file_handle=sys.stdout,
        ):
    """
    trim_paired_read_low_quality_bases trims off the low quality bases according to the
    subread_func and cutoff given.

    If one of a pair of reads is deleted, the remaining read is stored as an
    orphaned read and placed at the end of the left reads file.
    """
    phred_offset = determine_phred_offset(fname1)
    outname1 = os.path.splitext(fname1)[0] + '_%s_%s.fastq' % \
        (subread_func.__name__, cutoff)
    outname2 = os.path.splitext(fname2)[0] + '_%s_%s.fastq' % \
        (subread_func.__name__, cutoff)

    f1 = open(fname1)
    f2 = open(fname2)
    o1 = open(outname1, 'w')
    o2 = open(outname2, 'w')

    total_reads = 0
    trimmed_reads = 0
    deleted_reads = 0
    orphaned_reads = []
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

        subread_ends1 = subread_func(qualline1, phred_offset, cutoff, min_len)
        subread_ends2 = subread_func(qualline2, phred_offset, cutoff, min_len)

        if subread_ends1 and subread_ends2:
            # Do not delete either read. Output (trimmed) reads.
            start, end = subread_ends1
            if start != 0 or end != len(seqline1):
                trimmed_reads += 1
                seqline1 = seqline1[start:end]
                qualline1 = qualline1[start:end]

            start, end = subread_ends2
            if start != 0 or end != len(seqline2):
                trimmed_reads += 1
                seqline2 = seqline2[start:end]
                qualline2 = qualline2[start:end]

            o1.write('\n'.join([defline1, seqline1, plusline1, qualline1]) + '\n')
            o2.write('\n'.join([defline2, seqline2, plusline2, qualline2]) + '\n')

        elif subread_ends1:
            # Delete read 2, but keep read 1 as an orphaned read
            deleted_reads += 1
            start, end = subread_ends1
            if start != 0 or end != len(seqline1):
                trimmed_reads += 1
                seqline1 = seqline1[start:end]
                qualline1 = qualline1[start:end]
            orphaned_reads.append('\n'.join([defline1, seqline1, plusline1, qualline1]))

        elif subread_ends2:
            # Delete read 1, but keep read 2 as an orphaned read
            deleted_reads += 1
            start, end = subread_ends2
            if start != 0 or end != len(seqline2):
                trimmed_reads += 1
                seqline2 = seqline2[start:end]
                qualline2 = qualline2[start:end]
            orphaned_reads.append('\n'.join([defline2, seqline2, plusline2, qualline2]))

        else:
            # Delete both reads
            deleted_reads += 2

    # Output the orphaned reads to file 1
    o1.write('\n'.join(orphaned_reads))

    f1.close()
    f2.close()
    o1.close()
    o2.close()

    output_contaminant_removal_statistics(
        total_reads,
        trimmed_reads,
        deleted_reads,
        contaminants_found=None,
        log_file_handle=log_file_handle,
        contaminant_label='Low Quality Base',
        )
    log_file_handle.write('Orphaned reads: %d\n' % (len(orphaned_reads)))


def replace_bad_bases_with_N(filename, cutoff=20):
    """
    Replace all bases below the given cutoff (default 20) with 'N'. Outputs to new file.
    """
    phred_offset = determine_phred_offset(filename)
    coded_cutoff = cutoff + phred_offset
    outname = os.path.splitext(filename)[0] + '_Nified_%s.fastq' % cutoff
    with open(outname, 'w') as out:
        with open(filename) as infile:
            while True:
                defline = infile.readline().strip()
                if not defline:
                    break
                seqline = infile.readline().strip()
                plusline = infile.readline().strip()
                qualline = infile.readline().strip()

                seq = list(seqline)
                for i, c in enumerate(qualline):
                    if ord(c) < coded_cutoff:
                        seq[i] = 'N'

                out.write(defline + '\n')
                out.write(''.join(seq) + '\n')
                out.write(plusline + '\n')
                out.write(qualline + '\n')
