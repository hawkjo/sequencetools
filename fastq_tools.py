import sys, os, re
import numpy as np
from itertools import izip

def determine_phred_offset(filename):
    min_val = 126
    max_val = 0
    with open(filename) as f:
        while True:
            nameline = f.readline()
            if not nameline: sys.exit('Could not determine phred offset')
            seqline = f.readline()
            plusline = f.readline()
            qualline = f.readline()

            ascii_vals = map(ord, qualline.strip()) # Convert to numerical values
            min_val = min([min_val]+ascii_vals)
            max_val = max([max_val]+ascii_vals)
            if min_val < 59: return 33 # Illumina 1.8 and Sanger
            elif max_val > 80: return 64 # Illumina 1.3 - 1.7 and Solexa

def trim_low_quality_by_bookcase_good_reads(filename1, filename2=None):
    qual_cutoff = 30
    phred_offset = determine_phred_offset(filename1)
    coded_qual_cutoff = qual_cutoff + phred_offset
    if filename2 is None:
        out_filename = os.path.splitext(filename1)[0] + '_trimmed.fastq'
        with open(out_filename,'w') as out:
            with open(filename1) as infile:
                while True:
                    nameline = infile.readline()
                    if not nameline: break
                    seqline = infile.readline()
                    plusline = infile.readline()
                    qualline = infile.readline()
        
                    high_qual_inds = [i for i, c in enumerate(qualline.strip()) if ord(c) >= coded_qual_cutoff]

                    if not high_qual_inds: continue
                    start = high_qual_inds[0]
                    end = high_qual_inds[-1] + 1

                    if end - start < 25: continue
                    out.write(nameline)
                    out.write(seqline[start:end] + '\n')
                    out.write(plusline)
                    out.write(qualline[start:end] + '\n')
    else:
        out_filename1 = os.path.splitext(filename1)[0] + '_trimmed.fastq'
        out_filename2 = os.path.splitext(filename2)[0] + '_trimmed.fastq'

        out = [open(out_filename1,'w'), open(out_filename2,'w')]
        infile = [open(filename1), open(filename2)]

        while True:
            nameline = [infile[0].readline(), infile[1].readline()]
            if not (nameline[0] and nameline[1]): 
                if nameline[0] or nameline[1]:
                    sys.exit('Error: Unequal number of reads')
                else: break
            seqline = [infile[0].readline(), infile[1].readline()]
            plusline = [infile[0].readline(), infile[1].readline()]
            qualline = [infile[0].readline(), infile[1].readline()]
        
            # Accept part of read starting from the first position in which both reads are higher
            # than 30 and end on the last position where both are higher than 30.
            high_qual_inds = [i for i, (c1,c2) in enumerate( izip(qualline[0].strip(), 
                    qualline[1].strip()) ) if min(ord(c1),ord(c2)) >= coded_qual_cutoff]

            if not high_qual_inds: continue
            start = high_qual_inds[0]
            end = high_qual_inds[-1] + 1

            if end - start < 25: continue

            for i in [0,1]:
                out[i].write(nameline[i])
                out[i].write(seqline[i][start:end] + '\n')
                out[i].write(plusline[i])
                out[i].write(qualline[i][start:end] + '\n')

        for f in out + infile: f.close()

def bookend_qual_subread(qualline, phred_offset, qual_cutoff, min_len=25):
    """
    Accepts part of read starting from the first position in which read (or both reads) are higher
    than 30 and ending on the last position where both are higher than 30.
    """
    if type(qualline) is str:
        high_qual_inds = [i for i, c in enumerate(qualline.strip()) if ord(c) >= phred_offset + qual_cutoff]
    elif type(qualline) is list:
        if len(qualline) != 2: 
            sys.exit('Unexpected qualline input to bookend_qual_subreads: ' + str(qualline))
        high_qual_inds = [i for i, (c1,c2) in enumerate( izip(qualline[0].strip(), 
                qualline[1].strip()) ) if min(ord(c1),ord(c2)) >= phred_offset + qual_cutoff]
    else: sys.exit('Unexpected qualline input to bookend_qual_subreads: ' + str(qualline))

    if not high_qual_inds: return None
    start = high_qual_inds[0]
    end = high_qual_inds[-1] + 1
    if end - start < min_len: return None

    return (start, end)

def longest_subread_above_thresh(qualline, phred_offset, qual_cutoff, min_len=25):
    """
    Finds the longest subread where all bases are of quality above the cutoff, and returns the start
    and end indices. If shorter than min_len, returns None.
    """
    if type(qualline) is str:
        low_qual_inds = [i for i, c in enumerate(qualline.strip()) if ord(c) < phred_offset + qual_cutoff]
    elif type(qualline) is list:
        if len(qualline) != 2: 
            sys.exit('Unexpected qualline input to bookend_qual_subreads: ' + str(qualline))
        low_qual_inds = [i for i, (c1,c2) in enumerate( izip(qualline[0].strip(), 
                qualline[1].strip()) ) if min(ord(c1),ord(c2)) < phred_offset + qual_cutoff]

    # Add index before first and after last as 'bad bases' to consider as non-included end-points
    low_qual_inds.insert(0,-1)
    low_qual_inds.append( len(qualline) )

    max_len = 0
    for ind_i, ind_ip1 in izip(low_qual_inds[:-1], low_qual_inds[1:]):
        if ind_ip1 - ind_i > max_len:
            start = ind_i + 1 # Half-open, not including bad bases
            end = ind_ip1 
            max_len = end - start
    if max_len < min_len: return None
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
    # We then zero out the scores of all subreads along the first [min_len-1] diagonals, corresponding
    # to the reads shorter than the min length, thus enforcing the minimum length, since only
    # positive scores are accepted.

    if type(qualline) is str:
        base_scores = [ord(c) - phred_offset - score_zero for c in qualline]
    elif type(qualline) is list:
        if len(qualline) != 2: 
            sys.exit('Unexpected qualline input to max_score_subreads: ' + str(qualline))
        base_scores = [ord(c1) + ord(c2) - 2*(phred_offset + score_zero) \
                for c1,c2 in izip(qualline[0],qualline[1])]
    else: sys.exit('Unexpected qualline input to max_score_subreads: ' + str(qualline))

    score_mat = np.zeros((len(qualline),len(qualline)))
    for i in xrange(len(base_scores)): # Calculate scores
        score_mat[0:i+1, i:len(base_scores)] += base_scores[i]
    for i in xrange(len(base_scores)): # Enforce min length
        score_mat[i, i:i+min_len-1] = 0
    max_score = score_mat.max()
    if max_score <= 0: return None

    max_bases = 0
    start = 0
    end = 0
    max_inds = np.where(score_mat == max_score)
    for a,b in izip(max_inds[0], max_inds[1]):
        tmp_len = b-a+1
        if tmp_len > max_bases:
            max_bases = tmp_len
            start = a
            end = b+1
    return (start,end)

def trim_low_quality_bases(filename1, filename2=None, subread_func=bookend_qual_subread, cutoff=10, min_len=25):
    """
    trim_low_quality_by_score assigns a score to each phred score simply by subtracting 'cutoff'
    from the phred score. We then find the max-score length of read.
    """
    phred_offset = determine_phred_offset(filename1)
    if filename2 is None:
        out_filename = os.path.splitext(filename1)[0] + '_%s_%s.fastq' % \
                (subread_func.__name__, cutoff)
        with open(out_filename,'w') as out:
            with open(filename1) as infile:
                while True:
                    nameline = infile.readline().strip()
                    if not nameline: break
                    seqline = infile.readline().strip()
                    plusline = infile.readline().strip()
                    qualline = infile.readline().strip()
        
                    subread_ends = subread_func(qualline, phred_offset, cutoff, min_len)

                    if not subread_ends: continue
                    start, end = subread_ends 

                    out.write(nameline + '\n')
                    out.write(seqline[start:end] + '\n')
                    out.write(plusline + '\n')
                    out.write(qualline[start:end] + '\n')
    else:
        out_filename1 = os.path.splitext(filename1)[0] + '_%s_%s.fastq' % \
                (subread_func.__name__, cutoff)
        out_filename2 = os.path.splitext(filename2)[0] + '_%s_%s.fastq' % \
                (subread_func.__name__, cutoff)

        out = [open(out_filename1,'w'), open(out_filename2,'w')]
        infile = [open(filename1), open(filename2)]

        while True:
            nameline = [infile[0].readline(), infile[1].readline()]
            if not (nameline[0] and nameline[1]): 
                if nameline[0] or nameline[1]:
                    sys.exit('Error: Unequal number of reads')
                else: break
            seqline = [infile[0].readline(), infile[1].readline()]
            plusline = [infile[0].readline(), infile[1].readline()]
            qualline = [infile[0].readline(), infile[1].readline()]
        
            # Inputting list of two quallines combines the scores at each position
            subread_ends = subread_func(qualline, phred_offset, cutoff, min_len)

            if not subread_ends: continue
            start, end = subread_ends 

            for i in [0,1]:
                out[i].write(nameline[i])
                out[i].write(seqline[i][start:end] + '\n')
                out[i].write(plusline[i])
                out[i].write(qualline[i][start:end] + '\n')

        for f in out + infile: f.close()

def replace_bad_bases_with_N(filename, cutoff=20):
    """
    Replace all bases below the given cutoff (default 20) with 'N'. Outputs to new file.
    """
    phred_offset = determine_phred_offset(filename)
    coded_cutoff = cutoff + phred_offset
    out_filename = os.path.splitext(filename)[0] + '_Nified_%s.fastq' % cutoff
    with open(out_filename,'w') as out:
        with open(filename) as infile:
            while True:
                nameline = infile.readline().strip()
                if not nameline: break
                seqline = infile.readline().strip()
                plusline = infile.readline().strip()
                qualline = infile.readline().strip()
    
                seq = list(seqline)
                for i, c in enumerate(qualline):
                    if ord(c) < coded_cutoff:
                        seq[i] = 'N'

                out.write(nameline + '\n')
                out.write(''.join(seq) + '\n')
                out.write(plusline + '\n')
                out.write(qualline + '\n')

GSAF_fastq_id_line_re = re.compile(r"""
        ^@                              # Starts with @ sign
        (?P<instrument>[-A-Z0-9]+):     # Instrument name
        (?P<run>\d+):                   # Run id
        (?P<flowcell_id>[A-Z0-9]+):     # Flowcell id
        (?P<flowcell_lane>\d+):         # Flowcell lane
        (?P<tile>\d+):                  # tile number
        (?P<xcoord>\d+):                # x-coord
        (?P<ycoord>\d+)                 # y-coord
        \                               # space
        (?P<pair_member>[12]):          # member of pair (1 or 2)
        (?P<fail_status>[YN]):          # Y(es) if failed, N(o) otherwise
        (?P<control_bits>\d+):          # Control bits. 0 when all off
        (?P<barcode>[ACGT]+)            # barcode
    """, re.VERBOSE)

def get_GSAF_barcode(filename1,filename2=None):
    barcodes = []
    for fn in [filename1,filename2]:
        if fn == None: continue
        with open(fn) as f:
            line = f.readline().strip()
            m = GSAF_fastq_id_line_re.match(line)
            barcodes.append(m.group('barcode'))
    if len(barcodes) == 2 and barcodes[0] != barcodes[1]:
        sys.exit('Error: Different barcodes in read files')
    return barcodes[0]
