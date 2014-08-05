import sys, re, random, os
import numpy as np
from misc_tools import gzip_friendly_open

def determine_phred_offset(filename, num_reads_to_consider=1000):
    min_val = 126
    max_val = 0

    i = 0
    with gzip_friendly_open(filename) as f:
        while True:
            i += 1
            if i > num_reads_to_consider: break 

            nameline = f.readline()
            if not nameline: sys.exit('Could not determine phred offset')
            seqline = f.readline()
            plusline = f.readline()
            qualline = f.readline()

            ascii_vals = map(ord, qualline.strip()) # Convert to numerical values
            min_val = min([min_val]+ascii_vals)
            max_val = max([max_val]+ascii_vals)

    if min_val < 50:
        return 33 # Illumina 1.8 and Sanger
    elif max_val > 89: 
        return 64 # Illumina 1.3 - 1.7 and Solexa
    else:
        return determine_phred_offset(filename, num_reads_to_consider= 2*num_reads_to_consider)

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

def make_add_errors_to_seq_func( phred_offset ):
    p_given_coded_phred = np.zeros(150)
    for coded_Q in xrange(phred_offset, phred_offset + 42):
        Q = coded_Q - phred_offset
        p_given_coded_phred[coded_Q] = 10**(-0.1*Q)
        if p_given_coded_phred[coded_Q] > 1:
            sys.exit('Error probability greater than 1: p=%f, Q=%d, phred_offset=%d' % \
                        (p, Q, phred_offset) )

    def add_errors_to_seq( seq, qualline ):
        erroneous_seq = list(seq)
        for j in xrange(len(seq)):
            if random.random() < p_given_coded_phred[ord(qualline[j])]:
                erroneous_seq[j] = random.choice( 'ACGT'.replace( erroneous_seq[j], '') )
        return ''.join( erroneous_seq )

    return add_errors_to_seq

def make_add_errors_to_seq_func_given_example_file( fname ):
    phred_offset = determine_phred_offset( fname )
    return make_add_errors_to_seq_func( phred_offset )

def convert_phred_scores( fname, out_phred_offset ):
    if out_phred_offset not in [33,64]:
        sys.exit('Error: out_phred_offset must be 33 or 64. Received %s' % repr(out_phred_offset) )

    in_phred_offset = determine_phred_offset( fname )
    if in_phred_offset == out_phred_offset:
        print 'Cowardly refusing to convert %s from phred%d to phred%d' \
                % (fname, in_phred_offset, out_phred_offset) 
        return -1

    phred_diff = out_phred_offset - in_phred_offset

    fname_parts = fname.split('.')
    out_fname = fname_parts[0] + '_phred' + str(out_phred_offset) + '.' + '.'.join( fname_parts[1:] )

    from misc_tools import gzip_friendly_open
    with gzip_friendly_open( fname ) as f, gzip_friendly_open( out_fname, 'w' ) as out:
        while True:
            defline = f.readline().strip()
            if not defline: break
            seqline = f.readline().strip()
            plusline = f.readline().strip()
            qualline = f.readline().strip()

            out_qualline = ''.join( [chr(ord(c) + phred_diff) for c in qualline] )

            out.write( '\n'.join([defline, seqline, plusline, out_qualline]) + '\n' )
