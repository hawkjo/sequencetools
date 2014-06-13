import sys, re

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
