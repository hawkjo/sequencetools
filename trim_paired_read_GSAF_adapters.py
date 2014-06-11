import sys, re
from adapters import build_adapters, paired_end_adapter_position

fastq_id_line_re = re.compile(r"""
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

def get_barcode(filename1,filename2):
    barcodes = []
    for fn in [filename1,filename2]:
        with open(fn) as f:
            line = f.readline().strip()
            m = fastq_id_line_re.match(line)
            barcodes.append(m.group('barcode'))
    if barcodes[0] != barcodes[1]:
        sys.exit('Error: Different barcodes in read files')
    return barcodes[0]

if __name__ == "__main__":
    import os
    if len(sys.argv) != 3:
        sys.exit('Usage: trim_paired_read_adapters.py <fastq 1> <fastq 2>')
    
    filename1 = sys.argv[1]
    filename2 = sys.argv[2]
    outname1 = os.path.splitext(filename1)[0] + '_trimmed.fastq'
    outname2 = os.path.splitext(filename2)[0] + '_trimmed.fastq'

    barcode = get_barcode(filename1,filename2)

    len_adapter_match = 13 # Length of identical truseq prefix
    max_mismatches = 3 # Max hamming distance
    adapter_in_R1, adapter_in_R2 = build_adapters(index_sequence=barcode, max_length=len_adapter_match)

    f1 = open(filename1)
    f2 = open(filename2)
    o1 = open(outname1,'w')
    o2 = open(outname2,'w')

    total_reads = 0
    trimmed_reads = 0
    deleted_reads = 0
    while True:
        id_line1 = f1.readline().strip()
        if not id_line1: break # End of file
        total_reads += 1

        seq_line1 = f1.readline().strip()
        plus_line1 = f1.readline().strip()
        q_line1 = f1.readline().strip()

        id_line2 = f2.readline().strip()
        seq_line2 = f2.readline().strip()
        plus_line2 = f2.readline().strip()
        q_line2 = f2.readline().strip()

        adapter_position = position(seq_line1,
                                    seq_line2,
                                    adapter_in_R1,
                                    adapter_in_R2,
                                    len_adapter_match,
                                    max_mismatches)

        if adapter_position != None:
            if adapter_position == 0:
                deleted_reads += 1
                continue
            else:
                trimmed_reads += 1
                seq_line1 = seq_line1[:adapter_position]
                q_line1 = q_line1[:adapter_position]
                seq_line2 = seq_line2[:adapter_position]
                q_line2 = q_line2[:adapter_position]

        o1.write('\n'.join([id_line1, seq_line1, plus_line1, q_line1]) + '\n')
        o2.write('\n'.join([id_line2, seq_line2, plus_line2, q_line2]) + '\n')

    f1.close()
    f2.close()
    o1.close()
    o2.close()

    print 'Total reads:,',total_reads
    print 'Trimmed reads: %d, %f%%' % (trimmed_reads, 100*float(trimmed_reads)/total_reads)
    print 'Deleted reads: %d, %f%%' % (deleted_reads, 100*float(deleted_reads)/total_reads)
