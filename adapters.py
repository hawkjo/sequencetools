import sys, os
from adapters_cython import *
from fastq_tools import get_GSAF_barcode

tru_seq_R1_rc = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGA'
tru_seq_R2_rc = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'

# Previous generation of adapters
paired_end_R1_rc = tru_seq_R1_rc
paired_end_R2_rc = 'AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG'

P5_rc = 'TCTCGGTGGTCGCCGTATCATT'
P7_rc = 'ATCTCGTATGCCGTCTTCTGCTTG'

A_tail = 'A' * 10

def build_adapters(index_sequence='', max_length=None):
    adapter_in_R1 = tru_seq_R2_rc + index_sequence + P7_rc + A_tail
    adapter_in_R2 = tru_seq_R1_rc + P5_rc + A_tail
    truncated_slice = slice(None, max_length)
    return adapter_in_R1[truncated_slice], adapter_in_R2[truncated_slice]

def single_end_adapter_position(read_seq,
             adapter_seq,
             min_comparison_length,
             max_distance,
            ):
    read_positions = find_adapter_positions(read_seq, adapter_seq, min_comparison_length, max_distance)

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
             max_distance,
            ):
    R1_positions = find_adapter_positions(R1_seq, adapter_in_R1, min_comparison_length, max_distance)
    R2_positions = find_adapter_positions(R2_seq, adapter_in_R2, min_comparison_length, max_distance)

    R1_positions = set(R1_positions)
    R2_positions = set(R2_positions)
    common_positions = R1_positions & R2_positions
    if common_positions:
        return min(common_positions)
    else:
        return None

def get_contaminant_list(max_contaminant_length=None):
    contaminant_file = '/home/hawkjo/python_src/sequence_tools/contaminant_list.txt'
    contaminant_list = []
    for line in open(contaminant_file):
        if line[0] == '#': continue

        var = line.strip().split('\t')
        name = var[0]
        seq = var[-1]

        if not name: continue # Blank lines
        if 'A' not in seq and 'C' not in seq and 'G' not in seq and 'T' not in seq:
            sys.exit('Non-DNA contaminant, %s: %s' % (name, seq) )

        contaminant_list.append( (name, seq[:max_contaminant_length]) )

    if not contaminant_list:
        sys.exit('No contaminants found in ' + contaminant_file)

    return contaminant_list

def trim_GSAF_paired_read_adapters(filename1, filename2):
    outname1 = os.path.splitext(filename1)[0] + '_trimmed.fastq'
    outname2 = os.path.splitext(filename2)[0] + '_trimmed.fastq'

    barcode = get_GSAF_barcode(filename1,filename2)

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

        adapter_position = paired_end_adapter_position(seq_line1,
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

    return (outname1, outname2)
