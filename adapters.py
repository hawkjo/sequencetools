from adapters_cython import *

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
