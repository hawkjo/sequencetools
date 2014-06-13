import sys, os
from adapters_cython import *
from fastq_tools import get_GSAF_barcode
from contaminants import get_fastqc_contaminant_list
from collections import Counter

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

def get_fastqc_paired_read_adapters(fname1, fname2, max_comparison_length = 16 ):
    adapter_list = get_fastqc_contaminant_list( 'adapter_list.txt', 
            max_contaminant_length=max_comparison_length,
            include_rev_comp=False)
    adapter_tuples = [(name, adapter, adapter) for name, adapter in adapter_list]
    return adapter_tuples

def get_GSAF_paired_read_adapter_tuples(fname1, fname2, max_comparison_length=16):
    barcode = get_GSAF_barcode(fname1,fname2)
    adapter_in_R1, adapter_in_R2 = build_adapters(index_sequence=barcode,
            max_length=max_comparison_length)
    adapter_tuples = [('GSAF (TruSeq)', adapter_in_R1, adapter_in_R2)]
    return adapter_tuples

def trim_paired_read_adapters(fname1, fname2,  # Files to trim
        adapter_tuples_func = get_fastqc_paired_read_adapters,        
        min_read_len = 25,          # Minimum length of output read
        min_comparison_length = 9   # Minimum length to compare
        max_comparison_length = 16  # Maximum length to compare
        max_mismatches = 2,         # Max hamming distance
        log_file_handle = sys.stdout
        ):

    # Adapter_tuples_func generates a list of tuples:
    #       (adapter_name, adapter_in_R1, adapter_in_R2)
    adapter_tuples = adapter_tuples_func( fname1, fname2, max_comparison_length )

    outname1 = os.path.splitext(fname1)[0] + '_deadaptered.fastq'
    outname2 = os.path.splitext(fname2)[0] + '_deadaptered.fastq'

    f1 = open(fname1)
    f2 = open(fname2)
    o1 = open(outname1,'w')
    o2 = open(outname2,'w')

    total_reads = 0
    trimmed_reads = 0
    deleted_reads = 0
    adapters_found = Counter()
    while True:
        id_line1 = f1.readline().strip()
        id_line2 = f2.readline().strip()
        if not id_line1 or not id_line2: 
            if not id_line1 and not id_line2:
                break # End of file
            else:
                sys.exit('Unequal number of lines in paired files')

        total_reads += 1

        seq_line1 = f1.readline().strip()
        plus_line1 = f1.readline().strip()
        q_line1 = f1.readline().strip()

        seq_line2 = f2.readline().strip()
        plus_line2 = f2.readline().strip()
        q_line2 = f2.readline().strip()

        for adapter_name, adapter_in_R1, adapter_in_R2 in adapter_tuples:
            adapter_position = paired_end_adapter_position(seq_line1,
                                        seq_line2,
                                        adapter_in_R1,
                                        adapter_in_R2,
                                        min_comparison_length,
                                        max_mismatches)
    
            if adapter_position != None:
                adapters_found[ adapter_name ] += 1
                if adapter_position < min_read_len:
                    deleted_reads += 1
                    break
                else:
                    trimmed_reads += 1
                    seq_line1 = seq_line1[:adapter_position]
                    q_line1 = q_line1[:adapter_position]
                    seq_line2 = seq_line2[:adapter_position]
                    q_line2 = q_line2[:adapter_position]

        else:
            # Did not delete the read
            o1.write('\n'.join([id_line1, seq_line1, plus_line1, q_line1]) + '\n')
            o2.write('\n'.join([id_line2, seq_line2, plus_line2, q_line2]) + '\n')

    f1.close()
    f2.close()
    o1.close()
    o2.close()

    output_removal_statistics( 
        total_reads,
        deleted_reads,
        contaminants_found,
        log_file_handle,
        contaminant_label = 'Adapter')

    return (outname1, outname2)

def trim_single_read_adapters(fname, # File to trim
        min_read_len = 25,          # Minimum length of output read
        min_comparison_length = 12  # Minimum length to compare
        max_comparison_length = 16  # Maximum length to compare
        max_mismatches = 1,         # Max hamming distance
        log_file_handle = sys.stdout
        ):
    # Default min/max comparison len values of 12/16 have 2e-6 and 1e-8 probability of randomly
    # appearing within distance 1 of a given sequence.

    outname = os.path.splitext(fname)[0] + '_deadaptered.fastq'
    
    adapter_list = get_fastqc_contaminant_list( 'adapter_list.txt', 
            max_contaminant_length=max_comparison_length,
            include_rev_comp=False)
    
    f = open(fname)
    out = open(outname,'w')
    
    total_reads = 0
    trimmed_reads = 0
    deleted_reads = 0
    adapters_found = Counter()
    while True:
        id_line = f.readline().strip()
        if not id_line: break # End of file
        total_reads += 1
    
        seq_line = f.readline().strip()
        plus_line = f.readline().strip()
        q_line = f.readline().strip()
    
        for adapter_name, adapter_seq in adapter_list:
            adapter_position = single_end_adapter_position(seq_line,
                                    adapter_seq,
                                    min_comparison_length,
                                    max_mismatches)
    
            if adapter_position != None:
                adapters_found[adapter_name] += 1
                if adapter_position < min_read_len:
                    deleted_reads += 1
                    break
                else:
                    trimmed_reads += 1
                    seq_line = seq_line[:adapter_position]
                    q_line = q_line[:adapter_position]
        else:
            # Non-deleted read
            out.write('\n'.join([id_line, seq_line, plus_line, q_line]) + '\n')
    
    f.close()
    out.close()
    
    output_removal_statistics( 
        total_reads,
        deleted_reads,
        contaminants_found,
        log_file_handle,
        contaminant_label = 'Adapter')

    return outname

def output_removal_statistics( 
        total_reads,
        deleted_reads,
        contaminants_found,
        log_file_handle,
        contaminant_label = 'Contaminant',
        ):
    """
    Write a summary of the contaminants removed to log_file_handle.
    """
    contaminant_label = contaminant_label[0].upper() + contaminant_label[1:]
    log_file_handle.write( '%s removal statistics:\n' % contaminant_label)
    log_file_handle.write( 'Total reads: %d\n' % total_reads )
    log_file_handle.write( 'Trimmed reads: %d, %f%%\n' % \
            (trimmed_reads, 100*float(trimmed_reads)/total_reads) )
    log_file_handle.write( 'Deleted reads: %d, %f%%\n' % \
            (deleted_reads, 100*float(deleted_reads)/total_reads) )
    log_file_handle.write( '%ss found:\n' % contaminant_label)
    for name, quantity in reversed(sorted(contaminants_found.items(), key=lambda tup: tup[1])):
        log_file_handle.write( '  %12d  %s\n' % (quantity, name) )

def trim_paired_read_adapters(fname1, fname2, 
        min_read_len = 25,      # Minimum length of output read
        len_adapter_match = 13, # Length of identical truseq prefix
        max_mismatches = 2,     # Max hamming distance
        log_file_handle = sys.stdout
        ):
    outname1 = os.path.splitext(fname1)[0] + '_deadaptered.fastq'
    outname2 = os.path.splitext(fname2)[0] + '_deadaptered.fastq'
    
    max_mismatches = 1 # Max hamming distance
    min_comparison_length = 12 # 2e-6 chance of randomly appearing within distance 1
    max_comparison_length = 16 # 1e-8 chance of randomly appearing within distance 1
    contaminant_list = get_contaminant_list(max_comparison_length)
    
    f1 = open(fname1)
    f2 = open(fname2)
    out1 = open(outname1,'w')
    out2 = open(outname2,'w')
    
    total_reads = 0
    trimmed_reads = 0
    deleted_reads = 0
    contaminants_found = Counter()
    orphaned_reads = []
    while True:
        id_line1 = f1.readline().strip()
        id_line2 = f2.readline().strip()
        if not id_line1 or not id_line2: 
            if not id_line1 and not id_line2:
                break # End of file
            else:
                sys.exit('Unequal number of lines in paired files')
        total_reads += 1
    
        seq_line1 = f1.readline().strip()
        plus_line1 = f1.readline().strip()
        q_line1 = f1.readline().strip()
    
        seq_line2 = f2.readline().strip()
        plus_line2 = f2.readline().strip()
        q_line2 = f2.readline().strip()
    
        delete_read1 = False
        delete_read2 = False
        for contaminant_name, contaminant_seq in contaminant_list:
            # In this loop, we check for contaminants in both sequences, trimming or deleting as
            # appropriate. If both need deleted, that happens as normal. If only one needs deleted,
            # then the second read is added to the orphaned sequence list to be placed at the end of
            # f1.

            # Read 1
            if not delete_read1:
                adapter_position = single_end_adapter_position(seq_line1,
                                        contaminant_seq,
                                        min_comparison_length,
                                        max_mismatches)
        
                if adapter_position != None:
                    contaminants_found[contaminant_name] += 1
                    if adapter_position < min_read_len:
                        deleted_reads += 1
                        delete_read1 = True
                    else:
                        trimmed_reads += 1
                        seq_line1 = seq_line1[:adapter_position]
                        q_line1 = q_line1[:adapter_position]

            # Read 2
            if not delete_read2:
                adapter_position = single_end_adapter_position(seq_line2,
                                        contaminant_seq,
                                        min_comparison_length,
                                        max_mismatches)
        
                if adapter_position != None:
                    contaminants_found[contaminant_name] += 1
                    if adapter_position < min_read_len:
                        deleted_reads += 1
                        delete_read2 = True
                    else:
                        trimmed_reads += 1
                        seq_line2 = seq_line2[:adapter_position]
                        q_line2 = q_line2[:adapter_position]
        else:
            # Non-deleted read
            out1.write('\n'.join([id_line1, seq_line1, plus_line1, q_line1]) + '\n')
    
    f.close()
    out.close()
    
    log_file_handle.write( 'Adapter removal statistics:\n' )
    log_file_handle.write( 'Total reads: %d\n' % total_reads )
    log_file_handle.write( 'Trimmed reads: %d, %f%%\n' % \
            (trimmed_reads, 100*float(trimmed_reads)/total_reads) )
    log_file_handle.write( 'Deleted reads: %d, %f%%\n' % \
            (deleted_reads, 100*float(deleted_reads)/total_reads) )
    log_file_handle.write( 'Contaminants found:\n' )
    for name, quantity in reversed(sorted(contaminants_found.items(), key=lambda tup: tup[1])):
        log_file_handle.write( '  %12d  %s\n' % (quantity, name) )
    
    return outname1, outname2
