import sys, os
from contaminants_cython import *
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
        if line[0] == '#': continue
        var = line.strip().split('\t')
        name = var[0]
        seq = var[-1]
        if not name: continue # Blank lines
        if 'A' not in seq and 'C' not in seq and 'G' not in seq and 'T' not in seq:
            sys.exit('Non-DNA contaminant, %s: %s' % (name, seq) )

        contaminant_list.append( (name, seq[:max_contaminant_length]) )
        if include_rev_comp:
            contaminant_list.append( (name + ' RevComp', dna_rev_comp(seq[:max_contaminant_length])) )

    if not contaminant_list:
        sys.exit('No contaminants found in ' + contaminant_file)

    return contaminant_list

def output_contaminant_removal_statistics( 
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

def trim_paired_read_contaminants(fname1, fname2, 
        min_read_len = 25,          # Minimum length of output read
        min_comparison_length = 9   # Minimum length to compare
        max_comparison_length = 16  # Maximum length to compare
        max_mismatches = 2,         # Max hamming distance
        log_file_handle = sys.stdout
        ):
    outname1 = os.path.splitext(fname1)[0] + '_decontaminated.fastq'
    outname2 = os.path.splitext(fname2)[0] + '_decontaminated.fastq'
    
    contaminant_list = get_fastqc_contaminant_list( 'contaminant_list.txt' )
    
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
