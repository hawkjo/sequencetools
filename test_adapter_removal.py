import sys, random, re
from fastq_tools import make_add_errors_to_seq_func_given_example_file
import adapters

def generate_single_end_adapter_data( example_file,
        outfile_prefix, 
        adapter_seq,
        num_output_w_adapter_per_position,
        log_file_handle = sys.stdout):

    print 'Getting number of seqs and read length in',example_file
    num_seqs_in_example = sum( 1 for line in open( example_file ) ) / 4
    with open(example_file) as f:
        defline = f.readline().strip()
        seqline = f.readline().strip()
        read_len = len(seqline)

    # Output num_output_w_adapter_per_position reads with adapter in each position, and then an
    # equivalent number of reads without adapters
    num_output_seqs = 2 * num_output_w_adapter_per_position * read_len

    # Randomly sample the qualline information from the example file
    print 'Getting seq numbers to use'
    seqs_to_copy = set()
    while len(seqs_to_copy) < num_output_seqs:
        seqs_to_copy.add( random.randint(1,num_seqs_in_example) )

    print 'Getting quallines'
    read_num = 0
    quallines_to_use = []
    with open(example_file) as f:
        while True:
            defline = f.readline().strip()
            if not defline: break
            read_num += 1

            seqline = f.readline().strip()
            plusline = f.readline().strip()
            qualline = f.readline().strip()

            if read_num in seqs_to_copy:
                quallines_to_use.append( qualline )
    print 'Shuffling quallines'
    random.shuffle( quallines_to_use )

    print 'Generating fake data'
    outfile = outfile_prefix + '.fastq'
    add_errors_to_seq = make_add_errors_to_seq_func_given_example_file( example_file )
    with open(outfile,'w') as out:
        for i, qualline in enumerate( quallines_to_use ):
            adapter_position = i / num_output_w_adapter_per_position
            if adapter_position < read_len:
                defline = '@%d read_len=%d adapter_pos=%d' % ( i, read_len, adapter_position )
                seqline = ''.join([random.choice('ACGT') for _ in xrange(adapter_position)]) \
                            + adapter_seq[:read_len-adapter_position]
            else:
                defline = '@%d read_len=%d adapter_pos=%d' % (i, read_len, read_len)
                seqline = ''.join([random.choice('ACGT') for _ in xrange(read_len)])
            seqline = add_errors_to_seq( seqline, qualline )
            out.write( '\n'.join( [defline, seqline, '+', qualline] ) + '\n')

    log_file_handle.write( '---------- Fake data generation stats ----------\n')
    log_file_handle.write( 'Number of seqs in example file: %d\n' % (num_seqs_in_example) )
    log_file_handle.write( 'Read length: %d\n' % read_len )
    log_file_handle.write( 'Total reads generated: %d\n' % num_output_seqs )

    return outfile

def analyze_deadaptered_single_read_fake_data( fname, 
        num_output_w_adapter_per_position, 
        min_read_len, 
        log_file_handle = sys.stdout ):

    fake_data_defline_re = re.compile( r"""
        @(?P<read_num>\d+)                  # Read number
        \                                   # Space
        read_len=(?P<read_len>\d+)           # Read len
        \                                   # Space
        adapter_pos=(?P<adapter_pos>\d+)    # Adapter position
        """, re.VERBOSE)

    reads_correctly_unedited = 0
    reads_correctly_deleted = 0
    reads_correctly_trimmed = 0
    reads_not_trimmed_enough = 0
    reads_trimmed_too_much = 0
    reads_which_should_have_been_deleted_but_were_not = 0
    reads_which_should_not_have_been_deleted_but_were = 0
    bases_correctly_unedited = 0
    bases_correctly_removed = 0
    bases_which_should_have_been_removed_but_were_not = 0
    bases_which_should_not_have_been_removed_but_were = 0

    # We need the read length to known the last read which should be deleted
    with open(fname) as f:
        defline = f.readline().strip()
        read_num, read_len, adapter_position = map(int, fake_data_defline_re.match(defline).groups() )

    last_read_num = 2 * num_output_w_adapter_per_position * read_len - 1
    last_read_num_which_should_be_deleted = num_output_w_adapter_per_position * min_read_len - 1
    prev_read_num = -1

    # It turns out that dealing with reads which should not have been deleted but were is more
    # complicated due to the fact that they are not in the output file at all. We here encapsulate
    # the needed logic in function form for two reasons. First, it makes the decision tree below
    # more readable. Second, we will need this code twice below: once in the decision tree below,
    # but also again later at the end in case a few of the last reads were deleted.

    def stats_on_reads_which_should_not_have_been_deleted_but_were( read_num, prev_read_num ):
        # Here we test if the read should not have been deleted but was. Is so, there will
        # be an unexpected gap between prev_read_num and read_num.
        last_read_not_to_deal_with = max( prev_read_num, last_read_num_which_should_be_deleted )
        local_reads_which_should_not_have_been_deleted_but_were = read_num - last_read_not_to_deal_with - 1
        # Go through each deleted read, find the number of bases which should not have been deleted
        # in that read (equal to the adapter position, which is the read_len if no adapter), and add them up
        local_bases_which_should_not_have_been_removed_but_were = 0
        for i in xrange( last_read_not_to_deal_with + 1, read_num ):
            adapter_position = min( i / num_output_w_adapter_per_position, read_len )
            local_bases_which_should_not_have_been_removed_but_were += adapter_position
        # Some of the bases in the incorrectly deleted reads possibly should indeed have
        # been trimmed: exactly the total bases in those reads minus those which should not
        # have been trimmed.
        local_bases_correctly_removed = local_reads_which_should_not_have_been_deleted_but_were * read_len \
                - local_bases_which_should_not_have_been_removed_but_were
        return (local_reads_which_should_not_have_been_deleted_but_were,
                local_bases_which_should_not_have_been_removed_but_were,
                local_bases_correctly_removed)

    with open(fname) as f:
        while True:
            defline = f.readline().strip()
            if not defline: break

            seqline = f.readline().strip()
            plusline = f.readline().strip()
            qualline = f.readline().strip()

            if len(seqline) != len(qualline):
                sys.exit( 'Adapter code not trimming seqline and qualline equally' )

            read_num, read_len, adapter_position = map(int, fake_data_defline_re.match(defline).groups() )

            #----------------------------------------
            # Deal with reads which are not present
            #----------------------------------------
            if read_num != prev_read_num + 1:
                if prev_read_num < last_read_num_which_should_be_deleted:
                    local_correctly_deleted = min( read_num - prev_read_num - 1, 
                                            last_read_num_which_should_be_deleted - prev_read_num )
                    reads_correctly_deleted += local_correctly_deleted
                    bases_correctly_removed += local_correctly_deleted * read_len
    
                if read_num > last_read_num_which_should_be_deleted + 1:
                    local_reads_which_should_not_have_been_deleted_but_were, \
                        local_bases_which_should_not_have_been_removed_but_were, \
                        local_bases_correctly_removed = \
                        stats_on_reads_which_should_not_have_been_deleted_but_were( read_num, prev_read_num)
                    reads_which_should_not_have_been_deleted_but_were += \
                            local_reads_which_should_not_have_been_deleted_but_were
                    bases_which_should_not_have_been_removed_but_were += \
                            local_bases_which_should_not_have_been_removed_but_were
                    bases_correctly_removed += local_bases_correctly_removed

            #----------------------------------------
            # Deal with reads which are present
            #----------------------------------------
            if read_num <= last_read_num_which_should_be_deleted or adapter_position < min_read_len:
                # Should have been deleted but was not
                if not (read_num <= last_read_num_which_should_be_deleted and adapter_position < min_read_len):
                    sys.exit( 'Methods of determining deletion status inconsistent.')
                reads_which_should_have_been_deleted_but_were_not += 1
                bases_which_should_have_been_removed_but_were_not += len(seqline) 
            elif adapter_position < read_len and adapter_position == len(seqline):
                # Correctly trimmed
                reads_correctly_trimmed += 1
                bases_correctly_unedited += adapter_position
                bases_correctly_removed += read_len - adapter_position
            elif adapter_position == len(seqline):
                # Correctly unedited
                reads_correctly_unedited += 1
                bases_correctly_unedited += read_len 
            elif len(seqline) > adapter_position:
                # Not trimmed enough
                reads_not_trimmed_enough += 1
                bases_correctly_unedited += adapter_position
                bases_which_should_have_been_removed_but_were_not += len(seqline) - adapter_position
                bases_correctly_removed += read_len - len(seqline)
            elif len(seqline) < adapter_position:
                # Trimmed too much
                reads_trimmed_too_much += 1
                bases_correctly_unedited += len(seqline)
                bases_which_should_not_have_been_removed_but_were += adapter_position - len(seqline)
                bases_correctly_removed += read_len - adapter_position

            prev_read_num = read_num

    #--------------------------------------------------
    # Deal with reads posssibly deleted from the end
    #--------------------------------------------------
    if read_num != last_read_num:
        #if last_read_num != read_num + 1:
            # If more than one read was deleted, the afore-used function will deal with them.
            # We use last_read_num + 1 as the first argument because it is the next read number
            # which was not mistakenly deleted.
            local_reads_which_should_not_have_been_deleted_but_were, \
                local_bases_which_should_not_have_been_removed_but_were, \
                local_bases_correctly_removed = \
                stats_on_reads_which_should_not_have_been_deleted_but_were( last_read_num + 1, read_num)
            reads_which_should_not_have_been_deleted_but_were += \
                    local_reads_which_should_not_have_been_deleted_but_were
            bases_which_should_not_have_been_removed_but_were += \
                    local_bases_which_should_not_have_been_removed_but_were
            bases_correctly_removed += local_bases_correctly_removed
        # The last read needs dealt with separately, though.
        #reads_which_should_not_have_been_deleted_but_were += 1
        #bases_which_should_not_have_been_removed_but_were += read_len

    #--------------------------------------------------
    # Output results
    #--------------------------------------------------
    reads_to_count = last_read_num + 1                  # Total we should have found
    bases_to_count = (last_read_num + 1) * read_len     # Total we should have found

    log_file_handle.write( '---------- Adapter trimming stats ----------\n' )
    log_file_handle.write( 'Reads:\n' )
    log_file_handle.write( '    Total pre-editing: %d\n' % reads_to_count )
    log_file_handle.write( '    Correctly unedited: %d\n' % reads_correctly_unedited )
    log_file_handle.write( '    Correctly deleted: %d\n' % reads_correctly_deleted )
    log_file_handle.write( '    Correctly trimmed: %d\n' % reads_correctly_trimmed )
    log_file_handle.write( '    Not trimmed enough: %d\n' % reads_not_trimmed_enough )
    log_file_handle.write( '    Trimmed too much: %d\n' % reads_trimmed_too_much )
    log_file_handle.write( '    Should have been deleted but were not: %d\n' % \
            reads_which_should_have_been_deleted_but_were_not )
    log_file_handle.write( '    Should not have been deleted but were: %d\n' % \
            reads_which_should_not_have_been_deleted_but_were )
    log_file_handle.write( '\n' )
    log_file_handle.write( 'Bases:\n' )
    log_file_handle.write( '    Total pre-editing: %d\n' % bases_to_count )
    log_file_handle.write( '    Correctly unedited: %d\n' % bases_correctly_unedited )
    log_file_handle.write( '    Correctly removed: %d\n' % bases_correctly_removed )
    log_file_handle.write( '    Should have been removed but were not: %d\n' % \
            bases_which_should_have_been_removed_but_were_not )
    log_file_handle.write( '    Should not have been removed but were: %d\n' % \
            bases_which_should_not_have_been_removed_but_were )

    #--------------------------------------------------
    # Check that all reads and bases are accounted for
    #--------------------------------------------------
    reads_counted = \
            reads_correctly_unedited + \
            reads_correctly_deleted + \
            reads_correctly_trimmed + \
            reads_not_trimmed_enough + \
            reads_trimmed_too_much + \
            reads_which_should_have_been_deleted_but_were_not + \
            reads_which_should_not_have_been_deleted_but_were 

    bases_counted = \
            bases_correctly_unedited + \
            bases_correctly_removed + \
            bases_which_should_have_been_removed_but_were_not + \
            bases_which_should_not_have_been_removed_but_were

    if reads_counted != reads_to_count or bases_counted != bases_to_count:
        sys.exit( 
                """
    Error: Bases and/or reads not properly accounted for.
        Reads accounted for: %d
        Reads which should have been accounted for: %d
        Bases accounted for: %d
        Bases which should have been accounted for: %d
                """ % ( reads_counted, reads_to_count, bases_counted, bases_to_count )
                )

def test_single_read_adapter_removal( example_file, 
        outfile_prefix, 
        adapter_seq,
        num_output_w_adapter_per_position = 100,
        min_read_len = 25,          # Minimum length of output read
        min_min_comparison_length = 12, # Minimum length to compare
        max_max_comparison_length = 16, # Maximum length to compare
        max_max_mismatches = 1,         # Max hamming distance
        ):

    fake_data_log_file = outfile_prefix + '_fake_data_log.txt'
    with open(fake_data_log_file,'w') as data_gen_log:
        fake_data_file = generate_single_end_adapter_data( example_file,
            outfile_prefix, 
            adapter_in_R1,
            num_output_w_adapter_per_position, 
            log_file_handle=data_gen_log)

    for min_comparison_length in xrange(min_min_comparison_length, max_max_comparison_length + 1):
        for max_comparison_length in xrange(min_comparison_length, max_max_comparison_length + 1):
            for max_mismatches in xrange( max_max_mismatches ):

                log_file_name = '%s_%d_%d_%d_%d_%d_testing_log.txt' % \
                                    ( outfile_prefix,
                                    num_output_w_adapter_per_position ,
                                    min_read_len ,
                                    min_comparison_length ,
                                    max_comparison_length ,
                                    max_mismatches  )
            
                print 'Log file:',log_file_name
                with open(log_file_name, 'w') as log_file_handle:
                    log_file_handle.write( '-------------------- Adapter trimming test --------------------\n' )
                    log_file_handle.write( 'Parameters:\n' )
                    log_file_handle.write( 'Example File: %s\n' % example_file )
                    log_file_handle.write( 'Outfile prefix: %s\n' % outfile_prefix )
                    log_file_handle.write( 'Adapter sequence: %s\n' % adapter_seq )
                    log_file_handle.write( 'Number of reads per adapter position: %d\n' % num_output_w_adapter_per_position )
                    log_file_handle.write( 'Mininum read length: %d\n' % min_read_len )
                    log_file_handle.write( 'Minimum comparison length: %d\n' % min_comparison_length )
                    log_file_handle.write( 'Maximum comparison length: %d\n' % max_comparison_length )
                    log_file_handle.write( 'Maximum hamming distance: %d\n' % max_mismatches )
                    log_file_handle.write( '\n' )
                
                    log_file_handle.write( '---------- Adapter trimming stats ----------\n' )
                    deadaptered_file = adapters.trim_single_read_adapters( fake_data_file ,
                            min_read_len = min_read_len,
                            min_comparison_length = min_comparison_length,
                            max_comparison_length = max_comparison_length,
                            max_mismatches = max_mismatches,
                            log_file_handle=log_file_handle)
                
                    log_file_handle.write( '\n' )
                    analyze_deadaptered_single_read_fake_data( deadaptered_file , 
                        num_output_w_adapter_per_position, 
                        min_read_len, 
                        log_file_handle = log_file_handle )

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit( 'Usage: test_adapter_removal.py <example fastq> <output prefix>' )

    example_file = sys.argv[1]
    outfile_prefix = sys.argv[2]
    adapter_in_R1, adapter_in_R2 = adapters.build_adapters()
    num_output_w_adapter_per_position = 100

    test_single_read_adapter_removal( example_file, 
        outfile_prefix, 
        adapter_in_R1,
        num_output_w_adapter_per_position = num_output_w_adapter_per_position,
        min_read_len = 25,          
        min_min_comparison_length = 5,
        max_max_comparison_length = 20,
        max_max_mismatches = 3
        )
