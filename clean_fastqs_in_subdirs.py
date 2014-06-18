import sys, os
from adapters import trim_single_read_adapters, trim_paired_read_adapters
from contaminants import remove_single_read_contaminants, remove_paired_read_contaminants
from low_quality_bases import bookend_qual_subread
from low_quality_bases import trim_single_read_low_quality_bases, trim_paired_read_low_quality_bases

#---------------------------------------------------------------------------------------------------
# Gather paired- and single-end fastq file names 
#---------------------------------------------------------------------------------------------------
# Find all fastq files
fastq_files = []
for dpath, dnames, fnames in os.walk('.'):
    for fname in fnames:
        if fname[-6:] == '.fastq' and 'adapter' not in fname:
            fastq_files.append( os.path.join( dpath, fname ) )

# Separate fastq files into paired- and single-end reads
se_fastq_files = []
pe_fastq_files = []
for fname in fastq_files:
    if fname[-8:] == '_1.fastq':
        # If first of pair is found, add both
        fname2 = fname[:-8] + '_2.fastq' 
        if fname2 not in fastq_files: sys.exit( 'Unpaired _1 file: '+fname )
        pe_fastq_files.append( (fname, fname2) )
    elif fname[-8:] == '_2.fastq':
        # If second of pair is found, test for presence of first, but do nothing
        fname1 = fname[:-8] + '_1.fastq'
        if fname1 not in fastq_files: sys.exit( 'Unpaired _2 file: '+fname )
    else:
        se_fastq_files.append( fname )
if 2*len(pe_fastq_files) + len(se_fastq_files) != len(fastq_files):
    sys.exit('Filename lost when converting to se/pe lists')

#---------------------------------------------------------------------------------------------------
# For each single-end read file, trim adapters and low-quality bases
#---------------------------------------------------------------------------------------------------

for fname in se_fastq_files:
    print fname
    log_file = fname[:-6] + '_single_read_cleaning_log.txt'
    with open(log_file, 'w') as log:
        print 'Removing adapters'
        deadaptered_fname = trim_single_read_adapters(fname, log_file_handle=log)
        log.write( '\n' )
        print 'Removing contaminants'
        decontaminated_fname = remove_single_read_contaminants(deadaptered_fname, log_file_handle=log)
        log.write( '\n' )
        print 'Removing low quality bases'
        trim_single_read_low_quality_bases(decontaminated_fname, log_file_handle=log)

#---------------------------------------------------------------------------------------------------
# For each pair of paired-end read files, trim adapters and low-quality bases
#---------------------------------------------------------------------------------------------------

for fname1, fname2 in pe_fastq_files:
    print (fname1, fname2)
    log_file = fname1[:-8] + '_paired_read_cleaning_log.txt'
    with open(log_file, 'w') as log:
        print 'Removing adapters'
        deadaptered_fname1, deadaptered_fname2 \
                = trim_paired_read_adapters(fname1, fname2, log_file_handle=log)
        log.write( '\n' )
        print 'Removing contaminants'
        decontaminated_fname1, decontaminated_fname2 \
                = remove_paired_read_contaminants(deadaptered_fname1, deadaptered_fname2, log_file_handle=log)
        log.write( '\n' )
        print 'Removing low quality bases'
        trim_paired_read_low_quality_bases(decontaminated_fname1, decontaminated_fname2, log_file_handle=log)
