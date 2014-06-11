import sys
from adapters import trim_GSAF_paired_read_adapters

if len(sys.argv) != 3:
    sys.exit('Usage: trim_paired_read_adapters.py <fastq 1> <fastq 2>')

filename1 = sys.argv[1]
filename2 = sys.argv[2]

trim_GSAF_paired_read_adapters(filename1, filename2)
