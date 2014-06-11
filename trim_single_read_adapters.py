import sys, os
from adapters import trim_single_read_adapters
from collections import Counter

if len(sys.argv) != 2:
    sys.exit('Usage: trim_paired_read_adapters.py <fastq file>')

filename = sys.argv[1]

trim_single_read_adapters(filename)
