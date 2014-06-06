#!/usr/bin/env python

import sys, fastq_tools

if len(sys.argv) == 2:
    fastq_tools.trim_low_quality_bases( sys.argv[1] )
elif len(sys.argv) == 3:
    fastq_tools.trim_low_quality_bases( sys.argv[1], sys.argv[2] )
else:
    sys.exit('Usage: trim_low_quality_bases.py <fastq file> [<paired read fastq file>]')
