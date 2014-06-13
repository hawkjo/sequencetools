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
