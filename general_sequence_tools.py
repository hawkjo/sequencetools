import sys, string

dna_complements = string.maketrans('acgtACGT', 'tgcaTGCA')

def dna_rev_comp(dna):
     return dna.translate(dna_complements)[::-1]

