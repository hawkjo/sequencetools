import sys, string
from copy import deepcopy
from itertools import product

dna_complements = string.maketrans('acgtACGT', 'tgcaTGCA')

def dna_rev_comp(dna_string):
    # Note that this is string translation, not biopython
    return dna_string.translate(dna_complements)[::-1]

def transcribe(dna_string):
    return dna_string.replace( 'T', 'U' )

aa_given_codon = {
        'UUU': 'F',      'CUU': 'L',      'AUU': 'I',      'GUU': 'V',
        'UUC': 'F',      'CUC': 'L',      'AUC': 'I',      'GUC': 'V',
        'UUA': 'L',      'CUA': 'L',      'AUA': 'I',      'GUA': 'V',
        'UUG': 'L',      'CUG': 'L',      'AUG': 'M',      'GUG': 'V',
        'UCU': 'S',      'CCU': 'P',      'ACU': 'T',      'GCU': 'A',
        'UCC': 'S',      'CCC': 'P',      'ACC': 'T',      'GCC': 'A',
        'UCA': 'S',      'CCA': 'P',      'ACA': 'T',      'GCA': 'A',
        'UCG': 'S',      'CCG': 'P',      'ACG': 'T',      'GCG': 'A',
        'UAU': 'Y',      'CAU': 'H',      'AAU': 'N',      'GAU': 'D',
        'UAC': 'Y',      'CAC': 'H',      'AAC': 'N',      'GAC': 'D',
        'UAA': '*',      'CAA': 'Q',      'AAA': 'K',      'GAA': 'E',
        'UAG': '*',      'CAG': 'Q',      'AAG': 'K',      'GAG': 'E',
        'UGU': 'C',      'CGU': 'R',      'AGU': 'S',      'GGU': 'G',
        'UGC': 'C',      'CGC': 'R',      'AGC': 'S',      'GGC': 'G',
        'UGA': '*',      'CGA': 'R',      'AGA': 'R',      'GGA': 'G',
        'UGG': 'W',      'CGG': 'R',      'AGG': 'R',      'GGG': 'G',
        }

aa_or_X_given_codon = deepcopy( aa_given_codon )
for b in 'ACGU':
    for s in ['NN%s', 'N%sN', '%sNN']:
        aa_or_X_given_codon[ s % b ] = 'X'
for b1, b2 in product( 'ACGU', repeat=2 ):
    for s in ['N%s%s', '%sN%s', '%s%sN']:
        aa_or_X_given_codon[ s % (b1,b2) ] = 'X'
aa_or_X_given_codon['NNN'] = 'X'

def simple_translate(rna_string):
    return ''.join([aa_or_X_given_codon[codon] \
            for codon in ( rna_string[i:i+3] for i in range( 0, (len(rna_string)/3)*3, 3 ) )] )

def translate_with_warnings(rna_string, next_rna_codon=None):
    warnings = []
    if len(rna_string) % 3 != 0:
        warnings.append( 'non-triplet' )
    if rna_string[:3] != 'AUG':
        warnings.append( 'non-AUG first codon' )
    if 'N' in rna_string:
        warnings.append( 'Ns/Xs in sequence' )
    peptide = simple_translate(rna_string)
    if next_rna_codon is not None:
        next_aa = simple_translate(next_rna_codon)
    else:
        next_aa = 'X'
    if '*' in peptide and peptide.index('*') != len(peptide) - 1:
        warnings.append( 'stop codon in middle' )
    if peptide[-1] != '*' and next_aa == '*':
        warnings.append( 'non-stop last codon but next is' )
    elif peptide[-1] != '*':
        warnings.append( 'non-stop last codon' )
    return peptide, warnings
