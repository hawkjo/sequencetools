import sys
import os
import string
import re
from copy import deepcopy
from itertools import product
from collections import Counter
from Bio import SeqIO

dna_complements = string.maketrans('acgtnACGTN', 'tgcanTGCAN')


def dna_rev_comp(dna_string):
    # Note that this is string translation, not biopython
    return dna_string.translate(dna_complements)[::-1]


def transcribe(dna_string):
    dna_string = dna_string.upper()
    return dna_string.replace('T', 'U')

iupac_ambiguity_codes = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'T',
    'U': 'T',
    'M': 'AC',
    'R': 'AG',
    'W': 'AT',
    'S': 'CG',
    'Y': 'CT',
    'K': 'GT',
    'V': 'ACG',
    'H': 'ACT',
    'D': 'AGT',
    'B': 'CGT',
    'N': 'GATC',
    }

iupac_ambiguity_complements = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
    'U': 'A',
    'M': 'K',
    'R': 'Y',
    'W': 'W',
    'S': 'S',
    'Y': 'R',
    'K': 'M',
    'V': 'B',
    'H': 'D',
    'D': 'H',
    'B': 'V',
    'N': 'N',
    }

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

codon_set_given_aa = {}
for aa in set(aa_given_codon.values()):
    codon_set_given_aa[aa] = set([codon for codon, aaaa in aa_given_codon.items() if aaaa == aa])

# Now, build a general dict which accepts IUPAC ambiguity codes and returns the correct amino acid,
# possibly 'X'
aa_or_X_given_codon = {}
for (c1, nset1), (c2, nset2), (c3, nset3) in product(iupac_ambiguity_codes.items(), repeat=3):
    general_codon = c1 + c2 + c3
    possible_aas = set([aa_given_codon[transcribe(''.join(codon))] for codon in product(nset1, nset2, nset3)])
    if len(possible_aas) == 1:
        aa_or_X_given_codon[general_codon] = list(possible_aas)[0]
    else:
        aa_or_X_given_codon[general_codon] = 'X'
for c1, c2 in product(iupac_ambiguity_codes.keys(), repeat=2):
    aa_or_X_given_codon[c1+c2] = aa_or_X_given_codon[c1+c2+'N']  # Not necessarily X, e.g. CGN -> R
for c in iupac_ambiguity_codes.keys():
    aa_or_X_given_codon[c] = 'X'
for codon, aa in aa_given_codon.items():
    assert aa == aa_or_X_given_codon[codon]


def iterate_codons(na_string):
    assert isinstance(na_string, str)
    return (na_string[i:i+3] for i in xrange(0, len(na_string), 3))


def simple_translate(rna_string):
    return ''.join([aa_or_X_given_codon[codon] for codon in iterate_codons(rna_string)])


def translate_with_warnings(rna_string, next_rna_codon=None):
    rna_string = rna_string.upper()
    assert set(rna_string) <= set(iupac_ambiguity_codes.keys()), rna_string

    warnings = []
    if len(rna_string) == 0:
        return ('', ['empty sequence'])
    if len(rna_string) % 3 != 0:
        warnings.append('non-triplet')
    if rna_string[:3] not in ['AUG', 'ATG']:
        warnings.append('non-AUG first codon')
    if 'N' in rna_string:
        warnings.append('Ns/Xs in sequence')

    peptide = simple_translate(rna_string)
    if next_rna_codon is not None:
        next_aa = simple_translate(next_rna_codon)
    else:
        next_aa = 'X'
    try:
        if peptide.index('*') != len(peptide) - 1:
            warnings.append('stop codon in middle')
    except ValueError:
        pass
    if peptide[-1] != '*' and next_aa == '*':
        warnings.append('non-stop last codon but next is')
    elif peptide[-1] != '*':
        warnings.append('non-stop last codon')
    if not warnings:
        warnings = ['no warnings']
    return peptide, warnings

def records_to_strs(records):
    for record in records:
        if isinstance(record, str):
            yield record
        else:
            yield str(record.seq)

def translation_warnings_report(fname_or_seqiter):
    if (isinstance(fname_or_seqiter, str)
            and os.path.splitext(fname_or_seqiter)[-1] in ['.fa', '.fasta', '.fsa']):
        seq_iter = records_to_strs(SeqIO.parse(fname_or_seqiter, 'fasta'))
    elif hasattr(fname_or_seqiter, '__iter__'):
        seq_iter = records_to_strs(fname_or_seqiter)
    else:
        sys.exit('Cannot open or iterate over argument.')

    warnings = Counter()
    for seq in seq_iter:
        warnings.update(translate_with_warnings(seq)[1])
    for warning, count in warnings.items():
        print '%30s: %d' % (warning, count)
