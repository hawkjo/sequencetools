import os
import random
from Bio import SeqIO, pairwise2
from Bio.Align.Applications import MafftCommandline

def mafft_ginsi_alignment(*args):
    return mafft_alignment('mafft-ginsi', *args)

def mafft_linsi_alignment(s1, s2):
    return mafft_alignment('mafft-linsi', *args)

def mafft_alignment(mafft_cmd, *args):
    fa_fpath = '/dev/shm/tmp.fa'
    mafft_fpath = '/dev/shm/tmp.mafft'

    # Write seqs to fasta file
    with open(fa_fpath, 'w') as out:
        for i, s in enumerate(args):
            out.write('>{}\n'.format(i))
            out.write('{}\n'.format(s))

    # Align
    mf_cline = MafftCommandline(mafft_cmd, input=fa_fpath)
    stdout, stderr = mf_cline()
    with open(mafft_fpath, 'w') as out:
        out.write(stdout)
    #check_call([mafft_cmd, '--quiet', fa_fpath, '>', mafft_fpath], shell=True)

    # Read and order output
    alignment = [(i, str(rec.seq))
                 for i, rec in enumerate(SeqIO.parse(open(mafft_fpath), 'fasta'))]
    output = [s.upper() for i, s in sorted(alignment)]

    # Delete files
    os.remove(fa_fpath)
    os.remove(mafft_fpath)
    
    return output


def find_dna_synthesis_errors(ref_seq, observed_seq):
    alignments = pairwise2.align.globalms(ref_seq, observed_seq, 2, -1, -1, -0.9)
    ref_align, observed_align = random.choice(alignments)[:2]  # Randomly select one of equals
    goodsubdel, ins = [], []
    bases = 'ACGT'
    bases_or_deletion = bases + '-'
    i = 0
    for rc, oc in zip(ref_align, observed_align):
        if rc == oc or (rc in bases and oc in bases_or_deletion):
            # Matching bases or substitutions or deletions
            goodsubdel.append('{}{}{}'.format(rc, i, oc))
            i += 1
        elif rc == '-' and oc in bases:
            # Insertions
            ins.append('I{}{}'.format(i, oc))
        elif rc != '-':
            # Sequencing error
            i += 1
    return goodsubdel, ins


def find_dna_errors(ref_seq, observed_seq):
    alignments = pairwise2.align.globalms(ref_seq, observed_seq, 2, -1, -1, -0.9)
    ref_align, observed_align = random.choice(alignments)[:2]  # Randomly select one of equals
    goodsubdel, ins = [], []
    bases = 'ACGT'
    bases_or_deletion = bases + '-'
    i = 0
    for rc, oc in zip(ref_align, observed_align):
        if rc != '-':
            # Matching bases or substitutions or deletions
            goodsubdel.append('{}{}{}'.format(rc, i, oc))
            i += 1
        else:
            # Insertions
            ins.append('I{}{}'.format(i, oc))
    return goodsubdel, ins
