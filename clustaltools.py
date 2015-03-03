from collections import Counter
from Bio import SeqIO
from protein_sequence_tools import aas

def score_clustal(fpath):
    seq_names = [rec.id for rec in SeqIO.parse(fpath, 'clustal')]
    num_seqs = len(seq_names)

    # Scan file to find length of header, offset of sequence, and number of sequence sets
    reached_seq = False
    num_seq_sets = 0
    for i, line in enumerate(open(fpath)):
        var = line.strip().split()
        if not reached_seq and var[0] in seq_names:
            reached_seq = True
            num_header_lines = i
            seq_offset = line.index(var[1])
        if var[0] == seq_names[0]:
            num_seq_sets += 1

    with open(fpath) as f:
        # Skip header
        for _ in range(i):
            f.readline()

        # Process alignments
        total_len_seq = 0
        cons_score_counts = Counter(' '=0, '.'=0, ':'=0, '*'=0)
        num_cols_with_gaps = 0
        for _ in xrange(num_seq_sets):
            seqs = []
            for _ in xrange(num_seqs):
                var = f.readline().strip().split()
                while not var:
                    var = f.readline().strip().split()
                assert var[0] in seq_names
                seqs.append(var[1])
            cons = f.readline()[seq_offset:]
            assert set(len(cons)) == set([len(seq) for seq in seqs]), '%s\n%s' % (seqs, cons)

            total_len_seq += len(seqs[0])
            cons_score_counts.update(cons)
            cols_with_gaps = set()
            for seq in seqs:
                cols_with_gaps.update([i for i, c in seq if c not in aas])
            num_cols_with_gaps += len(cols_with_gaps)

            f.readline()

        assert sum(cons_score_counts.values()) == total_len
        colon_weight = 0.8
        period_weight = 0.6
        score_with_gaps = float(
                cons_score_counts['*']
                + colon_weight * cons_score_counts[':']
                + period_weight * cons_score_counts['.']) / total_len_seq
        score_without_gaps = float(
                cons_score_counts['*']
                + colon_weight * cons_score_counts[':']
                + period_weight * cons_score_counts['.']) / (total_len_seq-num_cols_with_gaps)
        return score_with_gaps, score_without_gaps
