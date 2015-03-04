from collections import Counter
from Bio import SeqIO
from protein_sequence_tools import aas

def score_clustal(fpath):
    seq_names = set([rec.id for rec in SeqIO.parse(fpath, 'clustal')])
    num_seqs = len(seq_names)

    # Scan file to find length of header and offset of sequence
    for i, line in enumerate(open(fpath)):
        var = line.strip().split()
        if var and var[0] in seq_names:
            num_header_lines = i
            seq_offset = line.index(var[1])
            break

    seqs = {seq_name: '' for seq_name in seq_names}
    cons = ''
    curr_seq_set = set()
    for i, line in enumerate(open(fpath)):
        if i < num_header_lines:
            continue
        elif curr_seq_set == seq_names:
            # All seqs in current set have been read. Next line in cons line.
            cons += line.replace('\n', '')[seq_offset:]
            curr_seq_set = set()
        else:
            var = line.strip().split()
            if var and var[0] in seq_names:
                seq_name = var[0]
                seqs[seq_name] += var[1]
                curr_seq_set.add(seq_name)
    # Assert that the cons string and all seq strings are the same length
    assert set([len(cons)]) == set([len(seq) for seq in seqs.values()]), '%s\n%s' % (seqs, cons)

    seq_len = len(cons)
    cons_score_counts = Counter({' ': 0, '.': 0, ':': 0, '*': 0})  # Must initialize properly
    cons_score_counts.update(cons)
    assert sum(cons_score_counts.values()) == seq_len

    cols_with_gaps = set()
    for seq in seqs.values():
        cols_with_gaps.update([i for i, c in enumerate(seq) if c in '.-'])
    num_head_gaps = next(i for i in xrange(seq_len) if i not in cols_with_gaps)
    num_tail_gaps = seq_len \
            - next(i for i in reversed(xrange(seq_len)) if i not in cols_with_gaps)

    colon_weight = 0.8
    period_weight = 0.6
    numerator = float(
            cons_score_counts['*']
            + colon_weight * cons_score_counts[':']
            + period_weight * cons_score_counts['.'])
    score_with_gaps = numerator / seq_len
    score_without_gaps = numerator / (seq_len - len(cols_with_gaps))
    score_without_head_tail_gaps = numerator / (seq_len - num_head_gaps - num_tail_gaps)
    return score_with_gaps, score_without_gaps, score_without_head_tail_gaps
