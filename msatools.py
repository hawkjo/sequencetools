import numpy as np
from collections import Counter
from protein_sequence_tools import is_strong_group, is_weak_group


def find_col_type(col_set):
    if col_set <= set('.-'):
        return None  # no amino acids
    elif len(col_set) == 1:
        return 'u'  # unanimous
    elif col_set & set('.-'):
        return 'g'  # gapped
    elif is_strong_group(col_set):
        return 's'  # strong conservation
    elif is_weak_group(col_set):
        return 'w'  # weak conservation
    else:
        return 'n'  # no conservation


def alignment_types_string(SeqIO_records):
    rec_strs = [str(rec.seq) for rec in SeqIO_records]
    assert len(set(len(rec_str) for rec_str in rec_strs)) == 1, 'MSA with unequal seq lens'
    msa_len = len(rec_strs[0])
    col_sets = [set(''.join(rs[i] for rs in rec_strs)) for i in xrange(msa_len)]
    col_types = [find_col_type(col_set) for col_set in col_sets]
    return ''.join([ct for ct in col_types if ct is not None])


def tally_alignment_types(SeqIO_records,
                          gap_treatment='include_gaps',
                          start=None,
                          end=None):
    rec_strs = [str(rec.seq[start:end]) for rec in SeqIO_records]
    assert len(set(len(rec_str) for rec_str in rec_strs)) == 1, 'MSA with unequal seq lens'
    msa_len = len(rec_strs[0])

    assert gap_treatment in ['include_gaps', 'no_gaps', 'no_head_or_tail_gaps'], \
            'Invalid gap_treatment method: %s' % gap_treatment

    col_sets = [set(''.join(rs[i] for rs in rec_strs)) for i in xrange(msa_len)]
    col_types = [find_col_type(col_set) for col_set in col_sets]

    output = Counter({t: 0 for t in 'uswng'})

    if gap_treatment == 'include_gaps':
        output.update([ct for ct in col_types if ct is not None])

    elif gap_treatment == 'no_gaps':
        cols_with_gaps = set(i for i, cs in enumerate(col_sets) if set('.-') & cs)
        output.update([ct for i, ct in enumerate(col_types)
                       if ct is not None and i not in cols_with_gaps])

    elif gap_treatment == 'no_head_or_tail_gaps':
        try:
            first_no_gap_col = next(i for i in xrange(msa_len) if not set('.-') & col_sets[i])
            last_no_gap_col = next(i for i in reversed(xrange(msa_len)) if not set('.-') & col_sets[i])
            output.update([ct for ct in col_types[first_no_gap_col:last_no_gap_col+1]
                           if ct is not None])
        except StopIteration:
            # Every column has a gap
            pass
    return [output[t] for t in 'uswng']


def fraction_alignment_types(SeqIO_records, gap_treatment='include_gaps'):
    tallies = tally_alignment_types(SeqIO_records, gap_treatment)
    denom = float(max(1, sum(tallies)))
    return [t / denom for t in tallies]


def score_aa_msa(SeqIO_records,
                 score_weights=[0.71252797, 0.02239721, 0.00166703, -0.03624533, -0.70034688],
                 gap_treatment='include_gaps'):
    col_type_tallies = tally_alignment_types(SeqIO_records, gap_treatment)
    denom = max(float(sum(col_type_tallies)), 1)
    return np.dot(score_weights, col_type_tallies) / denom


def clean_msa(SeqIO_records, thresh=0.13320488208):
    max_score = score_aa_msa(SeqIO_records)
    def omit_record_i(records, i):
        return [r for j, r in records if j != i]

    while len(SeqIO_records) > 1 and max_score < thresh:
        scores = []
        for i in range(len(SeqIO_records)):
            scores.append(score_aa_msa(omit_record_i(SeqIO_records, i)))
        max_score = max(scores)
        max_idx = max(range(len(SeqIO_records)), key=lambda x: scores[x])
        SeqIO_records = omit_record_i(SeqIO_records, max_idx)
    if len(SeqIO_records) > 1:
        return max_score, SeqIO_records
    else:
        return None
