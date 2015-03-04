from collections import Counter
from Bio import SeqIO
from protein_sequence_tools import aas

def score_clustal(fpath, period_weight=0.6, colon_weight=0.8, gap_treatment='all'):
    """score_clustal(fpath, period_weight=0.6, colon_weight=0.8, gap_treatment='all')
    
    score_clustal gives a very rough scoring of the quality of a clustal alignment based on a
    arbitrary weighted average of the conservation categories for each site:
    
    score = (1*num_asterisks + colon_weight*num_colons + period_weight*num_periods) / num_sites

    The gap_treatment parameter has the following options:
        'include_gaps': Include all sites regardless of gaps.
        'no_gaps': Ignore any site that has any gaps.
        'no_head_or_tail_gaps': Consider exactly the sites between the first site with no gaps and
            the last site with no gaps.
        'all': Return a 3-tuple with all three scores.
    """
    assert gap_treatment in ['include_gaps', 'no_gaps', 'no_head_or_tail_gaps', 'all'], \
            'Invalid gap_treatment method: %s' % gap_treatment

    seq_names = set([rec.id for rec in SeqIO.parse(fpath, 'clustal')])
    num_seqs = len(seq_names)

    # Scan file to find length of header and offset of sequence
    for i, line in enumerate(open(fpath)):
        var = line.strip().split()
        if var and var[0] in seq_names:
            num_header_lines = i
            seq_offset = line.index(var[1])
            break

    # Build sequences into dict given seq_name
    seq_given_seq_name = {seq_name: '' for seq_name in seq_names}
    cons = ''
    curr_seq_set = set()
    for i, line in enumerate(open(fpath)):
        if i < num_header_lines:
            continue
        elif curr_seq_set == seq_names:
            # All seqs in current set have been read. line is cons line.
            cons += line.replace('\n', '')[seq_offset:]
            curr_seq_set = set()
        else:
            var = line.strip().split()
            if var and var[0] in seq_names:
                seq_name = var[0]
                seq_given_seq_name[seq_name] += var[1]
                curr_seq_set.add(seq_name)
    # Assert that the cons string and all seq strings are the same length. If the program
    # generating the clustal file fails to pad the cons string with spaces, this will fail.
    assert set([len(cons)]) == set([len(seq) for seq in seq_given_seq_name.values()]), \
            '%s\n%s' % (seq_given_seq_name, cons)

    seq_len = len(cons)  # from previous assertion
    cons_class_counts = Counter({' ': 0, '.': 0, ':': 0, '*': 0})  # Must initialize properly
    cons_class_counts.update(cons)
    assert sum(cons_class_counts.values()) == seq_len

    cols_with_gaps = set()
    for seq in seqs.values():
        cols_with_gaps.update([i for i, c in enumerate(seq) if c in '.-'])
    num_head_gaps = next(i for i in xrange(seq_len) if i not in cols_with_gaps)
    num_tail_gaps = seq_len \
            - next(i for i in reversed(xrange(seq_len)) if i not in cols_with_gaps)

    numerator = float(
            cons_class_counts['*']
            + colon_weight * cons_class_counts[':']
            + period_weight * cons_class_counts['.'])
    score_with_gaps = numerator / seq_len
    score_without_gaps = numerator / (seq_len - len(cols_with_gaps))
    score_without_head_tail_gaps = numerator / (seq_len - num_head_gaps - num_tail_gaps)
    if gap_treatment == 'include_gaps':
        return score_with_gaps
    elif gap_treatment == 'no_gaps':
        return score_without_gaps
    elif gap_treatment == 'no_head_or_tail_gaps':
        return score_without_head_tail_gaps
    elif gap_treatment == 'all':
        return score_with_gaps, score_without_gaps, score_without_head_tail_gaps
