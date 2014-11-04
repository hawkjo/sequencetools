aas = 'ARNDCQEGHILKMFPSTWYVBZX'

strong_groups = [set('STA'),
                 set('NEQK'),
                 set('NHQK'),
                 set('NDEQ'),
                 set('QHRK'),
                 set('MILV'),
                 set('MILF'),
                 set('HY'),
                 set('FYW'),
                 ]

weak_groups = [set('CSA'),
               set('ATV'),
               set('SAG'),
               set('STNK'),
               set('STPA'),
               set('SGND'),
               set('SNDEQK'),
               set('NDEQHK'),
               set('NEQHRK'),
               set('FVLIM'),
               set('HFY'),
               ]


def is_strong_group(in_set):
    assert isinstance(in_set, set)
    for group in strong_groups:
        if in_set <= group:
            return True
    return False


def is_weak_group(in_set):
    assert isinstance(in_set, set)
    for group in weak_groups:
        if in_set <= group:
            return True
    return False


def clustal_conservation_char(in_aas):
    assert isinstance(in_aas, list) or isinstance(in_aas, str) or isinstance(in_aas, set), \
        'Input to clustal_conservation_char must be set, list, or string'

    aa_set = set([c.upper() for c in in_aas])

    if len(aa_set) == 1:
        return '*'
    elif is_strong_group(aa_set):
        return ':'
    elif is_weak_group(aa_set):
        return '.'
    else:
        return ' '
