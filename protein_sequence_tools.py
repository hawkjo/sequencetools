#--------------------------------------------------------------------------------------------------
# Initialize useful aa abr/char mappings and lists
#--------------------------------------------------------------------------------------------------
aa_char_given_abr = { 
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V' 
    }
aa_abr_given_char = {char: abr for abr, char in aa_char_given_abr.items()}
aas = aa_char_given_abr.values()
aa_abrs = aa_char_given_abr.keys()
aa_abrs.sort()
aa_ind_given_abr = {abr: aa_abrs.index(abr) for abr in aa_abrs}
aa_abr_given_ind = {aa_abrs.index(abr): abr for abr in aa_abrs}

special_aa_char_given_abr = {
    'MLY': 'K', 'M2L': 'K', 'M3L': 'K', 'ACK': 'K', 'ARM': 'R', 'SRM': 'R', 'DRM': 'R'
    }

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


def is_same_aa(aa1, aa2):
    aa1 = aa1.upper()
    aa2 = aa2.upper()
    return (aa1 == aa2
            or (aa1, aa2) in aa_char_given_abr.items()
            or (aa2, aa1) in aa_char_given_abr.items()
            or (aa1, aa2) in special_aa_char_given_abr.items()
            or (aa2, aa1) in special_aa_char_given_abr.items())
