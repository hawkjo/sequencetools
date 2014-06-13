import numpy as np
cimport numpy as np
cimport cython

DTYPEINT = np.int
ctypedef np.int_t DTYPEINT_t

cdef int base_compare_func( char b1, char b2 ):
    if b1 == b2: return 1
    else: return -2

cdef int max_alignment_score_above_cutoff(char *s1, char *s2, int cutoff):
    """
    Performs a local Smith-Waterman alignment until it finds a path with score at least as high as a
    given cutoff. If such a path is found, it immediately returns a value of 1. Otherwise, 0.
    """
    # s1 is first dimension, s2 second
    cdef int sigma = 3   # Penalty for gap opening/extension. Should be positive.

    cdef np.ndarray max_paths = np.zeros( (len(s1)+1,len(s2)+1) , dtype=DTYPEINT)

    cdef int i,j
    for i in range(1, len(s1)+1):
        for j in range(1, len(s2)+1):
            max_paths[i,j] = max(max_paths[i,j-1] - sigma, max_paths[i-1,j] - sigma, \
                            max_paths[i-1,j-1] + base_compare_func(s1[i-1], s2[j-1]) , 0)
            if max_paths[i,j] >= cutoff: return 1
    return 0
