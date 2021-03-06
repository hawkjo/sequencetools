import numpy as np
cimport numpy as np
cimport cython
from adapters_cython cimport hamming_distance
from general_sequence_tools import dna_rev_comp

DTYPEINT = np.int
ctypedef np.int_t DTYPEINT_t

cdef int base_compare_func( char b1, char b2 ):
    if b1 == b2: return 1
    else: return -2

def max_alignment_score_above_cutoff(char *s1, char *s2, int cutoff):
    """
    Performs a local Smith-Waterman alignment until it finds a path with score at least as high as a
    given cutoff. If such a path is found, it immediately returns True.  Else, False.
    """
    # s1 is first dimension, s2 second
    cdef int sigma = 3   # Penalty for gap opening/extension. Should be positive.

    cdef np.ndarray[DTYPEINT_t, ndim=2] max_paths = np.zeros( (len(s1)+1,len(s2)+1) , dtype=DTYPEINT)

    cdef int i,j
    for i in range(1, len(s1)+1):
        for j in range(1, len(s2)+1):
            max_paths[i,j] = max(max_paths[i,j-1] - sigma, max_paths[i-1,j] - sigma, \
                            max_paths[i-1,j-1] + base_compare_func(s1[i-1], s2[j-1]) , 0)
            if max_paths[i,j] >= cutoff: return True
    return False

def test_for_contaminant( char *seq, int alignment_score_cutoff, int k, contaminant_kmers_dict, contaminant_dict):
    """
    Function to determine if seq is a contaminant. If so, return name, else return None.
    """
    # This function is designed to be an efficient contaminant finder. First, it tests whether any
    # kmers in the given seq are present in any contaminants. If so, those contaminants are then
    # tested against the full contaminant via the smith-waterman algorithm as implemented in
    # max_alignment_score_above_cutoff. If the alignment score is above the cutoff, the contaminant
    # name is returned. Otherwise, None is returned.
    candidate_contaminants = set()
    cdef int i
    for i in xrange(len(seq)-k+1):
        if seq[i:i+k] in contaminant_kmers_dict:
            candidate_contaminants.update( contaminant_kmers_dict[ seq[i:i+k] ] )

    for contaminant_name in candidate_contaminants:
        if max_alignment_score_above_cutoff(seq,
                contaminant_dict[ contaminant_name ],
                alignment_score_cutoff):
            return contaminant_name
    return None


def simple_test_for_contaminant(char *seq, char *contaminant, int max_hamming, int min_len):
    cdef int i
    cdef int len_seq = len(seq)
    cdef int len_contaminant = len(contaminant)

    # Test for matches overhanging beginning or end of seq
    for i in xrange(1, max_hamming+1):
        # beginning of seq
        subseq = seq[:len_contaminant-i]
        subcont = contaminant[i:]
        if hamming_distance(subseq,
                            subcont,
                            len_contaminant-i,
                            len_contaminant-i,
                            0) <= max_hamming-i:
            return True
        # end of seq
        subseq = seq[len_seq-len_contaminant+i:]
        subcont = contaminant[:len_contaminant-i]
        if hamming_distance(subseq,
                            subcont,
                            len_contaminant-i,
                            len_contaminant-i,
                            0) <= max_hamming-i:
            return True
    for i in xrange(max_hamming+1, len_contaminant-min_len):
        if (seq[:len_contaminant-i] == contaminant[i:]
                or seq[len_seq-len_contaminant+i:] == contaminant[:len_contaminant-i]):
            return True

    # Test for matches entirely within seq
    for i in xrange(len_seq-len_contaminant+1):
        if hamming_distance(seq,
                            contaminant,
                            len_seq,
                            len_contaminant,
                            i) <= max_hamming:
            return True
    return False

