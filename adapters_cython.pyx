import numpy as np
cimport numpy as np
cimport cython

DTYPEINT = np.int
ctypedef np.int_t DTYPEINT_t

cdef int hamming_distance(char *read,
                          char *adapter,
                          int read_length,
                          int adapter_length,
                          int start,
                         ):
    cdef int compare_length = min(adapter_length, read_length - start)
    cdef int mismatches = 0
    cdef int i

    for i in range(compare_length):
        if read[start + i] != adapter[i]:
            mismatches += 1

    return mismatches

def simple_hamming_distance(first, second):
    return hamming_distance(first, second, len(first), len(second), 0)

cdef int cython_hamming_with_N(char *s1, char *s2, int compare_length):
    cdef int mismatches = 0
    cdef int i

    for i in range(compare_length):
        if 'N' != s1[i] != s2[i] != 'N':
            mismatches += 1

    return mismatches

def simple_hamming_with_N(first, second):
    return cython_hamming_with_N(first, second, min(len(first), len(second)))

def find_adapter_positions(read, adapter, int min_comparison_length, int max_distance):
    cdef int read_length = len(read)
    cdef int adapter_length = len(adapter)
    cdef int max_start = len(read) - min_comparison_length
    cdef int distance, start
        
    positions = [] 
    for start in range(max_start + 1):
        distance = hamming_distance(read, adapter, read_length, adapter_length, start)
        if distance <= max_distance:
            positions.append(start)
    return positions
