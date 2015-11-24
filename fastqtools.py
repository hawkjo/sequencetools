import sys
import glob
import os
import re
import random
import numpy as np
from misctools import gzip_friendly_open
from collections import defaultdict, Counter
from itertools import izip, islice
from general_sequence_tools import dna_rev_comp
from adapters_cython import simple_hamming_distance
from Bio import SeqIO
from Bio.Seq import Seq


def determine_phred_offset(filename, num_reads_to_consider=1000):
    min_val = 126
    max_val = 0

    i = 0
    with gzip_friendly_open(filename) as f:
        while True:
            i += 1
            if i > num_reads_to_consider:
                break
            nameline = f.readline()
            assert nameline, 'Could not determine phred offset'
            seqline = f.readline()
            plusline = f.readline()
            qualline = f.readline()

            ascii_vals = map(ord, qualline.strip())  # Convert to numerical values
            min_val = min([min_val]+ascii_vals)
            max_val = max([max_val]+ascii_vals)

    if min_val < 50:
        return 33  # Illumina 1.8 and Sanger
    elif max_val > 89:
        return 64  # Illumina 1.3 - 1.7 and Solexa
    else:
        return determine_phred_offset(filename, num_reads_to_consider=2*num_reads_to_consider)

GSAF_fastq_id_line_re = re.compile(r"""
        ^@                              # Starts with @ sign
        (?P<instrument>[-A-Z0-9]+):     # Instrument name
        (?P<run>\d+):                   # Run id
        (?P<flowcell_id>[A-Z0-9]+):     # Flowcell id
        (?P<flowcell_lane>\d+):         # Flowcell lane
        (?P<tile>\d+):                  # tile number
        (?P<xcoord>\d+):                # x-coord
        (?P<ycoord>\d+)                 # y-coord
        \                               # space
        (?P<pair_member>[12]):          # member of pair (1 or 2)
        (?P<fail_status>[YN]):          # Y(es) if failed, N(o) otherwise
        (?P<control_bits>\d+):          # Control bits. 0 when all off
        (?P<barcode>[ACGT]+)            # barcode
    """, re.VERBOSE)


def get_GSAF_barcode(filename1, filename2=None):
    barcodes = []
    for fn in [filename1, filename2]:
        if fn is None:
            continue
        with open(fn) as f:
            line = f.readline().strip()
            m = GSAF_fastq_id_line_re.match(line)
            barcodes.append(m.group('barcode'))
    if len(barcodes) == 2 and barcodes[0] != barcodes[1]:
        sys.exit('Error: Different barcodes in read files')
    return barcodes[0]


def make_add_errors_to_seq_func(phred_offset):
    p_given_coded_phred = np.zeros(150)
    for coded_Q in xrange(phred_offset, phred_offset + 42):
        Q = coded_Q - phred_offset
        p = 10**(-0.1*Q)
        p_given_coded_phred[coded_Q] = p
        assert p <= 1, \
            'Error probability greater than 1: p=%f, Q=%d, phred_offset=%d' % (p, Q, phred_offset)

    def add_errors_to_seq(seq, qualline):
        erroneous_seq = list(seq)
        for j in xrange(len(seq)):
            if random.random() < p_given_coded_phred[ord(qualline[j])]:
                erroneous_seq[j] = random.choice('ACGT'.replace(erroneous_seq[j], ''))
        return ''.join(erroneous_seq)

    return add_errors_to_seq


def merge_seqs_in_place(s1, s2):
    return ''.join([c1 == c2 and c1 or 'N' for c1, c2 in zip(s1, s2)])


def collapse_if_overlapped_pair(seq_list, min_overlap=10, frac_error=0.1):
    if len(seq_list) == 1:
        return seq_list[0]
    elif len(seq_list) == 2:
        R1 = seq_list[0]
        R2_rc = dna_rev_comp(seq_list[1])
        max_overlap = min(map(len, [R1, R2_rc]))
        for i in range(min_overlap, max_overlap):
            if simple_hamming_distance(R2_rc[:i], R1[-i:]) < frac_error * i:
                return R1[:-i] + merge_seqs_in_place(R1[-i:], R2_rc[:i]) + R2_rc[i:]
            elif simple_hamming_distance(R1[:i], R2_rc[-i:]) < frac_error * i:
                return R2_rc[:-i] + merge_seqs_in_place(R2_rc[-i:], R1[:i]) + R1[i:]
    return tuple(seq_list)  # If none of the above


def make_add_errors_to_seq_func_given_example_file(fname):
    phred_offset = determine_phred_offset(fname)
    return make_add_errors_to_seq_func(phred_offset)


def convert_phred_scores(fname, out_phred_offset):
    if out_phred_offset not in [33, 64]:
        sys.exit('Error: out_phred_offset must be 33 or 64. Received %s' % repr(out_phred_offset))

    in_phred_offset = determine_phred_offset(fname)
    if in_phred_offset == out_phred_offset:
        print 'Cowardly refusing to convert %s from phred%d to phred%d' \
            % (fname, in_phred_offset, out_phred_offset)
        return -1

    phred_diff = out_phred_offset - in_phred_offset

    fname_parts = fname.split('.')
    out_fname = fname_parts[0] + '_phred' + str(out_phred_offset) + '.' + '.'.join(fname_parts[1:])

    from misc_tools import gzip_friendly_open
    with gzip_friendly_open(fname) as f, gzip_friendly_open(out_fname, 'w') as out:
        while True:
            defline = f.readline().strip()
            if not defline:
                break
            seqline = f.readline().strip()
            plusline = f.readline().strip()
            qualline = f.readline().strip()

            out_qualline = ''.join([chr(ord(c) + phred_diff) for c in qualline])

            out.write('\n'.join([defline, seqline, plusline, out_qualline]) + '\n')


def iterate_seqs(fpath):
    with gzip_friendly_open(fpath) as f:
        while True:
            defline = f.readline().strip()
            if not defline:
                break
            seqline = f.readline().strip()
            plusline = f.readline().strip()
            qualline = f.readline().strip()
            yield defline, seqline, plusline, qualline

def write_seq(out, defline, seqline, plusline, qualline):
    out.write('\n'.join([defline, seqline, plusline, qualline]) + '\n')


def trim_fastq(input_fpath, output_fpath, length, from_start_or_end):
    assert from_start_or_end in ['start', 'end']
    if from_start_or_end == 'start':
        s = slice(length, None)
    else:
        s = slice(None, -length)
    with open(output_fpath, 'w') as out:
        for defline, seqline, plusline, qualline in iterate_seqs(input_fpath):
            seqline = seqline[s]
            qualline = qualline[s]
            write_seq(out, defline, seqline, plusline, qualline)

def extract_xy_data(input_fpaths, defline_format='illumina'):
    if defline_format == 'illumina':
        def tile_xy_pair_from_defline(defline):
            word1, word2 = defline.split()
            inst, run, flowcell_id, flowcell_lane, tile, x, y = word1[1:].split(':')
            pair, filtered, cntrl, idx = word2.split(':')
            run, tile, x, y, pair = map(int, [run, tile, x, y, pair])
            return tile, x, y, pair
    elif defline_format == 'illumina_old':
        def tile_xy_pair_from_defline(defline):
            inst, flowcell_lane, tile, x, last_bits = defline[1:].split(':')
            y, last_bits = last_bits.split('#')
            idx, pair = last_bits.split('/')
            tile, x, y, pair = map(int, [tile, x, y, pair])
            return tile, x, y, pair
    else:
        sys.exit('Bad defline format: %s' % defline_format)

    tmp_dict = defaultdict(list)
    for fpath in input_fpaths:
        for defline, _, _, _ in iterate_seqs(fpath):
            tile, x, y, pair = tile_xy_pair_from_defline(defline)
            if pair == 1:
                tmp_dict[tile].append((x,y))

    output_dict = {}
    for key, val in tmp_dict.items():
        output_dict[key] = np.array(val)
    return output_dict

def find_paired_and_unpaired_files(fq_dir, skip_I_files=True):
    se_fpaths = []
    pe_fpaths = []
    fq_fpaths = glob.glob(os.path.join(fq_dir, '*.fastq*'))
    fq_fpaths.sort()
    i = 0
    while i < len(fq_fpaths):
        if skip_I_files and '_I1_' in fq_fpaths[i] or '_I2_' in fq_fpaths[i]:
            i += 1
            continue
        if i + 1 == len(fq_fpaths):
            se_fpaths.append(fq_fpaths[i])
            break
        m1 = re.match('^(.*)_R1(.*)$', fq_fpaths[i])
        m2 = re.match('^(.*)_R2(.*)$', fq_fpaths[i+1])
        if m1 and m2 and m1.group(1) == m2.group(1) and m1.group(2) == m2.group(2):
            pe_fpaths.append((fq_fpaths[i], fq_fpaths[i+1]))
            i += 2
        else:
            se_fpaths.append(fq_fpaths[i])
            i += 1
    return pe_fpaths, se_fpaths

def insert_lengths(paired_fq_fpaths, short=True):
    fpath1, fpath2 = paired_fq_fpaths
    assert simple_hamming_distance(fpath1, fpath2) == 1, paired_fq_fpaths
    ins_len_counter = Counter()
    it = izip(SeqIO.parse(gzip_friendly_open(fpath1), 'fastq'),
              SeqIO.parse(gzip_friendly_open(fpath2), 'fastq'))
    if short:
        it = islice(it, None, 10000)
    for rec1, rec2 in it:
        assert rec1.id == rec2.id, '%s\n%s' % (rec1.id, rec2.id)
        seq_obj = collapse_if_overlapped_pair([str(rec1.seq), str(rec2.seq)])
        if isinstance(seq_obj, str):
            ins_len_counter[len(seq_obj)] += 1
        else:
            ins_len_counter[sum(len(s) for s in seq_obj)] += 1
    return ins_len_counter


def get_reads_by_name(fastq_dir_or_fpath, read_names):
    if os.path.isfile(fastq_dir_or_fpath):
        se_fpaths = [fastq_dir_or_fpath]
        pe_fpaths = []
    elif os.path.isdir(fastq_dir_or_fpath):
        pe_fpaths, se_fpaths = find_paired_and_unpaired_files(fastq_dir_or_fpath)
    else:
        assert False, 'Fastq dir or fpath required'

    if not isinstance(read_names, set):
        read_names = set(read_names)

    out_r1_records = []
    out_r2_records = []
    for i, (fpath1, fpath2) in enumerate(pe_fpaths):
        print '%d of %d: %s' % (i + 1, len(pe_fpaths) + len(se_fpaths), (fpath1, fpath2))
        for rec1, rec2 in izip(SeqIO.parse(gzip_friendly_open(fpath1), 'fastq'),
                               SeqIO.parse(gzip_friendly_open(fpath2), 'fastq')):
            assert rec1.id == rec2.id, (fpath1, fpath2, rec1.id, rec2.id)
            if str(rec1.id) in read_names:
                out_r1_records.append(rec1)
                out_r2_records.append(rec2)

    for i, fpath in enumerate(se_fpaths):
        print '%d of %d: %s' % (len(pe_fpaths) + i + 1, len(pe_fpaths) + len(se_fpaths), fpath)
        for rec in SeqIO.parse(gzip_friendly_open(fpath), 'fastq'):
            if str(rec.id) in read_names:
                out_r1_records.append(rec)

    return out_r1_records, out_r2_records
