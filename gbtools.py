from copy import deepcopy
from itertools import izip
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from misctools import gzip_friendly_open
from general_sequence_tools import aa_or_X_given_codon

def iterate_ncbi_rna_cds_and_tranlation(fpath):
    """From NCBI RNA annotation file, iterate over CDSs and tranlations."""
    standard_mismatches = set(['X*'])  # NCBI translates * as X
    for rec in SeqIO.parse(gzip_friendly_open(fpath), 'gb'):
        cds_feats = [feat for feat in rec.features if feat.type.upper() == 'CDS']
        if not cds_feats:
            continue
        assert len(cds_feats) == 1, rec.name
        cds_feat = cds_feats[0]
        assert len(cds_feat.qualifiers['translation']) == 1
        ncbi_translation = str(cds_feat.qualifiers['translation'][0])
        codon_start = int(cds_feat.qualifiers['codon_start'][0])
        cds, cds_translation = cds_extract(rec.seq, cds_feat.location, codon_start)
        mismatches = set(c1+c2 for c1, c2 in izip(ncbi_translation, cds_translation) if c1 != c2)
        assert mismatches <= standard_mismatches
        gene_name = cds_feat.qualifiers['gene']
        yield rec.name, gene_name, cds, ncbi_translation


def cds_extract(seq, location, codon_start=1):
    """Specialized wrapper for extract method of FeatureLocation.

    Returns cds and cds_translation as strings.
    Corrects BeforePosition behavior by specifying default codon_start.
    Correctly handles hanging incomplete codons in tail.
    """
    if location.start.__class__.__name__ == 'ExactPosition':
        cds_Seq = location.extract(seq)
    else:
        new_location = FeatureLocation(
                codon_start - 1,
                location.end,
                location.strand,
                location.ref,
                location.ref_db)
        cds_Seq = new_location.extract(seq)
    cds = str(cds_Seq)
    cds_translation = str(cds_Seq.translate())

    # Fix hanging codons (non-triplets at 3' end)
    # Final cds should be a triplet, either trimmed or with additional 'N'
    if len(cds) % 3 != 0:
        hanging_codon = cds[-(len(cds)%3):]
        hanging_aa = aa_or_X_given_codon[hanging_codon]  # E.g. 'CG(N)' -> 'R'
        if hanging_aa != 'X':
            assert len(hanging_codon) == 2
            cds += 'N'
            cds_translation += hanging_aa
        else:
            cds = cds[:-(len(cds)%3)]

    # Remove stop codons
    if cds_translation.endswith('*'):
        cds = cds[:-3]
        cds_translation = cds_translation[:-1]
    return cds, cds_translation
