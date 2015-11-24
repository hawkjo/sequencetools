import os
from copy import deepcopy
from itertools import izip
from Bio import SeqIO
from Bio import Seq
from Bio.SeqFeature import FeatureLocation
from misctools import gzip_friendly_open
from general_sequence_tools import aa_or_X_given_codon

def iterate_ncbi_rna_cds_and_tranlation(fpath):
    """From NCBI RNA annotation file, iterate over CDSs and tranlations."""
    standard_mismatches = set(['X*', 'U*'])  # NCBI often translates * as X

    def get_singular_qualifier(feature, qual_name):
        """Return qualifier known to have exactly one entry given feature and qualifier name."""
        quals = feature.qualifiers[qual_name]
        assert len(quals) == 1, str(quals)
        return quals[0]

    for rec in SeqIO.parse(gzip_friendly_open(fpath), 'gb'):
        cds_feats = [feat for feat in rec.features if feat.type.upper() == 'CDS']
        if not cds_feats:
            continue
        assert len(cds_feats) == 1, rec.name
        cds_feat = cds_feats[0]

        gene_name = get_singular_qualifier(cds_feat, 'gene')
        protein_id = get_singular_qualifier(cds_feat, 'protein_id')
        ncbi_translation = get_singular_qualifier(cds_feat, 'translation')
        codon_start = int(get_singular_qualifier(cds_feat, 'codon_start'))
        cds, cds_translation = cds_extract(rec.seq, cds_feat.location, codon_start)
        # Note that mismatching first aa is handled separately for alternative start codons
        mismatches = set(c1+c2 for c1, c2 in
                         izip(ncbi_translation[1:], cds_translation[1:]) if c1 != c2)
        if ncbi_translation[0] != 'M' and ncbi_translation[0] != cds_translation[0]:
            mismatches.add(ncbi_translation[0] + cds_translation[0])
        assert mismatches <= standard_mismatches, \
                '\n'.join([repr(mismatches), rec.name, ncbi_translation, cds_translation])
        yield rec, gene_name, rec.name, protein_id, cds, ncbi_translation


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

def extract_cds_and_protein_fasta_from_rna(fpath):
    basename = os.path.splitext(fpath)[0]
    cds_fpath = basename + '_extracted_cds.fa'
    prot_fpath = basename + '_extracted_protein.fa'
    with open(cds_fpath, 'w') as cds_out, open(prot_fpath, 'w') as prot_out:
        for rec, gene_id, rna_id, prot_id, cds, prot in iterate_ncbi_rna_cds_and_tranlation(fpath):
            cds_rec = SeqIO.SeqRecord(Seq.Seq(cds.upper()),
                    rna_id,
                    '',
                    'gene_id=%s prot_id=%s' % (gene_id, prot_id))
            prot_rec = SeqIO.SeqRecord(Seq.Seq(prot.upper()),
                    prot_id,
                    '',
                    'gene_id=%s rna_id=%s' % (gene_id, rna_id))
            SeqIO.write(cds_rec, cds_out, 'fasta')
            SeqIO.write(prot_rec, prot_out, 'fasta')
