import sys, random, fastq_tools
import numpy as np
from itertools import izip
from Bio import SeqIO
"""
This program takes in three files:
    A fasta file of oracle truth transcripts
    An RSEM.isoforms.results file for example expression levels
    A fastq file with quantity and quality of desired read output

Fake RNA Seq reads are then generated and output into desired output file.
"""

if len(sys.argv) != 5:
    sys.exit('Usage: generate_fake_rnaseq_data.py <transcripts fasta> ' + \
            '<example RSEM.isoforms.results file> <example fastq file> <output file>')

transcripts_file = sys.argv[1]
RSEM_file = sys.argv[2]
fastq_file = sys.argv[3]
outfile = sys.argv[4]

#-------------------------------------------------
# Read tpm from RSEM
#-------------------------------------------------
print 'Reading TPM...'
example_tpm = []
with open(RSEM_file) as f:
    next(f)
    for line in f:
        example_tpm.append( float(line.strip().split()[5]) )

#-------------------------------------------------
# Read transcripts from fasta
#-------------------------------------------------
print 'Reading transcripts...'
with open(fastq_file) as f:
    f.readline()
    read_len = len( f.readline().strip() )

transcript_list = [record for record in SeqIO.parse(open(transcripts_file), 'fasta') if len(record.seq) > read_len]

print 'Adding poly-A tails...'
for i in xrange(len(transcript_list)):
    transcript_list[i].seq += random.randint(10,200) * 'A'

#---------------------------------------------------------------------------------------------------
# Generate multinomial probabilities of drawing fragment from each transcript
#---------------------------------------------------------------------------------------------------
print 'Generating multinomial probabilities...'
# Draw unnormalized transcript expression levels
if len(transcript_list) < len(example_tpm):
    # Without replacement if possible
    sample_expressions = random.sample(example_tpm) 
else:
    # With replacement otherwise
    sample_expressions = [random.choice(example_tpm) for _ in xrange(len(transcript_list))]

# Scale by length of transcript and normalize to get probability of drawing fragment from transcript
sample_expressions = [len(record) * expression \
        for record, expression in izip(transcript_list, sample_expressions)]
sample_expressions_total = sum(sample_expressions)
sample_expressions = map(lambda x: x/sample_expressions_total, sample_expressions)

print 'Writing sample expressions...'
with open('sample_expressions.txt','w') as expression_out:
    for transcript, expression in izip(transcript_list, sample_expressions):
        expression_out.write( '%s\t%e\n' % (transcript.name, expression) )
print 'Sample expressions:',np.array(sample_expressions)

#-------------------------------------------------
# Drawing numbers of fragments per transcript
#-------------------------------------------------
print 'Counting sequences in fastq file...'
fastq_seqs = sum(1 for line in open(fastq_file)) / 4
print fastq_seqs, 'sequences in file'

print 'Drawing numbers of fragments per transcript...'
all_draws = np.random.multinomial(fastq_seqs, sample_expressions, size=1)[0]
print 'All draws:',all_draws
print 'Sum all draws:', sum(all_draws)

#---------------------------------------------------------------------------------------------------
# Draw fragments, and make errors according to quality of analogous transcript
#---------------------------------------------------------------------------------------------------

print 'Drawing fragments, making errors, outputting results...'
phred_offset = fastq_tools.determine_phred_offset(fastq_file)

p_given_coded_phred = np.zeros(150)
for coded_Q in xrange(phred_offset, phred_offset + 42):
    Q = coded_Q - phred_offset
    p_given_coded_phred[coded_Q] = 10**(-0.1*Q)
    if p_given_coded_phred[coded_Q] > 1:
        sys.exit('Error probability greater than 1: p=%f, Q=%d, phred_offset=%d' % \
                    (p, Q, phred_offset) )

out = open(outfile,'w')
fastq_in = open(fastq_file)

for i, draws in enumerate(all_draws):
    for _ in xrange(draws):
        nameline = fastq_in.readline().strip()
        seqline = fastq_in.readline().strip()
        plusline = fastq_in.readline().strip()
        qualline = fastq_in.readline().strip()

        start = random.randint(0,len(transcript_list[i].seq) - len(seqline) )

        # Extract fragment from transcript
        seq = transcript_list[i][start: start+len(seqline)]
        if random.random() < 0.5:
            seq = seq.reverse_complement()
        seq = str(seq.seq)

        # Add errors according to phred score
        erroneous_seq = list(seq)
        for j in xrange(len(seq)):
            if random.random() < p_given_coded_phred[ord(qualline[j])]:
                erroneous_seq[j] = random.choice( 'ACGT'.replace( erroneous_seq[j], '') )

        out.write(nameline + '\n')
        out.write( ''.join(erroneous_seq)  + '\n')
        out.write(plusline + '\n')
        out.write(qualline + '\n')

out.close()
fastq_in.close()
