import sys, os

# First we add a cs-to-basespace dict. The dict works as
#
#       base_given_cs_dict[num][prev_base] = current_base
#
# where num is in {0,1,2,3,'.'}.  We allow the keys to be strings or integers.
base_given_cs_dict = {
        0: {'A':'A', 'C':'C', 'G':'G', 'T':'T', 'N':'N'},
        1: {'A':'C', 'C':'A', 'G':'T', 'T':'G', 'N':'N'},
        2: {'A':'G', 'C':'T', 'G':'A', 'T':'C', 'N':'N'},
        3: {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'},
        '.': {'A':'N', 'C':'N', 'G':'N', 'N':'N', 'N':'N'},
        }
for i in range(4): 
    base_given_cs_dict[str(i)] = base_given_cs_dict[i]

def naive_csfastq_to_fastq(fname, num_prefix_bases_to_drop=1):
    outname = os.path.splitext(fname)[0] + '_basespace.fastq'
    with open(fname) as f, open(outname,'w') as out:
        while True:
            defline = f.readline().strip()
            if not defline: break

            seqline = f.readline().strip()
            plusline = f.readline().strip()
            qualline = f.readline().strip()

            current_base = seqline[0]
            if current_base not in 'ACGT':
                sys.exit('Unexpected first base: ' + current_base)

            outseq = [current_base]
            for c in seqline[1:]:
                current_base = base_given_cs_dict[c][current_base]
                outseq.append( current_base )

            outseqline = ''.join(outseq[num_prefix_bases_to_drop:])
            qualline = qualline[num_prefix_bases_to_drop:]
            
            out.write( '\n'.join( [defline, outseqline, plusline, qualline] ) + '\n')
