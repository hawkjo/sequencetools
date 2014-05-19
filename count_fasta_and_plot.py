import sys
import matplotlib.pyplot as plt
from Bio import SeqIO

if len(sys.argv) != 2:
    sys.exit("Usage: count_fasta_and_plot.py <fasta file>")

filename = sys.argv[1]
seq_lens = []
total_GC = 0
#for record in SeqIO.parse( open(filename), 'fasta'):
#    seq_lens.append(len(record))
#    total_GC += sum(1 for c in record.seq if c in 'GC')
l=0
with open(filename) as f:
    next(f)
    for line in f:
        if line[0] == '>':
            seq_lens.append(l)
            l=0
        else:
            line = line.strip()
            l += len(line)
            total_GC += sum(1 for c in line if c in 'GC')
    seq_lens.append(l)

seq_lens.sort(reverse=True)
total_bases = sum(seq_lens)

levels = [10,25,50,75,90] # Must be in ascending order
running_tally = 0
percent = 0
current_level = 0
n_seq_inds = []
ns = []
for i, l in enumerate(seq_lens):
    running_tally += l
    prev_percent = percent
    percent = float(running_tally)/total_bases*100
    if prev_percent < levels[current_level] < percent:
        n_seq_inds.append(i)
        ns.append(l)
        current_level += 1
        if current_level >= len(levels): break

w = 40
print 'Total length of sequence:'.ljust(w),'%d bp' % total_bases
print 'Total number of sequences:'.ljust(w), len(seq_lens)
for i in range(len(ns)):
    print 'N%d:'.ljust(w) % levels[i],
    print '%d bp, %d sequences' % (ns[i], n_seq_inds[i]+1)
print 'Total GC count:'.ljust(w), '%d bp' % total_GC
print 'GC %:'.ljust(w), '%.2f %%' % (float(total_GC)/total_bases*100)

plt.plot(seq_lens)
for i in range(len(n_seq_inds)):
    plt.plot([n_seq_inds[i],n_seq_inds[i]], [0,ns[i]], label='n%d'%(levels[i]))
plt.xlabel('Sequence Rank')
plt.ylabel('Sequence Length (bp)')
plt.title('Sequence Length Distribution')
plt.legend()
plt.show()
