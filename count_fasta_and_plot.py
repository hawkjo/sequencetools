import sys
import matplotlib.pyplot as plt

if len(sys.argv) != 2:
    sys.exit("Usage: count_fasta_and_plot.py <fasta file>")
filename = sys.argv[1]

#-------------------------------------------------------------------------------
# Read file, gather sequence lengths, count GC bases
#-------------------------------------------------------------------------------
seq_lens = []
total_GC = 0
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
total_bases = sum(seq_lens)

#-------------------------------------------------------------------------------
# Calculate N-Statistics
#-------------------------------------------------------------------------------
seq_lens.sort(reverse=True)
##################################################
# Enter here the desired N-Statistic values 
# !!!!  Must be in ascending order  !!!!
Ns = [10,25,50,75,90] 
##################################################
running_base_count = 0
running_percent = 0
current_N_ind = 0
N_seq_ind = []
N_bp = []
for i, l in enumerate(seq_lens):
    running_base_count += l
    prev_percent = running_percent
    running_percent = float(running_base_count)/total_bases*100
    if prev_percent < Ns[current_N_ind] <= running_percent:
        N_seq_ind.append(i)
        N_bp.append(l)
        current_N_ind += 1
        if current_N_ind >= len(Ns): break

#-------------------------------------------------------------------------------
# Print and plot results
#-------------------------------------------------------------------------------
w = 35
print 'Total length of sequence:'.ljust(w) + '%d bp' % total_bases
print 'Total number of sequences:'.ljust(w) + str(len(seq_lens))
for i in range(len(N_bp)):
    print 'N%d:'.ljust(w) % Ns[i] + '%d bp, %d sequences' % (N_bp[i], N_seq_ind[i]+1)
print 'Total GC count:'.ljust(w) + '%d bp' % total_GC
print 'GC %:'.ljust(w) + '%.2f %%' % (float(total_GC)/total_bases*100)

plt.plot(seq_lens, linewidth=2)
for i in range(len(N_seq_ind)):
    plt.plot([N_seq_ind[i],N_seq_ind[i]], [0,N_bp[i]], linewidth=2, label='N%d'%(Ns[i]), zorder=-1)
plt.xlabel('Sequence Length Rank')
plt.ylabel('Sequence Length (bp)')
plt.title('Sequence Length Distribution')
plt.legend()
plt.show()
