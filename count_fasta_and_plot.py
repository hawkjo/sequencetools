import sys
import matplotlib.pyplot as plt


def count_fasta_and_plot(fasta_fname,
                         show=True,
                         stats_fname=None,
                         fig_fname=None,
                         fig_title='Sequence Length Distribution'):
    if stats_fname is not None:
        out = open(stats_fname, 'w')
    else:
        out = sys.stdout
    # ------------------------------------------------------------------------------
    # Read file, gather sequence lengths, count GC bases
    # ------------------------------------------------------------------------------
    seq_lens = []
    total_GC = 0
    l = 0
    with open(fasta_fname) as f:
        next(f)
        for line in f:
            if line[0] == '>':
                seq_lens.append(l)
                l = 0
            else:
                line = line.strip()
                l += len(line)
                total_GC += sum(1 for c in line if c in 'GC')
        seq_lens.append(l)
    total_bases = sum(seq_lens)

    # ------------------------------------------------------------------------------
    # Calculate N-Statistics
    # ------------------------------------------------------------------------------
    seq_lens.sort(reverse=True)
    ##################################################
    # Enter here the desired N-Statistic values
    Ns = [10, 25, 50, 75, 90]
    Ns.sort()
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
            if current_N_ind >= len(Ns):
                break

    # ------------------------------------------------------------------------------
    # Print and plot results
    # ------------------------------------------------------------------------------
    w = 35
    out.write('Total length of sequence:'.ljust(w) + '%d bp' % total_bases + '\n')
    out.write('Total number of sequences:'.ljust(w) + str(len(seq_lens)) + '\n')
    for i in range(len(N_bp)):
        out.write(
            'N%d:'.ljust(w) % Ns[i] + '%d bp, %d sequences' % (N_bp[i], N_seq_ind[i]+1) + '\n')
    out.write('Total GC count:'.ljust(w) + '%d bp' % total_GC + '\n')
    out.write('GC %:'.ljust(w) + '%.2f %%' % (float(total_GC)/total_bases*100) + '\n')
    out.close()

    if not show and fig_fname is None:
        return
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(seq_lens, linewidth=2)
    for i in range(len(N_seq_ind)):
        ax.plot([N_seq_ind[i], N_seq_ind[i]], [0, N_bp[i]],
                linewidth=2, label='N%d=%d, %d seqs' % (Ns[i], N_bp[i], N_seq_ind[i]), zorder=-1)
    ax.set_xlabel('Sequence Length Rank')
    ax.set_ylabel('Sequence Length (bp)')
    ax.set_title(fig_title + '\nTotal bp: %d' % total_bases)
    plt.legend()
    if fig_fname is not None:
        fig.savefig(fig_fname)
    if show:
        plt.show()
    plt.close(fig)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: count_fasta_and_plot.py <fasta_file> <plot_fpath>")
    fasta_fname = sys.argv[1]
    plot_fpath = sys.argv[2]
    count_fasta_and_plot(fasta_fname, fig_fname=plot_fpath)
