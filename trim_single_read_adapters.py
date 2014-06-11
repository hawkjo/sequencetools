import sys, os
from adapters import build_adapters, single_end_adapter_position, get_contaminant_list
from collections import Counter

if len(sys.argv) != 2:
    sys.exit('Usage: trim_paired_read_adapters.py <fastq file>')

filename = sys.argv[1]
outname = os.path.splitext(filename)[0] + '_unadaptered.fastq'

max_mismatches = 1 # Max hamming distance
min_comparison_length = 12 # 2e-6 chance of randomly appearing within distance 1
max_comparison_length = 16 # 1e-8 chance of randomly appearing within distance 1
contaminant_list = get_contaminant_list(max_comparison_length)

f = open(filename)
out = open(outname,'w')

total_reads = 0
trimmed_reads = 0
deleted_reads = 0
contaminants_found = Counter()
while True:
    id_line = f.readline().strip()
    if not id_line: break # End of file
    total_reads += 1

    seq_line = f.readline().strip()
    plus_line = f.readline().strip()
    q_line = f.readline().strip()

    for contaminant in contaminant_list:
        adapter_position = single_end_adapter_position(seq_line,
                                contaminant[1], # Tuple (name, seq)
                                min_comparison_length,
                                max_mismatches)

        if adapter_position != None:
            contaminants_found[contaminant[0]] += 1
            if adapter_position == 0:
                deleted_reads += 1
                break
            else:
                trimmed_reads += 1
                seq_line = seq_line[:adapter_position]
                q_line = q_line[:adapter_position]
    else:
        # Non-deleted read
        out.write('\n'.join([id_line, seq_line, plus_line, q_line]) + '\n')

f.close()
out.close()

print 'Total reads:,',total_reads
print 'Trimmed reads: %d, %f%%' % (trimmed_reads, 100*float(trimmed_reads)/total_reads)
print 'Deleted reads: %d, %f%%' % (deleted_reads, 100*float(deleted_reads)/total_reads)
print 'Contaminants found:'
for name, quantity in reversed(sorted(contaminants_found.items(), key=lambda tup: tup[1])):
    print '  %12d  %s' % (quantity, name)
