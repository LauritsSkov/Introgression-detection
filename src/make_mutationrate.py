import numpy as np
from collections import defaultdict

from helper_functions import sortby, make_callability_from_bed, Make_folder_if_not_exists

def make_mutation_rate(freqfile, outfile, callablefile, window_size):

    snps_counts_window = defaultdict(lambda: defaultdict(int))
    with open(freqfile) as data:
        for line in data:
            if not line.startswith('chrom'):
                chrom, pos = line.strip().split()[0:2]
                pos = int(pos)
                window = pos - pos%window_size
                snps_counts_window[chrom][window] += 1


    mutations = []
    genome_positions = []
    for chrom in sorted(snps_counts_window, key=sortby):
        lastwindow = max(snps_counts_window[chrom]) + window_size

        for window in range(0, lastwindow, window_size):
            mutations.append(snps_counts_window[chrom][window])
            genome_positions.append([chrom, window, window + window_size])

    mutations = np.array(mutations)

    if callablefile is not None:
        callability = make_callability_from_bed(callablefile, window_size)
        callable_region = []
        for chrom in sorted(snps_counts_window, key=sortby):
            lastwindow = max(snps_counts_window[chrom]) + window_size
            for window in range(0, lastwindow, window_size):
                callable_region.append(callability[chrom][window]/window_size)
    else:
        callable_region = np.ones(len(mutations)) * window_size

    genome_mean = np.sum(mutations) / np.sum(callable_region)

    Make_folder_if_not_exists(outfile)
    with open(outfile,'w') as out:
        print('chrom', 'start', 'end', 'mutationrate', sep = '\t', file = out)
        for genome_pos, mut, call in zip(genome_positions, mutations, callable_region):
            chrom, start, end = genome_pos
            if mut * call == 0:
                ratio = 0
            else:
                ratio = round(mut/call/genome_mean, 2)

            print(chrom, start, end, ratio, sep = '\t', file = out)