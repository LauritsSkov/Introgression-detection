from collections import defaultdict
import sys

_, freqfile, window_size, callablefile, outfile = sys.argv
window_size = int(window_size)

snps_counts_window = defaultdict(int)

with open(freqfile) as data:
	for line in data:
		if 'N_ALLELES' not in line:
			chrom, pos, _, total, ref, alt = line.strip().split()

			window = int(pos) - int(pos)%window_size

			# we want to calculate the frequence of variants that are present in the outgroup
			if ref.split(':')[1] != '0' or alt.split(':')[1] != '0':
				snps_counts_window[window] += 1


genome_mean = sum(snps_counts_window.values()) / float(len(snps_counts_window))


with open(callablefile) as data, open(outfile,'w') as out:
	for line in data:
		chrom, start, _ = line.strip().split()

		window = int(start) - int(start)%window_size
		mutrate = snps_counts_window[window]/float(genome_mean)
		#print chrom, start, snps_counts_window[window], mutrate

		if mutrate > 5:
			mutrate = 5.0

		out.write('{}\t{}\t{}\n'.format(chrom, start, mutrate))
