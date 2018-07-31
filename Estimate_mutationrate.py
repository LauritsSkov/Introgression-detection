from collections import defaultdict
import sys

_, freqfile, window_size, smallwindow, callablefile, outfile = sys.argv
window_size = int(window_size)
smallwindow = int(smallwindow)

snps_counts_window = defaultdict(int)

with open(freqfile) as data:
	for line in data:
		if 'N_ALLELES' not in line:
			chrom, pos, _, total, ref, alt = line.strip().split()

			window = int(pos) - int(pos)%window_size

			# we want to calculate the frequence of variants that are present in the outgroup
			if ref.split(':')[1] not in ['0',total] and alt.split(':')[1] not in ['0', total]:
				snps_counts_window[window] += 1

genome_mean = 0
total_called = 0
callable_big_region = defaultdict(int)


with open(callablefile) as data:
	for line in data:
		chrom, start, callable_region = line.strip().split()
		
		callable_region = float(callable_region)
		window = int(start) - int(start)%window_size

		# How many bases can we call in each window?
		callable_big_region[window] += callable_region*smallwindow

		# How many bases can we call in total
		total_called += callable_region*smallwindow




genome_mean = sum(snps_counts_window.values())/(total_called/window_size)


print sum(snps_counts_window.values()), genome_mean, total_called

with open(callablefile) as data, open(outfile,'w') as out:
	for line in data:
		chrom, start, _ = line.strip().split()

		window = int(start) - int(start)%window_size
		mut_in_window = 0
		if callable_big_region[window] > 0:
			mut_in_window = snps_counts_window[window]/(callable_big_region[window]/float(window_size))

		mutrate = mut_in_window/genome_mean 

		out.write('{}\t{}\t{}\n'.format(chrom, start, mutrate))
