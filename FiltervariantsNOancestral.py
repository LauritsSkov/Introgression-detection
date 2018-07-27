import sys
from collections import defaultdict

_, outgroupfile, window_size, callablefile, outfile = sys.argv
window_size = int(window_size)

derived_found = defaultdict(int)

with open(outgroupfile) as data:
	for line in data:
		if 'N_ALLELES' not in line:
			chrom, pos, _, total, ref, alt = line.strip().split()

			ref_allele, ref_count = ref.split(':')
			alt_allele, alt_count = alt.split(':')

			if alt_count != '0' or ref_count != '0':
				derived_found[pos] = 1


private_variants_to_keep = defaultdict(list)

for line in sys.stdin:
    if 'N_ALLELES' not in line:
    	chrom, pos, _, total, ref, alt = line.strip().split()

    	window = int(pos) - int(pos)%window_size

    	if derived_found[pos] != 1:
			private_variants_to_keep[window].append(pos)



with open(callablefile) as data, open(outfile,'w') as out:
	for line in data:
		chrom, start, _ = line.strip().split()

		window = int(start) - int(start)%window_size
		private_obs = len(private_variants_to_keep)
		out.write('{}\t{}\t{}\t{}\n'.format(chrom, start, private_obs, ','.join(private_variants_to_keep)))
