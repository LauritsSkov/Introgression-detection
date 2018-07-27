import gzip
import sys
from collections import defaultdict


def readFasta(infile):
	sequence = ''
	if '.gz' in infile:
		with gzip.open(infile) as data:
			for line in data:
				if '>' in line:
					seqname = line.strip().replace('>','')
				else:
					sequence += line.strip().replace(' ','')

	else:
		with open(infile) as data:
			for line in data:
				if '>' in line:
					seqname = line.strip().replace('>','')
				else:
					sequence += line.strip().replace(' ','')

	return sequence




_, ancestral_file, outgroupfile, window_size, callablefile, outfile = sys.argv
window_size = int(window_size)

ancestral_allele = readFasta(ancestral_file)

derived_found = defaultdict(int)

with open(outgroupfile) as data:
	for line in data:
		if 'N_ALLELES' not in line:
			chrom, pos, _, total, ref, alt = line.strip().split()

			ref_allele, ref_count = ref.split(':')
			alt_allele, alt_count = alt.split(':')

			if (ancestral_allele[int(pos)-1].upper() == ref_allele.upper() and alt_count != '0') or (ancestral_allele[int(pos)-1].upper() == alt_allele.upper() and ref_count != '0'):
				derived_found[pos] = 1


private_variants_to_keep = defaultdict(list)

for line in sys.stdin:
	if 'N_ALLELES' not in line:
		chrom, pos, _, total, ref, alt = line.strip().split()

		window = int(pos) - int(pos)%window_size
		ref_allele, ref_count = ref.split(':')
		alt_allele, alt_count = alt.split(':')

		is_derived = False

		if (ancestral_allele[int(pos)-1].upper() == ref_allele.upper() and alt_count != '0') or (ancestral_allele[int(pos)-1].upper() == alt_allele.upper() and ref_count != '0'):
			is_derived = True

		if is_derived and derived_found[pos] != 1:

			private_variants_to_keep[window].append(pos)



with open(callablefile) as data, open(outfile,'w') as out:
	for line in data:
		chrom, start, _ = line.strip().split()

		window = int(start) - int(start)%window_size
		private_obs = len(private_variants_to_keep[window])
		out.write('{}\t{}\t{}\t{}\n'.format(chrom, start, private_obs, ','.join(private_variants_to_keep[window])))
