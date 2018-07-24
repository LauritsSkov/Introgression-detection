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




_, ancestral_file, outgroupfile, derived_snps_file = sys.argv
ancestral_allele = readFasta(ancestral_file)


with open(outgroupfile) as data, open(derived_snps_file,'w') as out:
	for line in data:
		if 'N_ALLELES' not in line:
			chrom, pos, _, total, ref, alt = line.strip().split()

			ref_allele, ref_count = ref.split(':')
			alt_allele, alt_count = alt.split(':')

			if (ancestral_allele[int(pos)-1].upper() == ref_allele.upper() and alt_count != '0') or (ancestral_allele[int(pos)-1].upper() == alt_allele.upper() and ref_count != '0'):
				out.write('{}\t{}\n'.format(chrom, pos))
