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


_, repeatmask_file, callable_mask_file, window_size, chrom, outprefix = sys.argv
window_size = int(window_size)

#repeatmask_file =  "helperfiles/RepeatMasks/chr{}.fa.masked"
#callable_mask_file = "helperfiles/AccessibilityMasks/20140520.chr{}.strict_mask.fasta.gz"


bases_called = 0

# Mask file for repetitative regions
repeatmask = readFasta(repeatmask_file)
callable_mask = readFasta(callable_mask_file)

with open(outprefix + '.bed','w') as outbed, open (outprefix + '.txt','w') as out:
	d = defaultdict(int)

	prev_base = 'Notcalled'
	start = 0

	for i in range(len(callable_mask)): 

		repeat_base = repeatmask[i]
		callable_base = callable_mask[i]


		# Round down to nearest window start
		window = i - i%window_size
		d[window] += 0
			
		if repeat_base != 'N' and callable_base == 'P':

			current_base = 'Called'
			d[window] += 1
		else:
			current_base = 'Notcalled'


		# extend
		if current_base == prev_base:
			end = i

		# Make a new one
		if current_base != prev_base:
			if prev_base == 'Called':
				outbed.write('{}\t{}\t{}\t{}\n'.format(chrom, start, end, prev_base)) 

			start = i 
			end = i 

		prev_base = current_base

	if prev_base == 'Called':
		outbed.write('{}\t{}\t{}\t{}\n'.format(chrom, start, end, prev_base)) 


	# Write output files
	for window in range(0, max(d)+window_size, window_size):
		out.write('{}\t{}\t{}\n'.format(chrom, window, d[window] / float(window_size)))


