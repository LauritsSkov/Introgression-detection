from templates import *

# Parameters (path to observations file, output file, model, weights file)
_, infile, outprefix, model, weights_file = sys.argv

# Load data
state_names, transitions, emissions, starting_probabilities, weights = make_hmm_from_file(model, weights_file) 
obs = read_observations_from_file(infile)

# Posterior decode the file
post_seq = Posterior_decoding(starting_probabilities, transitions, emissions, weights, obs)

with open(outprefix + '.All_posterior_probs.txt','w') as posterior_sequence, open(outprefix + '.mostLikely_probs.txt','w') as posterior_sequence_mostlikely, open(outprefix + '.Summary.txt','w') as summary: 
    
    previos_seg = ''
    previos_chrom = ''
    previos_state = ''

    counter = 1
    start = 1
    snp_counter = 0


    for i, (x,v) in enumerate(zip(obs, post_seq)):
        index, value = max(enumerate([float(y) for y in v]), key=operator.itemgetter(1))
        posterior_sequence.write('{0}\t{1}\n'.format(x, '\t'.join(v) ))
        posterior_sequence_mostlikely.write('{0}\t{1}\n'.format(x, state_names[index]))


        snps = x
        current_seg = state_names[index]

        if i > 0:

            # Extend block
            if current_seg == previos_seg:
                counter += 1
                snp_counter += int(snps)


            # Begin new block
            if current_seg != previos_seg:


                summary.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(infile, start, i, counter, previos_seg, snp_counter)) 
                
                counter = 1
                start = i + 1
                snp_counter = int(snps)


        previos_seg = current_seg


    summary.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(infile, start, i+1, counter, previos_seg, snp_counter))       
