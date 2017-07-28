from templates import *

# Parameters (path to observations file, output file, model, weights file)
_, infile, outprefix, model, weights_file, mutfile = sys.argv

# Load data
state_names, transitions, emissions, starting_probabilities, weights, mutrates = make_hmm_from_file(model, weights_file, mutfile) 
obs = read_observations_from_file(infile)

# Posterior decode the file
post_seq = Posterior_decoding(starting_probabilities, transitions, emissions, weights, obs, mutrates)

with open(outprefix + '.All_posterior_probs.txt','w') as posterior_sequence, open(outprefix + '.mostLikely_probs.txt','w') as posterior_sequence_mostlikely, open(outprefix + '.Summary.txt','w') as summary: 
    
    previos_seg = ''
    counter = 1
    start = 1
    snp_counter = 0
    total_prob = 0.0

    # Make headers
    summary.write('name\tstart\tend\tlength\tstate\tsnps\tmean_prob\n')
    posterior_sequence.write('observations\t{}\n'.format('\t'.join(state_names)))
    posterior_sequence_mostlikely.write('observations\tstate\n')



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
                total_prob += value

            # Begin new block
            if current_seg != previos_seg:

                print total_prob, counter
                mean_prob = total_prob / counter
                summary.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(infile, start, i, counter, previos_seg, snp_counter, mean_prob)) 
                
                counter = 1
                start = i + 1
                snp_counter = int(snps)
                total_prob = value


        previos_seg = current_seg

    mean_prob = total_prob / counter
    summary.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(infile, start, i+1, counter, previos_seg, snp_counter, mean_prob))       
