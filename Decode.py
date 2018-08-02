from templates import *

# Parameters (path to observations file, output file, model, weights file)
_, infile, outprefix, model, weights_file, mutfile, window_size = sys.argv
window_size = int(window_size)


# Load data
state_names, transitions, emissions, starting_probabilities, weights, mutrates = make_hmm_from_file(model, weights_file, mutfile) 
obs, chroms, starts, variants = read_observations_from_file(infile)




# Posterior decode the file
post_seq = Posterior_decoding(starting_probabilities, transitions, emissions, weights, obs, mutrates)

with open(outprefix + '.All_posterior_probs.txt','w') as posterior_sequence, open(outprefix + '.Summary.txt','w') as summary: 
    
    previos_seg = ''
    previous_chrom = ''

    counter = 1
    snp_counter = 0
    total_prob = 0.0
    
    start = 0
    end = 0

    # Make headers
    summary.write('name\tchrom\tstart\tend\tlength\tstate\tsnps\tmean_prob\n')
    posterior_sequence.write('chrom\tstart\tobservations\tvariants\tMostlikely\t{}\n'.format('\t'.join(state_names)))



    for i, (x,v, chrom, current_start, var) in enumerate(zip(obs, post_seq, chroms, starts, variants)):
        index, value = max(enumerate([float(y) for y in v]), key=operator.itemgetter(1))
        posterior_sequence.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(chrom, current_start, x, var, state_names[index], '\t'.join(v) ))

        snps = x
        current_seg = state_names[index]       

        if i > 0:

            # Extend block
            if current_seg == previos_seg and chrom == previous_chrom:
                counter += 1
                snp_counter += int(snps)
                total_prob += value
                end = current_start

            # Begin new block
            if current_seg != previos_seg or chrom != previous_chrom:

                mean_prob = total_prob / counter
                summary.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(outprefix, previous_chrom, start, end, end - start + window_size, previos_seg, snp_counter, mean_prob)) 
                
                counter = 1
                start = current_start
                snp_counter = int(snps)
                total_prob = value


        previos_seg = current_seg
        previous_chrom = chrom

    mean_prob = total_prob / counter
    summary.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(outprefix, previous_chrom, start, end, end - start + window_size, previos_seg, snp_counter, mean_prob))       
