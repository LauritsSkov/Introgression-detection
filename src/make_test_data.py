from wsgiref.simple_server import demo_app
import numpy as np
from collections import defaultdict

from hmm_functions import HMMParam, get_default_HMM_parameters, write_HMM_to_file

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Make test data
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
def create_test_data(data_set_length, write_out_files = False):
    '''Create test data set of size data_set_length. Also create uniform weights and uniform mutation rates'''
    
    # Config
    np.random.seed(42)
    window_size = 1000
    mutation_rate_window = 1000000
    mutation_matrix = {
        'A': [0, 0.16, 0.68, 0.16],
        'C': [0.16, 0,0.16, 0.68],
        'G': [0.68, 0.16, 0, 0.16],
        'T': [0.16, 0.68, 0.16, 0],
    }
    bases = ['A','C','G','T']
    base_composition = [0.31, 0.19, 0.19, 0.31]


    # Initialization HMM parameters, prob of staring in states and print parameters to user
    state_values = [0,1]
    hmm_parameters = get_default_HMM_parameters()

    print(f'creating 2 chromosomes each with {data_set_length} kb of test data with the following parameters..\n')
    print(hmm_parameters)  

    # Initialize data
    observations_for_obsfile = []
    obs_counter = defaultdict(lambda: defaultdict(int))
    variants_dict = defaultdict(lambda: defaultdict(list))   

    initial_guess = HMMParam(['Human', 'Archaic'], [0.5, 0.5], [[0.99,0.01],[0.02,0.98]], [0.03, 0.3]) 

    for chrom in ['chr1', 'chr2']:
        for index in range(data_set_length):
            
            # Use prior dist if starting window
            if index == 0:
                current_state = np.random.choice(state_values, p=hmm_parameters.starting_probabilities)
            else:
                current_state = np.random.choice(state_values, p=hmm_parameters.transitions[prevstate] )

            n_mutations = np.random.poisson(lam=hmm_parameters.emissions[current_state]) 
            for mutation in [int(x) for x in np.random.uniform(low=index*window_size, high=index*window_size + window_size, size=n_mutations)]: 
                ancestral_base = np.random.choice(bases, p=base_composition)
                derived_base = np.random.choice(bases, p=mutation_matrix[ancestral_base])
                observations_for_obsfile.append(f'{chrom}\t{mutation}\t{ancestral_base}\t{ancestral_base + derived_base}') 
                obs_counter[chrom][index*window_size] += 1
                variants_dict[chrom][index*window_size].append(str(mutation))

            prevstate = current_state

    observations = []
    chroms = []
    starts = []
    variants = []
    for chrom in ['chr1', 'chr2']:
        lastwindow = max(obs_counter[chrom]) + window_size

        for window in range(0, lastwindow, window_size):
            observations.append(obs_counter[chrom][window])
            chroms.append(chrom)
            starts.append(window)
            variants.append(','.join(variants_dict[chrom][window]))

    
    weights = np.ones(len(observations))
    mutrates = np.ones(len(observations))

    if write_out_files:
        # Make obs file
        with open('obs.txt','w') as obs_file:
            print('chrom', 'pos', 'ancestral_base', 'genotype', sep = '\t', file = obs_file)
            for line in observations_for_obsfile:
                print(line, file = obs_file)

        # Make mutation file
        with open('mutrates.bed','w') as mutrates_file:
            for chrom in ['chr1', 'chr2']:
                for start in range(int(data_set_length * window_size / mutation_rate_window)):
                    print(chrom, start * mutation_rate_window, (start + 1) * mutation_rate_window, 1, sep = '\t', file = mutrates_file)

        # Make weights file
        with open('weights.bed','w') as weights_file:
            for chrom in ['chr1', 'chr2']:
                print(chrom, 0, data_set_length * window_size, sep = '\t', file = weights_file)

        # Make initial guesses
        write_HMM_to_file(initial_guess, 'Initialguesses.json')

    return observations, chroms, starts, variants, weights, mutrates