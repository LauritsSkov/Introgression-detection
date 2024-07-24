import numpy as np

from hmm_functions import HMMParam, write_HMM_to_file, read_HMM_parameters_from_file
from helper_functions import find_runs

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Make test data
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
def create_test_data(data_set_length, n_chromosomes = 2, write_out_files = False, parameters_file = None):
    '''Create test data set of size data_set_length. Also create uniform weights and uniform mutation rates'''
    
    # Config
    np.random.seed(42)
    window_size = 1000
    mutation_matrix = {
        'A': [0, 0.16, 0.68, 0.16],
        'C': [0.16, 0,0.16, 0.68],
        'G': [0.68, 0.16, 0, 0.16],
        'T': [0.16, 0.68, 0.16, 0],
    }
    bases = ['A','C','G','T']
    base_composition = [0.31, 0.19, 0.19, 0.31]
    CHROMOSOMES = [f'chr{x + 1}' for x in range(n_chromosomes)] 


    # Initialization HMM parameters, prob of staring in states and print parameters to user
    hmm_parameters = read_HMM_parameters_from_file(parameters_file)
    state_values = [x for x in range(len(hmm_parameters.state_names))]
    
    print(f'> creating {len(CHROMOSOMES)} chromosomes each with {data_set_length} kb of test data with the following parameters..')
    print(f'> hmm parameters file: {parameters_file}')
    print(hmm_parameters)  

    # Initialize data
    observations_for_obsfile = []
    path = []
    observations = []
    chroms = []
    starts = []
    variants = []

    for chrom in CHROMOSOMES:
        for index in range(data_set_length):

            obs_counter = 0
            variants_list = []
            
            # Use prior dist if starting window
            if index == 0:
                current_state = np.random.choice(state_values, p=hmm_parameters.starting_probabilities)
            else:
                current_state = np.random.choice(state_values, p=hmm_parameters.transitions[prevstate] )

            path.append(current_state)

            n_mutations = np.random.poisson(lam=hmm_parameters.emissions[current_state]) 
            for mutation in [int(x) for x in np.random.uniform(low=index*window_size, high=index*window_size + window_size, size=n_mutations)]: 
                ancestral_base = np.random.choice(bases, p=base_composition)
                derived_base = np.random.choice(bases, p=mutation_matrix[ancestral_base])
                
                observations_for_obsfile.append(f'{chrom}\t{mutation}\t{ancestral_base}\t{ancestral_base + derived_base}') 
                
                obs_counter += 1
                variants_list.append(str(mutation))

            observations.append(obs_counter)
            chroms.append(chrom)
            starts.append(index * window_size)
            variants.append(','.join(variants_list))

            prevstate = current_state

    

    if write_out_files:
        # Make obs file
        with open('obs.txt','w') as obs_file:
            print('chrom', 'pos', 'ancestral_base', 'genotype', sep = '\t', file = obs_file)
            for line in observations_for_obsfile:
                print(line, file = obs_file)

        # Make weights file and mutation file
        with open('weights.bed','w') as weights_file, open('mutrates.bed','w') as mutrates_file:
            for chrom in CHROMOSOMES:
                print(chrom, 0, data_set_length * window_size, sep = '\t', file = weights_file)
                print(chrom, 0, data_set_length * window_size, 1, sep = '\t', file = mutrates_file)

        # Make initial guesses
        initial_guess = HMMParam(['Human', 'Archaic'], [0.5, 0.5], [[0.99,0.01],[0.02,0.98]], [0.03, 0.3]) 
        write_HMM_to_file(initial_guess, 'Initialguesses.json')

        # Write the "true" simulated segments
        with open('simulated_segments.txt', 'w') as out:
            print('chrom', 'start', 'end', 'length', 'state', sep = '\t', file = out)
            for (chrom, chrom_start_index, chrom_length_index) in find_runs(chroms):
                for (state_id, start_index, length_index) in find_runs(path[chrom_start_index:chrom_start_index + chrom_length_index]):
                    
                    state = hmm_parameters.state_names[state_id]
                    
                    start_index = start_index + chrom_start_index
                    genome_start = starts[start_index]
                    genome_length =  length_index * window_size
                    genome_end = genome_start + genome_length

                    print(chrom, genome_start, genome_end, genome_length, state, sep = '\t', file = out)


    weights = np.ones(len(observations))
    mutrates = np.ones(len(observations))                

    return np.array(observations).astype(int), chroms, starts, variants, np.array(mutrates).astype(float), np.array(weights).astype(float)
