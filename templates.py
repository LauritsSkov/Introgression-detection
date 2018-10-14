import msprime as msp
import numpy as np
from collections import defaultdict
import sys
import operator

#--------------------------------------------------------------------------------------
# Msprime functions
#--------------------------------------------------------------------------------------

def remove_outgroup(ts):
    sites = msp.SiteTable()
    mutations = msp.MutationTable()
    outgroup_IDs = ts.get_samples(0)

    for tree in ts.trees():
        for site in tree.sites():
            mut = site.mutations[0]
            pos = int(site.position)

            for leaf in tree.leaves(mut.node):
                if leaf in outgroup_IDs:
                    break
            else:
                site_id = sites.add_row(
                        position=site.position,
                        ancestral_state=site.ancestral_state)
                mutations.add_row(
                    site=site_id, node=mut.node, derived_state=mut.derived_state)
    tables = ts.dump_tables()
    new_ts = msp.load_tables(
        nodes=tables.nodes, edges=tables.edges, sites=sites, mutations=mutations)
    return new_ts

def combine_segs(segs, get_segs = False):
    merged = np.empty([0, 2])
    if len(segs) == 0:
        if get_segs:
            return([])
        else:
            return(0)
    sorted_segs = segs[np.argsort(segs[:, 0]), :]
    for higher in sorted_segs:
        if len(merged) == 0:
            merged = np.vstack([merged, higher])            
        else:
            lower = merged[-1, :]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1, :] = (lower[0], upper_bound) 
            else:
                merged = np.vstack([merged, higher])
    if get_segs:
        return(merged)
    else:
        return(np.sum(merged[:, 1] - merged[:, 0])/seq_len)

#--------------------------------------------------------------------------------------
# HMM functions
#--------------------------------------------------------------------------------------

def log_with_inf(x):
    if x == 0:
        return -np.inf
    else:
        return np.log(x)



def add_in_log_space(y):    

    if len(set(y)) == 1 and -np.inf in y: 
        return -np.inf

    else:
        x_star = max(y)  
        result = x_star + np.log(np.exp(y - x_star).sum())


    return result

def poisonprob(k, lamb):
    a = lamb**k
    b = np.exp(-lamb)
    c = np.math.factorial(k)
    d = a * b / c
    return d 


def log_with_inf_array(matrix):

    res = np.zeros((len(matrix), len(matrix[0])))
    for rows in range(len(matrix)):
        for col in range(len(matrix[0])):
            res[rows,col] = log_with_inf(matrix[rows,col])

    return res




def make_hmm_from_file(markov_param, weights_file, mut_file):
    
    with open(markov_param) as data:
       for line in data:
           exec(line.strip())


    # Load weights file
    weights = []
    with open(weights_file) as data:
        for line in data:
            weights.append(float(line.strip().split()[2]))


     # Load mutation rate file
    mutrates = []
    with open(mut_file) as data:
        for line in data:
            mutrates.append(float(line.strip().split()[2]))
            #mutrates.append(1.0)

    # Log transform the transitions
    for i, row in enumerate(transitions):
        for j, col in enumerate(row):
            transitions[i][j] = log_with_inf(col)

    # Log transform the starting_probabilities
    for i, start_prob in enumerate(starting_probabilities):
        starting_probabilities[i] = log_with_inf(start_prob)


    return (states, np.array(transitions), emissions, starting_probabilities, weights, mutrates)

def read_observations_from_file(f):
    obs = []
    chroms = []
    starts = []
    variants = []
    with open(f) as data:
        for line in data:
            chroms.append(line.strip().split()[0])
            starts.append(int(line.strip().split()[1]))
            obs.append(int(line.strip().split()[2]))
            
            if line.strip().split()[2] != '0':
                variants.append(line.strip().split()[3])
            else:
                variants.append('')

    return obs, chroms, starts, variants



def GET_forward_prob(init_start, transitions, emissions, weights, observations, mutrates):
    """
    Returns the probability of seeing the given `observations` sequence,
    using the Forward algorithm.
    """
    fractorials = {}
    frac_sum = 0

    for i in range(max(observations)+1):

        if i < 2:

            frac_sum = 0
        elif i == 2:
            frac_sum = 0
            frac_sum += log_with_inf(i)
        else:
            frac_sum += log_with_inf(i) 

        fractorials[i] = frac_sum
    

    state_nums = range(len(init_start))
    number_observations = len(observations)

    # Make and initialise forwards matrix
    forwards = np.zeros( (len(state_nums), number_observations) ) 


    for state in state_nums:
        forwards[state][0] = init_start[state]  + log_with_inf(poisonprob(observations[0],emissions[state]))


    # Fill out the matrix (we already filled out the first state)
    for t in range(1, number_observations): 
        for state in state_nums:            

            to_add = np.zeros(len(state_nums))
            calculate_poisson = - emissions[state] * weights[t] * mutrates[t] - fractorials[observations[t]] 
            if observations[t] != 0:                
                calculate_poisson += log_with_inf(emissions[state] * weights[t] * mutrates[t])*observations[t]

            for i,old_state in enumerate(state_nums):
                to_add[i] = transitions[old_state][state] + calculate_poisson + forwards[old_state][t-1]

            if len(set(to_add)) == 1 and -np.inf in to_add:
                forwards[state][t] =  -np.inf
            else:
                x_star = max(to_add)  
                forwards[state][t] = x_star + np.log(sum(np.exp(to_add - x_star))) 

    to_add = np.zeros(len(state_nums))
    for i,state in enumerate(state_nums):
        to_add[i] = forwards[state][-1] 


    final = add_in_log_space(to_add)

    return (final, forwards)

def GET_backward_prob(init_start, transitions, emissions, weights, observations, mutrates):
    """
    Returns the probability of seeing the given `observations` sequence,
    using the Backward algorithm.
    """    

    fractorials = {}
    frac_sum = 0
    for i in range(max(observations)+1):

        if i < 2:

            frac_sum = 0
        elif i == 2:
            frac_sum = 0
            frac_sum += log_with_inf(i)
        else:
            frac_sum += log_with_inf(i) 

        fractorials[i] = frac_sum

    state_nums = range(len(init_start))
    number_observations = len(observations)


    # Make and initialise backwards matrix
    backwards = np.zeros( (len(state_nums), number_observations) ) 

    # Fill out the matrix
    for t in reversed(range(1, number_observations)):
        for state in state_nums:
            to_add = np.zeros(len(state_nums))

            for i,next_state in enumerate(state_nums):
                calculate_poisson = - emissions[next_state] * weights[t] * mutrates[t] - fractorials[observations[t]] 
                if observations[t] != 0:
                    calculate_poisson += log_with_inf(emissions[next_state] * weights[t] * mutrates[t])*observations[t]

                to_add[i] = transitions[state, next_state] + calculate_poisson + backwards[next_state][t] 
  
            if len(set(to_add)) == 1 and -np.inf in to_add: 
                backwards[state][t-1] =  -np.inf

            else:
                x_star = max(to_add)
                backwards[state][t-1] = x_star + np.log(sum(np.exp(to_add - x_star))) 

    to_add = np.zeros(len(state_nums))
    for i,state in enumerate(state_nums):
        to_add[i] = init_start[state] + log_with_inf(poisonprob(observations[0],emissions[state])) + backwards[state][0]

    

    final = add_in_log_space(to_add)

    return (final, backwards) 




def Posterior_decoding(init_start, transitions, emissions, weights, observations, mutrates):
    """
    Posterior decoding, using the forward-backward algorithm. 
    """
    number_observations = len(observations)
    forward_prob,  forwards  = GET_forward_prob(init_start, transitions, emissions, weights, observations, mutrates)
    backward_prob, backwards = GET_backward_prob(init_start, transitions, emissions, weights, observations, mutrates)
    state_nums = range(len(init_start))

    results = []

    for t in range(number_observations):
        results.append([str(np.exp(forwards[state][t] + backwards[state][t] - forward_prob)) for state in state_nums])


    return results





def train_on_obs_pure_baum(init_start, transitions, emissions, weights, observations, mutrates):
    """
    Trains the model once, using the forward-backward algorithm. 
    """

    fractorials = {}
    frac_sum = 0
    for i in range(max(observations)+1):

        if i < 2:

            frac_sum = 0
        elif i == 2:
            frac_sum = 0
            frac_sum += log_with_inf(i)
        else:
            frac_sum += log_with_inf(i) 

        fractorials[i] = frac_sum 

    number_observations = len(observations)
    forward_prob,  forwards  = GET_forward_prob(init_start, transitions, emissions, weights, observations, mutrates)
    backward_prob, backwards = GET_backward_prob(init_start, transitions, emissions, weights, observations, mutrates)


    state_nums = range(len(init_start))
    
    posat = np.zeros( (len(state_nums), number_observations) )  
    for state in state_nums:
        for t in range(number_observations):
            posat[state][t] = forwards[state][t] + backwards[state][t] - forward_prob



    pot = np.zeros( (len(state_nums), len(state_nums), number_observations - 1)  )
    for state1 in state_nums:
        for state2 in state_nums:

            for t in range(number_observations-1):

                calculate_poisson = - emissions[state2] * weights[t+1] * mutrates[t+1] - fractorials[observations[t+1]] 

                if observations[t+1] != 0:
                    calculate_poisson += log_with_inf(emissions[state2] * weights[t+1] * mutrates[t+1])*observations[t+1]

                pot[state1][state2][t] = forwards[state1][t]  + transitions[state1, state2] + calculate_poisson + backwards[state2][t+1] - forward_prob



    # Initial starting probabilities
    start_prob = np.zeros((len(init_start)))

    for state in state_nums:
        start_prob[state] = np.exp(posat[state][0])



    # Transition probs
    trans = np.zeros((len(init_start), len(init_start)))            

    for state in state_nums:
        state_prob = add_in_log_space(posat[state])
        for oth in state_nums:
            trans[state][oth] = np.exp(add_in_log_space(pot[state][oth]) - state_prob) 

    for i,row in enumerate(trans):
        old_sum = sum(row)

        for j,col in enumerate(row):
            if old_sum == 0:
                trans[i][j] = 0
            else:
                trans[i][j] = col/old_sum


    # Emissions probs
    emit = np.zeros((len(init_start)))

    poissons = defaultdict(float)

    for state in state_nums:

        top = np.zeros(number_observations)
        bottom = np.zeros(number_observations)   

        top[0] = -np.inf
        bottom[0] = -np.inf        
        for t in range(number_observations):  
  
            top[t] = forwards[state][t] + backwards[state][t] + log_with_inf(observations[t])
            bottom[t] = forwards[state][t] + backwards[state][t] + log_with_inf(weights[t] * mutrates[t]) 

        top = add_in_log_space(top)
        bottom = add_in_log_space(bottom)
        emit[state] = np.exp(top - bottom)

    return (start_prob, trans, emit, forward_prob) 

def MakeHMMfile(state_names, starting_probabilities, transitions, emissions, outprefix):
    with open(outprefix + '.hmm','w') as out:
        out.write('# State names (only used for decoding)\n')
        out.write("states = [{states}]\n\n".format(states = ','.join(["'{}'".format(x) for x in state_names])))

        out.write('# Initialization parameters (prob of staring in states)\n')
        out.write("starting_probabilities = {values}\n\n".format(values = [x for x in starting_probabilities]))

        out.write('# transition matrix\n')
        out.write("transitions = [{values}]\n\n".format(values = ','.join(['[{}]'.format(','.join([str(y) for y in x])) for x in transitions])))

        out.write('# emission matrix (poisson parameter)\n')
        out.write("emissions = {values}\n".format(values = [x for x in emissions]))
        
        
def TrainModel(infile, outprefix, model, weights_file, mutfile):

    # Parameters (path to observations file, output file, model, weights file)

    # Load data
    state_names, transitions, emissions, starting_probabilities, weights, mutrates = make_hmm_from_file(model, weights_file, mutfile) 
    obs, _, _, _ = read_observations_from_file(infile)

    # Train model
    epsilon = 0.0001
    starting_probabilities, transitions, emissions, old_prob = train_on_obs_pure_baum(starting_probabilities, transitions, emissions, weights, obs, mutrates)

    with open(outprefix + '.log','w') as out:

        out.write('name\titeration\tstate\tvalue\tcomment\tmodel\n')

        for i in range(1000):

            transitions = log_with_inf_array(transitions)
            starting_probabilities, transitions, emissions, new_prob = train_on_obs_pure_baum(starting_probabilities, transitions, emissions, weights, obs, mutrates)
            
            print 'doing iteration {0} with old prob {1} and new prob {2}'.format(i, old_prob, new_prob)


            # Report emission values, transition values and likelihood of sequence
            out.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(infile, i, 1,new_prob, 'forward.probability', model))

            for state in range(len(state_names)):
                out.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(infile, i, state, emissions[state],'emission.state.{}'.format(state+1), model))
            
            for from_state in  range(len(state_names)):
                for to_state in  range(len(state_names)):
                    out.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(infile, i, state, transitions[from_state][to_state],'transition.state.{0}.to.{1}'.format(from_state+1,to_state+1), model))

            out.flush()

            if new_prob - old_prob < epsilon:        
                break   


            old_prob = new_prob


    # Write the optimal parameters
    with open(outprefix + '.hmm','w') as out:
        out.write('# State names (only used for decoding)\n')
        out.write("states = [{states}]\n\n".format(states = ','.join(["'{}'".format(x) for x in state_names])))

        out.write('# Initialization parameters (prob of staring in states)\n')
        out.write("starting_probabilities = {values}\n\n".format(values = [x for x in starting_probabilities]))

        out.write('# transition matrix\n')
        out.write("transitions = [{values}]\n\n".format(values = ','.join(['[{}]'.format(','.join([str(y) for y in x])) for x in transitions])))

        out.write('# emission matrix (poisson parameter)\n')
        out.write("emissions = {values}\n".format(values = [x for x in emissions]))

    return 0



def Decode(infile, outprefix, model, weights_file, mutfile, window_size, cutoff):
    
    # Parameters (path to observations file, output file, model, weights file)
    window_size = int(window_size)


    # Load data
    state_names, transitions, emissions, starting_probabilities, weights, mutrates = make_hmm_from_file(model, weights_file, mutfile) 
    obs, chroms, starts, variants = read_observations_from_file(infile)


    # Posterior decode the file
    post_seq = Posterior_decoding(starting_probabilities, transitions, emissions, weights, obs, mutrates)

    with open(outprefix + '.All_posterior_probs.txt','w') as posterior_sequence, open(outprefix + '.Summary.txt','w') as summary, open(outprefix + '.bed','w') as outbed: 
        
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
                        
                    if previos_seg == 'Archaic' and mean_prob > cutoff:
                        outbed.write('chr{0}\t{1}\t{2}\n'.format(previous_chrom, start, end))


                    counter = 1
                    start = current_start
                    snp_counter = int(snps)
                    total_prob = value


            previos_seg = current_seg
            previous_chrom = chrom

        mean_prob = total_prob / counter
        summary.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(outprefix, previous_chrom, start, end, end - start + window_size, previos_seg, snp_counter, mean_prob))       
        if previos_seg == 'Archaic' and mean_prob > cutoff:
            outbed.write('chr{0}\t{1}\t{2}\n'.format(previous_chrom, start, end))

        return 0       



