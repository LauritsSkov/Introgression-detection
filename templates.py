import numpy as np
import sys
import operator
import math
from numba import jit
from collections import defaultdict


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


@jit
def log_with_inf(x):
    if x == 0:
        return -np.inf
    else:
        return np.log(x)



@jit
def poisonprob(k, lamb):
    a = lamb**k
    b = np.exp(-lamb)
    c = np.math.factorial(k)
    d = a * b / c
    return d 


@jit(nopython=True)
def add_in_log_space(y):    

    if len(set(y)) == 1 and np.any(y == -np.inf): 
        result = -np.inf
    else:

        x_star = np.max(y)  
        result = x_star + np.log(np.exp(y - x_star).sum())

    return result

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


    return (states, np.array(transitions), np.array(emissions), np.array(starting_probabilities), np.array(weights), np.array(mutrates))


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

    return np.array(obs), chroms, starts, variants





@jit(nopython=True)
def Forward_prob(init_start, transitions, emissions, observations, probabilities, state_nums, number_observations, forwards_in):
    """
    Returns the probability of seeing the given `observations` sequence,
    using the Forward algorithm.
    """

    for t in range(1, number_observations): 
        for state in state_nums: 
            toadd = np.zeros(len(state_nums))
            for state2 in state_nums:
               toadd[state2] = transitions[state2,state] + probabilities[state,t] + forwards_in[state2,t-1]

            forwards_in[state,t] = add_in_log_space(toadd)

    toadd = np.zeros(len(init_start))
    for state in state_nums:
        toadd[state] = forwards_in[state,-1] 

    final = add_in_log_space(toadd)   

    return (final, forwards_in)





# weights, mutrates in prob also add factorials

@jit(nopython=True)
def Backward_prob(init_start, transitions, emissions, observations, probabilities, state_nums, number_observations, backwards, reversedlist):
    """
    Returns the probability of seeing the given `observations` sequence,
    using the Backward algorithm.
    """    

    # Fill out the matrix
    for t in reversedlist:
        for state in state_nums:

            toadd = np.zeros(len(state_nums))
            for state2 in state_nums:
               toadd[state2] = transitions[state,state2] + probabilities[state2,t] + backwards[state2,t]

            backwards[state,t-1] = add_in_log_space(toadd)

    toadd = np.zeros(len(init_start))
    for state in state_nums:
        toadd[state] = init_start[state] + probabilities[state,0] + backwards[state,0]

    final = add_in_log_space(toadd)
    
    return (final, backwards) 


def Forward_backward(init_start, transitions, emissions, weights, observations, mutrates):
    """
    Posterior decoding, using the forward-backward algorithm. 
    """

    fractorials = np.zeros(len(observations))
    for i, obs in enumerate(observations):
        fractorials[i] = np.log(math.factorial(obs))

    number_observations = len(observations)
    state_nums = range(len(init_start))

    probabilities = np.zeros( (len(state_nums), number_observations) ) 
    for state in state_nums: 
        probabilities[state,:] = - (emissions[state] * weights * mutrates) - fractorials +  np.log( (emissions[state] * weights * mutrates)**observations )

    forwards_in = np.zeros( (len(init_start), number_observations) ) 
    forwards_in[:,0] = init_start  + probabilities[:,0] 
    backwards_in = np.zeros( (len(init_start), number_observations) ) 

    forward_prob,  forwards  = Forward_prob(init_start, transitions, emissions, observations, probabilities, state_nums, number_observations, forwards_in)

    reversedlist = [x for x  in range(number_observations-1, 0 ,-1)]
    backward_prob, backwards = Backward_prob(init_start, transitions, emissions, observations, probabilities, state_nums, number_observations, backwards_in, reversedlist)


    results = np.exp(forwards + backwards - forward_prob)

    return results


@jit(nopython=True)
def makeprobability_of_transition_matrix(state_nums, number_observations,forwards,transitions,probabilities,backwards,forward_prob):

    pot = np.zeros( (len(state_nums), len(state_nums), number_observations - 1)  )
    for state1 in state_nums:
        for t in range(number_observations-1):
            pot[state1,:,t] = forwards[state1,t]  + transitions[state1, :] + probabilities[:,t+1] + backwards[:,t+1] - forward_prob

    return pot

#@profile
def TrainBaumWelsch(init_start, transitions, emissions, weights, observations, mutrates):
    """
    Trains the model once, using the forward-backward algorithm. 
    """

    fractorials = np.zeros(len(observations))
    for i, obs in enumerate(observations):
        fractorials[i] = np.log(math.factorial(obs))

    number_observations = len(observations)
    state_nums = range(len(init_start))

    probabilities = np.zeros( (len(state_nums), number_observations) ) 
    for state in state_nums: 
        probabilities[state,:] = - (emissions[state] * weights * mutrates) - fractorials +  np.log( (emissions[state] * weights * mutrates)**observations )


    # Make and initialise forwards matrix
    forwards_in = np.zeros( (len(init_start), number_observations) ) 
    forwards_in[:,0] = init_start  + probabilities[:,0] #np.log(poisonprob(observations[0],emissions))
    backwards_in = np.zeros( (len(init_start), number_observations) ) 

    
    forward_prob,  forwards  = Forward_prob(init_start, transitions, emissions, observations, probabilities, state_nums, number_observations, forwards_in)

    reversedlist = [x for x  in range(number_observations-1, 0 ,-1)]
    backward_prob, backwards = Backward_prob(init_start, transitions, emissions, observations, probabilities, state_nums, number_observations, backwards_in, reversedlist)

    posat = forwards + backwards - forward_prob 
    pot = makeprobability_of_transition_matrix(state_nums, number_observations,forwards,transitions,probabilities,backwards,forward_prob)


    # Initial starting probabilities
    normalize = np.exp(posat).sum()
    start_prob = np.zeros(len(state_nums))

    # Transition probs
    trans = np.zeros((len(init_start), len(init_start)))            

    for state in state_nums:
        state_prob = np.logaddexp.reduce(posat[state])
        start_prob[state] = np.exp(state_prob) / normalize

        for oth in state_nums:
            trans[state,oth] = np.exp(np.logaddexp.reduce(pot[state,oth]) - state_prob) 

    for i,row in enumerate(trans):
        old_sum = row.sum()

        for j,col in enumerate(row):
            if old_sum == 0:
                trans[i][j] = 0
            else:
                trans[i][j] = col/old_sum


    # Emissions probs
    emit = np.zeros((len(init_start)))

    for state in state_nums:
        top = np.logaddexp.reduce(forwards[state,:] + backwards[state,:] + np.log(observations))
        bottom = np.logaddexp.reduce(forwards[state,:] + backwards[state,:] + np.log(weights * mutrates) )
        emit[state] = np.exp(top - bottom)

    return (start_prob, trans, emit, forward_prob) 


def TrainModel(infile, outprefix, model, weights_file, mutfile):

    # Parameters (path to observations file, output file, model, weights file)

    # Load data
    state_names, transitions, emissions, starting_probabilities, weights, mutrates = make_hmm_from_file(model, weights_file, mutfile) 
    obs, _, _, _ = read_observations_from_file(infile)

    # Train model
    epsilon = 0.0001
    starting_probabilities, transitions, emissions, old_prob = TrainBaumWelsch(starting_probabilities, transitions, emissions, weights, obs, mutrates)

    with open(outprefix + '.log','w') as out:

        out.write('name\titeration\tstate\tvalue\tcomment\tmodel\n')

        for i in range(1000):

            transitions = log_with_inf_array(transitions)
            starting_probabilities, transitions, emissions, new_prob = TrainBaumWelsch(starting_probabilities, transitions, emissions, weights, obs, mutrates)
            
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
    post_seq = Forward_backward(starting_probabilities, transitions, emissions, weights, obs, mutrates)


    with open(outprefix + '.All_posterior_probs.txt','w') as posterior_sequence, open(outprefix + '.Summary.txt','w') as summary, open(outprefix + '.bed','w') as outbed: 
        
        previos_seg = ''
        previous_chrom = ''

        counter = 1
        snp_counter = 0
        total_prob = 0.0
        
        start = 0

        # Make headers
        summary.write('name\tchrom\tstart\tend\tlength\tstate\tsnps\tmean_prob\n')
        posterior_sequence.write('chrom\tstart\tobservations\tvariants\tMostlikely\t{}\n'.format('\t'.join(state_names)))



        for i, (x,chrom, currentstart, var) in enumerate(zip(obs, chroms, starts, variants)):

            v = post_seq[:,i]

            index, value = max(enumerate([float(y) for y in v]), key=operator.itemgetter(1))
            posterior_sequence.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(chrom, currentstart, x, var, state_names[index], '\t'.join([str(val) for val in v]) ))

            snps = x
            current_seg = state_names[index]       

           

            # Extend block
            if current_seg == previos_seg and chrom == previous_chrom:
                counter += 1
                snp_counter += int(snps)
                total_prob += value

            # Begin new block
            if current_seg != previos_seg or chrom != previous_chrom:
                
                # Write previous block to output
                if i > 0:
                    mean_prob = total_prob / counter
                    end = starts[i-1] + window_size
                    summary.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(outprefix, previous_chrom, start, end, end - start, previos_seg, snp_counter, mean_prob)) 
                
                    if previos_seg == 'Archaic' and mean_prob > cutoff:
                        outbed.write('chr{0}\t{1}\t{2}\n'.format(previous_chrom, start, end))


                # Keep track of length, number of snps and prob of new block
                counter = 1
                start = currentstart
                snp_counter = int(snps)
                total_prob = value


           
            previos_seg = current_seg
            previous_chrom = chrom


        # when the file is done write the last segment
        mean_prob = total_prob / counter
        end = currentstart + window_size
        summary.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(outprefix, previous_chrom, start, end, end - start, previos_seg, snp_counter, mean_prob))

        if previos_seg == 'Archaic' and mean_prob > cutoff:
            outbed.write('chr{0}\t{1}\t{2}\n'.format(previous_chrom, start, end))

    return 0 
