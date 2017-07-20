import numpy as np
from collections import defaultdict
import sys
import operator

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




def make_hmm_from_file(markov_param, weights_file):
    
    with open(markov_param) as data:
        for line in data:
            exec(line.strip())

    # Load weights file
    weights = []
    with open(weights_file) as data:
        for line in data:
            weights.append(float(line.strip()))

    # Log transform the transitions
    for i, row in enumerate(transitions):
        for j, col in enumerate(row):
            transitions[i][j] = log_with_inf(col)

    # Log transform the starting_probabilities
    for i, start_prob in enumerate(starting_probabilities):
        starting_probabilities[i] = log_with_inf(start_prob)


    return (states, np.array(transitions), emissions, starting_probabilities, weights)

def read_observations_from_file(f):
	obs = []
	with open(f) as data:
		for line in data:
			obs.append(int(line.strip()))

	return obs



def GET_forward_prob(init_start, transitions, emissions, weights, observations):
    """
    Returns the probability of seeing the given `observations` sequence,
    using the Forward algorithm.
    """
    fractorials = {}
    frac_sum = 0
    for i in range(100):

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
            calculate_poisson = - emissions[state] * weights[t] - fractorials[observations[t]] 
            if observations[t] != 0:                
                calculate_poisson += log_with_inf(emissions[state] * weights[t])*observations[t]

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

def GET_backward_prob(init_start, transitions, emissions, weights, observations):
    """
    Returns the probability of seeing the given `observations` sequence,
    using the Backward algorithm.
    """    

    fractorials = {}
    frac_sum = 0
    for i in range(100):

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
                calculate_poisson = - emissions[next_state] * weights[t] - fractorials[observations[t]] 
                if observations[t] != 0:
                    calculate_poisson += log_with_inf(emissions[next_state] * weights[t])*observations[t]

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




def Posterior_decoding(init_start, transitions, emissions, weights, observations):
    """
    Posterior decoding, using the forward-backward algorithm. 
    """
    number_observations = len(observations)
    forward_prob,  forwards  = GET_forward_prob(init_start, transitions, emissions, weights, observations)
    backward_prob, backwards = GET_backward_prob(init_start, transitions, emissions, weights, observations)
    state_nums = range(len(init_start))

    results = []

    for t in range(number_observations):
        results.append([str(np.exp(forwards[state][t] + backwards[state][t] - forward_prob)) for state in state_nums])


    return results




def train_on_obs_pure_baum(init_start, transitions, emissions, weights, observations):
    """
    Trains the model once, using the forward-backward algorithm. 
    """

    fractorials = {}
    frac_sum = 0
    for i in range(100):

        if i < 2:

            frac_sum = 0
        elif i == 2:
            frac_sum = 0
            frac_sum += log_with_inf(i)
        else:
            frac_sum += log_with_inf(i) 

        fractorials[i] = frac_sum 

    number_observations = len(observations)
    forward_prob,  forwards  = GET_forward_prob(init_start, transitions, emissions, weights, observations)
    backward_prob, backwards = GET_backward_prob(init_start, transitions, emissions, weights, observations)


    state_nums = range(len(init_start))
    
    posat = np.zeros( (len(state_nums), number_observations) )  
    for state in state_nums:
        for t in range(number_observations):
            posat[state][t] = forwards[state][t] + backwards[state][t] - forward_prob



    pot = np.zeros( (len(state_nums), len(state_nums), number_observations - 1)  )
    for state1 in state_nums:
        for state2 in state_nums:

            for t in range(number_observations-1):

                calculate_poisson = - emissions[state2] * weights[t+1] - fractorials[observations[t+1]] 

                if observations[t+1] != 0:
                    calculate_poisson += log_with_inf(emissions[state2] * weights[t+1])*observations[t+1]

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
            bottom[t] = forwards[state][t] + backwards[state][t] + log_with_inf(weights[t]) 

        top = add_in_log_space(top)
        bottom = add_in_log_space(bottom)
        emit[state] = np.exp(top - bottom)

    return (start_prob, trans, emit, forward_prob) 
