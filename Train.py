from templates import *

# Parameters (path to observations file, output file, model, weights file)
_, infile, outprefix, model, weights_file, mutfile = sys.argv

# Load data
state_names, transitions, emissions, starting_probabilities, weights, mutrates = make_hmm_from_file(model, weights_file, mutfile) 
obs = read_observations_from_file(infile)

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
