# Simulation studies

In this example I will be showing how to check the accuracy of the model given some simulated data. If you want to get an idea of how well this method will work on a given scenario then keep on reading. 

As I say in the paper (Skov et al 2018):
> We note that this method will likely only work in cases where the coalescence time distribution of the ingroup and archaic segments are 
> sufficiently different. This will work better in cases where the variation in the ingroup is a subset of variation in the outgroup so the 
> majority of variation in the common ancestor can be removed, as the case of Non-Africans and Africans.

The separation of coalescence times between ingroup/outgroup and ingroup/archaic depends on many factors, number of individuals in the outgroup you have, effective population size in the common ancestor of ingroup/outgroup, population structure and bottlenecks. 
Furthermore the model will work better for more recent introgression events, because recombination hasn't yet had time to break down the segments into pieces that are so small they cannot be detected. 


### Dependencies
To run the python script you will need numpy. I am using this version of python and numpy:

```
# To run the HMM model
python      version 2.7.12
numpy       version 1.11.1

# If you want to follow my simulation example 
msprime     version 0.5.0
pybedtools  version 0.7.10
```

I have a script that simulates one diploid individual and 50 diploid (100 haploid) individuals from an outgroup. The simulation is supposed to approximate the demographic history of a European test individual and 50 Sub-Saharan African individuals. 

```python
#!/usr/bin/python
import msprime as msp
import numpy as np
from collections import defaultdict
import pybedtools
from templates import *

#--------------------------------------------------------------------------------------
# Simulate
#--------------------------------------------------------------------------------------

chrom = '1'

# Generation time, mutation rate and recomination rate
gen_time = 29.0 
rec_rate = 1.2e-8 
mu = 1.2e-8

# Population sizes
N_DE = 5000
N_AF = 27122
N_ancestral_eurasians = 5000
N_ancestral_humans = 7000
N_neanderthal_human_ancestral = 18296
N_Europe = 3899
N_Asia = 5054

N_bottleneck_eurasians = 1305
N_bottleneck_nonAfricans = 250

# Split Times
T_human_archaic = 656908 / gen_time
T_bottleneck_nonAfricans = 62041 / gen_time
T_eurasians = 57100 / gen_time
T_europa_asia = 41997 / gen_time

# Bottleneck duration
T_bottlenecks = 100

# Admixture times and proportions from archaic humans
admixtureproportion = 5 #(in percent)
T_GF_DE = 54000 / gen_time


# Number of samples
n_outgroup = 100
n_ingroup = 2
samples = [msp.Sample(0, 0)] * n_outgroup + [msp.Sample(1, 0)]*n_ingroup

population_configurations = [
    msp.PopulationConfiguration(initial_size = N_AF),
    msp.PopulationConfiguration(initial_size = N_Europe),
    msp.PopulationConfiguration(initial_size = N_Asia),
    msp.PopulationConfiguration(initial_size = N_DE)   
]

demographic_events = [
    # European and Asian populations merge
    msp.MassMigration(
        time = T_europa_asia, source = 2, destination = 1, proportion = 1.0),
    msp.PopulationParametersChange(
        time = T_europa_asia, initial_size = N_ancestral_eurasians, 
        growth_rate = 0, population_id = 1),

    # Archaic admixture into Europeans
    msp.MassMigration(
        time = T_GF_DE, source = 1, destination = 3,
	proportion = admixtureproportion/100.0),

    # Eurasians experience bottleneck
    msp.PopulationParametersChange(
        time = T_eurasians - T_bottlenecks, 
        initial_size = N_bottleneck_eurasians, growth_rate = 0, population_id = 1),
    msp.PopulationParametersChange(
        time = T_eurasians, 
        initial_size = N_ancestral_eurasians, growth_rate = 0, population_id = 1),

    # Out of africa bottleneck in Eurasians and merge in Africans
    msp.PopulationParametersChange(
        time = T_bottleneck_nonAfricans - T_bottlenecks, 
        initial_size = N_bottleneck_nonAfricans, growth_rate = 0, population_id = 1),
    msp.MassMigration(
        time = T_bottleneck_nonAfricans, source = 1, destination = 0, proportion = 1.0),
    msp.PopulationParametersChange(
        time = T_bottleneck_nonAfricans, 
        initial_size = N_ancestral_humans, growth_rate = 0, population_id = 0),

    # Archiac and modern humans merge 650,000 years ago
    msp.MassMigration(
        time = T_human_archaic, source = 3, destination = 0, proportion = 1.0),
    msp.PopulationParametersChange(
         time = T_human_archaic, 
         initial_size = N_neanderthal_human_ancestral, 
         growth_rate = 0, population_id = 0), 
]


# Run simulations
#map_file = 'Recombination_map/genetic_map_GRCh37_chr{}.txt'.format(chrom)
#recomb_map = msp.RecombinationMap.read_hapmap(map_file)


print '-'*40
print 'Simulating sequences {} haplotypes from ingroup and {} haplotypes from outgroup...'.format(n_ingroup, n_outgroup) 
print '-'*40

ts = msp.simulate(
    samples = samples,
    population_configurations = population_configurations,
    demographic_events = demographic_events,

    # If you have a recombination map use:
    #recombination_map = recomb_map,

    # otherwise simulate with a uniform rate for a chromosome of some size
    length = 5000000, recombination_rate = rec_rate,
    random_seed=10,
    
    mutation_rate = mu,
    record_migrations = True
)

#--------------------------------------------------------------------------------------
# Get archaic segments
#--------------------------------------------------------------------------------------

print '-'*40
print 'Finding introgressed archaic sequences and making observations file...'
print '-'*40



# Keep track of which segments are actually introgressed (in this case from pop 3 into pop 1)
Testpopulation = ts.get_samples(1)
AF_ids = ts.get_samples(0)

de_seg = {i: [] for i in Testpopulation}

for mr in ts.migrations():
    if mr.source == 1 and mr.dest == 3:
        for tree in ts.trees(leaf_lists=True):
            if mr.left > tree.get_interval()[0]:
                continue
            if mr.right <= tree.get_interval()[0]:
                break
            for l in tree.leaves(mr.node):
                if l in Testpopulation:
                    #print l, mr
                    de_seg[l].append(tree.get_interval())

true_de_segs = [combine_segs(np.array(de_seg[i]), True) for i in sorted(de_seg.keys())]


with open('archaic_segments_chr{}.bed'.format(chrom),'w') as out:
    for haplotype, archaic_segments in enumerate(true_de_segs):
        for archaic_segment in archaic_segments:
            out.write('chr{}\t{}\t{}\t{}\n'.format(chrom, int(archaic_segment[0]), int(archaic_segment[1]), haplotype))

#--------------------------------------------------------------------------------------
# Get observation file, weights file and mutrates file
#--------------------------------------------------------------------------------------

ts = remove_outgroup(ts)

# Write genotype output
with open('haplo_chr{}.txt'.format(chrom),'w') as out:
    for variant in ts.variants():
        haplotypes_as_string = '\t'.join([str(x) for x in variant.genotypes[Testpopulation]])
        out.write('{}\t{}\n'.format(int(variant.site.position), haplotypes_as_string))

# Convert haplo file and make mutrates and weights
window_size = 1000
obs = defaultdict(list)
with open('haplo_chr{}.txt'.format(chrom)) as data:
    for line in data:
        pos, hap1, hap2 = line.strip().split()
        obs[int(pos) - int(pos)%window_size].append(pos)

with open('weights.txt','w') as w, open('mutrates.txt','w') as mut, open('observations.txt','w') as observations:
    for window in range(0, max(obs.keys())+window_size,window_size):
        
        w.write('{}\t{}\t1.0\n'.format(chrom, window))
        mut.write('{}\t{}\t1.0\n'.format(chrom, window))
        observations.write('{}\t{}\t{}\t{}\n'.format(chrom, window, len(obs[window]), ','.join(obs[window])))



#--------------------------------------------------------------------------------------
# Train and decode the HMM
#--------------------------------------------------------------------------------------

# Initialize your guesses
states = ['Human','Archaic']
starting_probabilities = [0.98, 0.02]
transitions = [[0.999,1-0.999],[0.02,0.980]]
emissions = [0.01, 0.1]

MakeHMMfile(state_names = states, starting_probabilities = starting_probabilities, transitions = transitions, emissions = emissions, outprefix = 'myhmm')

print '-'*40
print 'Started training the model...'
print '-'*40

# Train the HMM
TrainModel(infile = 'observations.txt', outprefix = 'trained', model = 'myhmm.hmm', weights_file = 'weights.txt', mutfile = 'mutrates.txt')

print '-'*40
print 'Started decoding the model...'
print '-'*40


# Decode the HMM
cutoff = 0.9
Decode(infile = 'observations.txt', outprefix = 'decoded', model = 'trained.hmm', weights_file = 'weights.txt', mutfile = 'mutrates.txt', window_size = window_size, cutoff = cutoff)


#--------------------------------------------------------------------------------------
# Checking accuracy
#--------------------------------------------------------------------------------------

print '-'*40
print 'Checking accuracy of the model with a posterior cutoff at {}...'.format(cutoff)
print '-'*40

truth = pybedtools.BedTool('archaic_segments_chr{}.txt'.format(chrom)).sort().merge()
HMM =  pybedtools.BedTool('decoded.bed').sort().merge()

total_HMM = sum([x.stop - x.start for x in (HMM)])
truth_seq =  sum([x.stop - x.start for x in (truth)])
true_positives = sum([x.stop - x.start for x in HMM.intersect(truth)])

print 'Precision =', true_positives/float(total_HMM)*100
print 'Sensitivity =', true_positives/float(truth_seq)*100
```
