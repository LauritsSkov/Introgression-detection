# Introgression-detection
These are the scripts needed to infere archaic introgression in modern human populations. 

## Dependencies
To run the python script you will need numpy. I am using this version of python and numpy:

```
python version 2.7.12
numpy  version 1.11.1
```

## Running the scripts

How to train the model
```
python Train.py <path to observations file> <output_prefix> <path to markov model> <path to callability file> <path to mutation rate file>
```

How to decode the model
```
python Decode.py <path to observations file> <output_prefix> <path to markov model> <path to callability file> <path to mutation rate file>
```


## Test data
In the folder Test is the files needed to run the analysis. You will need a file with observations, which is the number of private variants for each window. I have chosen a window size of 1000 bp - meaning that the first number in the file will be how many private variants that is observed in the first 1000 bp of the chromosome 17. The weights.txt files is the number of called basepairs for each window and the Mutrates.txt file is the local mutation rate in this window. If the local mutation rate or number of callable bases are not known one can create a mutation rate and callability file of the same length as the observations.txt file and fill it out with 1.0. This corresponds to having called all bases in a window and have the genome-wide average mutation rate. 

In the .hmm file you can insert your initial guesses of coalescene time for the different states into the outgroup and the time of admixure. 


To train the model on data from chromosome 17 for a Papuan individual write:

```
python Train.py observations.txt trained Model1.hmm weights.txt Mutrates.txt
```

This will return two files. One log file shows the emission and transition probabilities for each iteration and .hmm file is the final parameters found. This .hmm can be used for posterior decoding.

The trained.hmm file should look like this:
```
# State names (only used for decoding)
states = ['Human','Archaic']

# Initialization parameters (prob of staring in states)
starting_probabilities = [0.98369416906408647, 0.016305830840918301]

# transition matrix
transitions = [[0.999206173415,0.000793826585448],[0.0114703208336,0.988529679166]]

# emission matrix (poisson parameter)
emissions = [0.041525372286415022, 0.38744095667508627]


```

From the emission values it can be seen that the number of snps is around 10 times higher in the archaic state compare to the human state. If one a mutation rate of 1.25E-8 and a window size of 1000 bp an emission value of 0.38 corresponds to 30,000 generation until coalescene with the outgroup (0.38/(1000*1.25E-8)). The transition probability of staying within the archaic state is 0.98 meaning that the segments are around 50 kb on average. This correspond to the admixture event happening around 2000 generations ago. 



For decoding the sequence write:

```
python Decode.py observations.txt decoded trained.hmm weights.txt Mutrates.txt
```

This will return three files, the decoded.Summary.txt, decoded.All_posterior_probs.txt and the decoded.mostLikely_probs.txt. The decoded.All_posterior_probs.txt will have a column with the number of private snps plus a column for each state, stating what the probability of being in this state is. The decoded.mostLikely_probs.txt has two columns: One with the number of private snps and one with the most likely state for this position. The decoded.Summary.txt has 6 columns. One with the input file name, the starting of a segment, the end of the segment, the length of the segment, the name of the segment and the number of private variants in the segment.

The first few lines in decoded.Summary.txt are shown below:

```
name	start	end	length	state	snps
observations.txt	1	285	285	Human	5
observations.txt	286	300	15	Archaic	2
observations.txt	301	699	399	Human	8
observations.txt	700	760	61	Archaic	28
observations.txt	761	1812	1052	Human	18
observations.txt	1813	1843	31	Archaic	10
observations.txt	1844	3190	1347	Human	28
observations.txt	3191	3385	195	Archaic	26
observations.txt	3386	5416	2031	Human	36
```







