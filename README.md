# Introgression-detection
These are the scripts needed to infere archaic introgression in modern human populations. 

## Dependencies
python version 2.7.12
numpy  version 1.11.1

## Test data
To train the model on data from chromosome 17 for a Papuan individual write:

```
python Train.py observations.txt trained Model1.hmm weights.txt
```

This will return two files. One log file shows the emission and transition probabilities for each iteration and .hmm file is the final parameters found. This .hmm can be used for posterior decoding.

The trained.hmm file should look like this:
```
# State names (only used for decoding)
states = ['Human','Archaic']

# Initialization parameters (prob of staring in states)
starting_probabilities = [0.98263395309868673, 0.017366047186862847]

# transition matrix
transitions = [[0.999186089122,0.000813910878183],[0.0113823793075,0.988617620692]]

# emission matrix (poisson parameter)
emissions = [0.044167227445329589, 0.41773283862878263]

```



For decoding the sequence write:

```
python Decode.py observations.txt decoded trained.hmm weights.txt
```

This will return three files, the decoded.Summary.txt, decoded.All_posterior_probs.txt and the decoded.mostLikely_probs.txt. The decoded.All_posterior_probs.txt will have a column with the number of private snps plus a column for each state, stating what the probability of being in this state is. The decoded.mostLikely_probs.txt has two columns: One with the number of private snps and one with the most likely state for this position. The decoded.Summary.txt has 6 columns. One with the input file name, the starting of a segment, the end of the segment, the length of the segment, the name of the segment and the number of private variants in the segment.

The first few lines in decoded.Summary.txt are shown below:

```
observations.txt	1	697	697	Human	15
observations.txt	698	764	67	Archaic	28
observations.txt	765	1811	1047	Human	18
observations.txt	1812	1845	34	Archaic	10
observations.txt	1846	3186	1341	Human	28
observations.txt	3187	3386	200	Archaic	26
observations.txt	3387	5415	2029	Human	36
observations.txt	5416	5457	42	Archaic	7
observations.txt	5458	5961	504	Human	14
observations.txt	5962	5971	10	Archaic	4
```







