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
python Train.py {observationfile} {output_prefix} {parameterfile} {weightfile} {mutationratefile}
```

How to decode the model
```
python Decode.py {observationfile} {output_prefix} {parameterfile} {weightfile} {mutationratefile}
```


## Example (from VCF file to decoded segments)
I thought it would be nice to have an entire reproduceble example of how to use this model. From a common starting point such as a VCF file to the final output. In this example I will analyse an individual (HG00096) from the 1000 genomes project phase 3. 

First you will need to know which 1) bases can be called in the genome and 2) which variants are found in the outgroup. So I start out by downloading the files from the following directories.

### 1) Which bases could be called?

```
Callability file (hg37 - get all the files)
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/StrictMask/

Repeatmask file (hg37 - all files in one zip)
hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/chromFaMasked.tar.gz

Ancestral state file (hg37 alignment e71 - get all the files)
http://web.corral.tacc.utexas.edu/WGSAdownload/resources/human_ancestor_GRCh37_e71/

```

A site has to be callable in the 1000 genomes project (denoted P in the callability file) and not in a repetitive region (denoted N) in the repeatmask file. The small script below does this for chromosome 17:

```
python MakeMask.py {repeatmaskfile} {callabilityfile} {windowsize} {chromosomename} {outprefix}

So for chromosome 17
python MakeMask.py chr17.fa.masked 20140520.chr17.strict_mask.fasta.gz 1000 17 chr17
```

This will generate two files. One is chr17.txt which reports the fraction of called bases for each window and chr17.bed which is a bedfile of all the callable regions. You can see the first 10 lines of each file below. Notice how there was 68 (67 from 416 to 482 and 1 from 864 to 864) bases called between position 0 and position 1000.


```
head chr17.txt
17      0       0.068
17      1000    0.642
17      2000    0.662
17      3000    0.377
17      4000    0.058
17      5000    0.723
17      6000    0.528
17      7000    0.729
17      8000    0.494
17      9000    0.151

head chr17.bed
17      416     482     Called
17      864     864     Called
17      1086    1105    Called
17      1330    1951    Called
17      2024    2278    Called
17      2593    3212    Called
17      3508    3671    Called
17      4821    4878    Called
17      5082    5515    Called
17      5616    5769    Called
```

You can do it for all chromosomes (remember to change the outprefix name) and then concatenate them using cat. 

```
cat chr*.txt > weights.txt
cat chr*.bed > weights.bed
```




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
name	start	end	length	state	snps	mean_prob
observations.txt	1	285	285	Human	5	0.96009303904
observations.txt	286	300	15	Archaic	2	0.518621392803
observations.txt	301	699	399	Human	8	0.910527915201
observations.txt	700	760	61	Archaic	28	0.809617990419
observations.txt	761	1812	1052	Human	18	0.993446178829
observations.txt	1813	1843	31	Archaic	10	0.953483598978
observations.txt	1844	3190	1347	Human	28	0.986725347829
observations.txt	3191	3385	195	Archaic	26	0.951324931506
observations.txt	3386	5416	2031	Human	36	0.995510947589
```







