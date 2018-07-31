# Introgression-detection
These are the scripts needed to infere archaic introgression in modern human populations. 


### Dependencies
To run the python script you will need numpy. I am using this version of python and numpy:

```
# To run the HMM model
python    version 2.7.12
numpy     version 1.11.1

# If you want to follow my example (running on a 1000 genomes individual)
vcftools  version 0.1.14
tabix     version 0.2.5
```

### Running the scripts

How to train the model
```
python Train.py {observationfile} {output_prefix} {parameterfile} {weightfile} {mutationratefile}
```

How to decode the model
```
python Decode.py {observationfile} {output_prefix} {parameterfile} {weightfile} {mutationratefile} {windowsize}
```

Basically the three files observationfile, weightfile and mutationratefile says how many private snps, how many bases could be called and what the mutation rate is for each window in the genome. A toy example could look like this with column chromosome, window start and value:

```
head observationfile
1      0       0
1      1000    0
1      2000    1
1      3000    2
1      4000    1
1      5000    3
1      6000    0
1      7000    0
1      8000    0
1      9000    0

head weightfile
1      0       1.0
1      1000    1.0
1      2000    1.0
1      3000    1.0
1      4000    1.0
1      5000    1.0
1      6000    0.5
1      7000    0.1
1      8000    0.0
1      9000    1.0

head mutationratefile
1      0       1.5
1      1000    1.5
1      2000    1.5
1      3000    1.5
1      4000    1.5
1      5000    1.5
1      6000    1.5
1      7000    1.5
1      8000    1.5
1      9000    1.5
```

So looking at the observations file we a region from 2 kb to 6 kb with high number of variants. We can also see from the weightfile that we can call all bases in the region 0-6 kb but we could only call half the bases from 6-7 kb. From the mutationratefile we can see that the average mutation rate in this example is 1.5 times higher than the average for the chromosome.

The parameterfile contains the number and names of states and the transition matrix and emission probabilities.

```bash

head parameterfile
# State names (only used for decoding)
states = ['Human','Archaic']

# Initialization parameters (prob of staring in states)
starting_probabilities = [0.98, 0.02]

# transition matrix
transitions = [[0.9995,0.0001],[0.012,0.98]]

# emission matrix (poisson parameter)
emissions = [0.04, 0.4]
```

You dont need to know these parameters in advance because you train the model. The emission values are the average number of variants per window in each state. 




# Example (from VCF file to decoded segments)
I thought it would be nice to have an entire reproduceble example of how to use this model. From a common starting point such as a VCF file to the final output. In this example I will analyse an individual (HG00096) from the 1000 genomes project phase 3. 

First you will need to know which 1) bases can be called in the genome and 2) which variants are found in the outgroup. So I start out by downloading the files from the following directories.

### 1) Which bases could be called?

```bash
Callability file (hg37 - get all the files)
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/StrictMask/

Repeatmask file (hg37 - all files in one zip)
hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/chromFaMasked.tar.gz

Ancestral state file (hg37 alignment e71 - get all the files)
http://web.corral.tacc.utexas.edu/WGSAdownload/resources/human_ancestor_GRCh37_e71/

```

A site has to be callable in the 1000 genomes project (denoted P in the callability file) and not in a repetitive region (denoted N) in the repeatmask file. The small script below does this for chromosome 17:

```bash
python MakeMask.py {repeatmaskfile} {callabilityfile} {windowsize} {chromosomename} {outprefix}

# So for chromosome 17
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

```bash
for file in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
  do echo $file
  cat chr$file.txt >> weights.txt
  cat chr$file.bed >> weights.bed
  done
```

To check everything worked so far look at the first 10 lines:

```bash
head weights.*

==> weights.bed <==
1       756781  756870  Called
1       756879  756920  Called
1       756948  757069  Called
1       757071  757073  Called
1       757077  757246  Called
1       757364  757406  Called
1       757474  757476  Called
1       757765  757766  Called
1       757777  757787  Called
1       757920  757951  Called

==> weights.txt <==
1       0       0.0
1       1000    0.0
1       2000    0.0
1       3000    0.0
1       4000    0.0
1       5000    0.0
1       6000    0.0
1       7000    0.0
1       8000    0.0
1       9000    0.0
```


### 2) Which variants are found in the outgroup

Now we can download the 1000 genomes VCF files and remove all variants found in an outgroup (this case the YRI, ESN and MSL Subsaharan-Africans).

```bash
# VCF files are in this directory
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

# Metadata is here
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
```

The individual names in the outgroup should be in a file (could be called outgroups.txt) with one per line like shown below.

```
NA18486
NA18488
NA18489
NA18498
NA18499
NA18501
NA18502
NA18504
NA18505
```

We first want to get all derived alleles from the outgroup that fall within regions we can call (the weights.bed file that we just made above). Now is a good time to make sure that the "chromosomename" argument is indeed the same as in the vcf file i.e. if in the VCF file the first column is chr1 your "chromosomename" should also have been chr1 and not just 1. 

To do this I made a python script but for big vcffiles I found out tabix and vcftools works much faster. vcftools and tabix can be download from their websites:

https://vcftools.github.io/man_latest.html
http://www.htslib.org/doc/tabix.html


To get the frequency of bi-allelic snps in a given outgroup you can run (make sure to install vcftools and tabix!). I have tabix version 0.2.5 where one uses -B to keep regions from a bed file but the newest version have changed this to -R (http://www.htslib.org/doc/tabix.html). The example below is for the 1000 genomes where each chromosome is in a separate vcf file but you could have a chromosomes in one vcf file and use the weights.bed file instead of a chromosome specific bed file.

```bash
tabix -h ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -B chr17.bed | \
vcftools --vcf - --counts --stdout --keep outgroups.txt --remove-indels --min-alleles 2 --max-alleles 2 > chr17.freq
```

The first 10 lines of the file looks like this:

```bash
head chr17.freq
CHROM   POS     N_ALLELES       N_CHR   {ALLELE:COUNT}
17      439     2       584     C:584   A:0
17      460     2       584     G:583   A:1
17      1102    2       584     T:583   C:1
17      1352    2       584     A:584   T:0
17      1362    2       584     G:584   A:0
17      1382    2       584     G:584   A:0
17      1389    2       584     G:584   A:0
17      1397    2       584     C:469   T:115
17      1398    2       584     G:580   A:4
```

With this output file we can estimate the average mutation rate in a region (I recommend 1,000,000 bp or 100,000 for humans because of this paper http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007254) compared to the whole chromosome average. To do this run the script:

```bash
python Estimate_mutationrate.py {outgroup_frequencyfile} {windowsize_to_calculate_mutationrate_in} {windowsize} {maskfile} {outputfile}

# so for chromosome 17
python Estimate_mutationrate.py chr17.freq 1000000 1000 chr17.txt chr17.mut
```

The first 10 lines of the output file will look like this:

```bash
head chr17.mut
17      0       1.15732924021
17      1000    1.15732924021
17      2000    1.15732924021
17      3000    1.15732924021
17      4000    1.15732924021
17      5000    1.15732924021
17      6000    1.15732924021
17      7000    1.15732924021
17      8000    1.15732924021
17      9000    1.15732924021
```

You can do this for all chromosomes and concatonate like before:

```bash
for file in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
  do echo $file
  cat chr$file.mut >> mutationrates.txt
  done
```


### Getting variants from an individual for training and decoding

The last thing we need to run the scripts is the observation file for an individual. This will for each window show all the variants that individual has where the derived allele is not found in the outgroup and it is in a region we can call. For individual HG00096 we will prepare the observations like this:

```bash
tabix -fh 1000_genomes_phase3/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -B chr17.bed | \
vcftools --vcf - --indv HG00096 --remove-indels --min-alleles 2 --max-alleles 2 --stdout --counts | \
python Filtervariants.py homo_sapiens_ancestor_17.fa chr17.freq 1000 chr17.txt HG00096.chr17.observations.txt

# The python scripts takes the following arguments
python Filtervariants.py {ancestral} {outgroupfrequency} {windowsize} {weightsfile} {output}
```


If you dont have ancestral/derived allele information you can just make a file of all sites where the alternative allele is not found in the outgroup (assuming the reference allele is the ancestral): 

```bash
tabix -fh 1000_genomes_phase3/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -B chr17.bed | \
vcftools --vcf - --indv HG00096 --remove-indels --min-alleles 2 --max-alleles 2 --stdout --counts | \
python FiltervariantsNOancestral.py chr17.freq 1000 chr17.txt HG00096.chr17.observations_NOancestral.txt

# The python scripts takes the following arguments
python FiltervariantsNOancestral.py {outgroupfrequency} {windowsize} {weightsfile} {output}
```

If can now look at the observation file we have created for chromosome 17. This time we will look at the first 30 lines since nothing interestig is going on until then:

```bash
head -30 HG00096.chr17.observations.txt
17      0       0
17      1000    0
17      2000    0
17      3000    0
17      4000    0
17      5000    0
17      6000    0
17      7000    0
17      8000    0
17      9000    0
17      10000   0
17      11000   0
17      12000   0
17      13000   0
17      14000   0
17      15000   0
17      16000   0
17      17000   0
17      18000   0
17      19000   0
17      20000   0
17      21000   0
17      22000   0
17      23000   0
17      24000   0
17      25000   0
17      26000   0
17      27000   0
17      28000   1       28094
17      29000   0

```

We can see that in the window from 28,000 to 29,000 there is one private variant at position 28094. 

We can also look the effect of having ancestral information or not. If we count the number of windows with 1,2,3... variants in HG00096.chr17.observations.txt and HG00096.chr17.observations_NOancestral.txt we find the following:

```bash
for file in  HG00096.observations*
  do echo $file
  cut -f3 $file | sort -n | uniq -c 
  done
  
  
HG00096.chr17.observations_NOancestral.txt
  78986 0
   1871 1
    161 2
     87 3
     41 4
     15 5
     22 6
      5 7
      5 8
      2 9
      1 12
      
HG00096.chr17.observations.txt
  79003 0
   1959 1
    158 2
     51 3
     16 4
      4 5
      3 6
      1 7
      1 8

```

You can see that while the difference is not that great, having the ancestral information is still important. 

We can make observation files for all chromosomes and concatonate them using:
```bash
for file in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
  do echo $file
  cat HG00096.chr$file.observations.txt >> HG00096.observations.txt
  done
```

Now before we run the script lets check that all files needed for the analysis are right. We can check that all of them have the same lines:

```bash
wc -l weights.txt
3036315 weights.txt

wc -l mutationrates.txt
3036315 mutationrates.txt

wc -l HG00096.observations.txt
3036315 HG00096.observations.txt
```

We also have to make a file with some guesses for the starting parameters of the hidden markov model. I have made some in this file:

```bash
head StartingParameters.hmm
# State names (only used for decoding)
states = ['Human','Archaic']

# Initialization parameters (prob of staring in states)
starting_probabilities = [0.98, 0.02]

# transition matrix
transitions = [[0.9995,0.0005],[0.012,0.98]]

# emission matrix (poisson parameter)
emissions = [0.04, 0.1]
```

You can experiment with chosen different starting parameters but most of the time when you train on a whole genome the parameters will converge to some set. What these parameters mean (if we start with emissions) is that in the first state (that we call human) we expect a LOWER density of private variants than in the Archaic. In the human we expect 0.04 snps pr 1000 bp. If we for now just assume all variants are heterozygous then we find 0.02 variants per chromosome. This would mean that the coalescent to the outgroup is given by:
0.02 = Tcoal * L * mutrate = Tcoal * 1000 * 1.25 * 10^-8. If we solve for Tcoal this would mean a coalescence time of 1600 generation (48,000 years if you assume the generation time is 30 year/gen). The coalescence time for state 2 segments (Archaic) with the outgroup is 240,000 year (on the low end but hey we are just guessing now). 

From the transition matrix you can also see that we think probability of staying within an Archaic segment is 0.98. That means that the states will on average be 50 windows long (or 50 kb since the window size is 1000 bp). That would mean that the introgression event happend around 2000 generations ago. 

### Training the model

If you have made it until here then congratulations now it is (finally) time to train the model. We run the script:

```bash
python Train.py HG00096.observations.txt HG00096_trained StartingParameters.hmm weights.txt mutationrates.txt
```

The model will create two files. One is called HG00096_trained.log where it report the parameters and likelihood of the model for each iteration and HG00096_trained.hmm which is the same format as StartingParameters.hmm (just with the parameters that optimize the likelihood).


### Decoding the model

Now that we have a set of trained parameters that maximize the likelihood we can decode the model using the following script:

```bash
python Decode.py HG00096.observations.txt HG00096_decoded HG00096_trained.hmm weights.txt mutationrates.txt 1000




```

This will also produce two files. One is HG00096_decoded.Summary which is like a bedfile and tell you what part of the sequence belong to different states. It also tells you how many snps that are in each segment and what the average posterior probability for being in that segment is. The other is HG00096_decoded.All_posterior_probs.txt and this is a window by window assignment to each state. 

The files look like this:

```bash
head HG00096_decoded.Summary.txt
name    chrom   start   end     length  state   snps    mean_prob
HG00096_decoded 1       0       3424000 3425000 Human   114     0.975435021019
HG00096_decoded 1       3425000 3449000 25000   Archaic 19      0.956203499367
HG00096_decoded 1       3450000 4259000 810000  Human   24      0.984532458656
HG00096_decoded 1       4260000 4367000 108000  Archaic 20      0.792893841425
HG00096_decoded 1       4368000 6142000 1775000 Human   54      0.989351192407
HG00096_decoded 1       6143000 6150000 8000    Archaic 6       0.822766806209
HG00096_decoded 1       6151000 9319000 3169000 Human   81      0.990871350063
HG00096_decoded 1       9320000 9361000 42000   Archaic 10      0.853517013725
HG00096_decoded 1       9362000 12595000        3234000 Human   70      0.987032899434

head HG00096_decoded.All_posterior_probs.txt
chrom   start   observations    Mostlikely      Human   Archaic
1       0       0       Human   0.714740575355  0.285258832131
1       1000    0       Human   0.719097642964  0.280901764546
1       2000    0       Human   0.723371380417  0.276628027107
1       3000    0       Human   0.727563381381  0.27243602611
1       4000    0       Human   0.731675209186  0.268324198332
1       5000    0       Human   0.735708397086  0.264291010418
1       6000    0       Human   0.739664449158  0.26033495837
1       7000    0       Human   0.743544840571  0.256454566945
1       8000    0       Human   0.747351018397  0.252648389099
```





