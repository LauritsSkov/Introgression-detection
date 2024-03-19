
# Contact
https://sites.google.com/view/laurits-skov

---

# Helpful files
The outgroup files, mutation rate files and callability files are now premade! 

https://doi.org/10.5281/zenodo.10806733 (hg19 and hg38)

VCF file containing 4 high coverage archaic genomes (Altai, Vindija and Chagyrskaya Neanderthals and Denisovan) here: 

https://zenodo.org/records/7246376 (hg19) 

https://zenodo.org/records/10806726 (hg38)

---

# Introgression detection

If you are working with archaic introgression into present-day humans of non-African ancestry you can use these files and skip the following steps:
Find derived variants in outgroup and Estimate local mutation rate. 

These are the scripts needed to infere archaic introgression in modern populations using an unadmixed outgroup.

1. [Installation](#installation)
2. [Usage](#usage)
3. [Quick tutorial](#quick-tutorial)
4. [1000 genomes tutorial](#example-with-1000-genomes-data)
      - [Get data](#getting-data)
      - [Find derived variants in outgroup](#finding-snps-which-are-derived-in-the-outgroup)
      - [Estimate local mutation rate](#estimating-mutation-rate-across-genome)
      - [Find variants in ingroup](#find-a-set-of-variants-which-are-not-derived-in-the-outgroup)
      - [Train the HMM](#training)
      - [Decoding](#decoding)
      - [Phased data](#training-and-decoding-with-phased-data)
      - [Annotate](#annotate-with-known-admixing-population)
      - [Run in python](#annotate-with-known-admixing-population)

---

## Installation

Run the following command to install:

```bash
pip install hmmix 
```

If you want to work with bcf/vcf files you should also install vcftools and bcftools. You can either use conda or visit their websites.

```bash
conda install -c bioconda vcftools bcftools
```

![Overview of model](https://user-images.githubusercontent.com/30321818/43464826-4d11d46c-94dc-11e8-8f1a-6851aa5d9125.jpg)

The way the model works is by removing variation found in an outgroup population and then using the remaining variants to group the genome into regions of different variant density. If the model works well we would expect that introgressed regions have higher variant density than non-introgressed - because they have spend more time accumulation variation that is not found in the outgroup.

An example on simulated data is provided below:

![het_vs_archaic](https://user-images.githubusercontent.com/30321818/46877046-217eff80-ce40-11e8-9010-edb544e3e1ee.png)

In this example we zoom in on 1 Mb of simulated data for a haploid genome. The top panel shows the coalescence times with the outgroup across the region and the green segment is an archaic introgressed segment. Notice how much more deeper the coalescence time with the outgroup is. The second panel shows that probability of being in the archaic state. We can see that the probability is much higher in the archaic segment, demonstrating that in this toy example the model is working like we would hope. The next panel is the snp density if you dont remove all snps found in the outgroup. By looking at this one cant tell where the archaic segments begins and ends, or even if there is one. The bottom panel is the snp density when all variation in the outgroup is removed. Notice that now it is much clearer where the archaic segment begins and ends!

The method is now published in PlosGenetics and can be found here: [Detecting archaic introgression using an unadmixed outgroup](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007641) This paper is describing and evaluating the method.

---

## Usage

```note
Script for identifying introgressed archaic segments

> Turorial:
hmmix make_test_data 
hmmix train  -obs=obs.txt -weights=weights.bed -mutrates=mutrates.bed -param=Initialguesses.json -out=trained.json 
hmmix decode -obs=obs.txt -weights=weights.bed -mutrates=mutrates.bed -param=trained.json


Different modes (you can also see the options for each by writing hmmix make_test_data -h):
> make_test_data        
    -windows            Number of Kb windows to create (defaults to 50,000)
    -nooutfiles         Don't create obs.txt, mutrates.bed, weights.bed, Initialguesses.json (defaults to yes)

> mutation_rate         
    -outgroup           [required] path to variants found in outgroup
    -out                outputfile (defaults to mutationrate.bed)
    -weights            file with callability (defaults to all positions being called)
    -window_size        size of bins (defaults to 1 Mb)

> create_outgroup       
    -ind                [required] ingroup/outgrop list (json file) or comma-separated list e.g. ind1,ind2
    -vcf                [required] path to list of comma-separated vcf/bcf file(s) or wildcard characters e.g. chr*.bcf
    -weights            file with callability (defaults to all positions being called)
    -out                outputfile (defaults to stdout)
    -ancestral          fasta file with ancestral information - comma-separated list or wildcards like vcf argument (default none)
    -refgenome          fasta file with reference genome - comma-separated list or wildcards like vcf argument (default none)

> create_ingroup        
    -ind                [required] ingroup/outgrop list (json file) or comma-separated list e.g. ind1,ind2
    -vcf                [required] path to list of comma-separated vcf/bcf file(s) or wildcard characters e.g. chr*.bcf
    -outgroup           [required] path to variant found in outgroup
    -weights            file with callability (defaults to all positions being called)
    -out                outputfile prefix (default is a file named obs.<ind>.txt where ind is the name of individual in ingroup/outgrop list)
    -ancestral          fasta file with ancestral information - comma-separated list or wildcards like vcf argument (default none)

> train                 
    -obs                [required] file with observation data
    -weights            file with callability (defaults to all positions being called)
    -mutrates           file with mutation rates (default is mutation rate is uniform)
    -param              markov parameters file (default is human/neanderthal like parameters)
    -out                outputfile (default is a file named trained.json)
    -window_size        size of bins (default is 1000 bp)
    -haploid            Change from using diploid data to haploid data (default is diploid)

> decode                
    -obs                [required] file with observation data
    -weights            file with callability (defaults to all positions being called)
    -mutrates           file with mutation rates (default is mutation rate is uniform)
    -param              markov parameters file (default is human/neanderthal like parameters)
    -out                outputfile prefix <out>.hap1.txt and <out>.hap2.txt if -haploid option is used or <out>.diploid.txt (default is stdout)
    -window_size        size of bins (default is 1000 bp)
    -haploid            Change from using diploid data to haploid data (default is diploid)
    -admixpop ADMIXPOP  Annotate using vcffile with admixing population (default is none)
    -extrainfo          Add variant position for each SNP (default is off)
```

---

## Quick tutorial

Here is how we can simulate test data using hmmix. Lets make some test data and start using the program.

```note
> hmmix make_test_data
creating 2 chromosomes each with 50000 kb of test data with the following parameters..

> state_names = ['Human', 'Archaic']
> starting_probabilities = [0.98, 0.02]
> transitions = [[1.0, 0.0], [0.02, 0.98]]
> emissions = [0.04, 0.4]
```

This will generate 4 files, obs.txt, weights.bed, mutrates.bed and Initialguesses.json. obs.txt. These are the mutations that are left after removing variants which are found in the outgroup.

```note
chrom  pos     ancestral_base  genotype
chr1   5212    A               AG
chr1   32198   A               AG
chr1   65251   C               CG
chr1   117853  A               AG
chr1   122518  T               TC
chr1   142322  T               TC
chr1   144695  C               CG
chr1   206370  T               TG
chr1   218969  A               AT
```

weights.bed. This is the parts of the genome that we can accurately map to - in this case we have simulated the data and can accurately access the entire genome.

```note
chr1   0   50000000
chr2   0   50000000
```

mutrates.bed. This is the normalized mutation rate across the genome (in bins of 1 Mb).

```note
chr1  0        1000000   1
chr1  1000000  2000000   1
chr1  2000000  3000000   1
chr1  3000000  4000000   1
chr1  4000000  5000000   1
chr1  5000000  6000000   1
chr1  6000000  7000000   1
chr1  7000000  8000000   1
chr1  8000000  9000000   1
chr1  9000000  10000000  1
```

Initialguesses.json. This is our initial guesses when training the model - note these are different from those we simulated from.

```json
{
  "state_names": ["Human","Archaic"],
  "starting_probabilities": [0.5,0.5],
  "transitions": [[0.99,0.01],[0.02,0.98]],
  "emissions": [0.03,0.3]
}
```

We can find the best fitting parameters using BaumWelsch training. Here is how you use it: - note you can try to ommit the weights and mutrates arguments. Since this is simulated data the mutation is constant across the genome and we can asses the entire genome. Also notice how the parameters approach the parameters the data was generated from (jubii).

```note
> hmmix train  -obs=obs.txt -weights=weights.bed -mutrates=mutrates.bed -param=Initialguesses.json -out=trained.json
----------------------------------------
> state_names = ['Human', 'Archaic']
> starting_probabilities = [0.5, 0.5]
> transitions = [[0.99, 0.01], [0.02, 0.98]]
> emissions = [0.03, 0.3]
> number of windows: 99970 . Number of snps =  4230
> total callability: 1.0
> average mutation rate per bin: 1.0
> Output is trained.json
> Window size is 1000 bp
> Haploid False
----------------------------------------
iteration  loglikelihood  start1  start2  emis1   emis2   trans1_1  trans2_2
0          -18123.4432    0.5     0.5     0.03    0.3     0.99      0.98
1          -17506.017     0.96    0.04    0.035   0.2202  0.9969    0.9242
2          -17487.7921    0.971   0.029   0.0369  0.2235  0.9974    0.9141
...
16         -17401.3802    0.994   0.006   0.0398  0.4584  0.9999    0.9806
17         -17401.3786    0.994   0.006   0.0398  0.4586  0.9999    0.9807
18         -17401.3783    0.994   0.006   0.0398  0.4587  0.9999    0.9808


# run without mutrate and weights (only do this for simulated data)
> hmmix train  -obs=obs.txt -param=Initialguesses.json -out=trained.json
```

We can now decode the data with the best parameters that maximize the likelihood and find the archaic segments:

```note
> hmmix decode -obs=obs.txt -weights=weights.bed -mutrates=mutrates.bed -param=trained.json
----------------------------------------
> state_names = ['Human', 'Archaic']
> starting_probabilities = [0.994, 0.006]
> transitions = [[1.0, 0.0], [0.019, 0.981]]
> emissions = [0.04, 0.459]
> number of windows: 99970 . Number of snps =  4230
> total callability: 1.0
> average mutation rate per bin: 1.0
> Output is /dev/stdout
> Window size is 1000 bp
> Haploid False
----------------------------------------
chrom  start     end       length    state    mean_prob  snps
chr1   0         7233000   7234000   Human    0.9995     287
chr1   7234000   7246000   13000     Archaic  0.90427    9
chr1   7247000   21618000  14372000  Human    0.99946    610
chr1   21619000  21673000  55000     Archaic  0.9697     22
chr1   21674000  26859000  5186000   Human    0.99878    204
chr1   26860000  26941000  82000     Archaic  0.971      36
chr1   26942000  49989000  23048000  Human    0.99982    863
chr2   0         6793000   6794000   Human    0.99972    237
chr2   6794000   6822000   29000     Archaic  0.95461    14
chr2   6823000   12646000  5824000   Human    0.99927    244
chr2   12647000  12745000  99000     Archaic  0.97413    55
chr2   12746000  15461000  2716000   Human    0.99881    125
chr2   15462000  15547000  86000     Archaic  0.93728    38
chr2   15548000  32626000  17079000  Human    0.99951    709
chr2   32627000  32695000  69000     Archaic  0.98305    31
chr2   32696000  41087000  8392000   Human    0.9995     360
chr2   41088000  41178000  91000     Archaic  0.96092    43
chr2   41179000  49952000  8774000   Human    0.99789    328
chr2   49953000  49977000  25000     Archaic  0.98501    13

# Again here you could ommit weights and mutationrates. Actually one could also ommit trained.json because then the model defaults to using the parameters we used the generated the data
> hmmix decode -obs=obs.txt
```

---

## Example with 1000 genomes data

---

The whole pipeline we will run looks like this. In the following section we will go through all the steps on the way
NOTE: The outgroup files, mutation rate files and callability files are now premade! They can be downloaded in hg38 and hg19 here: https://doi.org/10.5281/zenodo.10806733
But keep reading along if you want to know HOW the files were generated!

```note
hmmix create_outgroup -ind=individuals.json -vcf=*.bcf -weights=strickmask.bed -out=outgroup.txt -ancestral=homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_*.fa -refgenome=referencegenome/*fa
hmmix mutation_rate -outgroup=outgroup.txt  -weights=strickmask.bed -window_size=1000000 -out mutationrate.bed
hmmix create_ingroup  -ind=individuals.json -vcf=*.bcf -weights=strickmask.bed -out=obs -outgroup=outgroup.txt -ancestral=homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_*.fa
hmmix train  -obs=obs.HG00096.txt -weights=strickmask.bed -mutrates=mutationrate.bed -out=trained.HG00096.json 
hmmix decode -obs=obs.HG00096.txt -weights=strickmask.bed -mutrates=mutationrate.bed -param=trained.HG00096.json 
```

### Getting data

I thought it would be nice to have an entire reproduceble example of how to use this model. From a common starting point such as a VCF file (well a BCF file in this case) to the final output.  The reason for using BCF files is because it is MUCH faster to extract data for each individual. You can convert a vcf file to a bcf file like this:

```note
bcftools view file.vcf -l 1 -O b > file.bcf
bcftools index file.bcf
```

In this example I will analyse an individual (HG00096) from the 1000 genomes project phase 3. All analysis are run on my lenovo thinkpad (8th gen) computer so it should run on yours too!

First we will need to know which 1) bases can be called in the genome and 2) which variants are found in the outgroup. So let's start out by downloading the files from the following directories.
To download callability regions, ancestral alleles information, ingroup outgroup information call this command:

```bash
# bcffiles (hg19)
ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/bcf_files/

# callability (remember to remove chr in the beginning of each line to make it compatible with hg19 e.g. chr1 > 1)
ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/20141020.strict_mask.whole_genome.bed
sed 's/^chr\|%$//g' 20141020.strict_mask.whole_genome.bed | awk '{print $1"\t"$2"\t"$3}' > strickmask.bed

# outgroup information
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

# Ancestral information
ftp://ftp.ensembl.org/pub/release-74/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71.tar.bz2

# Reference genome
wget 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz' -O chromFa.tar.gz

# Archaic variants (Altai, Vindija, Chagyrskaya and Denisova in hg19)
https://zenodo.org/record/7246376#.Y1cRBkrMJH4
```

For this example we will use all individuals from 'YRI','MSL' and 'ESN' as outgroup individuals. While we will only be decoding hG00096 in this example you can add as many individuals as you want to the ingroup.  

```json
{
  "ingroup": [
    "HG00096",
    "HG00097"
  ],
  "outgroup": [
    "HG02922",
    "HG02923",
    ...
    "HG02944",
    "HG02946"]
}
```

---

### Finding snps which are derived in the outgroup

First we need to find a set of variants found in the outgroup. We can use the wildcard character to loop through all bcf files. If you dont have ancestral information you can skip the ancestral argument.

```bash
(took an hour) > hmmix create_outgroup -ind=individuals.json -vcf=*.bcf -weights=strickmask.bed -out=outgroup.txt -ancestral=homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_*.fa

# Alternative usage (if you only have a few individual in the outgroup you can also provide a comma separated list)
> hmmix create_outgroup -ind=HG02922,HG02923,HG02938 -vcf=*.bcf -weights=strickmask.bed -out=outgroup.txt -ancestral=homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_*.fa

# Alternative usage (if you have no ancestral information)
> hmmix create_outgroup -ind=individuals.json -vcf=*.bcf -weights=strickmask.bed -out=outgroup.txt 

# Alternative usage (if you only want to run the model on a subset of chromosomes, with or without ancestral information)
> hmmix create_outgroup -ind=individuals.json -vcf=chr1.bcf,chr2.bcf -weights=strickmask.bed -out=outgroup.txt

> hmmix create_outgroup -ind=individuals.json -vcf=chr1.bcf,chr2.bcf -weights=strickmask.bed -out=outgroup.txt -ancestral=homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_1.fa,homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_2.fa
```

Something to note is that if you use an outgroup vcffile (like 1000 genomes) and an ingroup vcf file from a different dataset (like SGDP) there is an edge case which could occur. There could be recurrent mutations where every individual in 1000 genome has the derived variant and one individual in SGDP where the derived variant has mutated back to the ancestral allele. This means that this position will not be present in the outgroup file. However if a recurrent mutation occurs it will look like multiple individuals in the ingroup file have the mutation. This does not happen often but just in case you can create the outgroup file and adding the sites which are fixed derived in all individuals using the reference genome:

```bash
# Alternative usage (if you want to remove sites which are fixed derived in your outgroup/ingroup)
> hmmix create_outgroup -ind=individuals.json -vcf=*.bcf -weights=strickmask.bed -out=outgroup.txt -ancestral=homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_*.fa -refgenome=*fa
```

---

### Estimating mutation rate across genome

We can use the number of variants in the outgroup to estimate the substitution rate as a proxy for mutation rate.

```note
(took 30 sec) > hmmix mutation_rate -outgroup=outgroup.txt  -weights=strickmask.bed -window_size=1000000 -out mutationrate.bed
----------------------------------------
> Outgroupfile: outgroup.txt
> Outputfile is: mutationrate.bed
> Callability file is: strickmask.bed
> Window size: 1000000
----------------------------------------
```

---

### Find a set of variants which are not derived in the outgroup

Keep variants that are not found to be derived in the outgroup for each individual in ingroup. You can also speficy a single individual or a comma separated list of individuals.

```note
(took 20 min) > hmmix create_ingroup  -ind=individuals.json -vcf=*.bcf -weights=strickmask.bed -out=obs -outgroup=outgroup.txt -ancestral=homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_*.fa
----------------------------------------
> Ingroup individuals: 2
> Using vcf and ancestral files
vcffile: chr1.bcf ancestralfile: homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_1.fa
vcffile: chr2.bcf ancestralfile: homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_2.fa
vcffile: chr3.bcf ancestralfile: homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_3.fa
vcffile: chr4.bcf ancestralfile: homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_4.fa
vcffile: chr5.bcf ancestralfile: homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_5.fa
vcffile: chr6.bcf ancestralfile: homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_6.fa
vcffile: chr7.bcf ancestralfile: homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_7.fa
vcffile: chr8.bcf ancestralfile: homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_8.fa
vcffile: chr9.bcf ancestralfile: homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_9.fa
vcffile: chr10.bcf ancestralfile: homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_10.fa
vcffile: chr11.bcf ancestralfile: homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_11.fa
vcffile: chr12.bcf ancestralfile: homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_12.fa
vcffile: chr13.bcf ancestralfile: homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_13.fa
vcffile: chr14.bcf ancestralfile: homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_14.fa
vcffile: chr15.bcf ancestralfile: homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_15.fa
vcffile: chr16.bcf ancestralfile: homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_16.fa
vcffile: chr17.bcf ancestralfile: homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_17.fa
vcffile: chr18.bcf ancestralfile: homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_18.fa
vcffile: chr19.bcf ancestralfile: homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_19.fa
vcffile: chr20.bcf ancestralfile: homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_20.fa
vcffile: chr21.bcf ancestralfile: homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_21.fa
vcffile: chr22.bcf ancestralfile: homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_22.fa

> Using outgroup variants from: outgroup.txt 
> Callability file: strickmask.bed 
> Writing output to file with prefix: obs.<individual>.txt

----------------------------------------
Running command:
bcftools view -m2 -M2 -v snps -s HG00096 -T strickmask.bed chr1.bcf | vcftools --vcf - --exclude-positions outgroup.txt --recode --stdout 
...
bcftools view -m2 -M2 -v snps -s HG00097 -T strickmask.bed chr22.bcf | vcftools --vcf - --exclude-positions outgroup.txt --recode --stdout 


# Different way to define which individuals are in the ingroup
(took 20 min) > hmmix create_ingroup  -ind=HG00096,HG00097 -vcf=*.bcf -weights=strickmask.bed -out=obs -outgroup=outgroup.txt -ancestral=homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_*.fa
```

---

### Training

Now for training the HMM parameters and decoding

```note
(took 2 min) > hmmix train  -obs=obs.HG00096.txt -weights=strickmask.bed -mutrates=mutationrate.bed -out=trained.HG00096.json 
----------------------------------------
> state_names = ['Human', 'Archaic']
> starting_probabilities = [0.98, 0.02]
> transitions = [[1.0, 0.0], [0.02, 0.98]]
> emissions = [0.04, 0.4]
> number of windows: 2877010 . Number of snps =  129734
> total callability: 0.72
> average mutation rate per bin: 1.0
> Output is trained.HG00096_new.json
> Window size is 1000 bp
> Haploid False
----------------------------------------
iteration  loglikelihood  start1  start2  emis1   emis2   trans1_1  trans2_2
0          -492288.9165   0.98    0.02    0.04    0.4     0.9999    0.98
1          -488899.7271   0.961   0.039   0.0483  0.3913  0.9994    0.986
2          -488755.4068   0.958   0.042   0.0481  0.3911  0.9993    0.9835
...
19         -488642.2737   0.954   0.046   0.047   0.3862  0.9989    0.9771
20         -488642.2723   0.954   0.046   0.047   0.3862  0.9989    0.9771
21         -488642.2716   0.954   0.046   0.047   0.3862  0.9989    0.9771
```

---

### Decoding

```note
(took 30 sec) > hmmix decode -obs=obs.HG00096.txt -weights=strickmask.bed -mutrates=mutationrate.bed -param=trained.HG00096.json 
----------------------------------------
> state_names = ['Human', 'Archaic']
> starting_probabilities = [0.954, 0.046]
> transitions = [[0.999, 0.001], [0.023, 0.977]]
> emissions = [0.047, 0.386]
> number of windows: 2877010 . Number of snps =  129734
> total callability: 0.72
> average mutation rate per bin: 1.0
> Output prefix is /dev/stdout
> Window size is 1000 bp
> Haploid False
----------------------------------------
chrom  start    end      length   state    mean_prob  snps
1      0        2987000  2988000  Human    0.98484    91
1      2988000  2996000  9000     Archaic  0.71815    6
1      2997000  3424000  428000   Human    0.98944    30
1      3425000  3451000  27000    Archaic  0.95652    22
1      3452000  4301000  850000   Human    0.98203    36
1      4302000  4360000  59000    Archaic  0.84636    20
1      4361000  4499000  139000   Human    0.97136    4
1      4500000  4509000  10000    Archaic  0.84456    7
```

You can also save to an output file with the command:

```note
hmmix decode -obs=obs.HG00096.txt -weights=strickmask.bed -mutrates=mutationrate.bed -param=trained.HG00096.json -out=HG00096.decoded
```

This will create a file named HG00096.decoded.diploid.txt because the default option is treating the data as diploid (more on haploid decoding in next chapter)

---

### Training and decoding with phased data

It is also possible to tell the model that the data is phased with the -haploid parameter. For that we first need to train the parameters for haploid data and then decode. Training the model on phased data is done like this - and we also remember to change the name of the parameter file to include phased so future versions of ourselves don't forget. Another thing to note is that the number of snps is bigger than before 134799 vs 129149. This is because the program is counting snps on both haplotypes and homozygotes will be counted twice!

```note
(took 4 min) > hmmix train  -obs=obs.HG00096.txt -weights=strickmask.bed -mutrates=mutationrate.bed -out=trained.HG00096.phased.json -haploid
----------------------------------------
> state_names = ['Human', 'Archaic']
> starting_probabilities = [0.98, 0.02]
> transitions = [[1.0, 0.0], [0.02, 0.98]]
> emissions = [0.04, 0.4]
> number of windows: 5754020 . Number of snps =  135411
> total callability: 0.72
> average mutation rate per bin: 1.0
> Output is trained.HG00096.phased.json
> Window size is 1000 bp
> Haploid True
----------------------------------------
iteration  loglikelihood  start1  start2  emis1   emis2   trans1_1  trans2_2
0          -597753.9141   0.98    0.02    0.04    0.4     0.9999    0.98
1          -585032.3914   0.983   0.017   0.0261  0.4026  0.9998    0.9853
2          -584451.0968   0.979   0.021   0.0252  0.373   0.9996    0.9826
...
19         -584134.532    0.972   0.028   0.0241  0.3367  0.9993    0.9758
20         -584134.5305   0.972   0.028   0.0241  0.3366  0.9993    0.9758
21         -584134.5297   0.972   0.028   0.0241  0.3366  0.9993    0.975
```

Below I am only showing the first archaic segments on chromosome 1 for each haplotype (note you have to scroll down after chrom 22 before the new haplotype begins). The seem to fall more or less in the same places as when we used diploid data.

```note
(took 30 sec) > hmmix decode -obs=obs.HG00096.txt -weights=strickmask.bed -mutrates=mutationrate.bed -param=trained.HG00096.phased.json -haploid
----------------------------------------
> state_names = ['Human', 'Archaic']
> starting_probabilities = [0.972, 0.028]
> transitions = [[0.999, 0.001], [0.024, 0.976]]
> emissions = [0.024, 0.337]
> number of windows: 5754020 . Number of snps =  135411
> total callability: 0.72
> average mutation rate per bin: 1.0
> Output prefix is HG00096.decoded
> Window size is 1000 bp
> Haploid True
----------------------------------------
hap1
chrom  start    end      length  state    mean_prob  snps
1      2162000  2184000  23000   Archaic  0.61055    6
1      3425000  3451000  27000   Archaic  0.96595    22

...
hap2
1      2780000  2802000  23000   Archaic  0.61948    7
1      4302000  4336000  35000   Archaic  0.94008    13
1      4500000  4510000  11000   Archaic  0.87592    7
1      4989000  4999000  11000   Archaic  0.57897    4
```

You can also save to an output file with the command:

```note
hmmix decode -obs=obs.HG00096.txt -weights=strickmask.bed -mutrates=mutationrate.bed -param=trained.HG00096.phased.json -haploid -out=HG00096.decoded
```

This will create two files named HG00096.decoded.hap1.txt and HG00096.decoded.hap2.txt

---

### Annotate with known admixing population

Even though this method does not use archaic reference genomes for finding segments you can still use them to annotate your segments. 

I have uploaded a VCF file containing 4 high coverage archaic genomes (3 Neanderthals and 1 Denisovan) here:

https://zenodo.org/records/7246376 (hg19 - the one I use in this example)

https://zenodo.org/records/10806726 (hg38)

If you have a vcf from the population that admixed in VCF/BCF format you can write this:

```note
> hmmix decode -obs=obs.HG00096.txt -weights=strickmask.bed -mutrates=mutationrate.bed -param=trained.HG00096.json -admixpop=archaicvar/*bcf
----------------------------------------
> state_names = ['Human', 'Archaic']
> starting_probabilities = [0.954, 0.046]
> transitions = [[0.999, 0.001], [0.023, 0.977]]
> emissions = [0.047, 0.386]
> number of windows: 2877010 . Number of snps =  129734
> total callability: 0.72
> average mutation rate per bin: 1.0
> Output prefix is /dev/stdout
> Window size is 1000 bp
> Haploid False
----------------------------------------
bcftools view -a -R obs.HG00096.txttemp archaicvar/highcov_ind_9.bcf
bcftools view -a -R obs.HG00096.txttemp archaicvar/highcov_ind_19.bcf
bcftools view -a -R obs.HG00096.txttemp archaicvar/highcov_ind_7.bcf
bcftools view -a -R obs.HG00096.txttemp archaicvar/highcov_ind_21.bcf
bcftools view -a -R obs.HG00096.txttemp archaicvar/highcov_ind_20.bcf
bcftools view -a -R obs.HG00096.txttemp archaicvar/highcov_ind_15.bcf
bcftools view -a -R obs.HG00096.txttemp archaicvar/highcov_ind_10.bcf
bcftools view -a -R obs.HG00096.txttemp archaicvar/highcov_ind_3.bcf
bcftools view -a -R obs.HG00096.txttemp archaicvar/highcov_ind_17.bcf
bcftools view -a -R obs.HG00096.txttemp archaicvar/highcov_ind_6.bcf
bcftools view -a -R obs.HG00096.txttemp archaicvar/highcov_ind_X.bcf
bcftools view -a -R obs.HG00096.txttemp archaicvar/highcov_ind_16.bcf
bcftools view -a -R obs.HG00096.txttemp archaicvar/highcov_ind_1.bcf
bcftools view -a -R obs.HG00096.txttemp archaicvar/highcov_ind_18.bcf
bcftools view -a -R obs.HG00096.txttemp archaicvar/highcov_ind_14.bcf
bcftools view -a -R obs.HG00096.txttemp archaicvar/highcov_ind_4.bcf
bcftools view -a -R obs.HG00096.txttemp archaicvar/highcov_ind_2.bcf
bcftools view -a -R obs.HG00096.txttemp archaicvar/highcov_ind_22.bcf
bcftools view -a -R obs.HG00096.txttemp archaicvar/highcov_ind_5.bcf
bcftools view -a -R obs.HG00096.txttemp archaicvar/highcov_ind_8.bcf
bcftools view -a -R obs.HG00096.txttemp archaicvar/highcov_ind_11.bcf
bcftools view -a -R obs.HG00096.txttemp archaicvar/highcov_ind_12.bcf
bcftools view -a -R obs.HG00096.txttemp archaicvar/highcov_ind_13.bcf

chrom  start     end       length  state    mean_prob  snps  admixpopvariants  AltaiNeandertal  Vindija33.19  Denisova  Chagyrskaya-Phalanx
1      2988000   2996000   9000    Archaic  0.71815    6     4                 4                4             1         4
1      3425000   3451000   27000   Archaic  0.95652    22    17                17               15            3         17
1      4302000   4360000   59000   Archaic  0.84636    20    12                11               12            11        11
1      4500000   4509000   10000   Archaic  0.84456    7     5                 4                5             4         5
1      5339000   5346000   8000    Archaic  0.58909    4     3                 2                3             0         3
1      9322000   9354000   33000   Archaic  0.8475     9     0                 0                0             0         0
1      12599000  12653000  55000   Archaic  0.9142     18    11                4                11            0         10

```

For the first segment there are 6 derived snps. Of these snps 4 are shared with Altai,Vindija, Denisova and Chagyrskaya. Only 1 is shared with Denisova so this segment likeli introgressed from Neanderthals

---

### Run in python

You also have the choice to run the functions from within python. Can be handy if you are simulating data and don't want to generate a ton of outfiles.

```py
from make_test_data import create_test_data
from hmm_functions import TrainModel, DecodeModel, HMMParam, read_HMM_parameters_from_file
from helper_functions import Load_observations_weights_mutrates

# -----------------------------------------------------------------------------
# Test data from quick tutorial
# -----------------------------------------------------------------------------

# Initial HMM guess
initial_hmm_params = HMMParam(state_names = ['Human', 'Archaic'], 
                              starting_probabilities = [0.5, 0.5], 
                              transitions = [[0.99,0.01],[0.02,0.98]], 
                              emissions = [0.03, 0.3]) 

# Create test data
obs, chroms, starts, variants, weights, mutrates  = create_test_data(50000, write_out_files = False)

# Train model
hmm_parameters = TrainModel(obs, mutrates, weights, initial_hmm_params)

# Decode model
segments = DecodeModel(obs, chroms, starts, variants, mutrates, weights, hmm_parameters)

for segment_info in segments:
    chrom, genome_start, genome_end, genome_length, state, mean_prob, snp_counter, ploidity, variants = segment_info
    print(chrom, genome_start,  genome_end, genome_length, state, mean_prob, snp_counter, sep = '\t')



# -----------------------------------------------------------------------------
# Running on an individual from 1000 genomes
# -----------------------------------------------------------------------------

hmm_parameters = read_HMM_parameters_from_file('trained.HG00096.json')
obs, chroms, starts, variants, mutrates, weights = Load_observations_weights_mutrates(obs_file = 'obs.HG00096.txt', 
                                                                                      weights_file = 'strickmask.bed', 
                                                                                      mutrates_file = 'mutationrate.bed', 
                                                                                      window_size = 1000, 
                                                                                      haploid = False)

segments = DecodeModel(obs, chroms, starts, variants, mutrates, weights, hmm_parameters)

for segment_info in segments:  
    chrom, genome_start, genome_end, genome_length, state, mean_prob, snp_counter, ploidity, variants = segment_info
    print(chrom, genome_start,  genome_end, genome_length, state, mean_prob, snp_counter, sep = '\t')


```

And that is it! Now you have run the model and gotten a set of parameters that you can interpret biologically (see my paper) and you have a list of segments that belong to the human and Archaic state.

If you have any questions about the use of the scripts, if you find errors or if you have feedback you can contact my here (make an issue) or write to:
lauritsskov2@gmail.com



