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

And look at the number of lines:

```bash
wc -l weights.*
  3874613 weights.bed
  3036315 weights.txt
  6910928 total
```

### Which variants are found in the outgroup

Now we can download the 1000 genomes VCF files and remove all variants found in an outgroup (this case the YRI, ESN and MSL Subsaharan-Africans).

```
VCF files are in this directory
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

Metadata is here
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

To do this we run the script:




With this output file we can estimate the average mutation rate in a region compared to the 





