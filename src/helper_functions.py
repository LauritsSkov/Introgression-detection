import numpy as np
import json
from collections import defaultdict
import os, sys
from glob import glob

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Functions for handling observertions/bed files
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
def make_callability_from_bed(bedfile, window_size):
    callability = defaultdict(lambda: defaultdict(float))
    with open(bedfile) as data:
        for line in data:

            if line.startswith('chrom'):
                continue

            if len(line.strip().split('\t')) == 3:
                chrom, start, end = line.strip().split('\t')
                value = 1
            elif  len(line.strip().split('\t')) > 3:
                chrom, start, end, value = line.strip().split('\t')[0:4]
                value = float(value)

            start, end = int(start), int(end)

            firstwindow = start - start % window_size
            lastwindow = end - end % window_size
            
            # not spanning multiple windows (all is added to same window)
            if firstwindow == lastwindow:
                callability[chrom][firstwindow] += (end-start+1) * value

            # spanning multiple windows
            else:
                # add to end windows
                firstwindow_fill = window_size - start % window_size
                lastwindow_fill = end %window_size
            
                callability[chrom][firstwindow] += firstwindow_fill * value
                callability[chrom][lastwindow] += (lastwindow_fill+1) * value

                # fill in windows in the middle
                for window_tofil in range(firstwindow + window_size, lastwindow, window_size):
                    callability[chrom][window_tofil] += window_size * value

    return callability



def Load_observations_weights_mutrates(obs_file, weights_file, mutrates_file, window_size = 1000, haploid = False, chrom_to_look_for = 'All'):

    # get span of data
    chromosome_spans = defaultdict(int)
    with open(obs_file) as data:
        for line in data:
            if line.startswith('chrom'):
                continue
            chrom, pos, ancestral_base, genotype = line.strip().split()

            # convert 1-indexed position to 0-indexed position
            zero_based_pos = int(pos) - 1
            rounded_pos = zero_based_pos - zero_based_pos % window_size
            chromosome_spans[chrom] = rounded_pos


    # get span of callability
    if weights_file: 
        callability = make_callability_from_bed(weights_file, window_size)
        for chrom in callability:
            chromosome_spans[chrom] = max(callability[chrom])

    # which chromosomes are we interested in
    if chrom_to_look_for != 'All':
        chromosome_list = sorted(chrom_to_look_for.split(','), key=sortby)
    else:
        chromosome_list = sorted(list(chromosome_spans.keys()), key=sortby)

    obs_counter = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    haplotypes = defaultdict(int)
    chroms, starts, variants, obs = [], [], [], []

    # read observation data
    with open(obs_file) as data:
        for line in data:
            if line.startswith('chrom'):
                continue 

            chrom, pos, ancestral_base, genotype = line.strip().split()
            if not chrom in chromosome_list:
                continue

            # convert 1-indexed position to 0-indexed position
            zero_based_pos = int(pos) - 1
            rounded_pos = zero_based_pos - zero_based_pos % window_size

            if haploid:  
                for i, base in enumerate(genotype):
                    if base != ancestral_base:
                        obs_counter[chrom][rounded_pos][f'_hap{i+1}'].append(pos)
                        haplotypes[f'_hap{i+1}'] += 1
            else:
                obs_counter[chrom][rounded_pos][''].append(pos)
                haplotypes[''] += 1


    # take care of cases where there is 0 observations
    if len(haplotypes) == 0:
        if haploid:
            sys.exit('Could not determine haploidity because there is no data')
        else:
            haplotypes[''] += 1
        

    for haplotype in sorted(haplotypes, key=sortby_haplotype):
        
        for chrom in chromosome_list:
            lastwindow = chromosome_spans[chrom]

            for window in range(0, lastwindow, window_size):
                chroms.append(f'{chrom}{haplotype}')   
                starts.append(window)
                variants.append(','.join(obs_counter[chrom][window][haplotype]))  
                obs.append(len(obs_counter[chrom][window][haplotype])) 
            

    # Read weights file is it exists - else set all weights to 1
    if weights_file is None:
        weights = np.ones(len(obs)) 
    else:  
        callability = make_callability_from_bed(weights_file, window_size)
        weights = []
        for haplotype in sorted(haplotypes, key=sortby_haplotype):
            
            for chrom in chromosome_list:
                lastwindow = chromosome_spans[chrom]

                for window in range(0, lastwindow, window_size):
                    weights.append(callability[chrom][window] / float(window_size))


    # Read mutation rate file is it exists - else set all mutation rates to 1
    if mutrates_file is None:
        mutrates = np.ones(len(obs)) 
    else:  
        callability = make_callability_from_bed(mutrates_file, window_size)
        mutrates = []
        for haplotype in sorted(haplotypes, key=sortby_haplotype):
            
            for chrom in chromosome_list:
                lastwindow = chromosome_spans[chrom]

                for window in range(0, lastwindow, window_size):
                    mutrates.append(callability[chrom][window] / float(window_size))



    # Make sure there are no places with obs > 0 and 0 in mutation rate or weight
    for index, (observation, w, m) in enumerate(zip(obs, weights, mutrates)):
        if w*m == 0 and observation != 0:
            print(f'warning, you had {observation} observations but no called bases/no mutation rate at index:{index}. weights:{w}, mutrates:{m}')
            obs[index] = 0
            
    return np.array(obs).astype(int), chroms, starts, variants, np.array(mutrates).astype(float), np.array(weights).astype(float)


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# For decoding/training
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------

def find_runs(inarray):
    """ run length encoding. Partial credit to R rle function. 
        Multi datatype arrays catered for including non Numpy
        returns: tuple (runlengths, startpositions, values) """
    ia = np.asarray(inarray)                # force numpy
    n = len(ia)
    if n == 0: 
        return (None, None, None)
    else:
        y = ia[1:] != ia[:-1]               # pairwise unequal (string safe)
        i = np.append(np.where(y), n - 1)   # must include last element posi
        z = np.diff(np.append(-1, i))       # run lengths
        p = np.cumsum(np.append(0, z))[:-1] # positions

        for (a, b, c) in zip(ia[i], p, z):
            yield (a, b, c)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Various helper functions
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
def load_fasta(fasta_file):
    '''
    Read a fasta file with a single chromosome in and return the sequence as a string
    '''
    fasta_sequence = ''
    with open(fasta_file) as data:
        for line in data:
            if not line.startswith('>'):
                fasta_sequence += line.strip().upper()

    return fasta_sequence

def sortby(x):
    '''
    This function is used in the sorted() function. It will sort first by numeric values, then strings then other symbols

    Usage:
    mylist = ['1', '12', '2', 3, 'MT', 'Y']
    sortedlist = sorted(mylist, key=sortby)
    returns ['1', '2', 3, '12', 'MT', 'Y']
    '''

    lower_case_letters = 'abcdefghijklmnopqrstuvwxyz'
    if x.isnumeric():
        return int(x)
    elif type(x) == str and len(x) > 0:
        if x[0].lower() in lower_case_letters:
            return 1e6 + lower_case_letters.index(x[0].lower())
        else:
            return 2e6
    else:
        return 3e6
    
def sortby_haplotype(x):
    '''
    This function will sort haplotypes by number
    '''

    if '_hap' in x:
        return int(x.replace('_hap', ''))
    else:
        return x
    



def Make_folder_if_not_exists(path):
    '''
    Check if path exists - otherwise make it;
    '''
    path = os.path.dirname(path)
    if path != '':
        if not os.path.exists(path):
            os.makedirs(path)



def Annotate_with_ref_genome(vcffiles, obsfile):
    obs = defaultdict(list)
    shared_with = defaultdict(str)

    tempobsfile = obsfile + 'temp'
    with open(obsfile) as data, open(tempobsfile,'w') as out:
        for line in data:
            if not line.startswith('chrom'):
                out.write(line)
                chrom, pos, ancestral_base, genotype = line.strip().split()
                derived_variant = genotype.replace(ancestral_base, '')[0]
                ID = f'{chrom}_{pos}'
                obs[ID] = [ancestral_base, derived_variant]

    # handle case with 0 SNPs
    if len(obs) == 0:
        for vcffile in handle_infiles(vcffiles):
            command = f'bcftools view -h {vcffile}'
            print('You have no observations!')

            for line in os.popen(command):
                if line.startswith('#CHROM'):
                    individuals_in_vcffile = line.strip().split()[9:]

            # Clean log files generated by vcf and bcf tools
            clean_files('out.log')
            clean_files(tempobsfile)

            return shared_with, individuals_in_vcffile


    print('Loading in admixpop snp information')
    for vcffile in handle_infiles(vcffiles):
        command = f'bcftools view -a -R {tempobsfile} {vcffile}'
        print(command)

        for line in os.popen(command):
            if line.startswith('#CHROM'):
                individuals_in_vcffile = line.strip().split()[9:]

            if not line.startswith('#'):

                chrom, pos, _, ref_allele, alt_allele = line.strip().split()[0:5]
                ID =  f'{chrom}_{pos}'
                genotypes = [x.split(':')[0] for x in line.strip().split()[9:]]
                all_bases = [ref_allele] + alt_allele.split(',')

                ancestral_base, derived_base = obs[ID]
                found_in = []

                for original_genotype, individual in zip(genotypes, individuals_in_vcffile):

                    if '.' not in original_genotype:
                        genotype = convert_to_bases(original_genotype, all_bases)   

                        if genotype.count(derived_base) > 0:
                            found_in.append(individual)

                if len(found_in) > 0:
                    shared_with[ID] = '|'.join(found_in)


    # Clean log files generated by vcf and bcf tools
    clean_files('out.log')
    clean_files(tempobsfile)

    return shared_with, individuals_in_vcffile

def handle_individuals_input(argument, group_to_choose):
    if os.path.exists(argument):
        with open(argument) as json_file:
            data = json.load(json_file)
            return data[group_to_choose]
    else:
        return argument.split(',')




# Clean up
def clean_files(filename):
    if os.path.exists(filename):
        os.remove(filename)


# Find variants from admixed population
def flatten_list(variants_list):
    return ','.join([x for x in variants_list if x != ''])



def convert_to_bases(genotype, both_bases):

    return_genotype = 'NN'
    separator = None
    
    if '/' in genotype or '|' in genotype:
        separator = '|' if '|' in genotype else '/'

        base1, base2 = [x for x in genotype.split(separator)]
        if base1.isnumeric() and base2.isnumeric():
            base1, base2 = int(base1), int(base2)

            if both_bases[base1] in ['A','C','G','T'] and both_bases[base2] in ['A','C','G','T']:
                return_genotype = both_bases[base1] + both_bases[base2]

    return return_genotype





def get_consensus(infiles):
    '''
    Find consensus prefix, suffix and value that changes in set of files:

    myfiles = ['chr1.vcf', 'chr2.vcf', 'chr3.vcf']
    prefix, suffix, values = get_consensus(myfiles) 
    
    prefix  = chr
    suffix = .vcf
    values  = [1,2,3] 
    '''
    infiles = [str(x) for x in infiles]

    if len(infiles) <= 1:
        return None, None, None
    
    # find longest common prefix
    prefix = infiles[0]
    for s in infiles[1:]:
        while not s.startswith(prefix):
            prefix = prefix[:-1]

    # find longest common suffix
    suffix = infiles[0]
    for s in infiles[1:]:
        while not s.endswith(suffix):
            suffix = suffix[1:]

    values = [x[len(prefix):-len(suffix)] if suffix else x[len(prefix):] for x in infiles]

    return prefix, suffix, sorted(values, key=sortby)


# Check which type of input we are dealing with
def handle_infiles(input):

    file_list = glob(input)
    if len(file_list) > 0:
        return file_list
    else:
        return input.split(',')
    


    

# Check which type of input we are dealing with
def combined_files(ancestralfiles, vcffiles):

    if len(ancestralfiles) == len(vcffiles) and len(vcffiles) == 1 and ancestralfiles != ['']:
        return ancestralfiles, vcffiles

    # Get ancestral and vcf consensus
    prefix1, postfix1, values1 = get_consensus(vcffiles)
    prefix2, postfix2, values2 = get_consensus(ancestralfiles)


    # No ancestral files
    if ancestralfiles == ['']:

        ancestralfiles = [None for _ in vcffiles]
        vcffiles = [f'{prefix1}{x}{postfix1}' for x in values1]
        return ancestralfiles, vcffiles

    # Same length
    elif len(ancestralfiles) == len(vcffiles):

        vcffiles = [f'{prefix1}{x}{postfix1}' for x in values1]
        ancestralfiles = [f'{prefix2}{x}{postfix2}' for x in values2]
        return ancestralfiles, vcffiles

    # diff lengthts (both longer than 1)       
    elif len(ancestralfiles) > 1 and len(vcffiles) > 1:

        vcffiles = []
        ancestralfiles = []

        for joined in sorted(set(values1).intersection(set(values2)), key=sortby):
            vcffiles.append(''.join([prefix1, joined, postfix1]))
            ancestralfiles.append(''.join([prefix2, joined, postfix2]))
        return ancestralfiles, vcffiles

    # Many ancestral files only one vcf    
    elif len(ancestralfiles) > 1 and len(vcffiles) == 1:
        ancestralfiles = []
        
        for key in values2:
            if key in vcffiles[0]:
                ancestralfiles.append(''.join([prefix2, key, postfix2]))

        if len(vcffiles) != len(ancestralfiles):
            sys.exit('Could not resolve ancestral files and vcffiles (try comma separated values)')

        return ancestralfiles, vcffiles

    # only one ancestral file and many vcf files
    elif len(ancestralfiles) == 1 and len(vcffiles) > 1:
        vcffiles = []
       
        for key in values1:
            if key in ancestralfiles[0]:
                vcffiles.append(''.join([prefix1, key, postfix1]))

        if len(vcffiles) != len(ancestralfiles):
            sys.exit('Could not resolve ancestral files and vcffiles (try comma separated values)')


        return ancestralfiles, vcffiles
    else:
        sys.exit('Could not resolve ancestral files and vcffiles (try comma separated values)')
