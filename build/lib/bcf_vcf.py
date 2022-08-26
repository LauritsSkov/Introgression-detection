import os
import numpy as np
from collections import defaultdict
import sys 
from glob import glob

from helper_functions import get_consensus, sortby, Make_folder_if_not_exists, load_fasta, convert_to_bases, handle_individuals_input, handle_infiles, clean_files

# Check which type of input we are dealing with
def combined_files(ancestralfiles, vcffiles):

    # Get ancestral and vcf consensus
    prefix1, postfix1, values1 = get_consensus(vcffiles)
    prefix2, postfix2, values2 = get_consensus(ancestralfiles)

    
    # No ancestral files
    if ancestralfiles == ['']:
        ancestralfiles = [None for _ in vcffiles]
        return ancestralfiles, vcffiles

    # Same length
    elif len(ancestralfiles) == len(vcffiles):
        return ancestralfiles, vcffiles

    # diff lengthts (both longer than 1)       
    elif len(ancestralfiles) > 1 and len(vcffiles) > 1:
        vcffiles = []
        ancestralfiles = []

        for joined in sorted(values1.intersection(values2), key=sortby):
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


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Dealing with bcf/vcf files functions
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
def make_out_group(individuals_input, bedfile, vcffiles, outputfile, ancestralfiles, refgenomefiles):

    # Get list of outgroup individuals
    outgroup_individuals = handle_individuals_input(individuals_input, 'outgroup')

    # Get a list of vcffiles and ancestral files and intersect them
    vcffiles = handle_infiles(vcffiles)
    ancestralfiles = handle_infiles(ancestralfiles)
    refgenomefiles = handle_infiles(refgenomefiles)

    ancestralfiles, vcffiles = combined_files(ancestralfiles, vcffiles)
    refgenomefiles, vcffiles = combined_files(refgenomefiles, vcffiles)

    print('-' * 40)
    print('> Outgroup individuals:', len(outgroup_individuals))
    print('> Using vcf and ancestral files')
    for vcffile, ancestralfile, reffile in zip(vcffiles, ancestralfiles, refgenomefiles):
        print('vcffile:',vcffile, 'ancestralfile:',ancestralfile, 'reffile:', reffile)
    print()    
    print('> Callability file:', bedfile)
    print(f'> Writing output to:',outputfile)
    print('-' * 40)


    Make_folder_if_not_exists(outputfile)
    outgroup_individuals = ','.join(outgroup_individuals)

    with open(outputfile + '.unsorted', 'w') as out:

        print('chrom', 'pos', 'ref_allele_info', 'alt_allele_info', 'ancestral_base', sep = '\t', file = out)

        for vcffile, ancestralfile, reffile in zip(vcffiles, ancestralfiles, refgenomefiles):

            if ancestralfile is not None:
                ancestral_allele = load_fasta(ancestralfile)

            if bedfile is not None:
                command = f'bcftools view -m2 -M2 -v snps -s {outgroup_individuals} -T {bedfile} {vcffile} | vcftools --vcf - --counts --stdout'
            else:
                command = f'bcftools view -m2 -M2 -v snps -s {outgroup_individuals} {vcffile} | vcftools --vcf - --counts --stdout'
            
            print(f'Processing {vcffile}...')
            print('Running command:')
            print(command, '\n\n')

            variants_seen = defaultdict(int)
            for index, line in enumerate(os.popen(command)):
                if not line.startswith('CHROM'):

                    chrom, pos, _, _, ref_allele_info, alt_allele_info = line.strip().split()

                    ref_allele, ref_count = ref_allele_info.split(':')
                    alt_allele, alt_count = alt_allele_info.split(':')
                    pos, ref_count, alt_count = int(pos),  int(ref_count), int(alt_count)

                    # Always include polymorphic sites
                    if alt_count * ref_count > 0:
                        ancestral_base = ref_allele if ref_count > alt_count else alt_allele

                        # Use ancestral base info if available
                        if ancestralfile is not None:
                            ancestral_base_temp = ancestral_allele[pos-1]
                            if ancestral_base_temp in [ref_allele, alt_allele]:
                                 ancestral_base = ancestral_base_temp

                        print(chrom, pos, ref_allele_info, alt_allele_info, ancestral_base, sep = '\t', file = out)
                        variants_seen[pos-1] = 1

                    # Fixed sites
                    elif alt_count * ref_count == 0:
                        ancestral_base = ref_allele if ref_count > alt_count else alt_allele

                        # Use ancestral base info if available
                        if ancestralfile is not None:
                            ancestral_base_temp = ancestral_allele[pos-1]
                            if ancestral_base_temp in [ref_allele, alt_allele]:
                                 ancestral_base = ancestral_base_temp

                        if ancestral_base == alt_allele:
                            derived_count = ref_count
                        else:
                             derived_count = alt_count

                        if derived_count > 0:
                            print(chrom, pos, ref_allele_info, alt_allele_info, ancestral_base, sep = '\t', file = out)
                            variants_seen[pos-1] = 1


                    if index % 100000 == 0:
                        print(f'at line {index} at chrom {chrom} and position {pos}')

            # If reference genome is provided then remove positions where the reference and ancestral differ AND which is not found in the outgroup
            if reffile is not None and ancestralfile is not None:
                print('Find fixed derived sites')
                refgenome_allele = load_fasta(reffile)

                for index, (refbase, ancbase) in enumerate(zip(refgenome_allele, ancestral_allele)):
                    if ancbase in 'ACGT' and refbase in 'ACGT':
                        if refbase != ancbase and variants_seen[index] == 0:
                            print(chrom, index + 1, f'{refbase}:100', f'{ancbase}:0', ancbase, sep = '\t', file = out)

    # Sort outgroup file
    print('Sorting outgroup file')
    positions_to_sort = defaultdict(lambda: defaultdict(str))
    with open(outputfile + '.unsorted') as data, open(outputfile, 'w') as out:
        for line in data:
            if line.startswith('chrom'):
                out.write(line)
            else:
                chrom, pos = line.strip().split()[0:2]
                positions_to_sort[chrom][int(pos)] = line

        for chrom in sorted(positions_to_sort, key=sortby):
            for pos in sorted(positions_to_sort[chrom]):
                line =  positions_to_sort[chrom][pos]
                out.write(line)

    # Clean log files generated by vcf and bcf tools
    clean_files(outputfile + '.unsorted')
    clean_files('out.log')
    


def make_ingroup_obs(individuals_input, bedfile, vcffiles, outprefix, outgroupfile, ancestralfiles):

    # Get a list of ingroup individuals
    ingroup_individuals = handle_individuals_input(individuals_input, 'ingroup')

    # Get a list of vcffiles and ancestral files and intersect them
    vcffiles = handle_infiles(vcffiles)
    ancestralfiles = handle_infiles(ancestralfiles)
    ancestralfiles, vcffiles  = combined_files(ancestralfiles, vcffiles)

    print('-' * 40)
    print('> Ingroup individuals:', len(ingroup_individuals))
    print('> Using vcf and ancestral files')
    for vcffile, ancestralfile in zip(vcffiles, ancestralfiles):
        print('vcffile:',vcffile, 'ancestralfile:',ancestralfile)
    print()  
    print('> Using outgroup variants from:', outgroupfile)  
    print('> Callability file:', bedfile)
    print(f'> Writing output to file with prefix: {outprefix}.<individual>.txt')
    print('-' * 40)


    # handle output files
    Make_folder_if_not_exists(outprefix)
    outfile_handler = defaultdict(str)
    for individual in ingroup_individuals:
        outfile_handler[individual] = open(f'{outprefix}.{individual}.txt','w')
        print('chrom', 'pos', 'ancestral_base', 'genotype', sep = '\t', file = outfile_handler[individual])
        
    individuals_for_bcf = ','.join(ingroup_individuals)

    for vcffile, ancestralfile in zip(vcffiles, ancestralfiles):

        if ancestralfile is not None:
            ancestral_allele = load_fasta(ancestralfile)

        if bedfile is not None:
            command = f'bcftools view -m2 -M2 -v snps -s {individuals_for_bcf} -T {bedfile} {vcffile} | vcftools --vcf - --exclude-positions {outgroupfile} --recode --stdout'
        else:
            command = f'bcftools view -m2 -M2 -v snps -s {individuals_for_bcf} {vcffile} | vcftools --vcf - --exclude-positions {outgroupfile} --recode --stdout'

        print('Running command:')
        print(command, '\n\n')

        for index, line in enumerate(os.popen(command)):

            if line.startswith('#CHROM'):
                individuals_in_vcffile = line.strip().split()[9:]

            if not line.startswith('#'):

                chrom, pos, _, ref_allele, alt_allele = line.strip().split()[0:5]
                pos = int(pos)
                genotypes = [x.split(':')[0] for x in line.strip().split()[9:]]

                if ref_allele in 'ACGT' and alt_allele in 'ACGT':

                    for original_genotype, individual in zip(genotypes, individuals_in_vcffile):
                        ref_count = original_genotype.count('0')
                        alt_count = original_genotype.count('1')     
                        genotype = convert_to_bases(original_genotype, ref_allele, alt_allele)   

                        if ancestralfile is not None:
                            # With ancestral information look for derived alleles
                            ancestral_base = ancestral_allele[pos-1]
                            if ancestral_base in [ref_allele, alt_allele]:

                                derived_count = genotype.count(alt_allele) if ancestral_base == ref_allele else genotype.count(ref_allele)
                                if derived_count > 0:
                                    print(chrom, pos, ancestral_base, genotype, sep = '\t', file = outfile_handler[individual])

                        else:
                            # If no ancestral information is provided only include heterozygous variants
                            if alt_count * ref_count > 0:
                                ancestral_base = ref_allele
                                print(chrom, pos, ancestral_base, genotype, sep = '\t', file = outfile_handler[individual])
                

                if index % 100000 == 0:
                    print(f'at line {index} at chrom {chrom} and position {pos}')

    # Clean log files generated by vcf and bcf tools
    clean_files('out.log')

    for individual in ingroup_individuals:
        outfile_handler[individual].close()



