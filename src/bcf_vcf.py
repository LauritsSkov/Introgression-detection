import os
from collections import defaultdict

from helper_functions import sortby, Make_folder_if_not_exists, load_fasta, convert_to_bases, clean_files



# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Make Outgroup
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
def make_out_group(individuals_input, bedfile, vcffiles, outputfile, ancestralfiles, refgenomefiles):

    Make_folder_if_not_exists(outputfile)
    outgroup_individuals = ','.join(individuals_input)

    with open(outputfile + '.unsorted', 'w') as out:

        print('chrom', 'pos', 'ref_allele_info', 'alt_allele_info', 'ancestral_base', sep = '\t', file = out)

        for vcffile, ancestralfile, reffile in zip(vcffiles, ancestralfiles, refgenomefiles):

            if ancestralfile is not None:
                ancestral_allele = load_fasta(ancestralfile)

            if bedfile is not None:
                command = f'bcftools view -s {outgroup_individuals} -T {bedfile} {vcffile} | bcftools norm -m -any | bcftools view -v snps | vcftools --vcf - --counts --stdout'
            else:
                command = f'bcftools view -s {outgroup_individuals} {vcffile} | bcftools norm -m -any | bcftools view -v snps | vcftools --vcf - --counts --stdout'
            
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
                    if ancbase in ['A','C','G','T']  and refbase in ['A','C','G','T'] :
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
    

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Make ingroup
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
def make_ingroup_obs(ingroup_individuals, bedfile, vcffiles, outprefix, outgroupfile, ancestralfiles):

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
            command = f'bcftools view -v snps -s {individuals_for_bcf} -T {bedfile} {vcffile} | bcftools norm -m +any | vcftools --vcf - --exclude-positions {outgroupfile} --recode --stdout'
        else:
            command = f'bcftools view -v snps -s {individuals_for_bcf} {vcffile} | bcftools norm -m +any | vcftools --vcf - --exclude-positions {outgroupfile} --recode --stdout'

        print('Running command:')
        print(command, '\n\n')

        for index, line in enumerate(os.popen(command)):

            if line.startswith('#CHROM'):
                individuals_in_vcffile = line.strip().split()[9:]

            if not line.startswith('#'):

                chrom, pos, _, ref_allele, alt_allele = line.strip().split()[0:5]
                pos = int(pos)
                genotypes = [x.split(':')[0] for x in line.strip().split()[9:]]
                all_bases = [ref_allele] + alt_allele.split(',')

                if ref_allele in ['A','C','G','T']:

                    for original_genotype, individual in zip(genotypes, individuals_in_vcffile):  
                        genotype = convert_to_bases(original_genotype, all_bases)   

                        if ancestralfile is not None:
                            # With ancestral information look for derived alleles
                            ancestral_base = ancestral_allele[pos-1]
                            if ancestral_base in all_bases and genotype.count(ancestral_base) != 2 and genotype != 'NN':
                                print(chrom, pos, ancestral_base, genotype, sep = '\t', file = outfile_handler[individual])

                        else:
                            # If no ancestral information is provided only include heterozygous variants
                            if genotype[0] != genotype[1]:
                                print(chrom, pos, ref_allele, genotype, sep = '\t', file = outfile_handler[individual])
                

                if index % 100000 == 0:
                    print(f'at line {index} at chrom {chrom} and position {pos}')

    # Clean log files generated by vcf and bcf tools
    clean_files('out.log')

    for individual in ingroup_individuals:
        outfile_handler[individual].close()



