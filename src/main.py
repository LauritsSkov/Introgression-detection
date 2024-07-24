import argparse
import numpy as np

from hmm_functions import TrainModel, DecodeModel, write_HMM_to_file, read_HMM_parameters_from_file, Write_Decoded_output, inhomogeneous, DecodeModel_viterbi
from bcf_vcf import make_out_group, make_ingroup_obs
from make_test_data import create_test_data
from make_mutationrate import make_mutation_rate
from helper_functions import Load_observations_weights_mutrates, handle_individuals_input, handle_infiles, combined_files


VERSION = '0.7.3'


def print_script_usage():
    toprint = f'''
Script for identifying introgressed archaic segments (version: {VERSION})

> Turorial:
hmmix make_test_data 
hmmix train  -obs=obs.txt -weights=weights.bed -mutrates=mutrates.bed -param=Initialguesses.json -out=trained.json 
hmmix decode -obs=obs.txt -weights=weights.bed -mutrates=mutrates.bed -param=trained.json


Different modes (you can also see the options for each by writing hmmix make_test_data -h):
> make_test_data        
    -windows            Number of Kb windows to create (defaults to 50,000 per chromosome)
    -chromosomes        Number of chromosomes to simulate (defaults to 2)
    -nooutfiles         Don't create obs.txt, mutrates.bed, weights.bed, Initialguesses.json, simulated_segments.txt (defaults to yes)
    -param              markov parameters file (default is human/neanderthal like parameters)

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
    -chrom              Subset to chromosome or comma separated list of chromosomes e.g chr1 or chr1,chr2,chr3 (default is use all chromosomes)
    -weights            file with callability (defaults to all positions being called)
    -mutrates           file with mutation rates (default is mutation rate is uniform)
    -param              markov parameters file (default is human/neanderthal like parameters)
    -out                outputfile (default is a file named trained.json)
    -window_size        size of bins (default is 1000 bp)
    -haploid            Change from using diploid data to haploid data (default is diploid)

> decode                
    -obs                [required] file with observation data
    -chrom              Subset to chromosome or comma separated list of chromosomes e.g chr1 or chr1,chr2,chr3 (default is use all chromosomes)
    -weights            file with callability (defaults to all positions being called)
    -mutrates           file with mutation rates (default is mutation rate is uniform)
    -param              markov parameters file (default is human/neanderthal like parameters)
    -out                outputfile prefix <out>.hap1.txt and <out>.hap2.txt if -haploid option is used or <out>.diploid.txt (default is stdout)
    -window_size        size of bins (default is 1000 bp)
    -haploid            Change from using diploid data to haploid data (default is diploid)
    -admixpop           Annotate using vcffile with admixing population (default is none)
    -extrainfo          Add variant position for each SNP (default is off)
    -viterbi            decode using the viterbi algorithm (default is posterior decoding)

> inhomogeneous                
    -obs                [required] file with observation data
    -chrom              Subset to chromosome or comma separated list of chromosomes e.g chr1 or chr1,chr2,chr3 (default is use all chromosomes)
    -weights            file with callability (defaults to all positions being called)
    -mutrates           file with mutation rates (default is mutation rate is uniform)
    -param              markov parameters file (default is human/neanderthal like parameters)
    -out                outputfile prefix <out>.hap1_sim(0-n).txt and <out>.hap2_sim(0-n).txt if -haploid option is used or <out>.diploid_(0-n).txt (default is stdout)
    -window_size        size of bins (default is 1000 bp)
    -haploid            Change from using diploid data to haploid data (default is diploid)
    -samples            Number of simulated paths for the inhomogeneous markov chain (default is 100)
    -admixpop           Annotate using vcffile with admixing population (default is none)
    -extrainfo          Add variant position for each SNP (default is off)
    '''

    return toprint

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
def main():

    parser = argparse.ArgumentParser(description=print_script_usage(), formatter_class=argparse.RawTextHelpFormatter)

    subparser = parser.add_subparsers(dest = 'mode')

    # Run test
    test_subparser = subparser.add_parser('make_test_data', help='Create test data')
    test_subparser.add_argument("-windows", metavar='',help="Number of Kb windows to create (defaults to 50,000 per chromosome)", type=int, default = 50000)
    test_subparser.add_argument("-chromosomes", metavar='',help="Number of chromosomes to simulate (defaults to 2)", type=int, default = 2)
    test_subparser.add_argument("-nooutfiles",help="Don't create obs.txt, mutrates.bed, weights.bed, Initialguesses.json (defaults to yes)", action='store_false', default = True)
    test_subparser.add_argument("-param", metavar='',help="markov parameters file (default is human/neanderthal like parameters)", type=str)

    # Make outgroup
    outgroup_subparser = subparser.add_parser('create_outgroup', help='Create outgroup information')
    outgroup_subparser.add_argument("-ind",help="[required] ingroup/outgrop list (json file) or comma-separated list e.g. ind1,ind2", type=str, required = True)
    outgroup_subparser.add_argument("-vcf",help="[required] path to list of comma-separated vcf/bcf file(s) or wildcard characters e.g. chr*.bcf", type=str, required = True)
    outgroup_subparser.add_argument("-weights", metavar='',help="file with callability (defaults to all positions being called)")
    outgroup_subparser.add_argument("-out", metavar='',help="outputfile (defaults to stdout)", default = '/dev/stdout')
    outgroup_subparser.add_argument("-ancestral", metavar='',help="fasta file with ancestral information - comma-separated list or wildcards like vcf argument (default none)", default='')
    outgroup_subparser.add_argument("-refgenome", metavar='',help="fasta file with reference genome - comma-separated list or wildcards like vcf argument (default none)", default='')

    # Make mutation rates
    mutation_rate = subparser.add_parser('mutation_rate', help='Estimate mutation rate')
    mutation_rate.add_argument("-outgroup", help="[required] path to variants found in outgroup", type=str, required = True)
    mutation_rate.add_argument("-out", metavar='',help="outputfile (defaults to mutationrate.bed)", default = 'mutationrate.bed')
    mutation_rate.add_argument("-weights", metavar='',help="file with callability (defaults to all positions being called)")
    mutation_rate.add_argument("-window_size", metavar='',help="size of bins (defaults to 1 Mb)", type=int, default = 1000000)

    # Make ingroup observations
    create_obs_subparser = subparser.add_parser('create_ingroup', help='Create ingroup data')
    create_obs_subparser.add_argument("-ind", help="[required] ingroup/outgrop list (json file) or comma-separated list e.g. ind1,ind2", type=str, required = True)
    create_obs_subparser.add_argument("-vcf", help="[required] path to list of comma-separated vcf/bcf file(s) or wildcard characters e.g. chr*.bcf", type=str, required = True)
    create_obs_subparser.add_argument("-outgroup", help="[required] path to variant found in outgroup", type=str, required = True)
    create_obs_subparser.add_argument("-weights", metavar='',help="file with callability (defaults to all positions being called)")
    create_obs_subparser.add_argument("-out", metavar='',help="outputfile prefix (default is a file named obs.<ind>.txt where ind is the name of individual in ingroup/outgrop list)", default = 'obs')
    create_obs_subparser.add_argument("-ancestral", metavar='',help="fasta file with ancestral information - comma-separated list or wildcards like vcf argument (default none)", default='')

    # Train model
    train_subparser = subparser.add_parser('train', help='Train HMM')
    train_subparser.add_argument("-obs",help="[required] file with observation data", type=str, required = True)
    train_subparser.add_argument("-chrom",help="Subset to chromosome or comma separated list of chromosomes e.g chr1 or chr1,chr2,chr3", type=str, default='All')
    train_subparser.add_argument("-weights", metavar='',help="file with callability (defaults to all positions being called)")
    train_subparser.add_argument("-mutrates", metavar='',help="file with mutation rates (default is mutation rate is uniform)")
    train_subparser.add_argument("-param", metavar='',help="markov parameters file (default is human/neanderthal like parameters)", type=str)
    train_subparser.add_argument("-out", metavar='',help="outputfile (default is a file named trained.json)", default = 'trained.json')
    train_subparser.add_argument("-window_size", metavar='',help="size of bins (default is 1000 bp)", type=int, default = 1000)
    train_subparser.add_argument("-haploid",help="Change from using diploid data to haploid data (default is diploid)", action='store_true', default = False)

    # Decode model
    decode_subparser = subparser.add_parser('decode', help='Decode HMM')
    decode_subparser.add_argument("-obs",help="[required] file with observation data", type=str, required = True)
    decode_subparser.add_argument("-chrom",help="Subset to chromosome or comma separated list of chromosomes e.g chr1 or chr1,chr2,chr3", type=str, default='All')
    decode_subparser.add_argument("-weights", metavar='',help="file with callability (defaults to all positions being called)")
    decode_subparser.add_argument("-mutrates", metavar='',help="file with mutation rates (default is mutation rate is uniform)")
    decode_subparser.add_argument("-param", metavar='',help="markov parameters file (default is human/neanderthal like parameters)", type=str)
    decode_subparser.add_argument("-out", metavar='',help="outputfile prefix <out>.hap1.txt and <out>.hap2.txt if -haploid option is used or <out>.diploid.txt (default is stdout)", default = '/dev/stdout')
    decode_subparser.add_argument("-window_size", metavar='',help="size of bins (default is 1000 bp)", type=int, default = 1000)
    decode_subparser.add_argument("-haploid",help="Change from using diploid data to haploid data (default is diploid)", action='store_true', default = False)
    decode_subparser.add_argument("-admixpop",help="Annotate using vcffile with admixing population (default is none)")
    decode_subparser.add_argument("-extrainfo",help="Add archaic information on each SNP", action='store_true', default = False)
    decode_subparser.add_argument("-viterbi",help="Decode using the Viterbi algorithm", action='store_true', default = False)

    # inhomogeneous markov chain
    inhomogen_subparser = subparser.add_parser('inhomogeneous', help='Make inhomogen markov chain')
    inhomogen_subparser.add_argument("-obs",help="[required] file with observation data", type=str, required = True)
    inhomogen_subparser.add_argument("-chrom",help="Subset to chromosome or comma separated list of chromosomes e.g chr1 or chr1,chr2,chr3", type=str, default='All')
    inhomogen_subparser.add_argument("-weights", metavar='',help="file with callability (defaults to all positions being called)")
    inhomogen_subparser.add_argument("-mutrates", metavar='',help="file with mutation rates (default is mutation rate is uniform)")
    inhomogen_subparser.add_argument("-param", metavar='',help="markov parameters file (default is human/neanderthal like parameters)", type=str)
    inhomogen_subparser.add_argument("-out", metavar='',help="outputfile prefix <out>.hap1.txt and <out>.hap2.txt if -haploid option is used or <out>.diploid.txt (default is stdout)", default = '/dev/stdout')
    inhomogen_subparser.add_argument("-window_size", metavar='',help="size of bins (default is 1000 bp)", type=int, default = 1000)
    inhomogen_subparser.add_argument("-haploid",help="Change from using diploid data to haploid data (default is diploid)", action='store_true', default = False)
    inhomogen_subparser.add_argument("-samples",help="Number of paths to sample (default is 100)", type=int, default = 100)
    inhomogen_subparser.add_argument("-admixpop",help="Annotate using vcffile with admixing population (default is none)")
    inhomogen_subparser.add_argument("-extrainfo",help="Add archaic information on each SNP", action='store_true', default = False)






    args = parser.parse_args()

    # Make test data
    # ------------------------------------------------------------------------------------------------------------
    if args.mode == 'make_test_data':
        create_test_data(data_set_length = args.windows, n_chromosomes = args.chromosomes, write_out_files = args.nooutfiles, parameters_file = args.param)
        

    # Train parameters
    # ------------------------------------------------------------------------------------------------------------
    elif args.mode == 'train':

        hmm_parameters = read_HMM_parameters_from_file(args.param)
        obs, _, _, _, mutrates, weights = Load_observations_weights_mutrates(args.obs, args.weights, args.mutrates, args.window_size, args.haploid, args.chrom)
        
        print('-' * 40)
        print(hmm_parameters)
        print(f'> chromosomes to use: {args.chrom}')
        print(f'> number of windows: {len(obs)}. Number of snps = {sum(obs)}')
        print(f'> total callability: {np.sum(weights)} bp ({round(np.sum(weights) / len(obs) * 100,2)} %)')
        print('> average mutation rate per bin:', round(np.sum(mutrates * weights) / np.sum(weights), 2) )
        print('> Output is',args.out) 
        print('> Window size is',args.window_size, 'bp') 
        print('> Haploid',args.haploid) 
        print('-' * 40)

        hmm_parameters = TrainModel(obs, mutrates, weights, hmm_parameters)
        write_HMM_to_file(hmm_parameters, args.out)




    # Decode observations using parameters
    # ------------------------------------------------------------------------------------------------------------
    elif args.mode == 'decode':

        obs, chroms, starts, variants, mutrates, weights  = Load_observations_weights_mutrates(args.obs, args.weights, args.mutrates, args.window_size, args.haploid, args.chrom)
        hmm_parameters = read_HMM_parameters_from_file(args.param)

        print('-' * 40)
        print(hmm_parameters)  
        print(f'> chromosomes to use: {args.chrom}')
        print(f'> number of windows: {len(obs)}. Number of snps = {sum(obs)}')
        print(f'> total callability: {np.sum(weights)} bp ({round(np.sum(weights) / len(obs) * 100,2)} %)')
        print('> average mutation rate per bin:', round(np.sum(mutrates * weights) / np.sum(weights), 2) )
        print('> Output prefix is',args.out) 
        print('> Window size is',args.window_size, 'bp') 
        print('> Haploid',args.haploid)      
        if args.viterbi:
            print('> Decode using viterbi algorithm') 
        else:
            print('> Decode with posterior decoding') 
        print('-' * 40)

        # Find segments and write output
        if args.viterbi:
            segments = DecodeModel_viterbi(obs, chroms, starts, variants, mutrates, weights, hmm_parameters)
        else:
            segments = DecodeModel(obs, chroms, starts, variants, mutrates, weights, hmm_parameters)
        Write_Decoded_output(args.out, segments, args.obs, args.admixpop, args.extrainfo)

    # inhomogeneous markov chain
    # ------------------------------------------------------------------------------------------------------------
    elif args.mode == 'inhomogeneous':

        obs, chroms, starts, variants, mutrates, weights  = Load_observations_weights_mutrates(args.obs, args.weights, args.mutrates, args.window_size, args.haploid, args.chrom)
        hmm_parameters = read_HMM_parameters_from_file(args.param)
        
        print('-' * 40)
        print(hmm_parameters)  
        print(f'> chromosomes to use: {args.chrom}')
        print(f'> number of windows: {len(obs)}. Number of snps = {sum(obs)}')
        print(f'> total callability: {np.sum(weights)} bp ({round(np.sum(weights) / len(obs) * 100,2)} %)')
        print('> average mutation rate per bin:', round(np.sum(mutrates * weights) / np.sum(weights), 2) )
        print('> Output prefix is',args.out) 
        print('> Window size is',args.window_size, 'bp') 
        print('> Haploid',args.haploid) 
        print('-' * 40)

        # Find segments and write output
        segments = inhomogeneous(hmm_parameters, weights, obs, mutrates, args.samples, chroms, starts, variants)
        Write_Decoded_output(args.out, segments, args.obs, args.admixpop, args.extrainfo)

            



    # Create outgroup snps (set of snps to be removed)
    # ------------------------------------------------------------------------------------------------------------
    elif args.mode == 'create_outgroup':

        # Get list of outgroup individuals
        outgroup_individuals = handle_individuals_input(args.ind, 'outgroup')

        # Get a list of vcffiles and ancestral files and intersect them
        vcffiles = handle_infiles(args.vcf)
        ancestralfiles = handle_infiles( args.ancestral)
        refgenomefiles = handle_infiles(args.refgenome)

        ancestralfiles, vcffiles = combined_files(ancestralfiles, vcffiles)
        refgenomefiles, vcffiles = combined_files(refgenomefiles, vcffiles)

        print('-' * 40)
        print('> Outgroup individuals:', len(outgroup_individuals))
        print('> Using vcf and ancestral files')
        for vcffile, ancestralfile, reffile in zip(vcffiles, ancestralfiles, refgenomefiles):
            print('vcffile:',vcffile, 'ancestralfile:',ancestralfile, 'reffile:', reffile)
        print()    
        print('> Callability file:',  args.weights)
        print(f'> Writing output to:', args.out)
        print('-' * 40)

        
        make_out_group(outgroup_individuals, args.weights, vcffiles, args.out, ancestralfiles, refgenomefiles)




    # Create ingroup observations
    # ------------------------------------------------------------------------------------------------------------
    elif args.mode == 'create_ingroup':

        # Get a list of ingroup individuals
        ingroup_individuals = handle_individuals_input(args.ind,'ingroup')

        # Get a list of vcffiles and ancestral files and intersect them
        vcffiles = handle_infiles(args.vcf)
        ancestralfiles = handle_infiles(args.ancestral)

        ancestralfiles, vcffiles  = combined_files(ancestralfiles, vcffiles)

        print('-' * 40)
        print('> Ingroup individuals:', len(ingroup_individuals))
        print('> Using vcf and ancestral files')
        for vcffile, ancestralfile in zip(vcffiles, ancestralfiles):
            print('vcffile:',vcffile, 'ancestralfile:',ancestralfile)
        print()  
        print('> Using outgroup variants from:', args.outgroup)  
        print('> Callability file:', args.weights)
        print(f'> Writing output to file with prefix: {args.out}.<individual>.txt')
        print('-' * 40)


        make_ingroup_obs(ingroup_individuals, args.weights, vcffiles, args.out, args.outgroup, ancestralfiles)




    # Estimate mutation rate
    # ------------------------------------------------------------------------------------------------------------
    elif args.mode == 'mutation_rate':
        print('-' * 40)
        print(f'> Outgroupfile:', args.outgroup)
        print(f'> Outputfile is:', args.out)
        print(f'> Callability file is:', args.weights)
        print(f'> Window size:', args.window_size)
        print('-' * 40)

        make_mutation_rate(args.outgroup, args.out, args.weights, args.window_size)
    
    



    # Print usage
    # ------------------------------------------------------------------------------------------------------------
    else:
        print(print_script_usage())


if __name__ == "__main__":
    main()

