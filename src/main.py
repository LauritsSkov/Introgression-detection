import argparse
import numpy as np

from hmm_functions import TrainModel, DecodeModel, write_HMM_to_file, read_HMM_parameters_from_file, Write_Decoded_output
from bcf_vcf import make_out_group, make_ingroup_obs
from make_test_data import create_test_data
from make_mutationrate import make_mutation_rate
from helper_functions import Load_observations_weights_mutrates


def print_script_usage():
    toprint = '''
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
    -out                outputfile prefix (default is a file named trained.json)
    -window_size        size of bins (default is 1000 bp)
    -haploid            Change from using diploid data to haploid data (default is diploid)

> decode                
    -obs                [required] file with observation data
    -weights            file with callability (defaults to all positions being called)
    -mutrates           file with mutation rates (default is mutation rate is uniform)
    -param              markov parameters file (default is human/neanderthal like parameters)
    -out                outputfile prefix (default is stdout)
    -window_size        size of bins (default is 1000 bp)
    -haploid            Change from using diploid data to haploid data (default is diploid)
    -admixpop ADMIXPOP  Annotate using vcffile with admixing population (default is none)
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
    test_subparser.add_argument("-windows", metavar='',help="Number of Kb windows to create (defaults to 50,000)", type=int, default = 50000)
    test_subparser.add_argument("-nooutfiles",help="Don't create obs.txt, mutrates.bed, weights.bed, Initialguesses.json (defaults to yes)", action='store_false', default = True)

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
    train_subparser.add_argument("-weights", metavar='',help="file with callability (defaults to all positions being called)")
    train_subparser.add_argument("-mutrates", metavar='',help="file with mutation rates (default is mutation rate is uniform)")
    train_subparser.add_argument("-param", metavar='',help="markov parameters file (default is human/neanderthal like parameters)", type=str)
    train_subparser.add_argument("-out", metavar='',help="outputfile prefix (default is a file named trained.json)", default = 'trained.json')
    train_subparser.add_argument("-window_size", metavar='',help="size of bins (default is 1000 bp)", type=int, default = 1000)
    train_subparser.add_argument("-haploid",help="Change from using diploid data to haploid data (default is diploid)", action='store_true', default = False)

    # Decode model
    decode_subparser = subparser.add_parser('decode', help='Decode HMM')
    decode_subparser.add_argument("-obs",help="[required] file with observation data", type=str, required = True)
    decode_subparser.add_argument("-weights", metavar='',help="file with callability (defaults to all positions being called)")
    decode_subparser.add_argument("-mutrates", metavar='',help="file with mutation rates (default is mutation rate is uniform)")
    decode_subparser.add_argument("-param", metavar='',help="markov parameters file (default is human/neanderthal like parameters)", type=str)
    decode_subparser.add_argument("-out", metavar='',help="outputfile prefix (default is stdout)", default = '/dev/stdout')
    decode_subparser.add_argument("-window_size", metavar='',help="size of bins (default is 1000 bp)", type=int, default = 1000)
    decode_subparser.add_argument("-haploid",help="Change from using diploid data to haploid data (default is diploid)", action='store_true', default = False)
    decode_subparser.add_argument("-admixpop",help="Annotate using vcffile with admixing population (default is none)")


    args = parser.parse_args()


    if args.mode == 'make_test_data':
        create_test_data(data_set_length = args.windows, write_out_files = args.nooutfiles)


    elif args.mode == 'train':

        hmm_parameters = read_HMM_parameters_from_file(args.param)
        obs, _, _, _, mutrates, weights = Load_observations_weights_mutrates(args.obs, args.weights, args.mutrates, args.window_size, args.haploid)
        
        # Print parameters to screen
        print('-' * 40)
        print(hmm_parameters)
        print('> number of windows:', len(obs), '. Number of snps = ', sum(obs))
        print('> total callability:', round(np.sum(weights) / len(obs),2) )
        print('> average mutation rate per bin:', round(np.sum(mutrates * weights) / np.sum(weights), 2) )
        print('> Output is',args.out) 
        print('> Window size is',args.window_size, 'bp') 
        print('> Haploid',args.haploid) 
        print('-' * 40)

        hmm_parameters = TrainModel(obs, mutrates, weights, hmm_parameters)
        write_HMM_to_file(hmm_parameters, args.out)


    elif args.mode == 'decode':

        obs, chroms, starts, variants, mutrates, weights  = Load_observations_weights_mutrates(args.obs, args.weights, args.mutrates, args.window_size, args.haploid)
        hmm_parameters = read_HMM_parameters_from_file(args.param)
        
        # Print parameters to screen
        print('-' * 40)
        print(hmm_parameters)  
        print('> number of windows:', len(obs), '. Number of snps = ', sum(obs))
        print('> total callability:', round(np.sum(weights) / len(obs),2) )
        print('> average mutation rate per bin:', round(np.sum(mutrates * weights) / np.sum(weights), 2) )
        print('> Output is',args.out) 
        print('> Window size is',args.window_size, 'bp') 
        print('> Haploid',args.haploid) 
        print('-' * 40)

        # Find segments and write output
        segments = DecodeModel(obs, chroms, starts, variants, mutrates, weights, hmm_parameters)
        Write_Decoded_output(args.out, segments, args.obs, args.admixpop)



    elif args.mode == 'create_outgroup':
        make_out_group(args.ind, args.weights, args.vcf, args.out, args.ancestral, args.refgenome)


    elif args.mode == 'create_ingroup':
        make_ingroup_obs(args.ind, args.weights, args.vcf, args.out, args.outgroup, args.ancestral)


    elif args.mode == 'mutation_rate':
        print('-' * 40)
        print(f'> Outgroupfile:', args.outgroup)
        print(f'> Outputfile is:', args.out)
        print(f'> Callability file is:', args.weights)
        print(f'> Window size:', args.window_size)
        print('-' * 40)

        make_mutation_rate(args.outgroup, args.out, args.weights, args.window_size)
    else:
        print(print_script_usage())


if __name__ == "__main__":
    main()

