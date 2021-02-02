###simulate phylogentic tree and sequences 

def get_options():
    import argparse

    description = 'simulate bacterial genomes'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--num_species',
                        help='number of extant tips of the phylogenetic tree')
    parser.add_argument('--genome_length',
                        help='length of genome')

    parser.add_argument('--output',
                        default='simbac_out',
                        help='phylogeny filename [Default: simbac_sim]')
       
    parser.add_argument('--mutation_rate',
                default=0.01,
                help='mutaton rate [Default: 0.01]')
    
    parser.add_argument('--r_i',
            default=0.01,
            help='Sets the site-specific rate of internal (within species) recombination [Default: 0.01]')
    
    parser.add_argument('--r_e',
        default=0.0,
        help='Sets the site-specific rate of external (between species) recombination [Default: 0.0]')
    
    parser.add_argument('--random_seed',
        default=1487,
        help='given seed to initiate random number generation for reproducability [Default: 1487]')
    parser.add_argument('--simbac_path',
        help='path to simbac')

    return parser.parse_args()

options = get_options()
from subprocess import call

def genSim(num_species,genome_length,output,mutation_rate,r_i,r_e,random_seed,simbac_path):
    callString='%s -N %s -B %s -R %s -r %s -T %s -o %s -c %s -s %s' % (simbac_path,num_species,genome_length,
                                                              r_i,r_e,mutation_rate,output+'.fasta',output+'.nwk',random_seed)
    call(callString, shell=True)
    

genSim (options.num_species,options.genome_length,options.output,
        options.mutation_rate,options.r_i,options.r_e,
        options.random_seed,options.simbac_path)