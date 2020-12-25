#This can be done using birht-death modein dendropy which is available in conda as well

def get_options():
    import argparse

    description = 'simulate phylogenetic tree using a birth-death model'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--num_species',
                        help='number of extant tips of the phylogenetic tree')

    parser.add_argument('--output',
                        default='phylogeny.nwk',
                        help='phylogeny filename [Default: phylogeny.png]')
       
    parser.add_argument('--birth_rate',
                default=1.0,
                help='birth rate [Default: 1.0]')
    
    parser.add_argument('--death_rate',
            default=0.5,
            help='Death rate [Default: 0.5]')
    
    parser.add_argument('--birth_rate_sd',
        default=0.1,
        help='SD of birth rate [Default: 0.1]')
    
    parser.add_argument('--death_rate_sd',
        default=0.1,
        help='SD of death rate [Default: 0.1]')

    return parser.parse_args()

options = get_options()

def treesim (num_species,output,birth_rate=1.0,death_rate=0.5,birth_rate_sd=0.1,death_rate_sd=0.1):
    from dendropy.simulate import treesim
    t = treesim.birth_death_tree(birth_rate=float(birth_rate), death_rate=float(death_rate), num_extant_tips=int(num_species),
                                 birth_rate_sd=float(birth_rate_sd),death_rate_sd=float(death_rate_sd))
    t.write(path=output,schema='newick')

treesim (options.num_species,options.output,options.birth_rate,
         options.death_rate,options.birth_rate_sd,options.death_rate_sd)