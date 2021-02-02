
def get_options():
    import argparse

    description = 'Recieves list of all variants in plink.bim format and generates a list of randomly \
    chosen causal variants and OR values to be used by GCTA phenotype simulator'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--iterationsPhen',
                        help='Number of causal variants ')
    parser.add_argument('--replication',
                        help='Index of replicatin for phenotype simulation ')
    parser.add_argument('--PlinkFormatCausalBim',
                        help='List of all variants in plink.bim format from which causal variants are randomly chosen')
    parser.add_argument('--effectSize',
                        help='comma-separated OR values to be assigned to causal variants ')
    parser.add_argument('--out',
                        help='Complete path to output directory ')

    return parser.parse_args()

options = get_options()
   
import random, os

def phenPermutor(iterationsPhen,replication,PlinkFormatCausalBim,effectSize,out):
        effect_size_ODR=[int(item) for item in effectSize.split(',')]
        snpList=[]
        with open(PlinkFormatCausalBim,'r') as file:
            line=file.readline()
            while line:
                if "*" not in line.split()[1]:
                    snpList.append(line.split()[1])
                line=file.readline()
        variantNums= int(iterationsPhen)
        reps=replication
        causalNow=random.sample(snpList,variantNums)
        effectSizeNow=effect_size_ODR*int(variantNums/len(effect_size_ODR))
        causalVar=[item for item in zip(causalNow,effectSizeNow)]
        txt=open(out,'w')
        from math import log
        for variants in causalVar:
            txt.write('%s\t%s\n'%(variants[0],log(variants[1])))
        txt.close()
        
phenPermutor(options.iterationsPhen,options.replication,options.PlinkFormatCausalBim,options.effectSize,options.out)
