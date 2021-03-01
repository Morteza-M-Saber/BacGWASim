from math import log
import random
import os

# Snakemake variables
input_file = snakemake.input.causalPool_bim
output_file = snakemake.output.causal
iterationsPhen = snakemake.config["num_causal_var"]
effectSize = snakemake.config["effect_size_odr"]

snpList = []
variantNums = int(iterationsPhen)
effect_size_ODR = [int(item) for item in effectSize.split(',')]
with open(input_file, 'r') as file:
    line = file.readline()
    while line:
        if "*" not in line.split()[1]:
            snpList.append(line.split()[1])
        line = file.readline()
        
causalNow = random.sample(snpList, variantNums)
effectSizeNow = effect_size_ODR*int(variantNums/len(effect_size_ODR))
causalVar = [item for item in zip(causalNow, effectSizeNow)]
with open(output_file, 'w') as txt:
    for variants in causalVar:
        txt.write('%s\t%s\n' % (variants[0], log(variants[1])))
