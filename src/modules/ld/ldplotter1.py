"""
Haploview does not scale to lots of SNPs, to get it work you need to downsample snps (3000 is a reasonable choice),
you also need to input the data as male x chromosome where haploidy is accepted
"""

import os
import random
from operator import itemgetter
import numpy as np


def vcf2haploview(vcf, snplimit, out):
    txt = open(out, 'w')
    with open(vcf) as file:
        line = file.readline()
        intro_count = 0
        while line[0] == '#':
            txt.write(line)
            intro_count += 1
            line = file.readline()
    all_vcf = open(vcf).readlines()[intro_count:]
    if len(all_vcf) <= int(snplimit):
        for line in all_vcf:
            line = line.split('\t')
            line[0] = '25'
            txt.write('\t'.join(line))
    else:
        subset_vcf = np.random.choice(all_vcf, int(snplimit), replace=False)
        subset_vcf = [item.split('\t') for item in subset_vcf]
        subset = []
        for snp in subset_vcf:
            snp[1] = int(snp[1])
            subset.append(snp)
        subset = sorted(subset, key=itemgetter(1))
        for line in subset:
            line[0] = '25'
            line = list(map(str, line))
            txt.write('\t'.join(line))
    txt.close()


# Running vcf2haploview
vcf = snakemake.input.vcfmaf
snplimit = snakemake.config["snplimit"]
out = snakemake.output.subset
vcf2haploview(vcf, snplimit, out)
