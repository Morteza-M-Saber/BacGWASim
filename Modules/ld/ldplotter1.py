"""
Haploview does not scale to lots of SNPs, to get it work you need to downsample snps (3000 is a reasonable choice),
you also need to input the data as male x chromosome where haploidy is accepted
"""
def get_options():
    import argparse

    description = 'Change chromosome number to X chromosome and downsample snp size to 3000 for haploview'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--vcf',
                        help='path to vcf file ')
    parser.add_argument('--snplimit',
                        default=3000,
                        help='path to vcf file ')
    parser.add_argument('--out',
                        help='Complete path to output directory ')

    return parser.parse_args()

options = get_options()
   
import random, os
import numpy as np
from operator import itemgetter

def vcf2haploview(vcf,snplimit,out):
  txt=open(out,'w')
  with open (vcf) as file:
    line=file.readline()
    intro_count=0
    while line[0]=='#':
      txt.write(line)
      intro_count+=1
      line=file.readline()
  all_vcf=open(vcf).readlines()[intro_count:]
  if len(all_vcf) <= int(snplimit):
    for line in all_vcf:
      line=line.split('\t')
      line[0]='25'
      txt.write('\t'.join(line))
  else:
    subset_vcf=np.random.choice(all_vcf,int(snplimit),replace=False)
    subset_vcf=[item.split('\t') for item in subset_vcf]
    subset=[]
    for snp in subset_vcf:
      snp[1]=int(snp[1])
      subset.append(snp)
    subset=sorted(subset, key=itemgetter(1))
    for line in subset:
      line[0]='25'
      line=list(map(str,line))
      txt.write('\t'.join(line))
  txt.close()
        
vcf2haploview(options.vcf,options.snplimit,options.out)
