# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 14:23:47 2018

@author: morte
"""
#To filter the combined VCF file based on MAF and create smaller version of VCF file 
#which could be handled by python vcf parser in resonable time
from subprocess import call
vcfOri='/scratch/masih/output500/simulations/0/VarCall/BCFToolsRes'
vcfComp='/scratch/masih/output500/simulations/0/VarCall/BCFToolsRes.gz'


#compressing the original VCF file to be feed into bcftools
bgzipTools='/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/intel2016.4/htslib/1.5/bin/bgzip'
call('%s -c %s > %s'%(bgzipTools,vcfOri,vcfComp),shell=True)

#indexing the compressed vcf file
tabixTool=' /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/intel2016.4/htslib/1.5/bin/tabix'
call('%s -p vcf %s'%(tabixTool,vcfComp),shell=True)


#filtering the original VCFFiles using bcftools
bcfTools='/home/masih/BacterialSimulator/bcfTools/bcftools/bcftools'
vcfFilter='/scratch/masih/output500/simulations/0/VarCall/BCFToolsResFilterMAF0.05'
call('%s view -q 0.05:minor %s -o %s'%(bcfTools,vcfComp,vcfFilter ),shell=True)