# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 10:46:26 2018

@author: Masih

Construct LD plot for the result of each simulation
"""
import numpy as np
from subprocess import call

rule LDPlotter:
    input:
        vcf="{output_dir}/simulations/genSim/sims.vcf",
    output:
        vcfmaf=temp("{output_dir}/simulations/ld/sims_maf.vcf"),
        subset=temp("{output_dir}/simulations/ld/sims_subset_x.vcf"),
        vcfmap=temp("{output_dir}/simulations/ld/sims.map"),
        vcfped=temp("{output_dir}/simulations/ld/sims.ped"),
        haploview_inmap=temp("{output_dir}/simulations/ld/haploview_in.map"),
        haploview_inped=temp("{output_dir}/simulations/ld/haploview_in.ped"),
        haploview="{output_dir}/simulations/ld/ldPlot.LD.PNG",
        tmp1=temp("{output_dir}/simulations/ld/sims.nosex"),
        tmp2=temp("{output_dir}/simulations/ld/sims.log"),
    params:
        haploview=config['haploview'],
        snplimit=config['snplimit'],
        ldmaf=config['ldmaf'],
        heapSize=config['heapSize'],
        shellCallFile=os.path.join(config["outputDIR"],'BacGWASim.log'),
        ldplotter1=os.path.join('modules','ld','ldplotter1.py'),
        ldplotter2=os.path.join('modules','ld','ldplotter2.py'),
    log: 
        "{output_dir}/simulations/ld/ld.log",
    run:
        #Filter based on MAF to remove noise
        CallString='bcftools view -q %s %s -Ov -o %s'%(params.ldmaf,input.vcf,output.vcfmaf)
        call('echo %s >> %s' %(CallString,params.shellCallFile),shell=True)
        call(CallString,shell=True)
        #subseting snps to <3000 and converting chromosome to x
        CallString='python %s --vcf %s --snplimit %s --out %s'%(params.ldplotter1,output.vcfmaf,params.snplimit,output.subset)
        call('echo %s >> %s' %(CallString,params.shellCallFile),shell=True)
        call(CallString,shell=True)
        #converting vcf to plink ped/map file 
        CallString='plink --vcf %s --recode --out %s &> %s'%(output.subset,str(output.vcfped)[:-4],log)
        call('echo %s >> %s' %(CallString,params.shellCallFile),shell=True)
        call(CallString,shell=True)
        #changing the gender to male and creating info file out of plink ped file
        CallString='python %s --infile %s --out %s'%(params.ldplotter2,str(output.vcfped)[:-4],str(output.haploview_inped)[:-4])
        call('echo %s >> %s' %(CallString,params.shellCallFile),shell=True)
        call(CallString,shell=True)
        #generatng LD profile using Haploview
        CallString="java -jar %s -nogui -pedfile %s -info %s -dprime -compressedpng -chromosome X -maxDistance 10000 -memory %s -out %s" %(params.haploview,output.haploview_inped,output.haploview_inmap,params.heapSize,str(output.haploview)[:-7])
        call('echo %s >> %s' %(CallString,params.shellCallFile),shell=True)
        call(CallString,shell=True)
            
