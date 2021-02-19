#writing the snakemake file to call variants from simulation
# -*- coding: utf-8 -*-
"""
?? ?? ?????? ????? ?????, ???? ??? ?? ????? ???

"""

# Used for system calls.
from subprocess import call

# Used for timestamping the log files.
from time import localtime, strftime
import os
import numpy as np
#-----------------------------------------------------------------------------------
# SNAKEMAKE RULE #
   
rule snippy2sim:
    input: 
        vcf=expand("{outputDIR}/simulations/data/{name}.vcf",outputDIR=config["outputDIR"],name=np.array(range(config['num_species']))+1),
        ref=expand("{outputDIR}/simulations/data/reference.fa",outputDIR=config["outputDIR"]),
    output: 
        vcf=expand("{outputDIR}/simulations/data/sims.vcf",outputDIR=config["outputDIR"]),
        plink=expand("{outputDIR}/simulations/data/sims.bim",outputDIR=config["outputDIR"]),
        causalPool=expand("{outputDIR}/simulations/data/causalPool.bim",outputDIR=config["outputDIR"]),
    params:
        outDir=config['outputDIR'],
        cpus=config['cpus'],
        causalMAF=config['causalMAF'],
        causalMaxLD=config['causalMaxLD'],
        logNAME="Calling the variants from simulations." + strftime("%Y-%m-%d.%H-%M-%S", localtime()),
        snipp2sim=os.path.join('modules','snippy','snippy2sim.py'),
        shellCallFile=os.path.join(config["outputDIR"],'BacGWASim.log'),
    threads: 
        config['cpus'],
    run:  
        #merging vcf files
        vcflist=os.path.join(params.outDir,'simulations','data','vcflist.txt')
        txt=open(vcflist,'w')
        for vcf_ in input.vcf:
          call('bgzip %s'% vcf_,shell=True)
          call('tabix %s'% (vcf_+'.gz'),shell=True)
          txt.write('%s\n'% (vcf_+'.gz'))
        txt.close()
        callString="bcftools merge  -l %s --missing-to-ref --threads %s -Ob -m none >  %s" %(vcflist,threads, os.path.join(params.outDir,'simulations','data','core_badhead.bcf'))
        call('echo "' + str(params.logNAME) + ':\n ' + callString + '\n" >> ' + params.shellCallFile, shell=True)
        call(callString, shell=True)
        #correcting sample names
        vcflist=os.path.join(params.outDir,'simulations','data','vcflist.txt')
        txt=open(vcflist,'w')
        for vcf_ in input.vcf:
          txt.write('%s %s\n'%(str(vcf_)[:-4],os.path.split(str(vcf_))[1][:-4]))
        txt.close()
        callString="bcftools reheader --samples %s %s -o %s" %(vcflist, os.path.join(params.outDir,'simulations','data','core_badhead.bcf'),os.path.join(params.outDir,'simulations','data','core.bcf'))
        call('echo "' + str(params.logNAME) + ':\n ' + callString + '\n" >> ' + params.shellCallFile, shell=True)
        call(callString, shell=True)
        callString="bcftools view %s -Ov -o %s" %(os.path.join(params.outDir,'simulations','data','core.bcf'),os.path.join(params.outDir,'simulations','data','core.vcf'))
        call('echo "' + str(params.logNAME) + ':\n ' + callString + '\n" >> ' + params.shellCallFile, shell=True)
        call(callString, shell=True)

        #snippy2sim (normalize and annotate SNPs)
        callString='python3 %s  --vcf %s --ref %s  --out %s' % (params.snipp2sim,os.path.join(params.outDir,'simulations','data','core.vcf'),os.path.join(params.outDir,'simulations','data','reference.fa'),output.vcf)
        call('echo "' + str(params.logNAME) + ':\n ' + callString + '\n" >> ' + params.shellCallFile, shell=True)
        call(callString, shell=True)
        #sim2plink (GCTA accepted format)
        callString='plink  --vcf %s --maf 0.01 --allow-extra-chr --make-bed --out %s' % (output.vcf,str(output.plink)[:-4])
        call('echo "' + str(params.logNAME) + ':\n ' + callString + '\n" >> ' + params.shellCallFile, shell=True)
        call(callString, shell=True)
        #Pruning SNP pairs in highLD for phenotype simulation
        callString='bcftools +prune  -l %s -w 1000 %s -Ov -o %s' % (params.causalMaxLD,output.vcf,str(output.vcf)[:-4]+'_LDpruned.vcf')
        call('echo "' + str(params.logNAME) + ':\n ' + callString + '\n" >> ' + params.shellCallFile, shell=True)
        call(callString, shell=True)
        #causalMarkerPool
        callString='plink  --vcf %s --maf %s --allow-extra-chr --make-bed --out %s' % (str(output.vcf)[:-4]+'_LDpruned.vcf',params.causalMAF,str(output.causalPool)[:-4])
        call('echo "' + str(params.logNAME) + ':\n ' + callString + '\n" >> ' + params.shellCallFile, shell=True)
        call(callString, shell=True)