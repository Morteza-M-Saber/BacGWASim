#writing the snakemake file to call variants from simulation
# -*- coding: utf-8 -*-
"""
ای که خواهان تولدی دیگری, نخست مرگ را پذیرا باش

"""

# Used for system calls.
from subprocess import call

# Used for timestamping the log files.
from time import localtime, strftime
import os
import numpy as np
#-----------------------------------------------------------------------------------
# SNAKEMAKE RULE #
   
rule vcfRefiner:
    input: 
        vcf="{output_dir}/simulations/genSim/core.vcf",
    output: 
        vcf="{output_dir}/simulations/genSim/sims.vcf",
        core_rehead=temp("{output_dir}/simulations/genSim/core_rehead.vcf"),
        causalvcf=temp("{output_dir}/simulations/genSim/simsCausal.vcf"),
        plink=temp("{output_dir}/simulations/genSim/sims.bim"),
        plinkbed=temp("{output_dir}/simulations/genSim/sims.bed"),
        plinkfam=temp("{output_dir}/simulations/genSim/sims.fam"),
        plinklog=temp("{output_dir}/simulations/genSim/sims.log"),
        plinknosex=temp("{output_dir}/simulations/genSim/sims.nosex"),
        causalPool=temp("{output_dir}/simulations/genSim/causalPool.bim"),
        causalPoolbed=temp("{output_dir}/simulations/genSim/causalPool.bed"),
        causalPoolfam=temp("{output_dir}/simulations/genSim/causalPool.fam"),
        causalPoollog=temp("{output_dir}/simulations/genSim/causalPool.log"),
        causalPoolnosex=temp("{output_dir}/simulations/genSim/causalPool.nosex"),
    params:
        outDir=config['outputDIR'],
        causalMAF=config['causalMAF'],
        causalMaxMAF=config['causalMaxMAF'],
        causalMaxLD=config['causalMaxLD'],
        maf=config['maf'],
        varNumber=config['varNumber'],
        logNAME="Calling the variants from simulations." + strftime("%Y-%m-%d.%H-%M-%S", localtime()),
        vcfRefiner=os.path.join('modules','varCall','vcfRefiner.py'),
        rehead=os.path.join('modules','varCall','rehead'),
        shellCallFile=os.path.join(config["outputDIR"],'BacGWASim.log'),
    log: 
        temp("{output_dir}/simulations/genSim/genSim.log"),
    run:  
        #plink doesn't accept 0 as sample ID so we need to change it
        callString="bcftools reheader --samples %s %s -o %s" %(params.rehead, input.vcf,os.path.join(params.outDir,'simulations','genSim','core_rehead.vcf'))
        call('echo "' + str(params.logNAME) + ':\n ' + callString + '\n" >> ' + params.shellCallFile, shell=True)
        call(callString, shell=True)
        #vcfRefiner (convert multi-allele to single allele and annotate SNPs)
        callString='python3 %s  --vcf %s --out %s' % (params.vcfRefiner,os.path.join(params.outDir,'simulations','genSim','core_rehead.vcf'),os.path.join(params.outDir,'simulations','genSim','sims_no_selection.vcf'))
        call('echo "' + str(params.logNAME) + ':\n ' + callString + '\n" >> ' + params.shellCallFile, shell=True)
        call(callString, shell=True)
        # Simulating purifying selection by removing rare mutations 
        callString="bcftools view %s -q %s -Ov -o %s" %(os.path.join(params.outDir,'simulations','genSim','sims_no_selection.vcf'),params.maf,output.vcf)
        call('echo "' + str(params.logNAME) + ':\n ' + callString + '\n" >> ' + params.shellCallFile, shell=True)
        call(callString, shell=True)
        # Reducing number of markers per user request
        if params.varNumber != -1:
          info_lines=[]
          info_line_count=0
          with open(str(output.vcf),'r') as file:
            line=file.readline()
            while line[0] == '#':
              info_lines.append(line)
              info_line_count+=1
              line=file.readline()
          var_lines = open(output.vcf,'r').readlines()[info_line_count:]
          if len(var_lines) - info_line_count <= params.varNumber:
            print('Number of simulated markers is less than requested!\nConsider increasing mutation rate!')
          else:
            rand_lines=np.random.choice(range(len(var_lines)),params.varNumber,replace=False)
            rand_indx=np.sort(rand_lines)
            txt=open(output.vcf,'w')
            for info_ in info_lines:
              txt.write(info_)
            for var_choice in rand_indx:
              txt.write(var_lines[var_choice])
            txt.close()
        #Create matrix for machine/deep learning techniques
        info_line_c=0
        with open(str(output.vcf),'r') as file:
          line=file.readline()
          while line[:2] == '##':
            info_line_c+=1
            line=file.readline()
        import pandas as pd
        df=pd.read_csv(str(output.vcf),sep='\t',skiprows=info_line_c)
        df=df.drop(['#CHROM','POS','REF','ALT','QUAL','FILTER','INFO','FORMAT'],axis=1)
        df.set_index('ID',inplace=True)
        df=df.transpose()
        df.to_pickle(str(output.vcf)[:-4]+'.pickle')
        #sim2plink (GCTA accepted format), in Plink1.x if keep_allele_order not used, plink discards vcf's ref/alt allele order
                                           #and redefine it based on MAF which wreaks havoc in phenotype simulation!!
        callString='plink  --vcf %s --allow-extra-chr --keep-allele-order --make-bed --out %s &> %s' % (output.vcf,str(output.plink)[:-4],log)
        call('echo "' + str(params.logNAME) + ':\n ' + callString + '\n" >> ' + params.shellCallFile, shell=True)
        call(callString, shell=True)
        #Eliminate variants based on MAF for phenotype simulation (Plink1.9 --max-maf does have bug in it! use bcftools)
        callString="bcftools view %s -q %s -Q %s -Ov -o %s" %(output.vcf,params.causalMAF, params.causalMaxMAF, output.vcf[:-4] + 'Causal.vcf')
        call('echo "' + str(params.logNAME) + ':\n ' + callString + '\n" >> ' + params.shellCallFile, shell=True)
        call(callString, shell=True)
        #Pruning SNP pairs in highLD for phenotype simulation
        callString='bcftools +prune  -l %s -w 1000 %s -Ov -o %s' % (params.causalMaxLD,output.vcf[:-4]+'Causal.vcf',str(output.vcf)[:-4]+'_LDpruned.vcf')
        call('echo "' + str(params.logNAME) + ':\n ' + callString + '\n" >> ' + params.shellCallFile, shell=True)
        call(callString, shell=True)
        #causalMarkerPool
        callString='plink  --vcf %s  --allow-extra-chr --keep-allele-order --make-bed --out %s &> %s' % (str(output.vcf)[:-4]+'_LDpruned.vcf',str(output.causalPool)[:-4],log)
        call('echo "' + str(params.logNAME) + ':\n ' + callString + '\n" >> ' + params.shellCallFile, shell=True)
        call(callString, shell=True)