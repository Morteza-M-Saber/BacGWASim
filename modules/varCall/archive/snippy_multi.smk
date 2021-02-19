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
   
rule snippy:
    input: 
        fa=expand("{outputDIR}/simulations/data/{{name}}.fa",outputDIR=config["outputDIR"]),
        ref=expand("{outputDIR}/simulations/data/reference.fa",outputDIR=config["outputDIR"]),
    output: 
        fq1=temp(expand("{outputDIR}/simulations/data/{{name}}1.fq",outputDIR=config["outputDIR"])),
        fq2=temp(expand("{outputDIR}/simulations/data/{{name}}2.fq",outputDIR=config["outputDIR"])),
        scr=temp(expand("{outputDIR}/simulations/data/{{name}}_NGSsim.screen",outputDIR=config["outputDIR"])),
        vcf=expand("{outputDIR}/simulations/data/{{name}}/snps.decompose.vcf",outputDIR=config["outputDIR"]),
    params:
        outDir=config['outputDIR'],
        cpus=config['cpus'],
        logNAME="NGS simulation and Calling the variants from simulations." + strftime("%Y-%m-%d.%H-%M-%S", localtime()),
        shellCallFile=os.path.join(config["outputDIR"],'BacGWASim.log'),
    threads: 
        config['cpus'],
    log: expand("{outputDIR}/simulations/data/snippy{{name}}.log",outputDIR=config["outputDIR"]),
    run:  
        #simulate NGS data
        callString="art_illumina -ss HS25 -i %s -ir 0 -ir2 0 -dr 0 -dr2 0 -qs 93 -qs2 93 -na -p -l 150 -f 30 -m 200 -s 10 -o %s &> %s" % \
        (input.fa,str(output.fq1)[:-4],output.scr)
        call('echo "' + str(params.logNAME) + ':\n ' + callString + '\n" >> ' + params.shellCallFile, shell=True)
        call(callString,shell=True)
        #snippy
        callString='snippy --R1 %s --R2 %s --ref %s --cpus %s --outdir %s --force 2> %s' % \
        (output.fq1,output.fq2,input.ref,threads,str(output.fq1)[:-4],log)
        call('echo "' + str(params.logNAME) + ':\n ' + callString + '\n" >> ' + params.shellCallFile, shell=True)
        call(callString, shell=True)
        #convert MNPs to SNPs
        callString="vt decompose_blocksub %s -o %s" % \
        (os.path.join(os.path.split(str(output.vcf))[0],'snps.filt.vcf'),output.vcf)
        call('echo "' + str(params.logNAME) + ':\n ' + callString + '\n" >> ' + params.shellCallFile, shell=True)
        call(callString, shell=True)