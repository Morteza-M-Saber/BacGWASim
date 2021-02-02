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
        sam=temp(expand("{outputDIR}/simulations/data/{{name}}.sam",outputDIR=config["outputDIR"])),
        samef=temp(expand("{outputDIR}/simulations/data/{{name}}_errFree.sam",outputDIR=config["outputDIR"])),
        bamsorted=temp(expand("{outputDIR}/simulations/data/{{name}}",outputDIR=config["outputDIR"])),
        bamindex=temp(expand("{outputDIR}/simulations/data/{{name}}.bai",outputDIR=config["outputDIR"])),
        vcf=expand("{outputDIR}/simulations/data/{{name}}.vcf",outputDIR=config["outputDIR"]),
    params:
        outDir=config['outputDIR'],
        cpus=config['cpus'],
        logNAME="NGS simulation and Calling the variants from simulations." + strftime("%Y-%m-%d.%H-%M-%S", localtime()),
        shellCallFile=os.path.join(config["outputDIR"],'BacGWASim.log'),
    log: temp(expand("{outputDIR}/simulations/data/art{{name}}.log",outputDIR=config["outputDIR"])),
    run:  
        #simulate NGS data and error-fre sam alignment
        callString="art_illumina -ss HS25 -i %s -ir 0 -ir2 0 -dr 0 -dr2 0 -qs 93 -qs2 93 -ef -na -p -l 150 -f 30 -m 200 -s 10 -o %s &> %s" % \
        (input.fa,str(output.fq1)[:-4],log)
        call('echo "' + str(params.logNAME) + ':\n ' + callString + '\n" >> ' + params.shellCallFile, shell=True)
        call(callString,shell=True)
        #convert sam-file to bamfile and sort it 
        callString="samtools view -S -b %s | samtools sort - > %s" % \
        (output.samef,output.bamsorted)
        call('echo "' + str(params.logNAME) + ':\n ' + callString + '\n" >> ' + params.shellCallFile, shell=True)
        call(callString,shell=True)
        #index the bamfile
        callString="samtools index %s" % \
        (output.bamsorted)
        call('echo "' + str(params.logNAME) + ':\n ' + callString + '\n" >> ' + params.shellCallFile, shell=True)
        call(callString,shell=True)
        #call the variatn using bcftools mpileup,sample id will be whatever you insert as vcf so we need to correct it
        callString="bcftools mpileup -f %s %s | bcftools call -mv -Ov -o %s --ploidy 1 " % \
        (input.ref,output.bamsorted,output.vcf)
        call('echo "' + str(params.logNAME) + ':\n ' + callString + '\n" >> ' + params.shellCallFile, shell=True)
        call(callString,shell=True)
        