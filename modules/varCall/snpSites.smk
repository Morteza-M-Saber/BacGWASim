##use snp-sites to call the variants and retrieve pseuo-reference sequence

# -*- coding: utf-8 -*-
"""
@author: Masih
ای که خواهان تولدی دیگری, نخست مرگ را پذیرا باش

Call variants from genome simulations

"""

# Used for system calls.
from subprocess import call

# Used for timestamping the log files.
from time import localtime, strftime
import os
#-----------------------------------------------------------------------------------
# SNAKEMAKE RULE #
   
rule snpsites:
    input: 
        fasta=expand("{outputDIR}/simulations/genSim/genSim.fasta",outputDIR=config["outputDIR"]),
    output: 
        phylip=expand("{outputDIR}/simulations/genSim/core.phylip",outputDIR=config["outputDIR"]),
        vcf=temp(expand("{outputDIR}/simulations/genSim/core.vcf",outputDIR=config["outputDIR"])),
        aln=expand("{outputDIR}/simulations/genSim/core.snp_sites.aln",outputDIR=config["outputDIR"]),
    params:
        outDir=config['outputDIR'],
        logNAME="Calling variants from simulations." + strftime("%Y-%m-%d.%H-%M-%S", localtime()),
        shellCallFile=os.path.join(config["outputDIR"],'BacGWASim.log'),
    run:  
        callString='snp-sites -mvpr -o %s %s' % (str(output.vcf)[:-4],input.fasta)
        call('echo "' + str(params.logNAME) + ':\n ' + callString + '\n" >> ' + params.shellCallFile, shell=True)
        call(callString, shell=True)

