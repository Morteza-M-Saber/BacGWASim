#writing the snakemake file to simulate the tree
# -*- coding: utf-8 -*-
"""
Created on Fri june 11

@author: Masih
ای که خواهان تولدی دیگری, نخست مرگ را پذیرا باش

Simulate phylogentic tree

"""

# Used for system calls.
from subprocess import call

# Used for timestamping the log files.
from time import localtime, strftime

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# SNAKEMAKE RULE #
import os

    
rule treeSim:
    output: 
        phylogeny=expand("{outputDIR}/phylogeny/phylogeny.nwk",outputDIR=config["outputDIR"]),
    params:
        num_species=config["num_species"],
        birth_rate=config["birth_rate"],
        death_rate=config["death_rate"],
        birth_rate_sd=config["birth_rate_sd"],
        death_rate_sd=config["death_rate_sd"],
        logNAME="Simulating phylogenetic tree." + strftime("%Y-%m-%d.%H-%M-%S", localtime()),
        python=config['python'],
        treeSim=os.path.join('modules','treeSim','treeSim.py'),
        shellCallFile=os.path.join(config["outputDIR"],'BacGWASim.log'),
    run:  
        callString='python3 %s --num_species %s --output %s --birth_rate %s --death_rate %s \
                    --birth_rate_sd %s --death_rate_sd %s' % (params.treeSim,params.num_species,output.output,
                                                              params.birth_rate,params.death_rate,params.birth_rate_sd,params.death_rate_sd)
        call('echo "' + str(params.logNAME) + ':\n ' + callString + '\n" >> ' + params.shellCallFile, shell=True)
        call(callString, shell=True)


