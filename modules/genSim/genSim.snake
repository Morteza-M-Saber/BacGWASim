#writing the snakemake file to simulate the genome with simbac
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
import os
#-----------------------------------------------------------------------------------
# SNAKEMAKE RULE #
   
rule genSim:
    input: 
    output: 
        phylogeny=temp(expand("{outputDIR}/simulations/genSim/genSim.nwk",outputDIR=config["outputDIR"])),
        fasta=expand("{outputDIR}/simulations/genSim/genSim.fasta",outputDIR=config["outputDIR"]),
    params:
        num_species= config["num_species"],
        genome_length= config["genome_length"],
        mutation_rate= config["mutation_rate"],
        r_i= config["r_i"],
        random_seed= config['random_seed'],
        simbac_path=config['simbac_path'],
        outDir=config['outputDIR'],
        logNAME="Simulating genome and phylogeny." + strftime("%Y-%m-%d.%H-%M-%S", localtime()),
        genSim=os.path.join('modules','genSim','genSim.py'),
        shellCallFile=os.path.join(config["outputDIR"],'BacGWASim.log'),
    run:  
        callString='python3 %s --num_species %s --genome_length %s --output %s \
        --mutation_rate %s --r_i %s --r_e 0 --random_seed %s --simbac_path %s' % (params.genSim,params.num_species,params.genome_length,output.phylogeny[0][:-4],
                                                              params.mutation_rate,params.r_i,params.random_seed,params.simbac_path)
        call('echo "' + str(params.logNAME) + ':\n ' + callString + '\n" >> ' + params.shellCallFile, shell=True)
        call(callString, shell=True)

