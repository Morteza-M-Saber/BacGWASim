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
   
rule simbac2snippyPlugger:
    input: 
          expand("{outputDIR}/simulations/simbac_output.fasta",outputDIR=config["outputDIR"]),
    output: 
        fa=expand("{outputDIR}/simulations/data/{name}.fa",outputDIR=config["outputDIR"],name=np.array(range(config['num_species']))+1),
        ref=expand("{outputDIR}/simulations/data/reference.fa",outputDIR=config["outputDIR"]),
    params:
        outDir=config['outputDIR'],
        logNAME="Plugging simbac res to snippy." + strftime("%Y-%m-%d.%H-%M-%S", localtime()),
        simbac2snipp=os.path.join('modules','snippy','simbac2snippy.py'),
        shellCallFile=os.path.join(config["outputDIR"],'BacGWASim.log'),
    threads: 
        config['cpus'],
    run:  
        #simbac2snippyplugger
        callString='python3 %s  --infile %s --outdir %s' % \
        (params.simbac2snipp,input,os.path.join(params.outDir,'simulations','data'))
        call('echo "' + str(params.logNAME) + ':\n ' + callString + '\n" >> ' + params.shellCallFile, shell=True)
        call(callString, shell=True)