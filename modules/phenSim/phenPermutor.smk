#writing the snakemake file to call variants from simulation
# -*- coding: utf-8 -*-
"""
?? ?? ?????? ????? ?????, ???? ??? ?? ????? ???

"""

import os
   
def indexerRep(wildcards):
    return(int(wildcards.replication_index))



rule phenPermutor:
    input:
        causalPool_bim="{output_dir}/simulations/genSim/causalPool.bim",
    output:
        causal=temp("{output_dir}/simulations/phenSim/{replication_index}/causalLoci.snplist"),
    params:
        iterationsPhen= config['causal_variant_Num'],
        shellCallFile=os.path.join(config["outputDIR"],'BacGWASim.log'),
        phenpermutor=os.path.join('modules','phenSim','phenPermutor.py'),
        indexrep=indexerRep,
        effectSize=config['effect_size_ODR'],
        logNAME="Generate  causalvariant with OR values for phenotype simulation." + strftime("%Y-%m-%d.%H-%M-%S", localtime()),
    run:
        causals=input.causalPool_bim
        out=os.path.join('/'.join(str(output.causal).split('/')[:-2]),str(params.indexrep),'causalLoci.snplist')
        callString='python3 %s --iterationsPhen %s --replication %s --PlinkFormatCausalBim %s \
        --effectSize %s --out %s ' %(params.phenpermutor,params.iterationsPhen, \
                                     params.indexrep,input.causalPool_bim,params.effectSize,out)
        call('echo "' + str(params.logNAME) + ':\n ' + callString + '\n" >> ' + params.shellCallFile, shell=True)
        call(callString, shell=True)
        