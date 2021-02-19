# -*- coding: utf-8 -*-
"""
@author: Masih

Visualize population structure and distribution of causal variants
"""
import numpy as np
from subprocess import call

rule simViz:
    input:
        vcf=expand("{outputDIR}/simulations/genSim/sims.vcf",outputDIR=config["outputDIR"]),
        par=expand("{outputDIR}/simulations/phenSim/{{replication_index}}/phenSim.par",outputDIR=config["outputDIR"]),
        phen=expand("{outputDIR}/simulations/phenSim/{{replication_index}}/phenSim.phen",outputDIR=config["outputDIR"]),
        phylogeny=expand("{outputDIR}/simulations/genSim/genSim.nwk",outputDIR=config["outputDIR"]),
        haploview=expand("{outputDIR}/simulations/ld/ldPlot.LD.PNG",outputDIR=config["outputDIR"]),
    output:
        simViz=expand("{outputDIR}/simulations/phenSim/{{replication_index}}/sim{{replication_index}}.png",outputDIR=config["outputDIR"]),
    params:
        phenType=config['phenType'],
        simViz=os.path.join('modules','visualizers','simVis.py'),
        shellCallFile=os.path.join(config["outputDIR"],'BacGWASim.log'),
    run:
        #Replace sample_id:0 with sample_id:zero in the phylogenetic tree
        def phylo_tip_corrector (phylo, old_id,new_id):
                    pos1=phylo.find(','+old_id+':')
                    pos2=phylo.find('('+old_id+':')
                    if pos1 == -1 and pos2 == -1:
                      print('Sample %s does not exist in tree!'%old_id)
                      return (phylo)
                    else:
                      if pos1 != -1:
                        return phylo[:pos1]+','+new_id+':'+phylo[pos1+len(old_id)+2:]
                      else:
                        return phylo[:pos2]+'('+new_id+':'+phylo[pos2+len(old_id)+2:]
        txt= open(input.phylogeny[0],'r') 
        phylo=txt.readline().strip()
        old,new='0','zero'
        phylo_corrected=phylo_tip_corrector(phylo,old,new)
        path_=os.path.join(os.path.split(input.phylogeny[0])[0],'phylogeny.nwk')
        with open (path_,'w') as file:
          file.write(phylo_corrected)
        #run simViz on all simulations
        CallString='python %s --vcfIn %s --phen %s --phenType %s --par %s --phylo %s --out %s'%( \
                    params.simViz,input.vcf[0],input.phen[0],params.phenType,input.par[0],path_,output.simViz[0])
        call('echo %s >> %s' %(CallString,params.shellCallFile),shell=True)
        call(CallString,shell=True)
 