# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 13:37:32 2018

@author: morte
"""
from subprocess import call
import pandas as pd

rule PhenSim:
    input:
        plink=expand("{outputDIR}/simulations/genSim/sims.bim",outputDIR=config["outputDIR"]),
        plinkbed=expand("{outputDIR}/simulations/genSim/sims.bed",outputDIR=config["outputDIR"]),
        plinkfam=expand("{outputDIR}/simulations/genSim/sims.fam",outputDIR=config["outputDIR"]),
        causal=expand("{outputDIR}/simulations/phenSim/{{replication_index}}/causalLoci.snplist",outputDIR=config["outputDIR"]),
    output:
        gtcaPar=expand("{outputDIR}/simulations/phenSim/{{replication_index}}/phenSim.par",outputDIR=config["outputDIR"]),
        gtcaPhen=expand("{outputDIR}/simulations/phenSim/{{replication_index}}/phenSim.phen",outputDIR=config["outputDIR"]),
    params:
        phenType=config['phenType'],
        herit=config['heritability'],
        preval=config['diseasePrevalence'],
        case=config['case'],
        control=config['control'],
        gcta=config['gcta'],
        shellCallFile=os.path.join(config["outputDIR"],'BacGWASim.log'),
    log: temp(expand("{outputDIR}/simulations/phenSim/{{replication_index}}/phenSim.log",outputDIR=config["outputDIR"])),
    run:
        if params.phenType == 'cc':
          CallString="%s --bfile %s --simu-cc %s %s --simu-causal-loci %s --simu-hsq %s --simu-k %s  --out %s " %(params.gcta,str(input.plink)[:-4],params.case,params.control,input.causal,params.herit,params.preval,str(output.gtcaPar)[:-4])
        elif params.phenType == 'quant':
          CallString="%s --bfile %s --simu-qt --simu-causal-loci %s --simu-hsq %s --simu-k %s  --out %s " %(params.gcta,str(input.plink)[:-4],input.causal,params.herit,params.preval,str(output.gtcaPar)[:-4])
        call('echo %s >> %s' %(CallString,params.shellCallFile),shell=True)
        call(CallString,shell=True)
        #Convert phenotypes to matrix
        df=pd.read_csv(output.gtcaPhen[0],sep=' ',header=None,index_col=0)
        df.columns=['Sample','phenotype','nan']
        df.drop(['Sample','nan'],axis=1,inplace=True)
        df.to_pickle(output.gtcaPhen[0][:-5]+'.pickle')