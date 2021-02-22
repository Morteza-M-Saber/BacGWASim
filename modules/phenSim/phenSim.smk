# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 13:37:32 2018

@author: morte
"""
from subprocess import call
import pandas as pd

rule PhenSim:
    input:
        plink="{output_dir}/simulations/genSim/sims.bim",
        plinkbed="{output_dir}/simulations/genSim/sims.bed",
        plinkfam="{output_dir}/simulations/genSim/sims.fam",
        causal="{output_dir}/simulations/phenSim/{replication_index}/causalLoci.snplist",
    output:
        gtcaPar="{output_dir}/simulations/phenSim/{replication_index}/phenSim.par",
        gtcaPhen="{output_dir}/simulations/phenSim/{replication_index}/phenSim.phen",
        pickle="{output_dir}/simulations/phenSim/{replication_index}/phenSim.pickle",
    params:
        phenType=config['phenType'],
        herit=config['heritability'],
        preval=config['diseasePrevalence'],
        case=config['case'],
        control=config['control'],
        gcta=config['gcta'],
        shellCallFile=os.path.join(config["outputDIR"],'BacGWASim.log'),
    log: temp("{output_dir}/simulations/phenSim/{replication_index}/phenSim.log"),
    run:
        if params.phenType == 'cc':
          CallString="%s --bfile %s --simu-cc %s %s --simu-causal-loci %s --simu-hsq %s --simu-k %s  --out %s " %(params.gcta,str(input.plink)[:-4],params.case,params.control,input.causal,params.herit,params.preval,str(output.gtcaPar)[:-4])
        elif params.phenType == 'quant':
          CallString="%s --bfile %s --simu-qt --simu-causal-loci %s --simu-hsq %s --simu-k %s  --out %s " %(params.gcta,str(input.plink)[:-4],input.causal,params.herit,params.preval,str(output.gtcaPar)[:-4])
        call('echo %s >> %s' %(CallString,params.shellCallFile),shell=True)
        call(CallString,shell=True)
        #Convert phenotypes to matrix
        df=pd.read_csv(output.gtcaPhen,sep=' ',header=None,index_col=0)
        df.columns=['Sample','phenotype','nan']
        df.drop(['Sample','nan'],axis=1,inplace=True)
        df.to_pickle(output.pickle)