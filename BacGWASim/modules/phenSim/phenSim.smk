# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 13:37:32 2018

@author: morte
"""
from snakemake import shell


rule phenSim:
    input:
        plink = "{output_dir}/simulations/genSim/sims.bim",
        plinkbed = "{output_dir}/simulations/genSim/sims.bed",
        plinkfam = "{output_dir}/simulations/genSim/sims.fam",
        causal = "{output_dir}/simulations/phenSim/{replication_index}/causalLoci.snplist",
    output:
        gcta_par = "{output_dir}/simulations/phenSim/{replication_index}/phenSim.par",
        gcta_phen = "{output_dir}/simulations/phenSim/{replication_index}/phenSim.phen",
    params:
        plink_input = "{output_dir}/simulations/genSim/sims",
        gcta_output = "{output_dir}/simulations/phenSim/{replication_index}/phenSim",
    log:
        "{output_dir}/logs/phenSim.{replication_index}.log"
    run:
        if config["phen_type"] == "cc":
          simu = "--simu-cc {config[case]} {config[control]}"
        elif config["phen_type"] == "quant":
          simu = "--simu-qt"

        shell(
          "gcta64 "
          "--bfile {params.plink_input} "
          "--simu-causal-loci {input.causal} "
          "--simu-hsq {config[heritability]} "
          "--simu-k {config[disease_prevalence]} "
          "--out {params.gcta_output} " 
          + simu
          + " 1> {log}"
        )


rule phenSim_pickle:
    """ Convert phenotypes to matrix. """
    input:
        gtca_phen = rules.phenSim.output.gcta_phen,
    output:
        pickle="{output_dir}/simulations/phenSim/{replication_index}/phenSim.pickle",
    script:
        "phenSim_pickle.py"