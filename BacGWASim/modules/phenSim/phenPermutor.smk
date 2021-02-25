
rule phenPermutor:
    input:
        causalPool_bim = "{output_dir}/simulations/genSim/causalPool.bim",
    output:
        causal = temp("{output_dir}/simulations/phenSim/{replication_index}/causalLoci.snplist"),
    script:
        "phenPermutor.py"