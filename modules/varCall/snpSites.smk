##use snp-sites to call the variants and retrieve pseuo-reference sequence

# -*- coding: utf-8 -*-
"""
@author: Masih
ای که خواهان تولدی دیگری, نخست مرگ را پذیرا باش

Call variants from genome simulations

"""
   
rule snpsites:
    input: 
        fasta = rules.genSim.output.fasta
    output: 
        phylip = "{output_dir}/simulations/genSim/core.phylip",
        vcf = temp("{output_dir}/simulations/genSim/core.vcf"),
        aln = "{output_dir}/simulations/genSim/core.snp_sites.aln",
    params:
        output = "{output_dir}/simulations/genSim/core"
    shell:
        "snp-sites "
        "-mvpr "
        "-o {params.output} "
        "{input.fasta}" 

