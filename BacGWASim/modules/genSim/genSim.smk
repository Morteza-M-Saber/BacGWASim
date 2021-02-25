#writing the snakemake file to simulate the genome with simbac
# -*- coding: utf-8 -*-
"""
Created on Fri june 11

@author: Masih
ای که خواهان تولدی دیگری, نخست مرگ را پذیرا باش

Simulate phylogentic tree

"""

rule genSim:
    output: 
        fasta = "{output_dir}/simulations/genSim/genSim.fasta",
        phylogeny = temp("{output_dir}/simulations/genSim/genSim.nwk"),
    shell:
        "{config[simbac_path]} "
        "-N {config[num_species]} "
        "-B {config[genome_length]} "
        "-R {config[r_i]} "
        "-r 0 "
        "-T {config[mutation_rate]} "
        "-o {output.fasta} "
        "-c {output.phylogeny} "
        "-s {config[random_seed]}"
