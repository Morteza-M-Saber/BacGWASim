# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 10:46:26 2018

@author: Masih

Construct LD plot for the result of each simulation
"""


rule ld_plotter_bcf:
    """ Filter based on MAF to remove noise. """
    input:
        vcf = rules.vcfRefiner_markers.output.vcf,
    output:
        vcfmaf = temp("{output_dir}/simulations/ld/sims_maf.vcf"),
    shell:
        "bcftools view "
        "-q {config[ld_maf]} {input.vcf} "
        "-Ov -o {output.vcfmaf}"


rule ld_plotter_helper1:
    """ Subsetting SNPs to < 3000 and converting chromosome to X. """
    input:
        vcfmaf = rules.ld_plotter_bcf.output.vcfmaf,
    output:
        subset = "{output_dir}/simulations/ld/sims_subset.vcf",
    script:
        "ld_plotter1.py"


rule ld_plink_r2:
    """ Get linkage disequilibrium r2 from plink """
    input:
        vcf_subset = rules.ld_plotter_helper1.output.subset,
    output:
        ld = "{output_dir}/simulations/ld/sims.ld",
        nosex = temp("{output_dir}/simulations/ld/sims.nosex")
    params:
        output = "{output_dir}/simulations/ld/sims",
    threads: 8
    log:
        "{output_dir}/logs/ld_plink_r2.log"
    shell:
        "plink "
        "--vcf {input.vcf_subset} "
        "--r2 square0 "
        "--out {params.output} "
        "--threads {threads} "
        "&> {log}"


rule ld_plot:
    """ Plotting the LD from the vcf subset. """
    input:
        ld = rules.ld_plink_r2.output.ld,
        vcf_subset = rules.ld_plotter_helper1.output.subset,
    output:
        plot = "{output_dir}/simulations/ld/ld_plot.png"
    script:
        "ld_plot.py"


rule ld_stackplot:
    """ Generating LD stacked plot. """
    input:
        ld = rules.ld_plink_r2.output.ld,
        vcf_subset = rules.ld_plotter_helper1.output.subset,
    output:
        ld_stackplot = "{output_dir}/simulations/ld/ld_stackplot.png",
    script:
        "ld_stackplot.py"