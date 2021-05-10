# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 10:46:26 2018

@author: Masih

Construct LD plot for the result of each simulation
"""


rule LDPlotter_bcf:
    """ Filter based on MAF to remove noise. """
    input:
        vcf = rules.vcfRefiner_markers.output.vcf,
    output:
        vcfmaf = temp("{output_dir}/simulations/ld/sims_maf.vcf"),
    shell:
        "bcftools view "
        "-q {config[ld_maf]} {input.vcf} "
        "-Ov -o {output.vcfmaf}"


rule LDPlotter_helper1:
    """ Subsetting SNPs to < 3000 and converting chromosome to X. """
    input:
        vcfmaf = rules.LDPlotter_bcf.output.vcfmaf,
    output:
        subset = temp("{output_dir}/simulations/ld/sims_subset_x.vcf"),
    script:
        "ldplotter1.py"


rule LDPlotter_plink:
    """ Converting vcf to plink ped/map file. """
    input:
        subset = rules.LDPlotter_helper1.output.subset,
    output:
        vcfmap = temp("{output_dir}/simulations/ld/sims.map"),
        vcfped = temp("{output_dir}/simulations/ld/sims.ped"),
    params:
        output = "{output_dir}/simulations/ld/sims"
    log:
        "{output_dir}/logs/LDPlotter_plink.log"
    shell:
        "plink "
        "--vcf {input.subset} "
        "--recode "
        "--out {params.output} "
        "&> {log}"


rule LDPlotter_helper2:
    """ Changing gender to male and creating info file out of plink ped file. """
    input:
        vcfmap = rules.LDPlotter_plink.output.vcfmap,
        vcfped = rules.LDPlotter_plink.output.vcfped,
    output:
        haploview_inmap = temp("{output_dir}/simulations/ld/haploview_in.map"),
        haploview_inped = temp("{output_dir}/simulations/ld/haploview_in.ped"),
    script:
        "ldplotter2.py"


rule LDPlotter:
    """ Generating LD profile using Haploview. """
    input:
        haploview_inmap = rules.LDPlotter_helper2.output.haploview_inmap,
        haploview_inped = rules.LDPlotter_helper2.output.haploview_inped,
    output:
        haploview = "{output_dir}/simulations/ld/ldPlot.LD.PNG",
        haploview_ld = "{output_dir}/simulations/ld/ldPlot.LD",
    params:
        output = "{output_dir}/simulations/ld/ldPlot"
    log:
        "{output_dir}/logs/LDPlotter.log"
    shell:
        "haploview "
        "-nogui "
        "-pedfile {input.haploview_inped} "
        "-info {input.haploview_inmap} "
        "-dprime "
        "-compressedpng "
        "-chromosome X "
        "-maxDistance 10000 "
        "-memory {config[heap_size]} "
        "-out {params.output} "
        "&> {log}"            


rule ld_stackplot:
    """ Generating LD stacked plot. """
    input:
        haploview_ld = "{output_dir}/simulations/ld/ldPlot.LD",
    output:
        ld_stackplot = "{output_dir}/simulations/ld/ld_stackplot.png",
    script:
        "ld_stackplot.py"