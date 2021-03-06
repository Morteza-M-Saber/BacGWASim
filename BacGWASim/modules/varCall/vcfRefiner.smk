#writing the snakemake file to call variants from simulation
# -*- coding: utf-8 -*-
"""
ای که خواهان تولدی دیگری, نخست مرگ را پذیرا باش

"""

rule vcfRefiner_file:
    """ Creating the reheader file needed by bcftools. """
    output:
        file_rehead = temp("{output_dir}/rehead.tsv"),
    shell:
        "echo \"0 zero\" > {output.file_rehead}"


rule vcfRefiner_reheader:
    """ Changing sample 0 to sample zero as required by plink. """
    input:
        vcf = rules.snpsites.output.vcf,
        file_rehead = rules.vcfRefiner_file.output.file_rehead,
    output:
        core_rehead = temp("{output_dir}/simulations/genSim/core_rehead.vcf"),
    shell:
        "bcftools reheader "
        "--samples {input.file_rehead} "
        "--output {output.core_rehead} "
        "{input.vcf}"


rule vcfRefiner_norm:
    input:
        core_rehead = rules.vcfRefiner_reheader.output.core_rehead,
    output:
        core_norm = temp("{output_dir}/simulations/genSim/sims_norm.vcf"),
    log:
        "{output_dir}/logs/vcfRefiner_norm.log"
    shell:
        "bcftools norm "
        "-m "
        "-any "
        "-Ov -o {output.core_norm} "
        "{input.core_rehead} "
        "&> {log}"


rule vcfRefiner_annotate:
    """ Adding an ID tag with chrom, pos, ref and alt. """
    input:
        core_norm = rules.vcfRefiner_norm.output.core_norm,
    output:
        vcf_no_selec = "{output_dir}/simulations/genSim/sims_no_selection.vcf",
    log:
        "{output_dir}/logs/vcfRefiner_annotate.log"
    shell:
        "bcftools annotate "
        "-x ID "
        "-I +%CHROM:%POS:%REF:%ALT "
        "-Ov -o {output.vcf_no_selec} "
        "{input.core_norm} "
        "&> {log}"


rule vcfRefiner_purifying:
    """ Simulating purification selection by removing rare mutations. """
    input:
        vcf_no_selec = rules.vcfRefiner_annotate.output.vcf_no_selec,
    output:
        vcf = temp("{output_dir}/simulations/genSim/sims_raw.vcf"),
    shell:
        "bcftools view {input.vcf_no_selec} "
        "-q {config[maf]} "
        "-oV -o {output.vcf}"


rule vcfRefiner_markers:
    """ Reducing number of markers per user request. """
    input:
        vcf = rules.vcfRefiner_purifying.output.vcf,
    output:
        vcf = "{output_dir}/simulations/genSim/sims.vcf"
    script:
        "vcfRefiner_helper1.py"


rule vcfRefiner_sim2plink:
    """
        plink is the GCTA accepted format.
        In Plink1.x if keep_allel_order is not used, plink discards
        vcf's ref/alt allele order and redefine it based of MAF
        which wreaks havoc in phenotype simulations. 
    """
    input:
        vcf = rules.vcfRefiner_markers.output.vcf
    output:
        plink_bed = temp("{output_dir}/simulations/genSim/sims.bed"),
        plink_bim = temp("{output_dir}/simulations/genSim/sims.bim"),
        plink_fam = temp("{output_dir}/simulations/genSim/sims.fam"),
        plink_log = temp("{output_dir}/simulations/genSim/sims.log"),
        plink_nosex = temp("{output_dir}/simulations/genSim/sims.nosex"),
    params:
        output = "{output_dir}/simulations/genSim/sims",
    log:
        "{output_dir}/logs/vcfRefiner_sim2plink.log"
    shell:
        "plink "
        "--vcf {input.vcf} "
        "--allow-extra-chr "
        "--keep-allele-order "
        "--make-bed "
        "--out {params.output} "
        "&> {log}"


rule vcfRefiner_causal:
    """ Eliminating varians based of MAF for phenotype simulations. """
    input:
        vcf = rules.vcfRefiner_markers.output.vcf,
    output:
        vcf_causal = temp("{output_dir}/simulations/genSim/simsCausal.vcf"),
    shell:
        "bcftools view {input.vcf} "
        "-q {config[causal_maf_min]} "
        "-Q {config[causal_maf_max]} "
        "-Ov -o {output.vcf_causal}"


rule vcfRefiner_pruning:
    """ Pruning SNP pairs in highLD for phenotype simulation. """
    input:
        vcf_causal = rules.vcfRefiner_causal.output.vcf_causal,
    output:
        vcf_pruned = temp("{output_dir}/simulations/genSim/sims_LDpruned.vcf"),
    shell:
        "bcftools +prune "
        "-l {config[causal_ld_max]} "
        "-w 1000 {input.vcf_causal} "
        "-Ov -o {output.vcf_pruned} "


rule vcfRefiner_causalPool:
    """ Creating causal marker pool. """
    input:
        vcf_pruned = rules.vcfRefiner_pruning.output.vcf_pruned,
    output:
        causalPool_bed = temp("{output_dir}/simulations/genSim/causalPool.bed"),
        causalPool_bim = temp("{output_dir}/simulations/genSim/causalPool.bim"),
        causalPool_fam = temp("{output_dir}/simulations/genSim/causalPool.fam"),
        causalPool_log = temp("{output_dir}/simulations/genSim/causalPool.log"),
        causalPool_nosex = temp("{output_dir}/simulations/genSim/causalPool.nosex"),
    params:
        output = "{output_dir}/simulations/genSim/causalPool",
    log:
        "{output_dir}/logs/vcfRefiner_causalPool.log"
    shell:
        "plink "
        "--vcf {input.vcf_pruned} "
        "--allow-extra-chr "
        "--keep-allele-order "
        "--make-bed "
        "--out {params.output} "
        "&> {log}"


rule vcfRefiner_matrix:
    """ Creating a matrix for machine/deep leaning techniques. """
    input:
        vcf = rules.vcfRefiner_markers.output.vcf
    output:
        pickle = "{output_dir}/simulations/genSim/sims.pickle"
    script:
        "vcfRefiner_helper2.py"
