# -*- coding: utf-8 -*-
"""
@author: Masih

Visualize population structure and distribution of causal variants
"""

rule simVis_phyloCorrector:
    """ Replacing sample_id:0 with sample_id:zero in the phylogenetic tree. """
    input:
        phylo = "{output_dir}/simulations/genSim/genSim.nwk",
    output:
        phylo = "{output_dir}/simulations/genSim/phylogeny.nwk",
    script:
        "simVis_phyloCorrector.py"


rule simVis_extract_vcf:
    input:
        vcf = "{output_dir}/simulations/genSim/sims.vcf",
        par = "{output_dir}/simulations/phenSim/{replication_index}/phenSim.par",
    output:
        vcf = "{output_dir}/simulations/phenSim/{replication_index}/sims.vcf",
    shell:
        "(head -n 50 {input.vcf} | awk '/^#/ {{print $0}}' > {output.vcf}) && "
        "grep \"$(awk '{{if (NR!=1) {{print $1}}}}' {input.par})\" {input.vcf} >> {output.vcf}"


rule simVis:
    input:
        vcf = rules.simVis_extract_vcf.output.vcf,
        par = "{output_dir}/simulations/phenSim/{replication_index}/phenSim.par",
        phen = "{output_dir}/simulations/phenSim/{replication_index}/phenSim.phen",
        phylo = rules.simVis_phyloCorrector.output.phylo
    output:
        simVis = "{output_dir}/simulations/phenSim/{replication_index}/sim{replication_index}.png",
    script:
        "simVis.py"
