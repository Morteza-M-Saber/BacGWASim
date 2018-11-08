# A simulator for Bacterial Genome-wide Association studies (BacGWASim)
# Motivation
Identification of genomic elements underlying bacterial phenotypes such as resistance to antibiotics, virulence and fermentation is a fundamental step for a better understanding of the molecular mechanisms of these traits, and may be able to inform new clinical, industrial and agricultural interventions. So far several methods for Bacterial Genome-wide association study (GWAS) have been devised to tackle this challenge, however, due to unique characteristics of bacterial genome, each of the methods have been shown to have certain limitations. 

BacGWASim is designed to simulate whole-genomes and phenotypes of large bacterial populations focusing on unique characteristics of bacteria such as strong genome-wide linkage disequilibrium and homoplasy events under a range of realistc evolutionary regimens, which affects the power of bacterial GWAS tools. The data produced by BacGWASim provides a mean to 1) evaluate the power and limitations of existing bacterial GWAS methods, 2) develop and benchmark novel tools for bacterial GWAS.


# Prerequisites


Between parenthesis the versions the script was tested against:

**python3+** (3.6.6)    
**Snakemake** (5.3.0)     
**ALFsim** (1.0.0)    
**DAWG** (2.0.0)    
**numpy** (1.15.2)    
**scipy** (1.1.0)   
**pandas** (0.23.4)   
**ART** (2016.06.05)    
**Samtools** (1.9.0)    
**BWA** (0.7.17)    
**GATK** (3.8.0)    
**Plink** (1.9)   
**gcta** (1.26.0)

If visualization parameter is set to True:

**Matplotlib** (3.0.1)     
**ete3** (3.1.1)

# Installation

BacGWASim is build based on Snakemake and can be installed as following:

1)  clone workflow into working directory   
```    

git clone https://bitbucket.org/user/myworkflow.git path/to/workdir
cd path/to/workdir 

```
2) edit config file to include the path to prerequisites and workflow as needed
```
vim config.yaml
```

3) execute workflow, deploy software dependencies via conda
```
snakemake -n --use-conda
```
# Usage
BacGWASim accepts a bacterial genome in fasta format and its annotation in GFF format. The phylogenetic tree can be user-defined or simulated by BacGWASim using a birth-death model.     
The program outputs the simulated genotypes and phenotypes for each species of the phylogenetic tree.

To get the parameter information of the tool:
```
python BacGWASim.py -h
```
Output:
```
Usage: BacGWASim.py [options]

Options:
  -h, --help            show this help message and exit   

  -S SETTING, --Setting=SETTING
                        Absolute path to the setting file [required]    

  -J CORENUMBER, --CoreNumber=CORENUMBER
                        Use at most N cores in parallel (default: 1)    

  -D                    Print the directed acyclic
                        graph of jobs. (Default:False)    

  -R                    Print the dependency graph of rules.
                        *Caution*Rulegraph can not be set to True if
                        DirectedAcyclicGraph is set to True (default: False)    

  -Q                    Do not output any progress or rule information
                        (default: False)
```
# Examples
If visualization is set to True, BacGWASim graphs the simulated genotype-phenotypes as color-coded phylogenetic tree and causal variants presence-absence as heat-map.

1. Simulation of a population of size 40 with 50 segregating sites as causal variants. 
