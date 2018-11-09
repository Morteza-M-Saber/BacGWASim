# A simulator for Bacterial Genome-wide Association studies (BacGWASim)
# Motivation
Identification of genomic elements underlying bacterial phenotypes such as resistance to antibiotics, virulence and fermentation is a fundamental step for a better understanding of the molecular mechanisms of these traits, and may be able to inform new clinical, industrial and agricultural interventions. So far several methods for Bacterial Genome-wide association study (GWAS) have been devised to tackle this challenge, however, due to unique characteristics of bacterial genome, each of the methods have been shown to have certain limitations. 

BacGWASim is designed to simulate whole-genomes and phenotypes of large bacterial populations focusing on unique characteristics of bacteria such as strong genome-wide linkage disequilibrium and homoplasy events under a range of realistc evolutionary regimens, which affects the power of bacterial GWAS tools. The data produced by BacGWASim provides a mean to 1) evaluate the power and limitations of existing bacterial GWAS methods, 2) develop and benchmark novel tools for bacterial GWAS.

# Outline
BacGWASim is able to simulate whole-genomes of large bacterial populations along with the corresponding phenotypes based on the emergence and evolution of causal variants within the population. In the first step, using a forward-in-time birth-death model a phylogenetic tree is generated which is tunable in aspects of birth rate (λ), death rate(µ), population size and root-leaf distance. The resulting phylogenetic tree accurately simulates clonal expansion of bacterial populations. In the second step, BacGWASim accepts a bacterial genome and the corresponding genome annotations and divides the genome into four categories of DNA elements, namely, protein coding genes, intergenic regions, RNA genes and pseudogenes/repetitive sequences that may have variation in influence of evolutionary forces on them. The simulation of evolutionary events is performed independently for each DNA category at three levels: 1) Site level events such as codon-substitution model, nucleotide-substitution model, insertion/deletion events and rate heterogeneity across sites, 2) Gene level events such as gene duplication, gene loss and gene translocation and 3) Genome level events such as recombination, lateral gene transfer (HGT) and genome rearrangements. Evolutionary rates and parameters are tunable for each process and for each category of DNA elements. After simulation of genome sequences, in the final step, a set of causal variants are chosen based on the minor allele frequency (MAF) upon which phenotype of each simulated individual is then simulated using genetic additive model. The phenotype simulation can be tuned by a number of parameters which are relevant to bacterial phenotypes including 1) the number of causal variants, 2) Effect size range of causal variants, 3) Quantitative or binary phenotype, 4) Prevalence of the affected phenotype and 5) phenotype heritability. 
  
BacGWASim is written to be able to simulate various evolutionary scenarios by permuting over all possible combination of user-defined evolutionary parameters.  For this end, BacGWASim is written to maximize modularity and parallelizability so that the simulation of large number of genomes and their evolutionary events could be performed in reasonable amount of time provided enough computational resources. 

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
vim ConfigFile.yaml
```

3) execute workflow, deploy software dependencies via conda
```
python BacGWASim.py -S Snakefile -J 4
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
![alt text](https://github.com/Morteza-M-Saber/BacGWASim/blob/master/Img/mytree40_50.png)

2. Simulation of a population of size 100 with 50 segregating sites as causal variants.
![alt text](https://github.com/Morteza-M-Saber/BacGWASim/blob/master/Img/mytree100_50.png)

3. Simulation of a population of size 500 with 10 segregating sites as causal variants.
![alt text](https://github.com/Morteza-M-Saber/BacGWASim/blob/master/Img/mytree500_10.png)

# Verification
With implemented evolutionary models and parameters for simulation of phylogenetic tree, genome content and phenotype, BacGWASim seeks to accurately simulate the important characteristics of bacterial genomes and populations. The  main feature that distinguishes bacterial populations from eukaryotes is genome-wide linkage disequilibrium (LD).    

Genome-wide Linkage disequilibrium is one of the main features of bacterial population due to their clonal expansion. However, the genome-wide LD can be disrupted by homologous recombination which occur in bacteria in various degrees depending on the species. LD plots show that BacGWASim correctly captures the genome-wide LD while no recombination is simulated with decreasing LD levels in proportion with rate of homologous recombination simulated

![alt text](https://github.com/Morteza-M-Saber/BacGWASim/blob/master/Img/LDRange.jpg)

