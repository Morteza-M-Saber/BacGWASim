#BacGWASim1.2
A simulator for Bacterial Machine learning and Genome-wide Association studies (BacGWASim v1.2)

----------

Motivation
----------
Identification of genomic elements underlying bacterial phenotypes such as resistance to antibiotics, virulence and fermentation is a fundamental step for a better understanding of the molecular mechanisms of these traits, and may be able to inform new clinical, industrial and agricultural interventions.

BacGWASim is designed to simulate whole-genomes and phenotypes of large bacterial populations focusing on unique characteristics of bacteria such as strong genome-wide linkage disequilibrium and population stratification. The data produced by BacGWASim provides a mean to 1) evaluate the power and limitations of existing bacterial Machine learning and Genome-Wide assoctation study methods and 2) develop and benchmark novel tools for bacterial genotype to phenotype mapping tools.

Citations
--------

Original BacGWASim implementation paper: `Saber, Morteza M., Shapiro, B Jesse` [Benchmarking bacterial genome-wide association study methods using simulated genomes and phenotypes.](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000337#tab2) `Microbial Genomics https://doi.org/10.1099/mgen.0.000337`

Prerequisites
-------------

Between parenthesis the versions the script was tested against:

* `python` 3+ (3.6.12)
* `numpy` (1.15.2)
* `scipy` (1.1.0)
* `pandas` (1.1.5)
* `ete3` (3.1.2)
* `pyvcf` (0.6.8)
* `snakemake`(5.31.1)
* `bcftools` (1.10.2)
* `plink` (1.9)
* `GCTA` (1.93.2)
* `snp-sites` (2.5.1)
* `openjdk` (8.0.152)
* `Haploview` (4.2)
* `simbac` (commit `5015897`)


Installation
------------

BacGWASim is build based on Snakemake and can be installed as following:

1)  clone workflow into working directory   
```    
git clone https://github.com/Morteza-M-Saber/BacGWASim.
cd BacGWASim 
```
2) Install dependencies:
The easiest way to install BacGWASim dependencies is through [`mamba`](https://github.com/mamba-org/mamba) (a modern alternative to `conda`). **It is recommended to install dependencies inside an isolated environment to avoid package conflicts**:

```bash
mamba env create --file requirements.yaml
```

Two dependency packages of  [GCTA (1.93.2)](https://cnsgenomics.com/software/gcta/#Download) and [simbac (commit 5015897)](https://github.com/tbrown91/SimBac) which are not available in conda channels should be installed manually and path to their executables need to be defined in `configfile.yaml`


2) edit config file to include the path to prerequisites and other parameters as needed
```
vim configfile.yaml
```

3) execute workflow, determine number of available cpu cores for parallelization
```
snakemake --snakefile Snakefile --cores 5 --latency-wait 60
```
In case the filesystem suffer from latency, consider increasing the `latency-wait` accordingly.

Simulation parameters
------------
All the simulation parameters are included in `configfile.yaml` file and can be adjusted:

```
#output directory name
outputDIR: BACGWASIM_quant_s100_m0.05_r0.01


#BacGWASim dependencies not provided by conda channels
simbac_path: dependencies/sim/SimBac/SimBac
haploview: dependencies/ld/Haploview.jar
gcta: dependencies/gcta/gcta_1.93.2beta/gcta64


#Genome simulation parameters
num_species: 100                 #Number of samples in the simulated population
genome_length: 10000             #Length of the genome (.bp)
mutation_rate: 0.05              #Mutation rate
r_i: 0.01    #Recombination rate
random_seed: 1487                #Random seed for reproducibility of results
maf: 0.01                        #Minor allele frequency threshold of rare alleles to be discarded

#Phenotype simulation parameters
phenType: quant                  #Type of simulated phenotype,'cc':binary case-control, 'quant': quantitative
causalMAF: 0.1                   #Minimum MinorAlleleFrequency of causal markers
causalMaxMAF: 0.4                #Maximum MinorAlleleFrequency of causal markers
causalMaxLD: 0.6                 #Maximum permitted r2 score between pairs of causal markers in window size of 1000 candidate causal markers meeting causalMAF and causalMaxMaf thresholds
causal_variant_Num: 16           #Number of causal markers
effect_size_ODR: 2,3,4,7,10,11,15,20    #Effect sizes of causal markers (.odds ratios) (comma separated,Must be a multiple of causal_variant_Num)
phenReplication: 10              #Number of phenotype replication sets
heritability: 1    #Heritability of phenotype
diseasePrevalence: 0.5           #Prevalence of phenotype
#in case of case-control binary phenotype simulation, number of case and control samples must be defined by 'case' and 'control' parameters
case: 50
control: 50

#Linkage Disequilibrium plotting
snplimit: 3000                   #Number of SNPs randomly selected for plotting linkage map (Increasing this number will significatnly increase computation time and require increasing the java heap size
heapSize: 1000                   #java heap_size for ld plot visualization (.mb)
ldmaf: 0.1                       #Minimum MinorAlleleFrequency of markers for LD plotting (By reducing value, it gets harder to accurately estiamte statistical r2 values between pairs of markers leading to more noisy plot)
```
Outputs
------------
BacGWASim produces the following outputs:
```
genSim/
  sims.vcf               #Multi-sample variant call file excluding rare alleles in vcf format (used for phenotype simulation)
  sims.pickle            #Multi-sample variant call file excluding rare alleles in pandas pickle format (for machine learning analysis)
  sims_no_selection.vcf  #Multi-sample variant call file including rare alleles
  genSim.fasta           #Multiple-sequence alignment of simulated genomes in fasta format
  phylogeny.nwk         #phylogenetic tree of the simulated population in newick format

phenSim/
  phenSim.phen           #simulated phenotypes
  phenSim.par            #ID, MAF and effect size of causal variants
  phenSim.pickle         #simulated phenotpes in pandas pickle format
  simVis.png             #Distribution of phenotypes and casual variants in the population

ld/
  ldPlot.png             #Genome-wide linkage disequilibrium (LD) plot
  ldPlot.LD              #Estimated pairwise r2 scores between sites used for plotting genome-wide LD plot

```

Usage for Machine-learning models
------------
BacGWASim outputs the simulated genotype matrix and phenotypes in pickle format that could be used to train machine-learning models as below:

```
import pandas as pd
X=pd.read_pickle(sims.pickle)
Y=pd.read_pickle(phenSim.pickle)

```

Usage for GWAS tools
------------
BacGWASim also outputs the simulated genotype in vcf format and phenotypes in gcta-compliant format that could be used by various GWAS tools.

Examples
------------
BacGWASim graphs the simulated genotype-phenotypes as color-coded phylogenetic tree and causal variants presence-absence as heat-map.

1. Simulation of a population of size 40 with 50 segregating sites as causal variants. 
![alt text](https://github.com/Morteza-M-Saber/BacGWASim/blob/master/Img/mytree40_50.png)

2. Simulation of a population of size 100 with 50 segregating sites as causal variants.
![alt text](https://github.com/Morteza-M-Saber/BacGWASim/blob/master/Img/mytree100_50.png)

3. Simulation of a population of size 500 with 10 segregating sites as causal variants.
![alt text](https://github.com/Morteza-M-Saber/BacGWASim/blob/master/Img/mytree500_10.png)

# Verification
With implemented evolutionary models and parameters for simulation of phylogenetic tree, genome content and phenotype, BacGWASim seeks to accurately simulate the important characteristics of bacterial genomes and populations. The  main feature that distinguishes bacterial populations from eukaryotes is genome-wide linkage disequilibrium (LD).    

Genome-wide Linkage disequilibrium is one of the main features of bacterial population due to their clonal expansion. However, the genome-wide LD can be disrupted by homologous recombination which occur in bacteria in various degrees depending on the species. LD plots show that BacGWASim correctly captures the genome-wide LD while no recombination is simulated with decreasing LD levels in proportion with rate of homologous recombination simulated

![alt text](https://github.com/Morteza-M-Saber/BacGWASim/blob/master/Img/LDRangeComparison.png)

The distribution of linkage, measured as r2 scores binned into three categories, is shown as a function of distance in the genome, in units of kilobase pairs, within genomes of simulated (left panels) and real (right panels) bacteria.
