[![Anaconda-Server Badge](https://anaconda.org/bioconda/bacgwasim/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)
[![Bioconda](https://img.shields.io/conda/dn/bioconda/bacgwasim.svg?label=Bioconda)](https://bioconda.github.io/recipes/bacgwasim/README.html)

# BacGWASim

A simulator for Bacterial Machine learning and Genome-wide Association studies (BacGWASim v2.0.0)

---

## Motivation

Identification of genomic elements underlying bacterial phenotypes such as resistance to antibiotics, virulence and fermentation is a fundamental step for a better understanding of the molecular mechanisms of these traits, and may be able to inform new clinical, industrial and agricultural interventions.

BacGWASim is designed to simulate whole-genomes and phenotypes of large bacterial populations focusing on unique characteristics of bacteria such as strong genome-wide linkage disequilibrium and population stratification. The data produced by BacGWASim provides a mean to 1) evaluate the power and limitations of existing bacterial Machine learning and Genome-Wide assoctation study methods and 2) develop and benchmark novel tools for bacterial genotype to phenotype mapping tools.

## Citations

Original BacGWASim implementation paper: `Saber, Morteza M., Shapiro, B Jesse` [Benchmarking bacterial genome-wide association study methods using simulated genomes and phenotypes.](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000337#tab2) `Microbial Genomics https://doi.org/10.1099/mgen.0.000337`

## Requirements

Between parenthesis the versions the script was tested against:

- `python` 3+ (3.6.12)
- `numpy` (1.15.2)
- `scipy` (1.1.0)
- `pandas` (1.1.5)
- `ete3` (3.1.2)
- `pyvcf` (0.6.8)
- `snakemake`(5.31.1)
- `bcftools` (1.10.2)
- `plink` (1.9)
- `GCTA` (1.93.2beta)
- `snp-sites` (2.5.1)
- `openjdk` (8.0.152)
- `Haploview` (4.2)
- `simbac` (commit `5015897`)

## Installation - Using conda
BacGWASim is currently available on Bioconda. To install, run the following command. It is recommended to use `mamba` instead of `conda` for a quicker installation. 
```bash
conda create -n BacGWASim bacgwasim -c conda-forge -c bioconda
conda activate BacGWASim
```

## Installation - Using the repository

BacGWASim is build based on Snakemake and can be installed locally as following:

1.  Clone workflow into working directory

```
git clone https://github.com/Morteza-M-Saber/BacGWASim.
cd BacGWASim
```

2. Install dependencies:
   The easiest way to install BacGWASim dependencies is through [`mamba`](https://github.com/mamba-org/mamba) (a modern alternative to `conda`). **It is recommended to install dependencies inside an isolated environment to avoid package conflicts**:

```bash
mamba env create --file requirements.yaml
```

3. Build the python package.

```bash
python setup.py build install
```

4. Edit config file, or use the command-line to pass arguments to BacGWASim.

```
vim configfile.yaml
```

5. Execute workflow, determine number of available cpu cores for parallelization

```
BacGWASim --cores=5
```

In case the filesystem suffer from latency, consider increasing the `--latency-wait` accordingly.

## Simulation parameters

All the simulation parameters are included in `configfile.yaml` file and can be adjusted:

```
# Output directory name
output_dir: results_BacGWASim  # Path to the output directory

# Genome simulation parameters
num_species: 200      # Number of samples in the simulated population
genome_length: 20000  # Length of the genome (.bp)
mutation_rate: 0.06   # Mutation rate
recomb_rate: 0.01     # Recombination rate
maf: 0.01             # Minor allele frequency threshold of rare alleles to be discarded
num_var: -1           # Number of simulated variants, if kept '-1', variant number will be solely a function of mutation rate
random_seed: 1487     # Random seed for reproducibility of results

# Phenotype simulation parameters
phen_type: quant                      # Type of simulated phenotype,'cc':binary case-control, 'quant': quantitative
num_causal_var: 16                    # Number of causal markers
causal_maf_min: 0.1                   # Minimum MinorAlleleFrequency of causal markers
causal_maf_max: 0.4                   # Maximum MinorAlleleFrequency of causal markers
causal_ld_max: 0.6                    # Maximum permitted r2 score between pairs of causal markers in window size of 1000 candidate causal markers meeting causal_maf_min and causal_maf_max thresholds
effect_size_odr: 2,3,4,7,10,11,15,20  # Effect sizes of causal markers (.odds ratios) (comma separated, must be a multiple of num_causal_var)
phen_replication: 10                  # Number of phenotype replication sets
heritability: 1                       # Heritability of phenotype
disease_prevalence: 0.5               # Prevalence of phenotype
case: 50                              # Use when phen_type='cc', case + control must not be bigger than num_species
control: 50                           # Use when phen_type='cc', case + control must not be bigger than num_species

# Linkage Disequilibrium plotting
plot_ld: False   # Generate the LD plot
snp_limit: 3000  # Number of SNPs randomly selected for plotting linkage map (Increasing this number will significatnly increase computation time)
ld_maf: 0.1      # Minimum MinorAlleleFrequency of markers for LD plotting (Lower this value, it is more difficult to estimate accurate r2 values between pairs of markers leading to more noisy plot)

# Runtime parameters
cores: 1         # Number of cores available for computations
latency_wait: 3  # Time to wait (in sec) after a job to ensure all files are present
```

## Outputs

BacGWASim produces the following outputs:

```
genSim/
  sims.vcf               # Multi-sample variant call file excluding rare alleles in vcf format (used for phenotype simulation)
  sims.pickle            # Multi-sample variant call file excluding rare alleles in pandas pickle format (for machine learning analysis)
  sims_no_selection.vcf  # Multi-sample variant call file including rare alleles
  genSim.fasta           # Multiple-sequence alignment of simulated genomes in fasta format
  phylogeny.nwk          # phylogenetic tree of the simulated population in newick format

phenSim/
  phenSim.phen           # Simulated phenotypes
  phenSim.par            # ID, MAF and effect size of causal variants
  phenSim.pickle         # Simulated phenotpes in pandas pickle format
  simVis.png             # Distribution of phenotypes and casual variants in the population

ld/
  ld_plot.png            # Genome-wide linkage disequilibrium (LD) plot
  ld_stackplot.png       # Plot illustrating the relationship between R2 and distance of SNPs
  sims_subset.vcf        # Subset of the sims_no_selection.vcf file used for the plots

```

## Usage for Machine-learning models

BacGWASim outputs the simulated genotype matrix and phenotypes in pickle format that could be used to train machine-learning models as below:

```
import pandas as pd
X=pd.read_pickle(sims.pickle)
Y=pd.read_pickle(phenSim.pickle)

```

## Usage for GWAS tools

BacGWASim also outputs the simulated genotype in vcf format and phenotypes in gcta-compliant format that could be used by various GWAS tools.

## Examples

BacGWASim graphs the simulated genotype-phenotypes as color-coded phylogenetic tree and causal variants presence-absence as heat-map.

1. Simulation of a population of size 40 with 50 segregating sites as causal variants.
   ![alt text](https://github.com/Morteza-M-Saber/BacGWASim/blob/master/Img/mytree40_50.png)

2. Simulation of a population of size 100 with 50 segregating sites as causal variants.
   ![alt text](https://github.com/Morteza-M-Saber/BacGWASim/blob/master/Img/mytree100_50.png)

3. Simulation of a population of size 500 with 10 segregating sites as causal variants.
   ![alt text](https://github.com/Morteza-M-Saber/BacGWASim/blob/master/Img/mytree500_10.png)

# Verification

With implemented evolutionary models and parameters for simulation of phylogenetic tree, genome content and phenotype, BacGWASim seeks to accurately simulate the important characteristics of bacterial genomes and populations. The main feature that distinguishes bacterial populations from eukaryotes is genome-wide linkage disequilibrium (LD).

Genome-wide Linkage disequilibrium is one of the main features of bacterial population due to their clonal expansion. However, the genome-wide LD can be disrupted by homologous recombination which occur in bacteria in various degrees depending on the species. LD plots show that BacGWASim correctly captures the genome-wide LD while no recombination is simulated with decreasing LD levels in proportion with rate of homologous recombination simulated

![alt text](https://github.com/Morteza-M-Saber/BacGWASim/blob/master/Img/LDRangeComparison.png)

The distribution of linkage, measured as r2 scores binned into three categories, is shown as a function of distance in the genome, in units of kilobase pairs, within genomes of simulated (left panels) and real (right panels) bacteria.
