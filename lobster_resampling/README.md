# lobster_resampling

Resampling-based validation of a 79-SNP panel for hybrid detection across generations in European lobster, using empirical genotype data and snapclust (adegenet).

This repository contains scripts and summary outputs for:
- Generating resampled hybrid classes (parentals, F1, and backcross generations)
- Performing hybrid assignment using snapclust
- Summarising classification accuracy across sampling depths

## Contents

- lobster_resampling.R  
  Main analysis script implementing empirical resampling and snapclust-based hybrid assignment.

- snapclust_resampling_summary_table.tsv  
  Summary table reporting strict classification accuracy across hybrid classes and sample sizes.

- snapclust_resampling_runmeta.rds  
  Metadata for resampling runs, including replicate structure and random seeds.

## Contents

- lobster_resampling.R  
  Main analysis script implementing empirical resampling and snapclust-based hybrid assignment.

- snapclust_resampling_summary_table.tsv  
  Summary table reporting strict classification accuracy across hybrid classes and sample sizes.

- snapclust_resampling_runmeta.rds  
  Metadata for resampling runs, including replicate structure and random seeds.

- LDcheck.R  
  Script for evaluating linkage disequilibrium among SNPs within the parental reference groups. This script calculates pairwise r² values separately for the European and American references, exports LD summary tables, identifies SNP pairs with r² > 0.2, and generates the LD heatmap used for supplementary visualization.

- LDprun.R  
  Script for repeating the empirical resampling and snapclust classification analysis after LD pruning. SNPs involved in high-LD pairs identified by `LDcheck.R` are removed, and the full classification workflow is rerun on the LD-pruned marker panel.

- Downsampling.R  
  Script for the parental reference size sensitivity analysis. The European parental reference is randomly subsampled to match the American parental reference size (n = 38), after which the empirical resampling and snapclust classification workflow is repeated.

## Requirements

- R (tested with R 4.3.3)
- R packages:
  - adegenet
  - dplyr
  - tidyr
  - tibble
  - ggplot2
  - stringr

## Input data

The analysis relies on an empirical SNP genotype dataset stored as an adegenet genind object.
The dataset comprises 1,591 individuals genotyped at 79 SNP loci.

The SNP panel and empirical genotypes originate from:

Ellis, C. D., Jenkins, T. L., Svanberg, L., Eriksson, S. P. & Stevens, J. R.  
Crossing the pond: genetic assignment detects lobster hybridisation.  
Scientific Reports 10, 7781 (2020).

Empirical genotype data are not redistributed here. Users must ensure appropriate permission for data access and reuse.

## How to run

From the repository directory:

Rscript lobster_resampling.R

Input file paths are defined at the top of the script and may need to be adjusted locally.

## Reproducibility notes

- All resampling runs are controlled by explicit random seeds.
- snapclust assigns parental clusters as A and B without biological labels; downstream accuracy evaluation accounts for possible A/B label inversion.

