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

## Citation

If you use this code or results, please cite:

Ellis et al. (2020), Scientific Reports 10:7781.
