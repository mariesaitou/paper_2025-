# lobster_resampling

Empirical resampling workflow for validating a reduced SNP panel for hybrid detection across generations in European lobster.

This repository contains R scripts used to evaluate how a 79-SNP panel performs when assigning pure species, F1 hybrids, and later-generation backcrosses between European lobster (*Homarus gammarus*) and American lobster (*H. americanus*). The workflow generates synthetic hybrid and backcross genotypes by resampling from observed multilocus genotypes, rather than using forward-time demographic simulations, and evaluates hybrid assignment with `snapclust` from `adegenet`.

## Repository contents

- `01_empirical_resampling_snapclust.R`  
  Main empirical resampling and `snapclust` classification workflow. This script generates resampled parental, F1, and backcross genotype classes from empirical genotype data, runs hybrid assignment, and summarises classification performance across sampling depths.

- `02_parental_reference_downsampling.R`  
  Sensitivity analysis for parental reference size. The European parental reference is repeatedly downsampled to match the smaller American parental reference, and the resampling/classification workflow is rerun to evaluate whether unequal parental reference sizes affect assignment results.

- `03_balanced_vs_full_reference_comparison.R`  
  Direct comparison between the full-reference and balanced-reference analyses. This script compares classification outcomes from the full parental reference dataset and the downsampled balanced-reference dataset.

- `04_ld_descriptive_summary.R`  
  Descriptive linkage disequilibrium analysis of the SNP panel. This script summarises pairwise LD among loci and is used as a diagnostic description of marker independence. LD-pruned reclassification is not part of the current primary workflow.

## Analysis overview

The main workflow evaluates SNP-panel performance under empirically generated hybrid classes. The analysis focuses on:

1. Generating resampled multilocus genotypes for parental, F1, and backcross classes.
2. Running `snapclust`-based hybrid assignment.
3. Evaluating classification accuracy under strict and nearest-class criteria.
4. Summarising posterior probabilities as an indicator of assignment reliability.
5. Testing whether unequal parental reference sample sizes materially affect classification results.
6. Describing linkage disequilibrium among SNPs as a marker-panel diagnostic.

The main biological expectation is that pure species and F1 hybrids should be classified with high reliability, whereas uncertainty is expected to increase among adjacent later-generation backcross classes.

## Requirements

The scripts require R and the following R packages:

- `adegenet`
- `dplyr`
- `tidyr`
- `tibble`
- `ggplot2`
- `stringr`
- `purrr`
- `readr`

The analysis was developed for use with an empirical SNP genotype dataset stored as an `adegenet` `genind` object.

## Input data

The analysis uses empirical genotype data for 1,591 individuals genotyped at 79 SNP loci. The SNP panel and empirical genotypes originate from:

Ellis, C. D., Jenkins, T. L., Svanberg, L., Eriksson, S. P. & Stevens, J. R.  
Crossing the pond: genetic assignment detects lobster hybridisation.  
*Scientific Reports* 10, 7781 (2020).

The empirical genotype dataset is not redistributed in this repository. Users must ensure appropriate permission for data access and reuse.

Input file paths are defined near the top of each script and may need to be edited for local use.

## How to run

From the repository directory, run the scripts in the following order:

```bash
Rscript 01_empirical_resampling_snapclust.R
Rscript 02_parental_reference_downsampling.R
Rscript 03_balanced_vs_full_reference_comparison.R
Rscript 04_ld_descriptive_summary.R


The main empirical resampling script should be run first because later comparison scripts depend on its outputs.

Outputs

The scripts generate summary tables and diagnostic outputs describing:

strict classification accuracy
nearest-class classification accuracy
posterior probability distributions
classification performance across sampling depths
effects of balancing parental reference sample sizes
pairwise LD among SNP markers

Output file names and directories are defined inside the individual scripts.

Reproducibility notes

All resampling analyses use explicit random seeds. Because snapclust cluster labels are arbitrary, downstream classification summaries account for possible inversion of parental cluster labels.

Scope of the current workflow

The current repository focuses on empirical resampling, parental-reference sensitivity analysis, full-versus-balanced reference comparison, and descriptive LD assessment. Earlier exploratory scripts for LD pruning or older downsampling implementations are not part of the current documented workflow.
EOF
