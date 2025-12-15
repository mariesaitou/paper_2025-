Milk Transcriptome Analysis

Bulk RNA-seq, single-cell reference integration, WGCNA, evolutionary conservation, and GWAS enrichment

Overview

This repository contains scripts for the integrative analysis of mammary gland–related gene expression using bulk RNA-seq (GTEx), single-cell reference datasets, co-expression network analysis (WGCNA), evolutionary conservation across vertebrates, and GWAS enrichment (MAGMA).

The workflow is modular and organized so that each script performs a single, clearly defined analytical task.
All analyses are conducted in R and are intended to be reproducible using the provided environment specification.

Repository structure
.
├── environment.yml
├── milk_filtering.R
├── milk_DEG_bulk.R
├── milk_PCA.R
├── milk_WGCNA.R
├── milk_compevo.R
├── milk_single_cell.R
├── milk_MAGMA.R
├── milk_upset.R
└── README.md

File descriptions
environment.yml

Conda environment specification listing required R version and packages.
This file ensures computational reproducibility across systems.

milk_filtering.R

Sample and gene filtering for bulk RNA-seq data (GTEx)

Loads GTEx expression and sample metadata

Performs quality control filtering on samples

Applies expression-based gene filtering

Outputs a filtered DGEList object (y_filtered.rds) used downstream

milk_DEG_bulk.R

Differential expression analysis (bulk RNA-seq)

Uses edgeR + limma-voom

Performs:

Mammary vs non-glandular comparisons

Glandular vs non-glandular comparisons

Sex-biased expression analysis within mammary tissue

Outputs DEG tables used by downstream analyses

milk_PCA.R

Global transcriptome structure and tissue relationships

Normalizes expression data (TMM, logCPM)

Selects highly variable genes

Performs PCA (including IRLBA-based PCA for large matrices)

Generates PCA plots for tissue-level structure and quality control

milk_WGCNA.R

Weighted Gene Co-expression Network Analysis

Adjusts mammary expression data for covariates (sex, age, center, RNA quality, cell-type fractions)

Constructs WGCNA networks for mammary tissue

Computes module eigengenes and module membership

Compares mammary modules to other tissues using Jaccard similarity

Performs permutation tests for module overlap significance

Outputs per-tissue WGCNA results and Jaccard similarity tables

milk_compevo.R

Comparative evolutionary analysis

Maps mammary-associated genes to orthologs across vertebrate species

Computes ortholog presence/absence matrices

Estimates conservation rates and one-to-one orthology

Produces summary tables for evolutionary retention analyses

Focuses on a defined set of representative species

milk_single_cell.R

Single-cell reference integration

Loads single-cell RNA-seq reference datasets (H5AD)

Restricts analysis to bulk-derived DEG sets

Computes:

Mean expression per cell type (linear scale)

Fraction of cells expressing each gene

Uses on-disk / delayed computation for large matrices

Outputs cell-type–resolved reference tables for interpretation

milk_MAGMA.R

GWAS enrichment analysis

Integrates MAGMA gene-level association statistics

Links mammary gene categories to GWAS Atlas traits

Performs enrichment testing using Fisher’s exact tests

Applies multiple testing correction

Outputs trait–gene category enrichment tables

milk_upset.R

Visualization of cross-tissue module overlap (UpSet plots)

Uses WGCNA module assignments and Jaccard similarity results

Constructs gene set intersections between mammary modules and other glandular tissues

Restricts intersections to mammary genes for interpretability

Produces publication-ready UpSet plots (PDF/PNG)

Recommended execution order

milk_filtering.R

milk_DEG_bulk.R

milk_PCA.R

milk_WGCNA.R

milk_compevo.R

milk_single_cell.R

milk_MAGMA.R

milk_upset.R

Not all scripts are strictly dependent on all previous steps, but this order reflects the intended analytical flow.

Notes on reproducibility

All scripts assume consistent sample identifiers across input files

Internal tissue names are preserved throughout the analysis

Display labels are applied only at the visualization stage

Randomized procedures (e.g. permutations) use fixed seeds where applicable