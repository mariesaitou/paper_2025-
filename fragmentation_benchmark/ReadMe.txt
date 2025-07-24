
This repository contains scripts and environment files for running a structural variant (SV) mapping, filtering, and merging workflow based on Oxford Nanopore Technologies (ONT) data, designed for use in the BIO326 EUK course.

## Overview

The pipeline includes:

1. **Mapping and variant calling** using a SLURM job script
2. **VCF processing and merging** across samples from different years
3. A **conda environment YAML file** 
4.  **Downstream Analysis and Visualization** using an R script

---

## File Descriptions

### `ONT_workflow_mapping_and_variant.slurm`

- SLURM job script for alignment and structural variant calling
- Likely includes reference genome indexing, read mapping (e.g., using `minimap2`), and SV calling (e.g., using `Sniffles` or similar)
- Adjust SLURM parameters and input paths according to your HPC environment

### `VCF_processing_and_merging.sh`

- Bash script for post-calling processing and merging of VCF files
- Steps include:
  - Decompression of renamed VCFs
  - Optional reheadering via `bcftools`
  - Merging per year (2024 and 2025) using `SURVIVOR`
  - Annotation and field reduction via `bcftools annotate`
  - Splitting of merged VCFs into `TRA` and non-`TRA` subsets
  - Extraction of INFO columns for summary

**Dependencies:**
- `SURVIVOR`
- `bcftools`
- `awk`
- Conda environment assumed to be `sv` (`.conda/envs/sv`)

### `condaEUK_environment.yml`

- YAML file describing a reproducible conda environment
- Can be recreated via:

  ```bash
  conda env create -f condaEUK_environment.yml
