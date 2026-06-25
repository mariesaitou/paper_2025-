# Fixed diagnostic SNP hybrid sampling simulation

This repository contains a minimal R simulation for evaluating hybrid sampling and classification using fixed diagnostic SNPs. The simulation is designed for a two-reference-taxon system, here labelled **Salmon** and **Trout**, where every marker is assumed to be fully diagnostic in the absence of genotyping error or missing calls.

The script simulates sampled sets of individuals under user-defined population-level and individual-level hybridisation probabilities, records one row per simulation replicate, and produces one figure showing strict classification accuracy across sampling depth for selected hybrid classes.

## Script

```text
fixed_diagnostic_snp_hybrid_sampling_simulation_strict_accuracy.R
```

Run with:

```bash
Rscript fixed_diagnostic_snp_hybrid_sampling_simulation_strict_accuracy.R
```

## R packages

The script uses:

```r
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(readr)
library(ggplot2)
```

Install them if needed:

```r
install.packages(c("dplyr", "tidyr", "tibble", "purrr", "readr", "ggplot2"))
```

## Genetic model

Genotypes are encoded as Trout-allele dosage at fixed diagnostic SNPs:

| Dosage | Genotype interpretation |
|---:|---|
| `0` | Salmon / Salmon |
| `1` | Salmon / Trout |
| `2` | Trout / Trout |

The expected ancestry index used in the output is `q_trout_hat`, calculated as the mean Trout-allele dosage divided by 2 across observed SNPs for each individual.

## Output folder

The script writes all outputs to:

```text
fixed_diagnostic_snp_hybrid_sampling_outputs/
```

## Main outputs

```text
hybrid_class_definitions.csv
hybrid_sampling_replicate_results.csv
strict_classification_accuracy_by_sampling_depth.png
```

There is no scenario-level summary table and no individual-level genotype table. The main result file is `hybrid_sampling_replicate_results.csv`, where each row is one simulated sampled set, i.e. one parameter combination × one replicate.

## User-adjustable parameters

The main settings are defined near the top of the R script.

### General output setting

| Parameter | Meaning |
|---|---|
| `out_dir` | Output directory. |
| `seed` | Random seed used for reproducibility. |

### Sampling simulation settings

| Parameter | Meaning |
|---|---|
| `n_snps_grid` | Numbers of fixed diagnostic SNPs used in the sampling simulation. |
| `n_sample_grid` | Numbers of individuals sampled per simulated sampled set. |
| `n_replicates` | Number of replicate sampled sets per parameter combination. |
| `p_pop_hybridisation_grid` | Probability that hybridisation has occurred in the sampled population. This is population-level, not locus-level. |
| `p_individual_hybrid_given_pop_grid` | Probability that a sampled individual is hybrid, conditional on population-level hybridisation having occurred. |
| `p_nonhybrid_is_trout` | Among non-hybrid individuals, probability that the individual belongs to the Trout reference class. `0.5` gives an equal Salmon/Trout mixture among non-hybrids, `0` gives only nonhybrid Salmon, and `1` gives only nonhybrid Trout. |

The unconditional individual-level hybrid probability is calculated internally as:

```r
p_individual_hybrid_unconditional =
  p_pop_hybridisation * p_individual_hybrid_given_pop
```

The expected number of hybrid individuals in a sampled set is:

```r
expected_n_hybrid =
  n_sample * p_pop_hybridisation * p_individual_hybrid_given_pop
```

The theoretical probability that a sampled set contains at least one hybrid is:

```r
prob_sample_contains_at_least_one_hybrid_theory =
  p_pop_hybridisation * (1 - (1 - p_individual_hybrid_given_pop)^n_sample)
```

### Hybrid class scheme

Hybrid classes are defined in `hybrid_scheme`:

```r
hybrid_scheme <- tribble(
  ~class_id,        ~salmon_ancestry, ~trout_ancestry, ~genotype_model, ~bc_generation, ~weight,
  "F1",             0.5000,           0.5000,          "F1",           NA_integer_,    1,
  "F2",             0.5000,           0.5000,          "F2",           NA_integer_,    1,
  "BC1_to_Salmon",  0.7500,           0.2500,          "BC_to_Salmon", 1L,             1,
  "BC2_to_Salmon",  0.8750,           0.1250,          "BC_to_Salmon", 2L,             1,
  "BC3_to_Salmon",  0.9375,           0.0625,          "BC_to_Salmon", 3L,             1,
  "BC1_to_Trout",   0.2500,           0.7500,          "BC_to_Trout",  1L,             1,
  "BC2_to_Trout",   0.1250,           0.8750,          "BC_to_Trout",  2L,             1,
  "BC3_to_Trout",   0.0625,           0.9375,          "BC_to_Trout",  3L,             1
)
```

`weight` controls the relative frequency of each hybrid class among hybrid individuals. The weights are normalised internally. If all weights are `1`, all listed hybrid classes are sampled with equal probability among hybrid individuals. The ancestry columns should normally be treated as part of the class definition rather than as routine tuning parameters.

Non-hybrid reference classes are added internally as:

```text
nonhybrid_Salmon
nonhybrid_Trout
```

The term `pure` is not used in the output.

### Data imperfection settings

| Parameter | Meaning |
|---|---|
| `genotype_no_call_rate` | Probability that an individual × SNP genotype is missing. This represents genotype call failure or no-call rate. |
| `allele_call_error_rate` | Allele-level genotyping error rate. This represents incorrect allele calls, not missing genotypes. |

For example:

```r
genotype_no_call_rate <- 0.05
allele_call_error_rate <- 0.00
```

means that 5% of individual × SNP genotype calls are set to missing, while all non-missing genotype calls are assumed to be correct.

### Hybrid call threshold

The sampled-set output includes a simple hybrid call based on `q_trout_hat`:

```r
hybrid_q_trout_lower <- 0.0625
hybrid_q_trout_upper <- 0.9375
```

Individuals are called as hybrids when:

```r
q_trout_hat > hybrid_q_trout_lower &
q_trout_hat < hybrid_q_trout_upper
```

This threshold-based call is used for the detection counts in `hybrid_sampling_replicate_results.csv`. It is separate from the strict class-level classification used for the figure.

### Strict accuracy figure settings

The single figure is controlled by:

```r
n_per_class_grid <- c(5, 10, 20, 50, 100, 200)
n_snps_for_strict_accuracy_plot <- 20
strict_accuracy_classes <- c(
  "BC2_to_Salmon",
  "BC1_to_Salmon",
  "F1",
  "BC1_to_Trout",
  "BC2_to_Trout"
)
```

`n_per_class_grid` defines the x-axis of the figure. `n_snps_for_strict_accuracy_plot` fixes the number of diagnostic SNPs used for this figure. `strict_accuracy_classes` defines which true hybrid classes appear as panels.

## Output file: `hybrid_class_definitions.csv`

This file defines the simulated classes. It is the place to look up what each class means biologically and genetically.

| Column | Meaning |
|---|---|
| `class_id` | Class label used throughout the simulation output. Examples: `nonhybrid_Salmon`, `F1`, `BC2_to_Trout`. |
| `expected_ancestry_ratio` | Human-readable expected ancestry ratio, e.g. `Salmon 0.875 : Trout 0.125`. |
| `salmon_ancestry` | Expected Salmon ancestry proportion. |
| `trout_ancestry` | Expected Trout ancestry proportion. |
| `genotype_model` | Genotype-generating model used internally. |
| `bc_generation` | Backcross generation for backcross classes. `NA` for non-hybrids, F1, and F2. |
| `true_is_hybrid` | Whether the class is treated as a true hybrid class. |
| `expected_q_trout` | Expected Trout ancestry index. This is equal to `trout_ancestry`. |
| `expected_heterozygosity` | Expected heterozygosity used in strict class-level classification. |
| `weight` | Relative sampling weight among hybrid individuals. `NA` for non-hybrid classes. |
| `probability_among_hybrids` | Normalised probability of the class among hybrid individuals. `NA` for non-hybrid classes. |

Example rows:

```text
nonhybrid_Salmon   Salmon 1 : Trout 0
nonhybrid_Trout    Salmon 0 : Trout 1
F1                 Salmon 0.5 : Trout 0.5
F2                 Salmon 0.5 : Trout 0.5
BC1_to_Trout       Salmon 0.25 : Trout 0.75
BC2_to_Salmon      Salmon 0.875 : Trout 0.125
```

## Output file: `hybrid_sampling_replicate_results.csv`

Each row represents one simulated sampled set, i.e. one parameter combination × one replicate. This is not a scenario-level summary table.

| Column | Meaning |
|---|---|
| `simulation_id` | Unique row identifier for the full simulation output. |
| `scenario_id` | Human-readable label summarising the parameter combination used for the row. The same information is also stored in the numeric parameter columns. |
| `replicate` | Replicate number within a given scenario. |
| `n_snps` | Number of fixed diagnostic SNPs used for genotype simulation and hybrid detection. |
| `n_sample` | Number of individuals sampled in this simulated sampled set. |
| `p_pop_hybridisation` | Probability that hybridisation has occurred in the sampled population. |
| `p_individual_hybrid_given_pop` | Probability that a sampled individual is a hybrid, conditional on population-level hybridisation having occurred. |
| `p_individual_hybrid_unconditional` | Unconditional individual-level hybrid probability, calculated as `p_pop_hybridisation * p_individual_hybrid_given_pop`. |
| `expected_n_hybrid` | Expected number of hybrid individuals in the sampled set. |
| `prob_sample_contains_at_least_one_hybrid_theory` | Theoretical probability that the sampled set contains at least one hybrid individual. |
| `pop_hybridisation_occurred` | Whether population-level hybridisation occurred in this replicate. |
| `true_n_hybrid` | Number of true hybrid individuals actually present in the sampled set. |
| `true_n_nonhybrid` | Number of true non-hybrid individuals actually present in the sampled set. |
| `true_any_hybrid` | Whether the sampled set contains at least one true hybrid individual. |
| `detected_hybrid_n_by_q` | Number of individuals classified as hybrids using the `q_trout_hat` threshold rule. This may include true positives and false positives. |
| `detected_true_hybrid_n_by_q` | Number of true hybrid individuals correctly detected as hybrids using the `q_trout_hat` threshold rule. |
| `false_positive_n_by_q` | Number of non-hybrid individuals incorrectly classified as hybrids. |
| `false_negative_n_by_q` | Number of true hybrid individuals not classified as hybrids. |
| `detected_any_true_hybrid_by_q` | Whether at least one true hybrid individual was detected in the sampled set. |
| `any_false_positive_by_q` | Whether at least one false positive hybrid call occurred in the sampled set. |
| `mean_q_trout_hat` | Mean estimated Trout ancestry index across all sampled individuals. This is not a hybrid-prevalence estimate; it also reflects the mixture of nonhybrid Salmon and nonhybrid Trout individuals. |
| `min_q_trout_hat` | Minimum individual-level `q_trout_hat` value in the sampled set. |
| `max_q_trout_hat` | Maximum individual-level `q_trout_hat` value in the sampled set. |
| `class_counts` | Semicolon-separated counts of true simulated classes present in the sampled set. |

Example `class_counts` field:

```text
nonhybrid_Salmon=4; nonhybrid_Trout=6
```

This means that the sampled set contained 4 nonhybrid Salmon individuals and 6 nonhybrid Trout individuals. An example with hybrids might look like:

```text
nonhybrid_Salmon=7; BC1_to_Trout=2; F2=1
```

The ancestry ratios for these classes are not repeated in `class_counts`; they are defined once in `hybrid_class_definitions.csv`.

## Figure: `strict_classification_accuracy_by_sampling_depth.png`

The figure shows strict classification accuracy as a function of the number of individuals simulated per true hybrid class. Panels are faceted by true hybrid class. Points and connecting lines show the median strict accuracy across replicate runs. Shaded ribbons show the interquartile range, defined by the 25th and 75th percentiles across replicates. The x-axis uses a square-root scale.

Strict accuracy is calculated as:

```r
strict_accuracy = mean(predicted_class == true_class)
```

The class-level prediction used in the figure is based on the distance between each individual and the expected class values for both `q_trout_hat` and heterozygosity:

```r
distance = sqrt(
  (q_trout_hat - expected_q_trout)^2 +
  (heterozygosity_hat - expected_heterozygosity)^2
)
```

The predicted class is the candidate class with the smallest distance. This is used because F1 and F2 have the same expected ancestry ratio (`Salmon 0.5 : Trout 0.5`) but different expected heterozygosity.

Suggested caption:

**Figure 1. Strict classification accuracy across sampling depth for hybrid classes.** Strict classification accuracy is shown as a function of the number of individuals simulated per class (`N per class`) across replicate simulation runs. Panels are faceted by true hybrid class. Points and connecting lines indicate the median strict accuracy across replicates. Shaded ribbons represent the interquartile range, defined by the 25th and 75th percentiles of accuracy across replicates. The x-axis is shown on a square-root scale.

## Notes

This simulation assumes fixed diagnostic SNPs between the two reference taxa. Under the default settings, genotype calls are complete and error-free. Missing genotype calls can be introduced using `genotype_no_call_rate`, and allele-level call errors can be introduced using `allele_call_error_rate`.

The sampled-set output and the strict classification figure answer different questions. `hybrid_sampling_replicate_results.csv` records how many hybrids were present and detected in each simulated sampled set under varying prevalence and sampling designs. `strict_classification_accuracy_by_sampling_depth.png` evaluates how accurately selected hybrid classes can be classified when the true class is fixed and the number of individuals per class is varied.
