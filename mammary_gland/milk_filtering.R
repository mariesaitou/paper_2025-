# =============================================================================
# GTEx v10: Sample selection + QC filtering + phenotype merge + covariate table
# Purpose:
#   1) Read GTEx sample attributes and subject phenotypes
#   2) Select target tissues
#   3) Visualize QC metric distributions (optional)
#   4) Apply QC thresholds (as described in the manuscript)
#   5) Merge covariates and export a clean metadata table
# =============================================================================

# ---- 0) Setup ----------------------------------------------------------------

# Load packages
library(data.table)   # fread
library(dplyr)
library(stringr)
library(tidyr)        # pivot_longer / pivot_wider
library(ggplot2)
library(readr)        # read_tsv

# ---- 1) Input file paths -----------------------------------------------------
tpm_file <- "GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz"
sample_attr_file <- "GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"
subject_attr_file <- "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt"

# ---- 2) (Optional) Gene-level filtering using tissue median TPM --------------
# Note: This block loads the gene_median_tpm file (median TPM per tissue),
# and filters genes by max TPM across glandular tissues > 10.
# This is NOT sample-level QC; it is gene-level pre-filtering.

tpm <- fread(tpm_file)

glandular_cols <- c(
  "Breast_Mammary_Tissue",
  "Minor_Salivary_Gland",
  "Skin_Not_Sun_Exposed_Suprapubic",
  "Pancreas",
  "Stomach",
  "Small_Intestine_Terminal_Ileum"
)

non_glandular_cols <- c(
  "Muscle_Skeletal",
  "Heart_Left_Ventricle",
  "Artery_Tibial",
  "Esophagus_Muscularis",
  "Adipose_Subcutaneous",
  "Adipose_Visceral_Omentum"
)

tpm_subset <- tpm %>%
  select(
    gene_id = Name,
    all_of(glandular_cols),
    all_of(non_glandular_cols)
  )

tpm_filtered <- tpm_subset %>%
  filter(apply(select(., all_of(glandular_cols)), 1, max) > 10)

# ---- 3) Read sample attributes and select target tissues ---------------------
sample_info <- read_tsv(sample_attr_file)

target_tissues <- c(
  "Breast - Mammary Tissue",
  "Minor Salivary Gland",
  "Skin - Not Sun Exposed (Suprapubic)",
  "Skin - Sun Exposed (Lower Leg)",
  "Pancreas",
  "Stomach",
  "Small Intestine - Terminal Ileum",
  "Muscle - Skeletal",
  "Heart - Left Ventricle",
  "Artery - Tibial",
  "Esophagus - Muscularis",
  "Adipose - Subcutaneous",
  "Adipose - Visceral (Omentum)"
)

samples_selected <- sample_info %>%
  filter(SMTSD %in% target_tissues)

# ---- 4) QC visualization (optional) ------------------------------------------
# You can comment out this section when running on HPC / non-interactive mode.

qc_vars <- c(
  "SMRIN", "SMEXPEFF", "SMGNSDTC", "SMRDTTL", "SMMAPRT",
  "SMMPPDUN", "SMDPMPRT", "SMBSMMRT", "SMCHMRT", "SMUVCRT",
  "SM3PBMN", "SMESTLBS"
)

qc_long <- samples_selected %>%
  select(SAMPID, all_of(qc_vars)) %>%
  pivot_longer(cols = -SAMPID, names_to = "QC_metric", values_to = "value")

ggplot(qc_long, aes(x = value)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  facet_wrap(~QC_metric, scales = "free", ncol = 3) +
  theme_minimal(base_size = 12) +
  labs(title = "Distribution of QC Metrics for GTEx Samples",
       x = "Value", y = "Count")

# ---- 5) Apply QC thresholds ---------------------------------------------------
# Note: The manuscript text says "samples missing any key covariates were excluded".
# Here, we keep NA values at the QC filtering stage (i.e., NA passes),
# then later drop NA in key covariates explicitly.

sample_qc_filtered <- samples_selected %>%
  filter(
    is.na(SMRIN)     | SMRIN >= 6.0,
    is.na(SMMAPRT)   | SMMAPRT >= 0.8,
    is.na(SMEXPEFF)  | SMEXPEFF >= 0.6,
    is.na(SMGNSDTC)  | SMGNSDTC >= 10000,
    is.na(SMRDTTL)   | SMRDTTL >= 5e7,
    is.na(SMMPPDUN)  | SMMPPDUN >= 5e7,
    is.na(SMDPMPRT)  | SMDPMPRT <= 0.6,
    is.na(SMBSMMRT)  | SMBSMMRT <= 0.01,
    is.na(SMCHMRT)   | SMCHMRT <= 0.015,
    is.na(SMUVCRT)   | SMUVCRT <= 0.1,
    is.na(SM3PBMN)   | (SM3PBMN >= 0.4 & SM3PBMN <= 0.75),
    is.na(SMESTLBS)  | SMESTLBS >= 3e7
  )

# ---- 6) Merge with subject phenotypes ----------------------------------------
subject_phenotypes <- read_tsv(subject_attr_file)

# Extract subject ID (GTEX-XXXX) from sample ID
sample_qc_filtered <- sample_qc_filtered %>%
  mutate(SUBJID = str_extract(SAMPID, "GTEX-\\w+"))

sample_qc_with_pheno <- sample_qc_filtered %>%
  left_join(subject_phenotypes, by = "SUBJID")

# ---- 7) Build covariate table ------------------------------------------------
covariates <- sample_qc_with_pheno %>%
  select(
    SAMPID, SUBJID,
    SMTS, SMTSD,              # tissue (broad + detailed)
    SEX, AGE, DTHHRDY,        # subject-level
    SMCENTER, SMNABTCH, SMGEBTCH, ANALYTE_TYPE, SMAFRZE,  # technical
    SMRIN, SMTSISCH, SMTSPAX,  # RNA / collection
    SMRDTTL, SMMAPRT, SMMPPDUN, SMEXNCRT, SMEXPEFF,
    SMRRNART, SMCHMRT, SMUVCRT, SM3PBMN
  ) %>%
  # Drop samples missing key covariates required for downstream modeling
  drop_na(SEX, AGE, SMCENTER, SMNABTCH, SMGEBTCH, ANALYTE_TYPE, SMAFRZE, SMRIN, SMMPPDUN)

# ---- 8) Quick sanity check: counts by tissue and sex --------------------------
sex_counts <- covariates %>%
  mutate(SEX = case_when(
    SEX == 1 ~ "M",
    SEX == 2 ~ "F",
    TRUE ~ NA_character_
  )) %>%
  count(SMTS, SEX) %>%
  pivot_wider(names_from = SEX, values_from = n, values_fill = 0)

print(sex_counts)

# ---- 9) Add glandular classes ------------------------------------------------
covariates <- covariates %>%
  mutate(
    gland_class = case_when(
      SMTS %in% c("Breast", "Salivary Gland", "Skin", "Pancreas", "Stomach", "Small Intestine") ~ "Glandular",
      SMTS %in% c("Adipose Tissue", "Blood", "Blood Vessel", "Muscle", "Esophagus", "Heart") ~ "Non-glandular",
      TRUE ~ NA_character_
    ),
    gland_detail = case_when(
      SMTS == "Breast" ~ "Mammary",
      SMTS %in% c("Salivary Gland", "Skin", "Pancreas", "Stomach", "Small Intestine") ~ "Other Glandular",
      TRUE ~ NA_character_
    )
  )

# ---- 10) Export ----------------------------------------------------------------
write.csv(covariates, file = "sample_metadata_milk.csv", row.names = FALSE)
