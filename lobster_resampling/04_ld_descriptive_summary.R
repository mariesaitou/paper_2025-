#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(adegenet)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

## ---------------------------------------------------------------------------
## Descriptive LD analysis for the lobster 79-SNP panel
##
## Purpose:
##   1. Extract unphased genotype dosages from the genind object.
##   2. Define parental reference groups consistently with the manuscript.
##   3. Calculate genotype-dosage r^2 as cor(x_i, x_j)^2.
##   4. Output all 3,081 pairwise comparisons per reference, including
##      non-estimable pairs caused by fixed or insufficiently variable loci.
##   5. Produce:
##        - Figure S1 candidate: threshold-focused high-LD heatmap
##        - Figure S2 candidate: full pairwise r^2 distribution histogram
##
## This script does NOT perform LD-pruned reclassification and does NOT write
## LD_pruned_all_SNPs_in_high_LD_pairs.tsv.
## ---------------------------------------------------------------------------

ROOT <- "/Users/saitoumarie/Library/CloudStorage/Dropbox/Norway/grant/PhD2022/Erik/lob_erik/lobsim_redo"
RDS_FILE <- file.path(ROOT, "exeter_lobster_1591ind_79snps_52pop.rds")

setwd(ROOT)

obj <- readRDS(RDS_FILE)

## ---------------------------------------------------------------------------
## Population definitions
## ---------------------------------------------------------------------------

AMERICAN_POPS <- c("Americanus", "AmerCook")

MED_POPS <- c(
  "Adr", "Ale", "Chi", "Csa", "Ion", "Laz",
  "Sar", "Sky", "Spo", "The", "Tor"
)

HYBRID_POP <- "HybridX"

popv <- as.character(pop(obj))

idx_am <- popv %in% AMERICAN_POPS

idx_med <- popv %in% MED_POPS

idx_hybrid <- popv %in% HYBRID_POP

## Atlantic European H. gammarus parental reference:
## all non-Mediterranean H. gammarus samples, excluding H. americanus,
## Mediterranean populations, and HybridX.
##
## EuroCook is intentionally NOT excluded here, because it is non-Mediterranean
## H. gammarus and should be part of the Atlantic European reference if present.
idx_eu <- !(idx_am | idx_med | idx_hybrid)

reference_counts <- tibble(
  reference_group = c(
    "American reference",
    "Atlantic European reference",
    "Excluded Mediterranean H. gammarus",
    "HybridX",
    "Other / unclassified by these rules"
  ),
  n_individuals = c(
    sum(idx_am),
    sum(idx_eu),
    sum(idx_med),
    sum(idx_hybrid),
    sum(!(idx_am | idx_eu | idx_med | idx_hybrid))
  )
)

write_tsv(reference_counts, "LD_reference_group_counts.tsv")

population_counts <- as.data.frame(table(popv)) |>
  rename(population = popv, n_individuals = Freq) |>
  arrange(population)

write_tsv(population_counts, "population_counts.tsv")

print(reference_counts)

## ---------------------------------------------------------------------------
## Extract one genotype-dosage column per SNP locus
## ---------------------------------------------------------------------------

geno_all_alleles <- tab(
  obj,
  freq = FALSE,
  NA.method = "asis",
  NA.as.zero = FALSE
)

locfac <- locFac(obj)

## For a biallelic SNP represented by two allele columns, either allele column
## gives the same r^2 after squaring the Pearson correlation. We keep the first
## allele column per locus to match the existing SNP naming as closely as possible.
keep_one_allele_per_locus <- !duplicated(locfac)

geno012 <- geno_all_alleles[, keep_one_allele_per_locus, drop = FALSE]

if (ncol(geno012) != 79) {
  warning("Expected 79 SNP columns after keeping one allele per locus, but found ", ncol(geno012), ".")
}

## Write SNP x individual matrix for inspection / archiving
geno_df <- data.frame(
  ID = indNames(obj),
  pop = popv,
  geno012,
  check.names = FALSE
)

X_export <- t(geno_df[, -(1:2), drop = FALSE])
rownames(X_export) <- colnames(geno_df)[-(1:2)]
colnames(X_export) <- paste(geno_df$ID, geno_df$pop, sep = "|")

write.table(
  X_export,
  file = "lobster_genotypes.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

## ---------------------------------------------------------------------------
## Pairwise genotype-dosage r^2 calculation
## ---------------------------------------------------------------------------

calc_pairwise_r2 <- function(X_ref, reference_name) {
  loci <- colnames(X_ref)
  pair_index <- combn(seq_along(loci), 2)
  
  out <- lapply(seq_len(ncol(pair_index)), function(k) {
    i <- pair_index[1, k]
    j <- pair_index[2, k]
    
    locus1 <- loci[i]
    locus2 <- loci[j]
    
    x <- X_ref[, i]
    y <- X_ref[, j]
    
    complete <- complete.cases(x, y)
    xk <- x[complete]
    yk <- y[complete]
    
    n_complete <- length(xk)
    n_unique_locus1 <- length(unique(xk))
    n_unique_locus2 <- length(unique(yk))
    
    if (n_complete < 3) {
      r <- NA_real_
      r2 <- NA_real_
      status <- "not_estimable_fewer_than_3_complete_genotypes"
    } else if (var(xk) == 0 || var(yk) == 0) {
      r <- NA_real_
      r2 <- NA_real_
      status <- "not_estimable_zero_variance"
    } else {
      r <- suppressWarnings(cor(xk, yk, use = "complete.obs", method = "pearson"))
      r2 <- r^2
      status <- "estimable"
    }
    
    tibble(
      reference = reference_name,
      locus1 = locus1,
      locus2 = locus2,
      n_complete = n_complete,
      n_unique_locus1 = n_unique_locus1,
      n_unique_locus2 = n_unique_locus2,
      r = r,
      r2 = r2,
      status = status
    )
  })
  
  bind_rows(out)
}

X_am <- geno012[idx_am, , drop = FALSE]
X_eu <- geno012[idx_eu, , drop = FALSE]

ld_all <- bind_rows(
  calc_pairwise_r2(X_am, "American reference"),
  calc_pairwise_r2(X_eu, "European reference")
)

write_tsv(ld_all, "LD_pairwise_genotype_dosage_r2_all_pairs.tsv")

ld_all |>
  filter(reference == "American reference") |>
  write_tsv("LD_American_reference_pairs.tsv")

ld_all |>
  filter(reference == "European reference") |>
  write_tsv("LD_European_reference_pairs.tsv")

ld_all |>
  filter(reference == "American reference") |>
  select(r2) |>
  write_tsv("LD_American_reference_r2_values.tsv")

ld_all |>
  filter(reference == "European reference") |>
  select(r2) |>
  write_tsv("LD_European_reference_r2_values.tsv")

status_counts <- ld_all |>
  count(reference, status, name = "n_pairs") |>
  group_by(reference) |>
  mutate(prop_pairs = n_pairs / sum(n_pairs)) |>
  ungroup()

write_tsv(status_counts, "LD_pairwise_r2_status_counts.tsv")

ld_high <- ld_all |>
  filter(status == "estimable", r2 > 0.2) |>
  arrange(reference, desc(r2), locus1, locus2)

write_tsv(ld_high, "LD_high_r2_pairs.tsv")

ld_high |>
  filter(reference == "American reference") |>
  write_tsv("LD_American_reference_r2_gt_0.2.tsv")

ld_high |>
  filter(reference == "European reference") |>
  write_tsv("LD_European_reference_r2_gt_0.2.tsv")

ld_high_loci_by_ref <- ld_high |>
  select(reference, locus1, locus2) |>
  pivot_longer(
    cols = c(locus1, locus2),
    names_to = "locus_position",
    values_to = "locus"
  ) |>
  distinct(reference, locus) |>
  arrange(reference, locus)

write_tsv(ld_high_loci_by_ref, "LD_high_r2_loci_by_reference.tsv")

high_loci_union <- sort(unique(c(ld_high$locus1, ld_high$locus2)))

write_tsv(
  tibble(locus = high_loci_union),
  "LD_high_r2_loci_union.tsv"
)

summary_counts <- tibble(
  statistic = c(
    "n_total_loci",
    "n_pairwise_comparisons_per_reference",
    "n_high_r2_pairs_american",
    "n_high_r2_pairs_european",
    "n_high_r2_loci_american",
    "n_high_r2_loci_european",
    "n_high_r2_loci_union"
  ),
  value = c(
    ncol(geno012),
    choose(ncol(geno012), 2),
    sum(ld_high$reference == "American reference"),
    sum(ld_high$reference == "European reference"),
    ld_high_loci_by_ref |> filter(reference == "American reference") |> pull(locus) |> unique() |> length(),
    ld_high_loci_by_ref |> filter(reference == "European reference") |> pull(locus) |> unique() |> length(),
    length(high_loci_union)
  )
)

write_tsv(summary_counts, "LD_high_r2_summary_counts.tsv")

print(status_counts)
print(summary_counts)

## ---------------------------------------------------------------------------
## Figure S2 candidate:
## Full distribution of pairwise genotype-dosage r^2 values
## ---------------------------------------------------------------------------

ld_hist <- ld_all |>
  mutate(
    reference = factor(reference, levels = c("American reference", "European reference"))
  ) |>
  filter(status == "estimable")

hist_labels <- ld_all |>
  group_by(reference) |>
  summarise(
    n_all_pairs = n(),
    n_estimable = sum(status == "estimable"),
    n_not_estimable = sum(status != "estimable"),
    n_r2_gt_0.2 = sum(r2 > 0.2, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    reference = factor(reference, levels = c("American reference", "European reference")),
    label = paste0(
      "All pairs: ", n_all_pairs,
      "\nEstimable: ", n_estimable,
      "\nNot estimable: ", n_not_estimable,
      "\nr² > 0.2: ", n_r2_gt_0.2
    )
  )

p_hist <- ggplot(ld_hist, aes(x = r2)) +
  geom_histogram(binwidth = 0.02, boundary = 0, closed = "left") +
  geom_text(
    data = hist_labels,
    aes(x = 0.55, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1.2,
    size = 3
  ) +
  facet_wrap(~ reference, ncol = 1, scales = "free_y") +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    x = expression(paste("Genotype-dosage ", r^2)),
    y = "Number of estimable SNP pairs"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(color = "black")
  )

print(p_hist)

ggsave(
  "LD_pairwise_r2_distribution_all_pairs.pdf",
  p_hist,
  width = 7,
  height = 5
)

ggsave(
  "LD_pairwise_r2_distribution_all_pairs.png",
  p_hist,
  width = 7,
  height = 5,
  dpi = 300
)

## ---------------------------------------------------------------------------
## Figure S1 candidate:
## Heatmap of loci involved in at least one r^2 > 0.2 pair
## ---------------------------------------------------------------------------

make_ld_matrix <- function(ld_table, loci) {
  mat <- matrix(
    NA_real_,
    nrow = length(loci),
    ncol = length(loci),
    dimnames = list(loci, loci)
  )
  
  diag(mat) <- 1
  
  for (k in seq_len(nrow(ld_table))) {
    l1 <- ld_table$locus1[k]
    l2 <- ld_table$locus2[k]
    value <- ld_table$r2[k]
    
    if (l1 %in% loci && l2 %in% loci) {
      mat[l1, l2] <- value
      mat[l2, l1] <- value
    }
  }
  
  mat
}

if (length(high_loci_union) >= 2) {
  ld_am_mat <- make_ld_matrix(
    ld_all |> filter(reference == "American reference"),
    high_loci_union
  )
  
  ld_eu_mat <- make_ld_matrix(
    ld_all |> filter(reference == "European reference"),
    high_loci_union
  )
  
  ## Ordering is used only for plotting. NA values are treated as 0 for ordering,
  ## but remain NA in the plotted matrix.
  ld_am_for_order <- ld_am_mat
  ld_eu_for_order <- ld_eu_mat
  
  ld_am_for_order[is.na(ld_am_for_order)] <- 0
  ld_eu_for_order[is.na(ld_eu_for_order)] <- 0
  
  ld_mean_for_order <- (ld_am_for_order + ld_eu_for_order) / 2
  diag(ld_mean_for_order) <- 1
  
  hc <- hclust(as.dist(1 - ld_mean_for_order))
  snp_order <- hc$labels[hc$order]
  
  ld_am_df <- as.data.frame(as.table(ld_am_mat)) |>
    rename(locus_y = Var1, locus_x = Var2, r2 = Freq) |>
    mutate(reference = "American reference")
  
  ld_eu_df <- as.data.frame(as.table(ld_eu_mat)) |>
    rename(locus_y = Var1, locus_x = Var2, r2 = Freq) |>
    mutate(reference = "European reference")
  
  ld_heat_df <- bind_rows(ld_am_df, ld_eu_df) |>
    mutate(
      reference = factor(reference, levels = c("American reference", "European reference")),
      locus_x = factor(locus_x, levels = snp_order),
      locus_y = factor(locus_y, levels = rev(snp_order))
    )
  
  p_heat <- ggplot(ld_heat_df, aes(x = locus_x, y = locus_y, fill = r2)) +
    geom_tile() +
    facet_wrap(~ reference, nrow = 1) +
    scale_fill_gradient(
      low = "white",
      high = "blue",
      limits = c(0, 1),
      na.value = "grey70",
      name = expression(r^2)
    ) +
    coord_equal() +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
      axis.text.y = element_text(size = 7),
      strip.background = element_rect(fill = "white"),
      panel.grid = element_blank()
    )
  
  print(p_heat)
  
  ggsave(
    "LD_high_r2_heatmap_Europe_Americanus_facet.pdf",
    p_heat,
    width = 11,
    height = 5
  )
  
  ggsave(
    "LD_high_r2_heatmap_Europe_Americanus_facet.png",
    p_heat,
    width = 11,
    height = 5,
    dpi = 300
  )
} else {
  warning("Fewer than two loci were involved in r2 > 0.2 pairs; heatmap was not generated.")
}

cat("\nDone.\n")
cat("Key outputs:\n")
cat("  LD_pairwise_genotype_dosage_r2_all_pairs.tsv\n")
cat("  LD_pairwise_r2_status_counts.tsv\n")
cat("  LD_high_r2_pairs.tsv\n")
cat("  LD_high_r2_loci_by_reference.tsv\n")
cat("  LD_high_r2_loci_union.tsv\n")
cat("  LD_pairwise_r2_distribution_all_pairs.pdf/png\n")
cat("  LD_high_r2_heatmap_Europe_Americanus_facet.pdf/png\n")