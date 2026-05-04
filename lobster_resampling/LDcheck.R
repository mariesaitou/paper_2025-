#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
})

obj <- readRDS("/Users/saitoumarie/Library/CloudStorage/Dropbox/Norway/grant/PhD2022/Erik/lob_erik/lobsim_redo/exeter_lobster_1591ind_79snps_52pop.rds")

setwd("/Users/saitoumarie/Library/CloudStorage/Dropbox/Norway/grant/PhD2022/Erik/lob_erik/lobsim_redo")

library(adegenet)

# Convert genind object to a genotype matrix (0/1/2)
geno <- tab(obj, NA.method = "mean")

# Add individual metadata
df <- data.frame(
  ID = indNames(obj),
  pop = pop(obj),
  geno
)

library(adegenet)

geno <- tab(obj, freq = FALSE, NA.method = "asis", NA.as.zero = FALSE)

locfac <- locFac(obj)
keep <- !duplicated(locfac)
geno012 <- geno[, keep]

df <- data.frame(
  ID = indNames(obj),
  pop = pop(obj),
  geno012,
  check.names = FALSE
)

# Transpose to SNP x individual format and export
X <- t(df[, -(1:2)])

# Add row and column names
rownames(X) <- colnames(df)[-(1:2)]

# Add individual ID and population to column names
colnames(X) <- paste(df$ID, df$pop, sep = "|")

# Export
write.table(X, file = "lobster_genotypes.tsv",
            sep = "\t", quote = FALSE, col.names = NA)

# Load data
df <- read.table("lobster_genotypes.tsv", header = TRUE, sep = "\t", check.names = FALSE)

# Number of individuals per population
tab <- as.data.frame(table(pop(obj)))
colnames(tab) <- c("population", "n_individuals")

# Export
write.table(tab, file = "population_counts.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

library(adegenet)

# Create a 0/1/2 genotype matrix
geno <- tab(obj, freq = FALSE, NA.method = "asis", NA.as.zero = FALSE)
locfac <- locFac(obj)
keep <- !duplicated(locfac)
geno012 <- geno[, keep, drop = FALSE]

popv <- as.character(pop(obj))

# Define population groups
idx_am <- popv == "Americanus"
idx_eu <- !(popv %in% c("Americanus", "AmerCook", "EuroCook", "HybridX"))

# Extract genotype matrices
X_am <- geno012[idx_am, , drop = FALSE]
X_eu <- geno012[idx_eu, , drop = FALSE]

# Calculate LD as r2
ld_am <- cor(X_am, use = "pairwise.complete.obs")^2
ld_eu <- cor(X_eu, use = "pairwise.complete.obs")^2

# Extract the upper triangle to create vectors for the LD distribution
ld_am_vec <- ld_am[upper.tri(ld_am)]
ld_eu_vec <- ld_eu[upper.tri(ld_eu)]

# Create tables for SNP pairs
make_ld_table <- function(ldmat) {
  ut <- upper.tri(ldmat)
  data.frame(
    SNP1 = rownames(ldmat)[row(ldmat)[ut]],
    SNP2 = colnames(ldmat)[col(ldmat)[ut]],
    r2   = ldmat[ut],
    row.names = NULL,
    check.names = FALSE
  )
}

ld_am_tab <- make_ld_table(ld_am)
ld_eu_tab <- make_ld_table(ld_eu)

# Export
write.table(ld_am_tab, "LD_Americanus_pairs.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ld_eu_tab, "LD_Europe_pairs.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

write.table(data.frame(r2 = ld_am_vec), "LD_Americanus_r2_values.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(r2 = ld_eu_vec), "LD_Europe_r2_values.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Extract Americanus samples
X_am <- geno012[idx_am, , drop = FALSE]

# Exclude loci with fewer than three non-missing genotypes or zero variance
keep_am <- apply(X_am, 2, function(x) {
  x <- x[!is.na(x)]
  length(x) >= 3 && var(x) > 0
})

X_am2 <- X_am[, keep_am, drop = FALSE]

# Calculate LD as r2
ld_am <- cor(X_am2, use = "pairwise.complete.obs")^2

# Create SNP-pair table
make_ld_table <- function(ldmat) {
  ut <- upper.tri(ldmat)
  data.frame(
    SNP1 = rownames(ldmat)[row(ldmat)[ut]],
    SNP2 = colnames(ldmat)[col(ldmat)[ut]],
    r2   = ldmat[ut],
    row.names = NULL,
    check.names = FALSE
  )
}

ld_am_tab <- make_ld_table(ld_am)
write.table(ld_am_tab, "LD_Americanus_pairs.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

X_eu <- geno012[idx_eu, , drop = FALSE]

keep_eu <- apply(X_eu, 2, function(x) {
  x <- x[!is.na(x)]
  length(x) >= 3 && var(x) > 0
})

X_eu2 <- X_eu[, keep_eu, drop = FALSE]
ld_eu <- cor(X_eu2, use = "pairwise.complete.obs")^2
ld_eu_tab <- make_ld_table(ld_eu)
write.table(ld_eu_tab, "LD_Europe_pairs.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Extract r2 values
ld_eu_vec <- ld_eu[upper.tri(ld_eu)]

# Plot the LD distribution
hist(ld_eu_vec, breaks = 50, main = "LD (r²) distribution - Europe",
     xlab = "r²")

# Summarize LD values
eu_stats <- c(
  n_snps = ncol(X_eu2),
  n_pairs = length(ld_eu_vec),
  mean = mean(ld_eu_vec, na.rm = TRUE),
  median = median(ld_eu_vec, na.rm = TRUE),
  max = max(ld_eu_vec, na.rm = TRUE),
  n_r2_gt_02 = sum(ld_eu_vec > 0.2, na.rm = TRUE)
)

sum(ld_eu_vec < 2, na.rm = TRUE)

sum(ld_eu_vec > 0.2, na.rm = TRUE)

# Extract only SNP pairs with r2 > 0.2 in the European reference
ld_eu_high <- ld_eu_tab %>%
  filter(r2 > 0.2)

write.table(
  ld_eu_high,
  "LD_Europe_r2_gt_0.2.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

ld_eu_high

# Extract only SNP pairs with r2 > 0.2 in the American reference
ld_am_high <- ld_am_tab %>%
  filter(r2 > 0.2)

write.table(
  ld_am_high,
  "LD_America_r2_gt_0.2.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

library(dplyr)
library(tidyr)
library(ggplot2)

# Combine high-LD SNPs from the European and American references
high_ld_snps_all <- sort(unique(c(
  ld_eu_high$SNP1, ld_eu_high$SNP2,
  ld_am_high$SNP1, ld_am_high$SNP2
)))

# Align matrices using the same SNP order
ld_eu_full <- ld_eu[high_ld_snps_all, high_ld_snps_all]
ld_am_full <- ld_am[high_ld_snps_all, high_ld_snps_all]

# Treat NA values as zero, if required
ld_eu_full[is.na(ld_eu_full)] <- 0
ld_am_full[is.na(ld_am_full)] <- 0

# Cluster SNPs based on the mean LD across the European and American references
ld_mean <- (ld_eu_full + ld_am_full) / 2
diag(ld_mean) <- 1

hc <- hclust(as.dist(1 - ld_mean))
snp_order <- hc$labels[hc$order]

# Convert to long format
ld_eu_df <- as.data.frame(as.table(ld_eu_full)) %>%
  rename(SNP1 = Var1, SNP2 = Var2, r2 = Freq) %>%
  mutate(reference = "European reference")

ld_am_df <- as.data.frame(as.table(ld_am_full)) %>%
  rename(SNP1 = Var1, SNP2 = Var2, r2 = Freq) %>%
  mutate(reference = "American reference")

ld_facet_df <- bind_rows(ld_eu_df, ld_am_df) %>%
  mutate(
    SNP1 = factor(SNP1, levels = snp_order),
    SNP2 = factor(SNP2, levels = snp_order)
  )

p <- ggplot(ld_facet_df, aes(x = SNP1, y = SNP2, fill = r2)) +
  geom_tile() +
  facet_wrap(~ reference, nrow = 1) +
  scale_fill_gradient(
    low = "white",
    high = "black",
    limits = c(0, 1),
    name = expression(r^2)
  ) +
  coord_equal() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title = element_blank(),
    strip.background = element_rect(fill = "white"),
    panel.grid = element_blank()
  )

p

ggsave("LD_high_r2_heatmap_Europe_Americanus_facet.pdf", p, width = 11, height = 5)
ggsave("LD_high_r2_heatmap_Europe_Americanus_facet.png", p, width = 11, height = 5, dpi = 300)

snps_to_drop_all <- sort(unique(c(
  ld_eu_high$SNP1, ld_eu_high$SNP2,
  ld_am_high$SNP1, ld_am_high$SNP2
)))

snps_to_drop_all
length(snps_to_drop_all)

write.table(
  data.frame(SNP = snps_to_drop_all),
  "LD_pruned_all_SNPs_in_high_LD_pairs.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)