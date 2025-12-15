#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(edgeR)
  library(irlba)
})


# ---- Inputs ----
y_rds <- "y_filtered.rds"   
exclude_blood <- TRUE       
top_n <- 2000
k <- 10
seed <- 1

# ---- Outputs ----
out_y   <- if (exclude_blood) "y_noBlood.rds" else "y_used_for_pca.rds"
out_pca <- if (exclude_blood) "pca_noBlood_top2000_k10.rds" else "pca_top2000_k10.rds"

# ---- Load ----
y <- readRDS(y_rds)

# ---- Optional: exclude blood samples ----
if (exclude_blood) {
  keep <- y$samples$SMTS != "Whole_Blood"
  y$counts  <- y$counts[, keep, drop = FALSE]
  y$samples <- y$samples[keep, , drop = FALSE]
}

# ---- Normalize + logCPM ----
y <- calcNormFactors(y, method = "TMM")
logCPM <- cpm(y, log = TRUE, prior.count = 1)

# ---- Top variable genes ----
vars <- apply(logCPM, 1, var)
top_idx <- order(vars, decreasing = TRUE)[seq_len(min(top_n, length(vars)))]
X <- t(logCPM[top_idx, , drop = FALSE])   # samples x genes

# ---- PCA (IRLBA) ----
set.seed(seed)
pca <- prcomp_irlba(X, n = k, center = TRUE, scale. = TRUE)

# ---- Save ----
saveRDS(y, out_y)
saveRDS(pca, out_pca)

message("Saved: ", out_y)
message("Saved: ", out_pca)





#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})


# ---- Inputs ----
pca_rds <- "pca_noBlood_top2000_k10.rds"
y_rds   <- "y_noBlood.rds"

# ---- Output ----
out_pdf <- "PCA_top2000_noBlood_PC1_PC2.pdf"
out_png <- "PCA_top2000_noBlood_PC1_PC2.png"

# ---- Load ----
pca <- readRDS(pca_rds)
y   <- readRDS(y_rds)

pca_df <- data.frame(pca$x[, 1:5], y$samples)

expl_var <- round(100 * (pca$sdev^2) / sum(pca$sdev^2), 1)


stopifnot("SMTS" %in% names(pca_df))

centroids12 <- aggregate(cbind(PC1, PC2) ~ SMTS, pca_df, median)

# ---- A fixed palette (keep it here only) ----
cols <- c(
  "Adipose Tissue"  = "#CCBB44",
  "Blood Vessel"    = "#999999",
  "Skin"            = "#999999",
  "Salivary Gland"  = "#CCBB44",
  "Blood"           = "#999999",
  "Muscle"          = "#66C2A5",
  "Heart"           = "#EE99AA",
  "Esophagus"       = "#88CCEE",
  "Pancreas"        = "#EE99AA",
  "Stomach"         = "#66C2A5",
  "Small Intestine" = "#88CCEE",
  "Breast"          = "#EE8866"
)

shapes <- c(
  "Adipose Tissue"  = 16,
  "Blood Vessel"    = 17,
  "Skin"            = 2,
  "Salivary Gland"  = 0,
  "Blood"           = 1,
  "Muscle"          = 15,
  "Heart"           = 18,
  "Esophagus"       = 19,
  "Pancreas"        = 19,
  "Stomach"         = 18,
  "Small Intestine" = 15,
  "Breast"          = 17
)

# ---- Plot ----
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = SMTS, shape = SMTS)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_label_repel(
    data = centroids12,
    aes(label = SMTS, fill = SMTS),
    color = "black", label.size = 0.3,
    size = 3.5, box.padding = 0.5,
    show.legend = FALSE
  ) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_shape_manual(values = shapes) +
  labs(
    x = paste0("PC1 (", expl_var[1], "%)"),
    y = paste0("PC2 (", expl_var[2], "%)"),
    title = "PCA of top 2000 most variable genes"
  ) +
  theme_bw() +
  theme(legend.position = "right")

ggsave(out_pdf, p, width = 8, height = 6)
ggsave(out_png, p, width = 8, height = 6, dpi = 300)


