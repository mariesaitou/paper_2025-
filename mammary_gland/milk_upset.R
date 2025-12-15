#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(UpSetR)
})


# ============================================================
# Required inputs
# ============================================================
# results:
#   A list containing WGCNA results per tissue, where
#   results[["per_tissue"]][[tissue]]$module_labels
#   is a named vector mapping genes to module IDs (e.g. "Module1").
#
# tissue_annotations (optional):
#   A data.frame with columns:
#     - tissue: tissue name (must match results$per_tissue keys)
#     - gland_class: e.g. "Glandular" / "Non-glandular"
#
# Example (if not already loaded):
# results <- readRDS("WGCNA_per_tissue_results.rds")
# tissue_annotations <- read_tsv("tissue_annotations.tsv", show_col_types = FALSE)

jaccard_path   <- "jaccard_best_per_tissue_adipoReg.csv"
target_modules <- c("Module1")   # Mammary modules to visualise
threshold      <- NA_real_       # Set numeric value (e.g. 0.2) to filter by Jaccard, or NA to disable

# ============================================================
# Display-only label shortening
# (Internal tissue keys must NEVER be modified)
# ============================================================
label_map <- c(
  "Breast - Mammary Tissue"            = "Mammary",
  "Small Intestine - Terminal Ileum"   = "Intestine",
  "Minor Salivary Gland"               = "Salivary G",
  "Skin - Not Sun Exposed (Suprapubic)"= "Skin",
  "Skin - Sun Exposed (Lower Leg)"     = "Skin (Leg)"
)

short_label <- function(x) {
  ifelse(x %in% names(label_map), label_map[[x]], x)
}

# Output directory for figures
out_dir <- "upset_plots"
dir.create(out_dir, showWarnings = FALSE)

# ============================================================
# 1. Identify Mammary tissue key used internally in results
# ============================================================
pt_names <- names(results[["per_tissue"]])
mammary_tissue <- pt_names[grepl("Breast", pt_names) & grepl("Mammary", pt_names)][1]

if (is.na(mammary_tissue) || mammary_tissue == "") {
  stop("Mammary tissue not found in results$per_tissue.")
}

# ============================================================
# 2. Load Jaccard similarity table (best partner per tissue)
# ============================================================
jaccard_df <- read_tsv(jaccard_path, show_col_types = FALSE)

required_cols <- c("tissue", "module_mammary", "module_other", "jaccard")
if (!all(required_cols %in% colnames(jaccard_df))) {
  stop("Jaccard table is missing required columns.")
}

# ============================================================
# 3. Determine which partner tissues to include
#    (Restrict to glandular tissues if annotation is available)
# ============================================================
if (exists("tissue_annotations")) {
  stopifnot(all(c("tissue", "gland_class") %in% colnames(tissue_annotations)))
  glandular_tissues <- tissue_annotations %>%
    filter(gland_class == "Glandular") %>%
    pull(tissue) %>%
    intersect(names(results[["per_tissue"]]))
} else {
  # Fallback: allow all tissues present in results
  glandular_tissues <- names(results[["per_tissue"]])
}

# ============================================================
# 4. Construct gene sets for a single Mammary module
#    Optionally restrict partner sets to Mammary genes
# ============================================================
build_gene_sets <- function(module_mam, restrict_to_mammary = TRUE) {
  
  # Genes belonging to the Mammary module
  mam_labels <- results[["per_tissue"]][[mammary_tissue]]$module_labels
  mam_genes  <- names(mam_labels[mam_labels == module_mam])
  if (length(mam_genes) == 0) return(NULL)
  
  # Candidate partner modules from the Jaccard table
  sub_df <- jaccard_df %>%
    filter(module_mammary == module_mam,
           tissue %in% glandular_tissues)
  
  if (!is.na(threshold)) {
    sub_df <- sub_df %>% filter(jaccard >= threshold)
  }
  
  if (nrow(sub_df) == 0) return(NULL)
  
  gene_sets <- list()
  
  # Mammary reference set (display label only)
  mam_set_name <- paste0(short_label(mammary_tissue), "_", module_mam)
  gene_sets[[mam_set_name]] <- mam_genes
  
  # Partner tissue modules
  for (i in seq_len(nrow(sub_df))) {
    tissue_i <- sub_df$tissue[i]      # internal tissue key
    module_i <- sub_df$module_other[i]
    
    if (!tissue_i %in% names(results[["per_tissue"]])) next
    
    labels_i <- results[["per_tissue"]][[tissue_i]]$module_labels
    genes_i  <- names(labels_i[labels_i == module_i])
    if (length(genes_i) == 0) next
    
    # Restrict partner genes to Mammary genes if requested
    if (restrict_to_mammary) {
      genes_i <- intersect(genes_i, mam_genes)
    }
    if (length(genes_i) == 0) next
    
    set_name <- paste0(short_label(tissue_i), "_", module_i)
    set_name <- gsub("[^A-Za-z0-9._-]+", "_", set_name)
    gene_sets[[set_name]] <- genes_i
  }
  
  # Skip if Mammary-only
  if (length(gene_sets) < 2) return(NULL)
  
  gene_sets
}

# ============================================================
# 5. Generate UpSet plots for each target Mammary module
# ============================================================
for (m in target_modules) {
  
  gene_sets <- build_gene_sets(m, restrict_to_mammary = TRUE)
  
  if (is.null(gene_sets)) {
    message("Skipping ", m, ": no comparable gene sets.")
    next
  }
  
  mtx <- UpSetR::fromList(gene_sets)
  
  # Force Mammary set to appear at the top
  mam_col <- paste0(short_label(mammary_tissue), "_", m)
  if (mam_col %in% colnames(mtx)) {
    desired_order <- c(mam_col, setdiff(colnames(mtx), mam_col))
  } else {
    desired_order <- colnames(mtx)
  }
  upset_sets <- rev(desired_order)  # UpSetR plotting order quirk
  
  out_pdf <- file.path(out_dir, paste0("UpSet_", m, "_restrictedToMammary.pdf"))
  out_png <- file.path(out_dir, paste0("UpSet_", m, "_restrictedToMammary.png"))
  
  pdf(out_pdf, width = 9, height = 5)
  UpSetR::upset(
    mtx,
    nintersects = 20,
    order.by = "freq",
    sets = upset_sets,
    keep.order = TRUE,
    mainbar.y.label = paste0("Intersection size (within Mammary ", m, ")"),
    sets.x.label   = "Set size (Mammary genes only)"
  )
  dev.off()
  
  png(out_png, width = 1800, height = 1000, res = 200)
  UpSetR::upset(
    mtx,
    nintersects = 20,
    order.by = "freq",
    sets = upset_sets,
    keep.order = TRUE,
    mainbar.y.label = paste0("Intersection size (within Mammary ", m, ")"),
    sets.x.label   = "Set size (Mammary genes only)"
  )
  dev.off()
  
  message("Saved UpSet plots for ", m)
}
