
# ===== 0. Setup =====

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(WGCNA)
  library(stats)
  library(data.table)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# ===== 1. Parameters =====
y_obj      <- y_filtered            # DGEList object
tissue_var <- "SMTSD"               # tissue label
cov_rhs    <- ~ SEX + AGE + SMCENTER + SMEXPEFF
min_prop_expressed <- 0.20
logCPM_threshold    <- 0
power_fallback      <- 8
minModuleSize       <- 120
mergeCutHeight      <- 0.35

# ===== 2. Normalization and logCPM =====
y_norm <- calcNormFactors(y_obj)
logCPM_all <- cpm(y_norm, log = TRUE, prior.count = 1)

stopifnot(tissue_var %in% colnames(y_norm$samples))
tissue_cols <- split(seq_len(ncol(logCPM_all)), y_norm$samples[[tissue_var]])

# ===== 3. Per-tissue WGCNA =====
results <- list(per_tissue = list(), summary = data.frame())

for (tt in names(tissue_cols)) {
  message(sprintf("[ %s ] WGCNA...", tt))
  
  idx <- tissue_cols[[tt]]
  expr_mat <- logCPM_all[, idx, drop = FALSE]
  meta     <- y_norm$samples[idx, , drop = FALSE]
  
  cov_mat <- model.matrix(cov_rhs, data = meta)[, -1, drop = FALSE]
  
  # Covariate regression
  expr_adj <- tryCatch(
    removeBatchEffect(expr_mat, covariates = cov_mat),
    error = function(e) expr_mat
  )
  
  # Filter genes
  min_n <- ceiling(ncol(expr_adj) * min_prop_expressed)
  keep  <- rowSums(expr_adj > logCPM_threshold) >= min_n
  expr_adj <- expr_adj[keep, , drop = FALSE]
  
  # Transpose for WGCNA (samples × genes)
  datExpr <- t(expr_adj)
  gsg <- goodSamplesGenes(datExpr, verbose = 0)
  if(!gsg$allOK) datExpr <- datExpr[, gsg$goodGenes, drop = FALSE]
  
  # Soft-thresholding power
  powers <- 1:20
  sft <- suppressWarnings(pickSoftThreshold(datExpr, powerVector = powers, verbose = 0))
  power_est <- ifelse(is.na(sft$powerEstimate), power_fallback, sft$powerEstimate)
  
  # Network construction
  net <- blockwiseModules(
    datExpr,
    power = power_est,
    TOMType = "unsigned",
    networkType = "signed",
    corType = "bicor",
    minModuleSize = minModuleSize,
    deepSplit = 1,
    mergeCutHeight = mergeCutHeight,
    reassignThreshold = 0,
    pamRespectsDendro = TRUE,
    numericLabels = TRUE,
    verbose = 0
  )
  
  # Module labeling and eigengenes
  tab <- sort(table(net$colors), decreasing = TRUE)
  label_map <- setNames(paste0("Module", seq_along(tab)), names(tab))
  module_labels <- label_map[as.character(net$colors)]
  names(module_labels) <- colnames(datExpr)
  
  MEs <- orderMEs(moduleEigengenes(datExpr, colors = net$colors)$eigengenes)
  kME <- signedKME(datExpr, MEs)
  
  results$per_tissue[[tt]] <- list(
    datExpr = datExpr,
    module_labels = module_labels,
    MEs = MEs,
    kME = kME,
    power = power_est,
    n_genes = ncol(datExpr),
    n_samples = nrow(datExpr),
    module_sizes = tab
  )
  
  results$summary <- rbind(
    results$summary,
    data.frame(
      tissue = tt,
      n_samples = nrow(datExpr),
      n_genes = ncol(datExpr),
      n_modules = length(unique(module_labels)),
      power = power_est,
      stringsAsFactors = FALSE
    )
  )
}

saveRDS(results, "WGCNA_per_tissue_results.rds")

# ===== 4. Cross-tissue overlap (Mammary vs others) =====
results <- readRDS("WGCNA_per_tissue_results.rds")

mammary <- "Breast - Mammary Tissue"
ancestors <- setdiff(names(results$per_tissue), mammary)

module_overlap_map <- list()
mix_index <- list()

genes_mam <- names(results$per_tissue[[mammary]]$module_labels)
mods_mam  <- unique(results$per_tissue[[mammary]]$module_labels)

for (mm in mods_mam) {
  gset_m <- genes_mam[results$per_tissue[[mammary]]$module_labels == mm]
  score_tbl <- data.frame()
  
  for (tt in ancestors) {
    genes_t <- names(results$per_tissue[[tt]]$module_labels)
    mods_t  <- unique(results$per_tissue[[tt]]$module_labels)
    common  <- intersect(gset_m, genes_t)
    if (length(common) < 10) next
    
    cols_t <- results$per_tissue[[tt]]$module_labels[common]
    for (mt in mods_t) {
      a <- sum(cols_t == mt)
      b <- length(common) - a
      c <- sum(results$per_tissue[[tt]]$module_labels[genes_t %in% genes_mam] == mt) - a
      d <- length(genes_t) - a - b - c
      if (min(a,b,c,d) < 0) next
      p <- fisher.test(matrix(c(a,b,c,d), nrow=2))$p.value
      score_tbl <- rbind(score_tbl, data.frame(tissue=tt, module=mt, pval=p))
    }
  }
  
  if (nrow(score_tbl) > 0) {
    score_tbl$FDR <- p.adjust(score_tbl$pval, method="BH")
    score_tbl$log10FDR <- -log10(score_tbl$FDR)
    best_by_tissue <- aggregate(log10FDR ~ tissue, data=score_tbl, FUN=max)
    w <- pmax(best_by_tissue$log10FDR, 0)
    w <- if (sum(w) > 0) w / sum(w) else w
    mix_index[[mm]] <- setNames(w, best_by_tissue$tissue)
    module_overlap_map[[mm]] <- head(score_tbl[order(-score_tbl$log10FDR), ], 5)
  }
}

# ===== 5. Output summaries =====
label_tbl <- data.frame(
  SMTSD = c(
    "Pancreas","Minor Salivary Gland","Skin - Not Sun Exposed (Suprapubic)",
    "Small Intestine - Terminal Ileum","Stomach","Adipose - Subcutaneous",
    "Adipose - Visceral (Omentum)","Artery - Tibial","Esophagus - Muscularis",
    "Heart - Left Ventricle","Muscle - Skeletal","Breast - Mammary Tissue"
  ),
  gland_class = c(rep("Glandular",5), rep("Non-glandular",6), "Glandular"),
  stringsAsFactors = FALSE
)

mix_df <- do.call(rbind, lapply(names(mix_index), function(mm){
  if (length(mix_index[[mm]]) == 0) return(NULL)
  data.frame(mammary_module = mm, tissue = names(mix_index[[mm]]),
             weight = as.numeric(mix_index[[mm]]), stringsAsFactors = FALSE)
}))
mix_df <- merge(mix_df, label_tbl, by.x="tissue", by.y="SMTSD", all.x=TRUE)
write.table(mix_df, "mix_index.tsv", sep="\t", quote=FALSE, row.names=FALSE)

overlap_df <- do.call(rbind, lapply(names(module_overlap_map), function(mm){
  df <- module_overlap_map[[mm]]
  if (is.null(df) || nrow(df)==0) return(NULL)
  data.frame(mammary_module=mm, df, stringsAsFactors=FALSE)
}))
overlap_df <- merge(overlap_df, label_tbl, by.x="tissue", by.y="SMTSD", all.x=TRUE)
write.table(overlap_df, "module_overlap_map.tsv", sep="\t", quote=FALSE, row.names=FALSE)

message("Analysis complete. Results saved to: WGCNA_per_tissue_results.rds, mix_index.tsv, module_overlap_map.tsv")





######## Module preservation analysis


library(WGCNA)
library(dplyr)
library(tibble)
library(tidyr)
library(pheatmap)
library(viridis)
library(gprofiler2)
library(purrr)
library(ggplot2)
library(stringr)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()



#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(WGCNA)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

args <- commandArgs(trailingOnly = TRUE)
in_dir  <- if (length(args) >= 1) args[[1]] else "."
out_dir <- if (length(args) >= 2) args[[2]] else "."

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Input dir: ", normalizePath(in_dir, mustWork = FALSE))
message("Output dir: ", normalizePath(out_dir, mustWork = FALSE))





# =========================
rds_file <- file.path(in_dir, "WGCNA_per_tissue_results_adipoReg_breastOnly.rds")

if (!file.exists(rds_file)) {
  stop("Input RDS file not found: ", rds_file)
}

message("Reading RDS: ", rds_file)
obj <- readRDS(rds_file)

# =========================
# Tissues
# =========================
tissues_use <- c(
  "Breast - Mammary Tissue",
  "Minor Salivary Gland",
  "Pancreas",
  "Stomach",
  "Skin - Not Sun Exposed (Suprapubic)",
  "Small Intestine - Terminal Ileum"
)

missing_tissues <- setdiff(tissues_use, names(obj$per_tissue))
if (length(missing_tissues) > 0) {
  stop("These tissues are missing in obj$per_tissue: ",
       paste(missing_tissues, collapse = ", "))
}

message("Tissues used:")
message(paste0("  - ", tissues_use, collapse = "\n"))

# =========================
# Common genes
# =========================
common_genes <- Reduce(
  intersect,
  lapply(obj$per_tissue[tissues_use], function(x) colnames(x$datExpr))
)

if (length(common_genes) < 1000) {
  stop("Too few common genes found: ", length(common_genes))
}

message("Number of common genes: ", length(common_genes))

# =========================
# multiExpr 
# =========================
multiExpr <- lapply(obj$per_tissue[tissues_use], function(x) {
  dat <- x$datExpr[, common_genes, drop = FALSE]
  list(data = as.data.frame(dat))
})
names(multiExpr) <- tissues_use

# =========================
# Breast reference modules
# =========================
breast_name <- "Breast - Mammary Tissue"
breast_labels <- obj$per_tissue[[breast_name]]$module_labels[common_genes]

if (any(is.na(breast_labels))) {
  stop("NA found in breast_labels after gene matching.")
}

message("Breast module sizes on common genes:")
print(table(breast_labels))

multiColor <- lapply(tissues_use, function(x) breast_labels)
names(multiColor) <- tissues_use


n_perm <- as.integer(Sys.getenv("NPERM", unset = "200"))
if (is.na(n_perm) || n_perm < 1) {
  stop("Invalid NPERM value.")
}
message("nPermutations = ", n_perm)

# =========================
# module preservation
# =========================
set.seed(1)

message("Starting modulePreservation at: ", Sys.time())

mp <- modulePreservation(
  multiData = multiExpr,
  multiColor = multiColor,
  referenceNetworks = which(names(multiExpr) == breast_name),
  networkType = "signed",
  nPermutations = n_perm,
  verbose = 3
)

message("Finished modulePreservation at: ", Sys.time())

saveRDS(
  mp,
  file = file.path(out_dir, paste0("module_preservation_breast_reference_", n_perm, "perm.rds"))
)


extract_preservation <- function(mp_obj, ref_name) {
  ref_idx <- which(names(mp_obj$preservation$Z) == ref_name)
  
  out <- lapply(seq_along(mp_obj$preservation$Z[[ref_idx]]), function(i) {
    test_name <- names(mp_obj$preservation$Z[[ref_idx]])[i]
    
    zdf <- mp_obj$preservation$Z[[ref_idx]][[i]] %>%
      as.data.frame() %>%
      rownames_to_column("module")
    
    odf <- mp_obj$preservation$observed[[ref_idx]][[i]] %>%
      as.data.frame() %>%
      rownames_to_column("module")
    
    z_col <- grep("Zsummary", colnames(zdf), value = TRUE)[1]
    mr_col <- grep("medianRank", colnames(odf), value = TRUE)[1]
    size_col <- grep("moduleSize|size", colnames(odf), value = TRUE)[1]
    
    if (is.na(z_col)) stop("Zsummary column not found for ", test_name)
    if (is.na(mr_col)) stop("medianRank column not found for ", test_name)
    if (is.na(size_col)) stop("moduleSize column not found for ", test_name)
    
    zdf %>%
      select(module, Zsummary = all_of(z_col)) %>%
      left_join(
        odf %>% select(module, medianRank = all_of(mr_col), moduleSize = all_of(size_col)),
        by = "module"
      ) %>%
      mutate(test_tissue = test_name)
  }) %>%
    bind_rows()
  
  out
}

pres_df <- extract_preservation(mp, ref_name = breast_name) %>%
  filter(!module %in% c("gold", "grey")) %>%
  mutate(
    preservation_class = case_when(
      Zsummary > 10 ~ "Strong",
      Zsummary > 2  ~ "Moderate",
      TRUE          ~ "Weak_or_none"
    )
  ) %>%
  arrange(test_tissue, desc(Zsummary))

write.csv(
  pres_df,
  file = file.path(out_dir, paste0("module_preservation_summary_", n_perm, "perm.csv")),
  row.names = FALSE
)

# =========================
p_all <- ggplot(pres_df, aes(x = test_tissue, y = Zsummary)) +
  geom_hline(yintercept = c(2, 10), linetype = 2) +
  geom_point(size = 2) +
  facet_wrap(~ module, scales = "free_x") +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = "Module preservation Zsummary"
  )

ggsave(
  filename = file.path(out_dir, paste0("module_preservation_all_modules_", n_perm, "perm.pdf")),
  plot = p_all,
  width = 14,
  height = 10
)

if (nrow(module2_df) > 0) {
  p_m2 <- ggplot(module2_df, aes(x = reorder(test_tissue, Zsummary), y = Zsummary)) +
    geom_hline(yintercept = c(2, 10), linetype = 2) +
    geom_point(size = 3) +
    coord_flip() +
    theme_bw() +
    labs(
      x = NULL,
      y = "Module preservation Zsummary",
      title = "Preservation of Breast Module2"
    )
  
  ggsave(
    filename = file.path(out_dir, paste0("module2_preservation_", n_perm, "perm.pdf")),
    plot = p_m2,
    width = 7,
    height = 4.5
  )
}


# =========================
sink(file.path(out_dir, paste0("module_preservation_log_", n_perm, "perm.txt")))
cat("Run completed at:", as.character(Sys.time()), "\n\n")
cat("Input RDS:", rds_file, "\n")
cat("Common genes:", length(common_genes), "\n")
cat("nPermutations:", n_perm, "\n\n")
cat("Module counts in breast on common genes:\n")
print(table(breast_labels))
cat("\nSummary of preservation classes:\n")
print(table(pres_df$preservation_class, useNA = "ifany"))
cat("\nTop rows of pres_df:\n")
print(utils::head(pres_df, 20))
sink()

message("All done.")



# Downstream visualization and GO

# --------------------------------------------------
# 1. Load module preservation results
# --------------------------------------------------
mp <- readRDS("module_preservation_breast_reference_200perm.rds")

# --------------------------------------------------
# 2. Extract module preservation statistics
# --------------------------------------------------
extract_preservation <- function(mp_obj) {
  ref_name <- names(mp_obj$preservation$Z)[1]
  test_names <- names(mp_obj$preservation$Z[[ref_name]])
  
  bind_rows(lapply(test_names, function(test_name) {
    x <- mp_obj$preservation$Z[[ref_name]][[test_name]]
    
    if (is.logical(x) && length(x) == 1 && is.na(x)) {
      return(NULL)
    }
    
    tibble(
      module = rownames(x),
      Zsummary = x$Zsummary.pres,
      test_tissue = test_name
    )
  }))
}

pres_df <- extract_preservation(mp) %>%
  filter(!module %in% c("gold", "grey")) %>%
  mutate(
    test_tissue = recode(
      sub("^inColumnsAlsoPresentIn\\.", "", test_tissue),
      "Minor Salivary Gland" = "Salivary",
      "Pancreas" = "Pancreas",
      "Skin - Not Sun Exposed (Suprapubic)" = "Skin",
      "Stomach" = "Stomach",
      "Small Intestine - Terminal Ileum" = "Intestine"
    )
  )

# --------------------------------------------------
# 3. Prepare heatmap matrix
# --------------------------------------------------
mat <- pres_df %>%
  pivot_wider(names_from = test_tissue, values_from = Zsummary) %>%
  column_to_rownames("module") %>%
  as.matrix()

num_mat <- ifelse(mat >= 10, sprintf("%.1f", mat), "")
rownames(num_mat) <- rownames(mat)
colnames(num_mat) <- colnames(mat)

# --------------------------------------------------
# 4. Draw Zsummary heatmap
# --------------------------------------------------
pheatmap(
  mat,
  color = viridis(100),
  display_numbers = num_mat,
  number_color = "black",
  cluster_rows = TRUE,
  cluster_cols = TRUE
)

# --------------------------------------------------
# 5. Extract gene lists for selected modules
# --------------------------------------------------
breast_labels <- obj$per_tissue[["Breast - Mammary Tissue"]]$module_labels

modules_to_extract <- c("Module11", "Module7", "Module14", "Module6")

module_gene_lists <- lapply(modules_to_extract, function(m) {
  names(breast_labels[breast_labels == m])
})
names(module_gene_lists) <- modules_to_extract

for (m in modules_to_extract) {
  write.table(
    module_gene_lists[[m]],
    file = paste0(m, "_genes.txt"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

# --------------------------------------------------
# 6. Run g:Profiler enrichment
# --------------------------------------------------
gprof_results <- lapply(names(module_gene_lists), function(m) {
  gost(
    query = module_gene_lists[[m]],
    organism = "hsapiens",
    correction_method = "fdr",
    sources = c("GO:BP", "GO:CC", "GO:MF", "REAC")
  )
})
names(gprof_results) <- names(module_gene_lists)

# --------------------------------------------------
# 7. Save enrichment tables
# --------------------------------------------------
for (m in names(gprof_results)) {
  res <- gprof_results[[m]]$result
  
  if (!is.null(res) && nrow(res) > 0) {
    res_out <- res %>%
      select(
        query, significant, p_value, term_size, query_size,
        intersection_size, precision, recall,
        term_id, source, term_name,
        effective_domain_size, source_order
      )
    
    write.csv(res_out, paste0(m, "_gprofiler.csv"), row.names = FALSE)
  }
}

# --------------------------------------------------
# 8. Prepare enrichment plot data
# --------------------------------------------------
gprof_df <- map_dfr(names(gprof_results), function(m) {
  res <- gprof_results[[m]]$result
  if (is.null(res) || nrow(res) == 0) return(NULL)
  
  res %>%
    mutate(module = m)
})

plot_df <- gprof_df %>%
  arrange(module, p_value) %>%
  group_by(module) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(
    module = factor(module, levels = c("Module11", "Module7", "Module14", "Module6")),
    term_name = str_wrap(term_name, width = 24),
    source = factor(source, levels = c("GO:BP", "GO:CC", "GO:MF", "REAC"))
  )

shape_map <- c(
  "GO:BP" = 16,
  "GO:CC" = 15,
  "GO:MF" = 17,
  "REAC"  = 18
)

# --------------------------------------------------
# 9. Draw enrichment bubble plot
# --------------------------------------------------
p <- ggplot(
  plot_df,
  aes(
    x = -log10(p_value),
    y = reorder(term_name, -log10(p_value))
  )
) +
  geom_point(
    aes(size = intersection_size, shape = source),
    color = "#0a2a92"
  ) +
  scale_shape_manual(values = shape_map) +
  facet_wrap(
    ~ module,
    scales = "free_y",
    ncol = 2
  ) +
  theme_bw() +
  labs(
    x = expression(-log[10](p)),
    y = NULL,
    size = "Intersection size",
    shape = "Source"
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "plain", size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.3),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

p