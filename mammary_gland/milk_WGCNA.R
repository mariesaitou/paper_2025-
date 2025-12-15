
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
