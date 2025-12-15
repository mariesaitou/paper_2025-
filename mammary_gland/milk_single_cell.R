# =============================================================================
# Single-cell reference aggregation (cell type mean + detected fraction)
# - Input: h5ad (CELLxGENE), DEG.csv (column "gene")
# - Output (per dataset):
#     <prefix>_mean_linear_by_celltype.tsv.gz
#     <prefix>_detected_fraction_by_celltype.tsv.gz
# Notes:
#   - Works on-disk with DelayedArray/HDF5Array to avoid OOM
#   - Gene IDs are matched after stripping Ensembl version suffix
# =============================================================================

options(stringsAsFactors = FALSE)
options(error = function() {
  cat("ERROR at ", format(Sys.time()), "\n", file = stderr())
  traceback(max.lines = 20)
  quit(status = 1)
})

msg <- function(...) message(format(Sys.time()), " :: ", paste0(..., collapse = " "))

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(zellkonverter)   # readH5AD
  library(Matrix)
  library(DelayedArray)
  library(HDF5Array)
  library(dplyr)
  library(readr)
  library(edgeR)
  library(pbapply)
})


strip_ver <- function(x) sub("\\.\\d+$", "", x)

# ---------------------------------------------------------------------
# Helper 1: aggregate gene x cell matrix into gene x type sums, on-disk
#   M: DelayedMatrix (genes x cells)
#   G: sparse model matrix (cells x types)
# ---------------------------------------------------------------------
aggregate_cols_to_hdf5 <- function(M, G, filepath, name,
                                   msglabel = "aggregate",
                                   chunk_cols = 5L,
                                   row_chunk_hint = 4096L) {
  
  stopifnot(ncol(M) == nrow(G))
  
  k  <- ncol(G)   # types
  nr <- nrow(M)   # genes
  msg(msglabel, ": types = ", k, ", genes = ", nr)
  
  # Pre-allocate an on-disk HDF5 array (genes x types)
  S0 <- DelayedArray(matrix(0, nrow = nr, ncol = k,
                            dimnames = list(rownames(M), colnames(G))))
  S  <- writeHDF5Array(S0, filepath = filepath, name = name,
                       chunkdim = c(min(row_chunk_hint, nr), 1L))
  
  idx <- split(seq_len(k), ceiling(seq_len(k) / chunk_cols))
  done <- 0L
  pb <- txtProgressBar(min = 0, max = k, style = 3)
  
  for (ii in seq_along(idx)) {
    jj <- idx[[ii]]
    Gj <- G[, jj, drop = FALSE]     # cells x |jj|
    Sj <- M %*% Gj                  # genes x |jj| (Delayed)
    S[, jj] <- as.matrix(Sj)        # materialize only this small block
    done <- done + length(jj)
    setTxtProgressBar(pb, done)
    if ((ii %% 2L) == 0L) flush.console()
  }
  
  close(pb)
  S
}

# ---------------------------------------------------------------------
# Helper 2: write DelayedMatrix to TSV.GZ in row blocks
# ---------------------------------------------------------------------
write_delayedmatrix_tsv_gz <- function(H, out_path, gene_colname = "gene",
                                       row_block = 5000L) {
  
  con <- gzfile(out_path, open = "wt")
  on.exit(close(con), add = TRUE)
  
  # header
  hdr <- c(gene_colname, colnames(H))
  writeLines(paste(hdr, collapse = "\t"), con = con)
  
  nr <- nrow(H)
  rb <- split(seq_len(nr), ceiling(seq_len(nr) / row_block))
  
  for (ii in seq_along(rb)) {
    rr <- rb[[ii]]
    block <- as.matrix(H[rr, , drop = FALSE])  # small row block
    rn_clean <- strip_ver(rownames(H)[rr])
    out <- cbind(rn_clean, block)
    write.table(out, file = con, sep = "\t", row.names = FALSE,
                col.names = FALSE, quote = FALSE)
    if ((ii %% 10L) == 0L) flush(con)
  }
  
  invisible(TRUE)
}

# ---------------------------------------------------------------------
# Main runner for one h5ad
# ---------------------------------------------------------------------
run_one <- function(h5ad_path, deg_path,
                    celltype_col = "cell_type",
                    out_prefix = "ref") {
  
  msg("read H5AD: ", h5ad_path)
  sce <- readH5AD(h5ad_path, use_hdf5 = TRUE, reader = "R")
  
  stopifnot(celltype_col %in% colnames(colData(sce)))
  grp <- droplevels(as.factor(colData(sce)[[celltype_col]]))
  names(grp) <- colnames(sce)
  
  # Design matrix: cells x types
  G <- Matrix::sparse.model.matrix(~ grp - 1)
  colnames(G) <- levels(grp)
  n_by_type <- colSums(G)
  
  # DEG list
  msg("read DEG: ", deg_path)
  deg <- readr::read_csv(deg_path, show_col_types = FALSE)
  stopifnot("gene" %in% names(deg))
  deg_genes <- unique(strip_ver(deg$gene))
  
  # Subset genes early (speed + memory)
  rn_clean <- strip_ver(rownames(sce))
  genes_keep <- rownames(sce)[rn_clean %in% deg_genes]
  if (length(genes_keep) == 0) stop("No DEG genes found in this h5ad.")
  msg("subset genes: ", length(genes_keep), " / ", nrow(sce))
  sce <- sce[genes_keep, , drop = FALSE]
  
  # Choose assay
  assays_avail <- assayNames(sce)
  assay_to_use <- if ("counts" %in% assays_avail) "counts" else
    if ("X" %in% assays_avail) "X" else
      if ("logcounts" %in% assays_avail) "logcounts" else assays_avail[1]
  msg("using assay: ", assay_to_use)
  mat <- assay(sce, assay_to_use)  # genes x cells (Delayed)
  
  # ------------------------------------------------------------
  # 1) Detected fraction per cell type
  # ------------------------------------------------------------
  msg("compute detected fraction (on-disk)")
  det_tmp <- paste0(out_prefix, "_tmp_det.h5")
  if (file.exists(det_tmp)) unlink(det_tmp)
  
  det_hits <- (DelayedArray(mat) > 0)
  S_det <- aggregate_cols_to_hdf5(det_hits, G, filepath = det_tmp, name = "det_sum",
                                  msglabel = "sum detected")
  det_by_type <- sweep(S_det, 2, n_by_type, "/")
  
  # ------------------------------------------------------------
  # 2) Mean expression on linear scale per cell type
  #    - If counts available: CPM then mean
  #    - Else: expm1(log1p) then mean (assumes log1p scale, not z-scored)
  # ------------------------------------------------------------
  msg("compute mean expression on linear scale (on-disk)")
  mean_tmp <- paste0(out_prefix, "_tmp_mean.h5")
  if (file.exists(mean_tmp)) unlink(mean_tmp)
  
  if ("counts" %in% assays_avail) {
    cts <- assay(sce, "counts")
    # CPM via library-size scaling as DelayedArray
    libsize <- DelayedArray::colSums2(cts)
    scale <- 1e6 / pmax(libsize, 1)
    cpm_D <- sweep(DelayedArray(cts), 2, scale, `*`)
    S_mean <- aggregate_cols_to_hdf5(cpm_D, G, filepath = mean_tmp, name = "mean_sum",
                                     msglabel = "sum CPM")
    mean_by_type <- sweep(S_mean, 2, n_by_type, "/")
  } else {
    # WARNING: only valid if mat is log1p-like (not scaled / centered)
    lin_mat <- DelayedArray(mat)
    lin_mat <- DelayedArray::apply(lin_mat, 1, function(x) pmax(expm1(x), 0))
    S_mean <- aggregate_cols_to_hdf5(lin_mat, G, filepath = mean_tmp, name = "mean_sum",
                                     msglabel = "sum linear")
    mean_by_type <- sweep(S_mean, 2, n_by_type, "/")
  }
  
  # ------------------------------------------------------------
  # Outputs
  # ------------------------------------------------------------
  mean_out <- paste0(out_prefix, "_mean_linear_by_celltype.tsv.gz")
  det_out  <- paste0(out_prefix, "_detected_fraction_by_celltype.tsv.gz")
  
  msg("write: ", mean_out)
  write_delayedmatrix_tsv_gz(mean_by_type, mean_out, gene_colname = "gene")
  
  msg("write: ", det_out)
  write_delayedmatrix_tsv_gz(det_by_type, det_out, gene_colname = "gene")
  
  msg("done: ", out_prefix)
  invisible(list(mean = mean_out, det = det_out))
}

# =============================================================================
# Run (two datasets)
# =============================================================================

deg_path <- "DEG.csv"

run_one(
  h5ad_path    = "/Downloads/2c0072ed-656e-4b08-a2dc-f8537945ff54.h5ad",
  deg_path     = deg_path,
  celltype_col = "cell_type",
  out_prefix   = "mammary_ref"
)

run_one(
  h5ad_path    = "/Downloads/1c5acd5a-a0ef-4b89-a238-2ef153ad1129.h5ad",
  deg_path     = deg_path,
  celltype_col = "cell_type",
  out_prefix   = "lipid_ref"
)
