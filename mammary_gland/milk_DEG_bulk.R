# =============================================================================
# GTEx v10: Differential expression (voom+limma) with a reusable function
# Analyses included:
#   A) Mammary tissue: Female vs Male (within Mammary)
#   B) All tissues:   Glandular vs Non-glandular
# =============================================================================


library(data.table)
library(dplyr)
library(stringr)
library(readr)
library(edgeR)
library(limma)

# ---- Inputs ------------------------------------------------------------------
metadata_csv <- "sample_metadata_milk.csv"

# If you already have 'counts' in memory from the previous script, skip loading.
# Otherwise, load a saved object (recommended).
# Example: saveRDS(counts, "counts_matrix.rds") once, then read it here.
# counts <- readRDS("counts_matrix.rds")

meta <- fread(metadata_csv, na.strings = c("", "NA", "NaN"))
meta[, SAMPID := trimws(SAMPID)]
meta_df <- as.data.frame(meta)
rownames(meta_df) <- meta_df$SAMPID

strip_ver <- function(x) sub("\\.\\d+$", "", x)

# ---- Utility function ---------------------------------------------------------
run_voom_de <- function(counts_mat, meta_df, group_var, covars,
                        lfc_cut = 1, fdr_cut = 1e-3,
                        coef_pattern = NULL, coef_name = NULL,
                        up_label = NULL, out_csv = NULL) {
  
  stopifnot(all(rownames(meta_df) %in% colnames(counts_mat)))
  counts_use <- counts_mat[, rownames(meta_df), drop = FALSE]
  
  y <- edgeR::DGEList(counts = counts_use, samples = meta_df)
  
  # filterByExpr uses the group variable
  keep <- edgeR::filterByExpr(y, group = y$samples[[group_var]])
  y <- y[keep, , keep.lib.sizes = FALSE]
  
  y <- edgeR::calcNormFactors(y, method = "TMM")
  
  # model formula
  fml <- as.formula(paste("~", paste(c(group_var, covars), collapse = " + ")))
  design <- model.matrix(fml, data = y$samples)
  
  v <- limma::voom(y, design, plot = FALSE)
  fit <- limma::lmFit(v, design)
  fit <- limma::eBayes(fit)
  
  # coefficient selection
  if (!is.null(coef_name)) {
    coef_use <- coef_name
  } else {
    stopifnot(!is.null(coef_pattern))
    cn <- colnames(coef(fit))
    hits <- cn[grepl(coef_pattern, cn)]
    stopifnot(length(hits) == 1)
    coef_use <- hits
  }
  
  tt <- limma::topTable(fit, coef = coef_use, number = Inf, sort.by = "P")
  
  if (is.null(up_label)) up_label <- coef_use
  
  RES <- tt %>%
    dplyr::mutate(
      Regulation = dplyr::case_when(
        adj.P.Val < fdr_cut & logFC >=  lfc_cut ~ paste0("Upregulated_in_",   up_label),
        adj.P.Val < fdr_cut & logFC <= -lfc_cut ~ paste0("Downregulated_in_", up_label),
        TRUE ~ "Not_significant"
      )
    )
  
  if (!is.null(out_csv)) {
    utils::write.csv(RES, file = out_csv, row.names = TRUE)
  }
  
  return(RES)
}

# =============================================================================
# A) Mammary tissue: Female vs Male
# =============================================================================

meta_mam <- meta_df %>% dplyr::filter(gland_detail == "Mammary")
meta_mam <- droplevels(meta_mam)

# Ensure SEX factor baseline order (Male first, Female second)
# Your metadata currently uses 1/2; this keeps the same behavior as your script.
if (!is.factor(meta_mam$SEX)) meta_mam$SEX <- factor(meta_mam$SEX)
if (all(levels(meta_mam$SEX) %in% c("1","2"))) {
  meta_mam$SEX <- factor(meta_mam$SEX, levels = c("1","2"))
  up_label_sex <- "Female"
} else if (all(levels(meta_mam$SEX) %in% c("Male","Female"))) {
  meta_mam$SEX <- factor(meta_mam$SEX, levels = c("Male","Female"))
  up_label_sex <- "Female"
} else {
  # fallback: use the last level as "up"
  up_label_sex <- tail(levels(meta_mam$SEX), 1)
}

# Covariate types (match manuscript)
meta_mam$AGE      <- factor(meta_mam$AGE)
meta_mam$SMCENTER <- factor(meta_mam$SMCENTER)
meta_mam$SMEXPEFF <- as.numeric(meta_mam$SMEXPEFF)

RES_mam_sex <- run_voom_de(
  counts_mat   = counts,
  meta_df      = meta_mam,
  group_var    = "SEX",
  covars       = c("AGE", "SMCENTER", "SMEXPEFF"),
  coef_pattern = "^SEX",                # automatically picks the SEX coefficient
  up_label     = up_label_sex,
  out_csv      = "DEG_results_with_regulation_Mammary_SEX.csv"
)

print(table(RES_mam_sex$Regulation))

# =============================================================================
# B) All tissues: Glandular vs Non-glandular
# =============================================================================

meta_all <- droplevels(meta_df)

# Fix baseline order: Non-glandular as reference, Glandular as comparison
meta_all$gland_class <- factor(meta_all$gland_class,
                               levels = c("Non-glandular", "Glandular"))

meta_all$AGE      <- factor(meta_all$AGE)
meta_all$SMCENTER <- factor(meta_all$SMCENTER)
meta_all$SMEXPEFF <- as.numeric(meta_all$SMEXPEFF)
meta_all$SEX      <- factor(meta_all$SEX)

RES_gland_vs_non <- run_voom_de(
  counts_mat   = counts,
  meta_df      = meta_all,
  group_var    = "gland_class",
  covars       = c("AGE", "SMCENTER", "SMEXPEFF", "SEX"),
  coef_pattern = "^gland_class",        # picks gland_classGlandular
  up_label     = "Glandular",
  out_csv      = "DEG_results_with_regulation_gland_class_Glandular_vs_Non.csv"
)

print(table(RES_gland_vs_non$Regulation))




# =============================================================================
# C) Three-level contrast DE:
#    gland3 = {Other Glandular, Mammary, Non-glandular}
#    Contrasts:
#      - Other_vs_NonG
#      - Mammary_vs_NonG
# =============================================================================

run_voom_contrasts <- function(counts_mat, meta_df, group_var, group_levels,
                               covars,
                               contrasts_named,   # named list: name -> contrast expression string
                               lfc_cut = 1, fdr_cut = 1e-3,
                               out_prefix = "DEG_with_regulation_") {
  
  # Ensure sample alignment
  stopifnot(all(rownames(meta_df) %in% colnames(counts_mat)))
  counts_use <- counts_mat[, rownames(meta_df), drop = FALSE]
  
  # Factorize grouping
  meta_df[[group_var]] <- factor(meta_df[[group_var]], levels = group_levels)
  
  # Prepare DGEList
  y <- edgeR::DGEList(counts = counts_use, samples = meta_df)
  
  keep <- edgeR::filterByExpr(y, group = y$samples[[group_var]])
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y, method = "TMM")
  
  # Design without intercept to make contrasts explicit
  fml <- as.formula(paste("~ 0 +", paste(c(group_var, covars), collapse = " + ")))
  design <- model.matrix(fml, data = y$samples)
  colnames(design) <- make.names(colnames(design))
  
  # Voom + lmFit
  v <- limma::voom(y, design, plot = FALSE)
  fit <- limma::lmFit(v, design)
  
  # Build contrasts
  # Convert contrast strings to use make.names()-compatible column names
  # Example: group level "Other Glandular" becomes "groupVarOther.Glandular"
  make_group_col <- function(level) make.names(paste0(group_var, level))
  group_cols <- sapply(group_levels, make_group_col, USE.NAMES = TRUE)
  
  # Replace raw group names in contrast expressions with column-safe names
  contrast_exprs <- lapply(contrasts_named, function(expr) {
    out <- expr
    for (lvl in names(group_cols)) {
      out <- gsub(paste0("`", lvl, "`"), group_cols[[lvl]], out, fixed = TRUE)
    }
    out
  })
  
  cont_mat <- limma::makeContrasts(contrasts = contrast_exprs, levels = design)
  
  fitc <- limma::contrasts.fit(fit, cont_mat)
  fitc <- limma::eBayes(fitc)
  
  extract_res <- function(coef_name, up_label, file_stub) {
    tt <- limma::topTable(fitc, coef = coef_name, number = Inf, sort.by = "P")
    RES <- tt %>%
      dplyr::mutate(
        Regulation = dplyr::case_when(
          adj.P.Val < fdr_cut & logFC >=  lfc_cut ~ paste0("Upregulated_in_",   up_label),
          adj.P.Val < fdr_cut & logFC <= -lfc_cut ~ paste0("Downregulated_in_", up_label),
          TRUE ~ "Not_significant"
        )
      )
    utils::write.csv(RES,
                     file = paste0(out_prefix, file_stub, ".csv"),
                     row.names = TRUE)
    message("Saved: ", paste0(out_prefix, file_stub, ".csv"))
    print(table(RES$Regulation))
    invisible(RES)
  }
  
  out_list <- list()
  for (nm in names(contrasts_named)) {
    # Decide "up_label" based on the contrast name convention
    # Here we assume "X_vs_NonG" means "Upregulated_in_X"
    up_label <- sub("_vs_.*$", "", nm)
    out_list[[nm]] <- extract_res(
      coef_name = nm,
      up_label  = up_label,
      file_stub = nm
    )
  }
  
  return(out_list)
}

# ---- Build gland3 in metadata ------------------------------------------------
meta3 <- meta_df
meta3$gland3 <- ifelse(is.na(meta3$gland_detail),
                       "Non-glandular",
                       as.character(meta3$gland_detail))

meta3 <- meta3[meta3$gland3 %in% c("Other Glandular", "Mammary", "Non-glandular"), , drop = FALSE]
meta3 <- droplevels(meta3)

# Covariate types
meta3$AGE      <- factor(meta3$AGE)
meta3$SMCENTER <- factor(meta3$SMCENTER)
meta3$SMEXPEFF <- as.numeric(meta3$SMEXPEFF)
meta3$SEX      <- factor(meta3$SEX)

# ---- Run contrasts ------------------------------------------------------------
contrast_results <- run_voom_contrasts(
  counts_mat    = counts,
  meta_df       = meta3,
  group_var     = "gland3",
  group_levels  = c("Other Glandular", "Mammary", "Non-glandular"),
  covars        = c("AGE", "SMCENTER", "SMEXPEFF", "SEX"),
  contrasts_named = list(
    Other_vs_NonG   = "`Other Glandular` - `Non-glandular`",
    Mammary_vs_NonG = "`Mammary` - `Non-glandular`"
  ),
  out_prefix = "DEG_with_regulation_"
)



# =============================================================================
# D) Alluvial/Sankey-style plot from the three comparisons
#   - Other_vs_NonG
#   - Mammary_vs_NonG
#   - Mammary_SEX
# =============================================================================

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggalluvial)
library(viridisLite)

# ---- Input files -------------------------------------------------------------
f_other_vs_nonG   <- "DEG_with_regulation_Other_vs_NonG.csv"
f_mammary_vs_nonG <- "DEG_with_regulation_Mammary_vs_NonG.csv"
f_mammary_sex     <- "DEG_results_with_regulation_Mammary_SEX.csv"

# ---- Helper: read a DE result file and keep Regulation ------------------------
read_reg <- function(path, newname) {
  read.csv(path, row.names = 1, check.names = FALSE) %>%
    dplyr::select(Regulation) %>%
    dplyr::rename(!!newname := Regulation) %>%
    tibble::rownames_to_column("gene")
}

reg1 <- read_reg(f_other_vs_nonG,   "Other_vs_NonG")
reg2 <- read_reg(f_mammary_vs_nonG, "Mammary_vs_NonG")
reg3 <- read_reg(f_mammary_sex,     "Mammary_SEX")

# ---- Combine into a single table (wide) --------------------------------------
all_df <- reg1 %>%
  full_join(reg2, by = "gene") %>%
  full_join(reg3, by = "gene")

# Treat missing as NS
all_df <- all_df %>%
  mutate(
    Other_vs_NonG   = coalesce(Other_vs_NonG,   "Not_significant"),
    Mammary_vs_NonG = coalesce(Mammary_vs_NonG, "Not_significant"),
    Mammary_SEX     = coalesce(Mammary_SEX,     "Not_significant")
  )

# Drop genes that are NS in all three comparisons
all_df2 <- all_df %>%
  filter(!(Other_vs_NonG == "Not_significant" &
             Mammary_vs_NonG == "Not_significant" &
             Mammary_SEX == "Not_significant"))

# Map Regulation -> RegGroup (Up / Down / NS)
to_reggroup <- function(x) {
  dplyr::case_when(
    grepl("^Upregulated_in_", x)   ~ "Up",
    grepl("^Downregulated_in_", x) ~ "Down",
    TRUE                          ~ "NS"
  )
}

long_df <- all_df2 %>%
  pivot_longer(cols = c(Other_vs_NonG, Mammary_vs_NonG, Mammary_SEX),
               names_to = "Comparison", values_to = "Regulation") %>%
  mutate(RegGroup = to_reggroup(Regulation))

# Stage order (x-axis)
stages <- c("Other_vs_NonG", "Mammary_vs_NonG", "Mammary_SEX")
long_df$Comparison <- factor(long_df$Comparison, levels = stages)

# Create a per-gene flow ID
long_df <- long_df %>%
  group_by(gene) %>%
  mutate(combo_id = paste(RegGroup[Comparison == stages[1]],
                          RegGroup[Comparison == stages[2]],
                          RegGroup[Comparison == stages[3]], sep = "->")) %>%
  ungroup()

# Aggregate counts for plotting (reduces size and speeds up ggalluvial)
agg_long <- long_df %>%
  count(Comparison, RegGroup, combo_id, name = "n")

# Factor order within each bar (bottom->top)
reg_levels <- c("Up", "NS", "Down")
agg_long$RegGroup <- factor(agg_long$RegGroup, levels = reg_levels)

# Color scheme: three different "Up" colors (by stage), and grey for Down/NS
mcols <- viridisLite::magma(7)
fill_cols <- c(
  "Up_Other"   = mcols[5],
  "Up_Mammary" = mcols[6],
  "Up_Sex"     = mcols[7],
  "Down"       = "#C8C8C8",
  "NS"         = "#F0F0F0"
)

agg_long <- agg_long %>%
  mutate(
    FillKey = case_when(
      RegGroup == "Up"   & Comparison == "Other_vs_NonG"   ~ "Up_Other",
      RegGroup == "Up"   & Comparison == "Mammary_vs_NonG" ~ "Up_Mammary",
      RegGroup == "Up"   & Comparison == "Mammary_SEX"     ~ "Up_Sex",
      RegGroup == "Down"                                   ~ "Down",
      TRUE                                                 ~ "NS"
    )
  )

alpha_vals <- c("NS" = 0.80, "Down" = 0.55, "Up" = 0.75)

p <- ggplot(
  agg_long,
  aes(x = Comparison, stratum = RegGroup, alluvium = combo_id, y = n, fill = FillKey)
) +
  geom_flow(aes(alpha = RegGroup),
            stat = "alluvium", lode.guidance = "forward") +
  geom_stratum(size = 0.2, color = "grey30") +
  geom_text(
    stat = "stratum",
    aes(label = scales::comma(after_stat(count))),
    size = 3.2, color = "black"
  ) +
  scale_fill_manual(
    values = fill_cols, name = "Status",
    breaks = c("Up_Other","Up_Mammary","Up_Sex","Down","NS"),
    labels = c("Up (Other vs Non)", "Up (Mammary vs Non)", "Up (Female vs Male in Mammary)",
               "Down", "NS")
  ) +
  scale_alpha_manual(values = alpha_vals, guide = "none") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right",
        panel.grid.minor = element_blank()) +
  labs(
    title = "Differentially Expressed Genes Across Comparisons",
    x = NULL, y = "Number of genes"
  )

print(p)

ggsave("DEG_alluvial_plot.png", p, width = 9, height = 5, dpi = 300)
write.csv(all_df2, file = "DEG_combined_three_comparisons.csv", row.names = FALSE)

