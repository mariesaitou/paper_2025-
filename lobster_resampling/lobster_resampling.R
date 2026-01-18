#!/usr/bin/env Rscript

## ---------------------------------------------------------------------------
## Lobster hybrid detection via empirical resampling + snapclust
##
## Outputs:
##   Fig 1: Accuracy plots from resampling runs (median + IQR)
##   Fig 2: PCA (learn on real data, project resampled classes)
##   Fig 3: snapclust confusion heatmap (single run)
##   Fig 4: snapclust calibration plot (single run)
##
## Required inputs:
##   - exeter_lobster_1591ind_79snps_52pop.rds (genind)
##   - snapclust_resampling_all_runs.tsv (accuracy runs)
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(adegenet)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(stringr)
  library(ggrepel)
  library(ggnewscale)
  library(RColorBrewer)
  library(viridisLite)
  library(patchwork)
})

theme_set(theme_bw())

## =========================== User settings =================================



LOBSTERS_RDS <- "exeter_lobster_1591ind_79snps_52pop.rds"
OUT_PREFIX   <- "snapclust_resampling"
RUNS_TSV     <- paste0(OUT_PREFIX, "_all_runs.tsv")

Ns <- c(20, 50, 100, 200, 500)
R  <- 100

HYBRID_COEF_BASE <- c(0.125, 0.25, 0.5)
HYBRID_COEF <- sort(unique(c(HYBRID_COEF_BASE, 1 - HYBRID_COEF_BASE)))

N_use    <- 50
seed_use <- 1

N_heat    <- 500
seed_heat <- 1

MED_POPS <- c("Adr", "Ale", "Chi", "Csa", "Ion", "Laz", "Sar", "Sky", "Spo", "The", "Tor")

## =========================== Helper functions ===============================

make_resampled_classes <- function(lobsters, n_per_class = 200, seed = 1) {
  set.seed(seed)
  
  tab_full    <- lobsters@tab
  allele_cols <- colnames(tab_full)
  
  excluded <- c("AmerCook", "Americanus", "HybridX", MED_POPS)
  
  g_ids <- which(pop(lobsters) %in% setdiff(levels(pop(lobsters)), excluded))
  a_ids <- which(pop(lobsters) %in% c("Americanus", "AmerCook"))
  
  loc_fac    <- lobsters@loc.fac
  loc_levels <- levels(loc_fac)
  loc2idx    <- lapply(loc_levels, function(L) which(loc_fac == L))
  
  make_gamete_from_tabrow <- function(tab_row) {
    gam <- integer(length(tab_row))
    for (idx in loc2idx) {
      counts <- tab_row[idx]
      if (any(is.na(counts)) || sum(counts) == 0) {
        gam[idx] <- NA_integer_
        next
      }
      probs  <- counts / sum(counts)
      chosen <- sample(idx, 1, prob = probs)
      gam[chosen] <- 1L
    }
    gam
  }
  
  make_gamete_from_ids <- function(ids) {
    ind     <- sample(ids, 1)
    tab_row <- as.numeric(tab_full[ind, ])
    make_gamete_from_tabrow(tab_row)
  }
  
  make_child_tab <- function(gam1, gam2) {
    out <- gam1 + gam2
    storage.mode(out) <- "integer"
    out
  }
  
  resample_pure_tab <- function(ids, n) {
    sampled <- sample(ids, n, replace = TRUE)
    mat <- tab_full[sampled, , drop = FALSE]
    storage.mode(mat) <- "integer"
    colnames(mat) <- allele_cols
    mat
  }
  
  make_F1_tab <- function(n) {
    children <- vector("list", n)
    for (i in seq_len(n)) {
      g_gam <- make_gamete_from_ids(g_ids)
      a_gam <- make_gamete_from_ids(a_ids)
      children[[i]] <- make_child_tab(g_gam, a_gam)
    }
    mat <- do.call(rbind, children)
    storage.mode(mat) <- "integer"
    colnames(mat) <- allele_cols
    mat
  }
  
  make_backcross_tab <- function(hybrid_tab_mat, parent_ids, n) {
    children <- vector("list", n)
    for (i in seq_len(n)) {
      h_row <- as.numeric(hybrid_tab_mat[sample(seq_len(nrow(hybrid_tab_mat)), 1), ])
      h_gam <- make_gamete_from_tabrow(h_row)
      p_gam <- make_gamete_from_ids(parent_ids)
      children[[i]] <- make_child_tab(h_gam, p_gam)
    }
    mat <- do.call(rbind, children)
    storage.mode(mat) <- "integer"
    colnames(mat) <- allele_cols
    mat
  }
  
  G0_tab <- resample_pure_tab(g_ids, n_per_class)               # 0
  A1_tab <- resample_pure_tab(a_ids, n_per_class)               # 1
  F1_tab <- make_F1_tab(n_per_class)                            # 0.5
  
  BC1A_tab <- make_backcross_tab(F1_tab,   a_ids, n_per_class)  # 0.75
  BC2A_tab <- make_backcross_tab(BC1A_tab, a_ids, n_per_class)  # 0.875
  
  BC1G_tab <- make_backcross_tab(F1_tab,   g_ids, n_per_class)  # 0.25
  BC2G_tab <- make_backcross_tab(BC1G_tab, g_ids, n_per_class)  # 0.125
  
  list(
    `0`     = G0_tab,
    `0.125` = BC2G_tab,
    `0.25`  = BC1G_tab,
    `0.5`   = F1_tab,
    `0.75`  = BC1A_tab,
    `0.875` = BC2A_tab,
    `1`     = A1_tab
  )
}

rebuild_genind_from_tab_simple <- function(template_genind, tab_mat, pop_name, prefix) {
  tab_mat <- as.matrix(tab_mat)
  storage.mode(tab_mat) <- "integer"
  if (is.null(colnames(tab_mat))) {
    colnames(tab_mat) <- colnames(template_genind@tab)
  }
  
  base <- template_genind[rep(1, nrow(tab_mat)), ]
  base@tab <- tab_mat
  base@pop <- factor(rep(pop_name, nInd(base)), levels = pop_name)
  indNames(base) <- paste0(prefix, "_", seq_len(nInd(base)))
  base
}

classes_to_genind <- function(lobsters, class_tabs) {
  gens <- Map(
    f = function(tab, nm) {
      rebuild_genind_from_tab_simple(
        template_genind = lobsters,
        tab_mat         = tab,
        pop_name        = nm,
        prefix          = paste0("C", gsub("\\.", "", nm))
      )
    },
    tab = class_tabs,
    nm  = names(class_tabs)
  )
  do.call(adegenet::repool, gens)
}

apply_snapclust <- function(genind_data, hybrid_coef) {
  adegenet::snapclust(
    x           = genind_data,
    k           = 2,
    hybrids     = TRUE,
    hybrid.coef = hybrid_coef
  )
}

## ---------------------------- PCA helpers -----------------------------------

impute_col_mean <- function(X) {
  if (!anyNA(X)) return(X)
  for (j in seq_len(ncol(X))) {
    m <- mean(X[, j], na.rm = TRUE)
    X[is.na(X[, j]), j] <- m
  }
  X
}

pca_learn_prcomp <- function(genind_obj) {
  X <- as.matrix(genind_obj@tab)
  X <- impute_col_mean(X)
  prcomp(X, center = TRUE, scale. = TRUE)
}

pca_project <- function(pca_fit, X_new, pcs = 1:2) {
  X_new <- X_new[, names(pca_fit$center), drop = FALSE]
  X_new <- impute_col_mean(X_new)
  scale(X_new, center = pca_fit$center, scale = pca_fit$scale) %*% pca_fit$rotation[, pcs, drop = FALSE]
}

## ---------------------------- Labels ----------------------------------------

deg_to_class <- function(x) {
  dplyr::case_when(
    x == 0     ~ "Europe",
    x == 1     ~ "America",
    x == 0.5   ~ "F1",
    x == 0.125 ~ "BC2 (Europe)",
    x == 0.25  ~ "BC1 (Europe)",
    x == 0.75  ~ "BC1 (America)",
    x == 0.875 ~ "BC2 (America)",
    TRUE       ~ as.character(x)
  )
}

## FIX: Always map A/B strings after possible relabeling checks downstream.
expected_label_A <- c(
  `0`     = "A",
  `0.125` = "0.875_A-0.125_B",
  `0.25`  = "0.75_A-0.25_B",
  `0.5`   = "0.5_A-0.5_B",
  `0.75`  = "0.25_A-0.75_B",
  `0.875` = "0.125_A-0.875_B",
  `1`     = "B"
)
expected_label_B <- c(
  `0`     = "B",
  `0.125` = "0.875_B-0.125_A",
  `0.25`  = "0.75_B-0.25_A",
  `0.5`   = "0.5_A-0.5_B",      # symmetric
  `0.75`  = "0.25_B-0.75_A",
  `0.875` = "0.125_B-0.875_A",
  `1`     = "A"
)

## =============================== Main =======================================

lobsters <- readRDS(LOBSTERS_RDS)

## ============================================================
## Fig 1) Accuracy plots (from TSV)
## ============================================================

df <- read.table(RUNS_TSV, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

if (!("class_label" %in% colnames(df))) {
  if ("class" %in% colnames(df)) {
    df$class_label <- df$class
  } else {
    stop("Input TSV must contain either 'class_label' or 'class'.")
  }
}

df_strict <- df %>%
  filter(class_label != "overall") %>%
  mutate(class_label = as.character(class_label)) %>%
  filter(!class_label %in% c("0", "1"))

facet_levels <- c("0.125", "0.25", "0.5", "0.75", "0.875")
facet_labels <- c(
  "0.125" = "BC2 (Europe)",
  "0.25"  = "BC1 (Europe)",
  "0.5"   = "F1",
  "0.75"  = "BC1 (America)",
  "0.875" = "BC2 (America)"
)

df_sum <- df_strict %>%
  mutate(class_label = factor(class_label, levels = facet_levels),
         N = as.integer(N)) %>%
  group_by(class_label, N) %>%
  summarise(
    med = median(acc_strict, na.rm = TRUE),
    q25 = quantile(acc_strict, 0.25, na.rm = TRUE),
    q75 = quantile(acc_strict, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

p_strict_hyb <- ggplot(df_sum, aes(x = N, y = med, group = 1)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.25) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  facet_wrap(~ class_label, ncol = 5, labeller = as_labeller(facet_labels)) +
  scale_x_continuous(
    trans  = "sqrt",
    breaks = c(20, 50, 100, 200, 500),
    labels = c("20", "50", "100", "200", "500")
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "N per class (sqrt scale)", y = "Strict accuracy (median; IQR)") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
    strip.text = element_text(color = "black"),
    panel.grid.minor = element_blank()
  )

print(p_strict_hyb)

## ============================================================
## Fig 2) PCA: learn on real, project simulated samples
## ============================================================

pca_real <- pca_learn_prcomp(lobsters)
var_pct  <- (pca_real$sdev^2) / sum(pca_real$sdev^2) * 100
xlab_pca <- sprintf("PC1 (%.1f%%)", var_pct[1])
ylab_pca <- sprintf("PC2 (%.1f%%)", var_pct[2])

df_real <- tibble(
  Ind   = indNames(lobsters),
  Group = as.character(pop(lobsters)),
  PC1   = pca_real$x[, 1],
  PC2   = pca_real$x[, 2]
) %>%
  mutate(
    Group4 = case_when(
      Group %in% c("Americanus", "AmerCook") ~ "America",
      Group == "HybridX"                    ~ "HybridX",
      Group %in% MED_POPS                   ~ "Mediterranean",
      TRUE                                  ~ "Atlantic Europe"
    ),
    Group4 = factor(Group4, levels = c("Atlantic Europe", "Mediterranean", "HybridX", "America"))
  )

tabs <- make_resampled_classes(lobsters, n_per_class = N_use, seed = seed_use)
gi   <- classes_to_genind(lobsters, tabs)

Z_sim <- pca_project(pca_real, as.matrix(gi@tab), pcs = 1:2)

df_sim <- tibble(
  Ind   = indNames(gi),
  Group = as.character(pop(gi)),
  PC1   = Z_sim[, 1],
  PC2   = Z_sim[, 2],
  hyb   = suppressWarnings(as.numeric(as.character(pop(gi))))
) %>%
  mutate(
    Class_proj = case_when(
      Group == "1"   ~ "America",
      Group == "0"   ~ "Europe",
      Group == "0.5" ~ "F1",
      TRUE           ~ "Backcross"
    ),
    Class_proj = factor(Class_proj, levels = c("Europe", "Backcross", "F1", "America"))
  )

## FIX: define shared axis limits for patchwork
xlim_all <- range(c(df_real$PC1, df_sim$PC1), finite = TRUE)
ylim_all <- range(c(df_real$PC2, df_sim$PC2), finite = TRUE)

vir_cols <- viridis(256, option = "viridis")
col_eu   <- vir_cols[1]
col_hx   <- vir_cols[128]
col_am   <- vir_cols[256]

label_df_real <- df_real %>%
  group_by(Group4) %>%
  summarise(PC1 = median(PC1), PC2 = median(PC2), .groups = "drop")

p_real <- ggplot(df_real) +
  geom_point(aes(PC1, PC2, shape = Group4, fill = Group4),
             size = 3.0, alpha = 0.8, color = "black", stroke = 0.8) +
  scale_shape_manual(values = c("Atlantic Europe"=21,"Mediterranean"=4,"HybridX"=23,"America"=24),
                     name = "Group (real)") +
  scale_fill_manual(values = c("Atlantic Europe"=col_eu,"Mediterranean"=NA,"HybridX"=col_hx,"America"=col_am),
                    name = "Group (real)", na.translate = FALSE) +
  geom_label_repel(data = label_df_real,
                   aes(PC1, PC2, label = Group4),
                   inherit.aes = FALSE,
                   size = 4, fontface = "bold", color = "black",
                   fill = "white", label.size = 0.25,
                   box.padding = 0.35, point.padding = 0.25,
                   segment.color = "grey50", min.segment.length = 0) +
  coord_cartesian(xlim = xlim_all, ylim = ylim_all) +
  labs(title = "PCA from Real SNP Genotype Data", x = xlab_pca, y = ylab_pca) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

p_proj <- ggplot(df_sim) +
  geom_point(aes(PC1, PC2, fill = hyb, shape = Class_proj),
             size = 3.0, alpha = 0.8, color = "black", stroke = 0.8) +
  scale_shape_manual(values = c("Europe"=21,"Backcross"=22,"F1"=23,"America"=24),
                     name = "Class (Resampled)") +
  scale_fill_viridis_c(option = "viridis", limits = c(0,1),
                       breaks = c(0,0.25,0.5,0.75,1),
                       name = "Hybrid degree") +
  coord_cartesian(xlim = xlim_all, ylim = ylim_all) +
  labs(title = sprintf("Resampled data projected on the real PCA space"),
       x = xlab_pca, y = ylab_pca) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), legend.box = "vertical")

p_combined <- (p_real + theme(legend.position = "none")) | p_proj
print(p_combined)

## ============================================================
## Fig 3/4) snapclust confusion heatmap + calibration
## ============================================================

tabs_h <- make_resampled_classes(lobsters, n_per_class = N_heat, seed = seed_heat)
gi_h   <- classes_to_genind(lobsters, tabs_h)

sc <- apply_snapclust(gi_h, HYBRID_COEF)
proba <- as.data.frame(sc$proba)

snap_df <- tibble(
  Ind = seq_len(nrow(proba)),
  TrueGroup = as.character(pop(gi_h))
) %>%
  bind_cols(proba) %>%
  pivot_longer(cols = -c(Ind, TrueGroup),
               names_to = "AssignedGroup",
               values_to = "Probability")

## Confusion (max posterior per individual)
conf_df <- snap_df %>%
  group_by(Ind) %>%
  slice_max(Probability, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    TrueDeg = suppressWarnings(as.numeric(TrueGroup)),
    TrueClass = deg_to_class(TrueDeg)
  ) %>%
  count(TrueClass, AssignedGroup, name = "n") %>%
  group_by(TrueClass) %>%
  mutate(prop = n / sum(n), label = sprintf("%.2f", prop)) %>%
  ungroup()

lvl_order <- c("Europe","BC2 (Europe)","BC1 (Europe)","F1","BC1 (America)","BC2 (America)","America")

conf_df <- conf_df %>%
  mutate(TrueClass = factor(TrueClass, levels = lvl_order))

## NOTE: AssignedGroup labels are snapclust-internal (A/B + hybrids).
## Here we keep them as-is for the heatmap x-axis unless you want forced relabeling.

p_heat <- ggplot(conf_df, aes(x = AssignedGroup, y = TrueClass, fill = prop)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient(low = "white", high = "#1da2d8", limits = c(0, 1)) +
  labs(x = "Assigned class (snapclust label)", y = "True class", fill = "Classification fraction") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid = element_blank())

print(p_heat)

## Calibration (Fig 4): handle possible A/B inversion
calib_base <- snap_df %>%
  group_by(Ind) %>%
  slice_max(Probability, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    strict_ok_A = AssignedGroup == expected_label_A[TrueGroup],
    strict_ok_B = AssignedGroup == expected_label_B[TrueGroup],
    strict_ok   = pmax(as.integer(strict_ok_A), as.integer(strict_ok_B)) == 1L,
    prob_bin    = cut(Probability, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)
  ) %>%
  group_by(prob_bin) %>%
  summarise(n = n(), accuracy = mean(strict_ok), .groups = "drop")

p_calib <- ggplot(calib_base, aes(x = prob_bin, y = accuracy, group = 1)) +
  geom_point(size = 2) +
  geom_line() +
  labs(x = "Posterior probability of assigned class", y = "Fraction correctly classified") +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1),
        panel.grid.minor = element_blank())

print(p_calib)
