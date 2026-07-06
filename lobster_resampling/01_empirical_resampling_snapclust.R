#!/usr/bin/env Rscript

## ---------------------------------------------------------------------------
## Full-reference empirical resampling + snapclust analysis
##
## Outputs:
##   snapclust_resampling_all_runs.tsv
##   snapclust_resampling_summary_table.tsv
##
## Notes:
##   - American parental source: Americanus + AmerCook
##   - European parental source: non-Mediterranean H. gammarus, excluding
##     Americanus, AmerCook, HybridX, and Mediterranean populations
##   - Adds maximum-posterior summaries for direct comparison with the
##     size-balanced analysis.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(adegenet)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(purrr)
})

## =========================== User settings =================================

LOBSTERS_RDS <- "exeter_lobster_1591ind_79snps_52pop.rds"
OUT_PREFIX   <- "snapclust_resampling"

Ns <- c(20, 50, 100, 200, 500)
R  <- 100

HYBRID_COEF_BASE <- c(0.125, 0.25, 0.5)
HYBRID_COEF <- sort(unique(c(HYBRID_COEF_BASE, 1 - HYBRID_COEF_BASE)))

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
  `0.5`   = "0.5_A-0.5_B",
  `0.75`  = "0.25_B-0.75_A",
  `0.875` = "0.125_B-0.875_A",
  `1`     = "A"
)

class_order <- c("0", "0.125", "0.25", "0.5", "0.75", "0.875", "1")

orient_predictions <- function(proba, true_group) {
  proba <- as.data.frame(proba, check.names = FALSE)

  top <- tibble(
    Ind = seq_len(nrow(proba)),
    TrueGroup = as.character(true_group),
    TrueDeg = as.numeric(as.character(true_group)),
    AssignedGroup = colnames(proba)[max.col(proba, ties.method = "first")],
    Probability = apply(proba, 1, max, na.rm = TRUE)
  )

  degree_by_label_A <- setNames(as.numeric(names(expected_label_A)), unname(expected_label_A))
  degree_by_label_B <- setNames(as.numeric(names(expected_label_B)), unname(expected_label_B))

  eval_orientation <- function(degree_by_label, orientation) {
    assigned_deg <- unname(degree_by_label[top$AssignedGroup])
    assigned_chr <- as.character(assigned_deg)
    true_chr <- as.character(top$TrueGroup)

    true_index <- match(true_chr, class_order)
    assigned_index <- match(assigned_chr, class_order)

    top |>
      mutate(
        orientation = orientation,
        AssignedDeg = assigned_deg,
        AssignedClass = assigned_chr,
        strict_ok = !is.na(AssignedDeg) & AssignedDeg == TrueDeg,
        nearest_ok = !is.na(assigned_index) & !is.na(true_index) & abs(assigned_index - true_index) <= 1
      )
  }

  pred_A <- eval_orientation(degree_by_label_A, "A_Europe_B_America")
  pred_B <- eval_orientation(degree_by_label_B, "B_Europe_A_America")

  score_A <- mean(pred_A$strict_ok, na.rm = TRUE)
  score_B <- mean(pred_B$strict_ok, na.rm = TRUE)

  if (is.na(score_A)) score_A <- -Inf
  if (is.na(score_B)) score_B <- -Inf

  if (score_B > score_A) pred_B else pred_A
}

summarise_predictions <- function(pred, n_per_class, seed) {
  qfun <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE))

  overall <- pred |>
    summarise(
      N = n_per_class,
      n_per_class = n_per_class,
      seed = seed,
      class_label = "overall",
      metric = "overall",
      acc_strict = mean(strict_ok, na.rm = TRUE),
      strict_acc = mean(strict_ok, na.rm = TRUE),
      acc_nearest = mean(nearest_ok, na.rm = TRUE),
      nearest_acc = mean(nearest_ok, na.rm = TRUE),
      mean_max_post = mean(Probability, na.rm = TRUE),
      median_max_post = median(Probability, na.rm = TRUE),
      q25_max_post = qfun(Probability, 0.25),
      q75_max_post = qfun(Probability, 0.75),
      orientation = first(orientation)
    )

  per_class <- pred |>
    group_by(TrueGroup) |>
    summarise(
      N = n_per_class,
      n_per_class = n_per_class,
      seed = seed,
      class_label = first(TrueGroup),
      metric = first(TrueGroup),
      acc_strict = mean(strict_ok, na.rm = TRUE),
      strict_acc = mean(strict_ok, na.rm = TRUE),
      acc_nearest = mean(nearest_ok, na.rm = TRUE),
      nearest_acc = mean(nearest_ok, na.rm = TRUE),
      mean_max_post = mean(Probability, na.rm = TRUE),
      median_max_post = median(Probability, na.rm = TRUE),
      q25_max_post = qfun(Probability, 0.25),
      q75_max_post = qfun(Probability, 0.75),
      orientation = first(orientation),
      .groups = "drop"
    )

  bind_rows(overall, per_class)
}

run_one <- function(lobsters, n_per_class, seed) {
  message("Full reference: N = ", n_per_class, ", seed = ", seed)
  tabs <- make_resampled_classes(lobsters, n_per_class = n_per_class, seed = seed)
  gi <- classes_to_genind(lobsters, tabs)
  sc <- apply_snapclust(gi, HYBRID_COEF)
  pred <- orient_predictions(sc$proba, pop(gi))
  summarise_predictions(pred, n_per_class = n_per_class, seed = seed)
}

## =============================== Main =======================================

lobsters <- readRDS(LOBSTERS_RDS)

all_runs <- map_dfr(Ns, function(N) {
  map_dfr(seq_len(R), function(seed) {
    run_one(lobsters, n_per_class = N, seed = seed)
  })
})

write_tsv(all_runs, paste0(OUT_PREFIX, "_all_runs.tsv"))

summary_table <- all_runs |>
  group_by(N, n_per_class, class_label, metric) |>
  summarise(
    strict_mean = mean(acc_strict, na.rm = TRUE),
    strict_sd = sd(acc_strict, na.rm = TRUE),
    strict_median = median(acc_strict, na.rm = TRUE),
    strict_q25 = as.numeric(quantile(acc_strict, 0.25, na.rm = TRUE)),
    strict_q75 = as.numeric(quantile(acc_strict, 0.75, na.rm = TRUE)),
    nearest_mean = mean(acc_nearest, na.rm = TRUE),
    nearest_sd = sd(acc_nearest, na.rm = TRUE),
    mean_max_post_mean = mean(mean_max_post, na.rm = TRUE),
    mean_max_post_sd = sd(mean_max_post, na.rm = TRUE),
    median_max_post_median = median(median_max_post, na.rm = TRUE),
    median_max_post_q25 = as.numeric(quantile(median_max_post, 0.25, na.rm = TRUE)),
    median_max_post_q75 = as.numeric(quantile(median_max_post, 0.75, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  arrange(N, factor(class_label, levels = c("overall", class_order)))

write_tsv(summary_table, paste0(OUT_PREFIX, "_summary_table.tsv"))

message("Done. Wrote:")
message("  ", paste0(OUT_PREFIX, "_all_runs.tsv"))
message("  ", paste0(OUT_PREFIX, "_summary_table.tsv"))
