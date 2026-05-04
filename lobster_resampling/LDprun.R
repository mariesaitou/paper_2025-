#!/usr/bin/env Rscript

## ---------------------------------------------------------------------------
## Lobster hybrid detection via empirical resampling + snapclust
## LD-pruned version
##
## Inputs:
##   - exeter_lobster_1591ind_79snps_52pop.rds
##   - LD_pruned_all_SNPs_in_high_LD_pairs.tsv
##
## Outputs:
##   - snapclust_resampling_LDpruned_all_runs.tsv
##   - snapclust_resampling_LDpruned_summary_table.tsv
##   - snapclust_resampling_LDpruned_runmeta.rds
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(adegenet)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

## --------------------------- User settings ----------------------------------

LOBSTERS_RDS <- "exeter_lobster_1591ind_79snps_52pop.rds"
LD_DROP_TSV  <- "LD_pruned_all_SNPs_in_high_LD_pairs.tsv"

lobsters <- readRDS(LOBSTERS_RDS)

snps_to_drop <- scan(
  "LD_pruned_all_SNPs_in_high_LD_pairs.tsv",
  what = character(),
  skip = 1
)

lobsters <- lobsters[, !locNames(lobsters) %in% snps_to_drop]



## ------------------------------- Main ---------------------------------------

lobsters <- readRDS(LOBSTERS_RDS)

snps_to_drop <- scan(
  "LD_pruned_all_SNPs_in_high_LD_pairs.tsv",
  what = character(),
  skip = 1
)

lobsters <- lobsters[, !locNames(lobsters) %in% snps_to_drop]

all_runs <- do.call(
  bind_rows,
  lapply(Ns, function(N) {
    lapply(seq_len(R), function(r) {
      
      cat("Running N =", N, "replicate =", r, "\n")
      
      summarise_one_run(lobsters, n_per_class = N, seed = r)
    })
  }) |> unlist(recursive = FALSE)
)

summary_table <- all_runs |>
  group_by(n_per_class, metric) |>
  summarise(
    strict_mean  = mean(strict_acc),
    strict_sd    = sd(strict_acc),
    nearest_mean = mean(nearest_acc),
    nearest_sd   = sd(nearest_acc),
    .groups = "drop"
  ) |>
  arrange(n_per_class, factor(metric, levels = c("overall","0","0.125","0.25","0.5","0.75","0.875","1")))

print(summary_table)

out_raw <- file.path(getwd(), paste0(OUT_PREFIX, "_prunned_runs.tsv"))
write.table(all_runs, out_raw, sep = "\t", quote = FALSE, row.names = FALSE)
out_raw



meta <- list(
  Ns = Ns,
  R = R,
  hybrid_coef = HYBRID_COEF,
  lobsters_rds = LOBSTERS_RDS,
  timestamp = Sys.time()
)
saveRDS(meta, file.path(getwd(), paste0(OUT_PREFIX, "_runmeta.rds")))

cat("\nWrote:\n", out_raw, "\n", sep = "")


lobsters <- readRDS(LOBSTERS_RDS)

## ================= LD-pruned panel =================

LD_DROP_TSV <- "LD_pruned_all_SNPs_in_high_LD_pairs.tsv"

snps_to_drop <- read.table(
  LD_DROP_TSV,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)$SNP

snps_to_drop <- intersect(snps_to_drop, locNames(lobsters))

cat("Original SNP number:", nLoc(lobsters), "\n")
cat("Dropping SNP number:", length(snps_to_drop), "\n")
cat("Dropped SNPs:\n")
print(snps_to_drop)

lobsters <- lobsters[, loc = setdiff(locNames(lobsters), snps_to_drop)]

cat("Pruned SNP number:", nLoc(lobsters), "\n")

Ns <- c(20, 50, 100, 200, 500)
R  <- 100

HYBRID_COEF <- c(0.125, 0.25, 0.5)

OUT_PREFIX <- "snapclust_resampling_LDpruned"

## --------------------------- Helper functions -------------------------------

make_resampled_classes <- function(lobsters, n_per_class = 200, seed = 1) {
  set.seed(seed)
  
  tab_full <- lobsters@tab
  allele_cols <- colnames(tab_full)
  
  med_pops <- c("Adr","Ale","Chi","Csa","Ion","Laz","Sar","Sky","Spo","The","Tor")
  excluded <- c("AmerCook","Americanus","HybridX", med_pops)
  
  g_ids <- which(pop(lobsters) %in% setdiff(levels(pop(lobsters)), excluded))
  a_ids <- which(pop(lobsters) %in% c("Americanus","AmerCook"))
  
  loc_fac <- lobsters@loc.fac
  loc_levels <- levels(loc_fac)
  loc2idx <- lapply(loc_levels, function(L) which(loc_fac == L))
  
  make_gamete_from_tabrow <- function(tab_row) {
    gam <- integer(length(tab_row))
    for (idx in loc2idx) {
      counts <- tab_row[idx]
      if (any(is.na(counts)) || sum(counts) == 0) {
        gam[idx] <- NA_integer_
        next
      }
      probs <- counts / sum(counts)
      chosen <- sample(idx, 1, prob = probs)
      gam[chosen] <- 1L
    }
    gam
  }
  
  make_gamete_from_ids <- function(ids) {
    ind <- sample(ids, 1)
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
  
  G0_tab <- resample_pure_tab(g_ids, n_per_class)
  A1_tab <- resample_pure_tab(a_ids, n_per_class)
  F1_tab <- make_F1_tab(n_per_class)
  
  BC1A_tab <- make_backcross_tab(F1_tab,   a_ids, n_per_class)
  BC2A_tab <- make_backcross_tab(BC1A_tab, a_ids, n_per_class)
  
  BC1G_tab <- make_backcross_tab(F1_tab,   g_ids, n_per_class)
  BC2G_tab <- make_backcross_tab(BC1G_tab, g_ids, n_per_class)
  
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
  if (is.null(colnames(tab_mat))) colnames(tab_mat) <- colnames(template_genind@tab)
  
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
        tab_mat = tab,
        pop_name = nm,
        prefix = paste0("C", gsub("\\.", "", nm))
      )
    },
    tab = class_tabs,
    nm  = names(class_tabs)
  )
  do.call(adegenet::repool, gens)
}

apply_exeter_snapclust_pipeline <- function(genind_data, hybrid_coef = c(0.5), true_group_list) {
  gc()
  
  snapclust_result <- adegenet::snapclust(
    x = genind_data,
    k = 2,
    hybrids = TRUE,
    hybrid.coef = hybrid_coef
  )
  
  pop_labels_chr <- as.character(adegenet::pop(genind_data))
  
  group_names <- names(true_group_list)
  if (is.null(group_names) || any(group_names == "")) {
    stop("true_group_list must be a named list.")
  }
  
  pop_to_group <- function(p) {
    hit <- which(vapply(true_group_list, function(v) p %in% v, logical(1)))
    if (length(hit) == 0) return(NA_character_)
    group_names[hit[1]]
  }
  
  true_groups_chr <- vapply(pop_labels_chr, pop_to_group, character(1))
  true_groups_chr[is.na(true_groups_chr)] <- pop_labels_chr[is.na(true_groups_chr)]
  
  proba <- as.data.frame(snapclust_result$proba)
  proba$Ind <- seq_len(nrow(proba))
  proba$TrueGroup <- factor(true_groups_chr, levels = unique(group_names))
  
  out <- tibble::as_tibble(proba) |>
    dplyr::relocate(Ind, TrueGroup) |>
    tidyr::pivot_longer(
      cols = -c(Ind, TrueGroup),
      names_to = "AssignedGroup",
      values_to = "Probability"
    ) |>
    dplyr::mutate(
      Ind = as.factor(Ind),
      AssignedGroup = factor(AssignedGroup, levels = unique(AssignedGroup))
    )
  
  gc()
  out
}

expected_label <- c(
  `0`     = "A",
  `0.125` = "0.875_A-0.125_B",
  `0.25`  = "0.75_A-0.25_B",
  `0.5`   = "0.5_A-0.5_B",
  `0.75`  = "0.25_A-0.75_B",
  `0.875` = "0.125_A-0.875_B",
  `1`     = "B"
)

allowed_nearest <- list(
  `0`     = c("A", "0.875_A-0.125_B"),
  `0.125` = c("A", "0.875_A-0.125_B", "0.75_A-0.25_B"),
  `0.25`  = c("0.875_A-0.125_B", "0.75_A-0.25_B", "0.5_A-0.5_B"),
  `0.5`   = c("0.75_A-0.25_B", "0.5_A-0.5_B", "0.25_A-0.75_B"),
  `0.75`  = c("0.5_A-0.5_B", "0.25_A-0.75_B", "0.125_A-0.875_B"),
  `0.875` = c("0.25_A-0.75_B", "0.125_A-0.875_B", "B"),
  `1`     = c("0.125_A-0.875_B", "B")
)

summarise_one_run <- function(lobsters, n_per_class, seed) {
  tabs <- make_resampled_classes(lobsters, n_per_class = n_per_class, seed = seed)
  gi <- classes_to_genind(lobsters, tabs)
  
  true_group_list <- as.list(names(tabs))
  names(true_group_list) <- names(tabs)
  
  snap <- apply_exeter_snapclust_pipeline(
    genind_data = gi,
    hybrid_coef = HYBRID_COEF,
    true_group_list = true_group_list
  )
  
  pred <- snap |>
    group_by(Ind) |>
    slice_max(order_by = Probability, n = 1, with_ties = FALSE) |>
    ungroup() |>
    mutate(
      TrueGroup = as.character(TrueGroup),
      AssignedGroup = as.character(AssignedGroup),
      strict_ok = (AssignedGroup == unname(expected_label[TrueGroup])),
      nearest_ok = mapply(
        function(tg, ag) ag %in% allowed_nearest[[tg]],
        TrueGroup,
        AssignedGroup
      )
    )
  
  overall <- pred |>
    reframe(
      n_per_class = n_per_class,
      seed = seed,
      metric = "overall",
      strict_acc = mean(strict_ok),
      nearest_acc = mean(nearest_ok)
    )
  
  per_class <- pred |>
    group_by(TrueGroup) |>
    reframe(
      n_per_class = n_per_class,
      seed = seed,
      metric = TrueGroup,
      strict_acc = mean(strict_ok),
      nearest_acc = mean(nearest_ok)
    )
  
  bind_rows(overall, per_class)
}

## ------------------------------- Main ---------------------------------------

lobsters <- readRDS(LOBSTERS_RDS)

## LD pruning
snps_to_drop <- read.table(
  LD_DROP_TSV,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)$SNP

snps_to_drop <- intersect(snps_to_drop, locNames(lobsters))

cat("Original SNP number:", nLoc(lobsters), "\n")
cat("Dropping SNP number:", length(snps_to_drop), "\n")
cat("Dropped SNPs:\n")
print(snps_to_drop)

lobsters <- lobsters[, loc = setdiff(locNames(lobsters), snps_to_drop)]

cat("Pruned SNP number:", nLoc(lobsters), "\n")

all_runs <- do.call(
  bind_rows,
  lapply(Ns, function(N) {
    lapply(seq_len(R), function(r) {
      summarise_one_run(lobsters, n_per_class = N, seed = r)
    })
  }) |> unlist(recursive = FALSE)
)

summary_table <- all_runs |>
  group_by(n_per_class, metric) |>
  summarise(
    strict_mean  = mean(strict_acc),
    strict_sd    = sd(strict_acc),
    nearest_mean = mean(nearest_acc),
    nearest_sd   = sd(nearest_acc),
    .groups = "drop"
  ) |>
  arrange(
    n_per_class,
    factor(metric, levels = c("overall","0","0.125","0.25","0.5","0.75","0.875","1"))
  )




results_list <- list()

for (N in Ns) {
  cat("=== Running N =", N, "===\n")
  
  res_N <- bind_rows(lapply(seq_len(R), function(r) {
    cat(sprintf("N=%d replicate=%d\n", N, r))
    summarise_one_run(lobsters, n_per_class = N, seed = r)
  }))
  
  # Nごとに保存
  out_file <- paste0(OUT_PREFIX, "_N", N, ".tsv")
  write.table(res_N, out_file, sep="\t", quote=FALSE, row.names=FALSE)
  
  cat("Saved:", out_file, "\n")
  
  # 必要なら集計だけ保持
  results_list[[as.character(N)]] <- res_N
}

print(summary_table)

out_raw <- file.path(getwd(), paste0(OUT_PREFIX, "_all_runs.tsv"))
write.table(all_runs, out_raw, sep = "\t", quote = FALSE, row.names = FALSE)

out_summary <- file.path(getwd(), paste0(OUT_PREFIX, "_summary_table.tsv"))
write.table(summary_table, out_summary, sep = "\t", quote = FALSE, row.names = FALSE)

meta <- list(
  Ns = Ns,
  R = R,
  hybrid_coef = HYBRID_COEF,
  lobsters_rds = LOBSTERS_RDS,
  ld_drop_tsv = LD_DROP_TSV,
  dropped_snps = snps_to_drop,
  n_loci_after_pruning = nLoc(lobsters),
  timestamp = Sys.time()
)

out_meta <- file.path(getwd(), paste0(OUT_PREFIX, "_runmeta.rds"))
saveRDS(meta, out_meta)

cat(
  "\nWrote:\n",
  out_raw, "\n",
  out_summary, "\n",
  out_meta, "\n",
  sep = ""
)