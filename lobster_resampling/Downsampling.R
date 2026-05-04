#!/usr/bin/env Rscript

## ---------------------------------------------------------------------------
## Sensitivity analysis:
## Effect of parental reference size imbalance
##
## European parental pool is randomly subsampled to match the American parental
## reference size (n = 38) before empirical resampling.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(adegenet)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

## --------------------------- User settings ----------------------------------

LOBSTERS_RDS <- "exeter_lobster_1591ind_79snps_52pop.rds"

Ns <- c(20, 50, 100, 200)
R  <- 50

HYBRID_COEF <- c(0.125, 0.25, 0.5)

OUT_PREFIX <- "snapclust_resampling_balanced_parental_reference"

PARENTAL_POOL_N <- 38

## --------------------------- Helper functions -------------------------------

make_resampled_classes <- function(
    lobsters,
    n_per_class = 200,
    seed = 1,
    balance_parental_pool = TRUE,
    parental_pool_n = 38
) {
  set.seed(seed)
  
  tab_full <- lobsters@tab
  allele_cols <- colnames(tab_full)
  
  med_pops <- c("Adr","Ale","Chi","Csa","Ion","Laz","Sar","Sky","Spo","The","Tor")
  excluded <- c("AmerCook","Americanus","HybridX", med_pops)
  
  g_ids <- which(pop(lobsters) %in% setdiff(levels(pop(lobsters)), excluded))
  a_ids <- which(pop(lobsters) %in% c("Americanus","AmerCook"))
  
  if (balance_parental_pool) {
    g_ids <- sample(g_ids, parental_pool_n, replace = FALSE)
    a_ids <- sample(a_ids, parental_pool_n, replace = FALSE)
  }
  
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
        lobsters,
        tab,
        nm,
        paste0("C", gsub("\\.", "", nm))
      )
    },
    tab = class_tabs,
    nm  = names(class_tabs)
  )
  do.call(adegenet::repool, gens)
}

apply_exeter_snapclust_pipeline <- function(genind_data, hybrid_coef, true_group_list) {
  gc()
  
  snapclust_result <- adegenet::snapclust(
    x = genind_data,
    k = 2,
    hybrids = TRUE,
    hybrid.coef = hybrid_coef
  )
  
  pop_labels_chr <- as.character(adegenet::pop(genind_data))
  group_names <- names(true_group_list)
  
  pop_to_group <- function(p) {
    hit <- which(vapply(true_group_list, function(v) p %in% v, logical(1)))
    if (length(hit) == 0) return(NA_character_)
    group_names[hit[1]]
  }
  
  true_groups_chr <- vapply(pop_labels_chr, pop_to_group, character(1))
  
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
  tabs <- make_resampled_classes(
    lobsters,
    n_per_class = n_per_class,
    seed = seed,
    balance_parental_pool = TRUE,
    parental_pool_n = PARENTAL_POOL_N
  )
  
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
      strict_ok = AssignedGroup == unname(expected_label[TrueGroup]),
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

all_runs <- do.call(
  bind_rows,
  lapply(Ns, function(N) {
    lapply(seq_len(R), function(r) {
      cat(sprintf(
        "Running balanced parental reference: N = %d, replicate = %d / %d\n",
        N, r, R
      ))
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

print(summary_table)

out_raw <- file.path(getwd(), paste0(OUT_PREFIX, "_all_runs.tsv"))
out_summary <- file.path(getwd(), paste0(OUT_PREFIX, "_summary_table.tsv"))
out_meta <- file.path(getwd(), paste0(OUT_PREFIX, "_runmeta.rds"))

write.table(all_runs, out_raw, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary_table, out_summary, sep = "\t", quote = FALSE, row.names = FALSE)

meta <- list(
  Ns = Ns,
  R = R,
  hybrid_coef = HYBRID_COEF,
  lobsters_rds = LOBSTERS_RDS,
  parental_pool_n = PARENTAL_POOL_N,
  analysis = "European parental reference subsampled to match American parental reference size",
  timestamp = Sys.time()
)

saveRDS(meta, out_meta)

cat("\nWrote:\n", out_raw, "\n", out_summary, "\n", out_meta, "\n", sep = "")






library(readr)
library(dplyr)
library(ggplot2)

full <- read_tsv("snapclust_resampling_summary_table.tsv") %>%
  mutate(analysis = "Full panel")

ldpruned <- read_tsv("snapclust_resampling_LDpruned_summary_table.tsv") %>%
  mutate(analysis = "LD pruned")

balanced <- read_tsv("snapclust_resampling_balanced_parental_reference_summary_table.tsv") %>%
  mutate(analysis = "Balanced parental")

df <- bind_rows(full, ldpruned, balanced)

df_plot <- df %>%
  filter(metric != "overall")

p <- ggplot(df_plot, aes(x = n_per_class, y = strict_mean, color = analysis)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ metric, nrow = 2) +
  scale_x_continuous(breaks = sort(unique(df_plot$n_per_class))) +
  labs(
    x = "N per class",
    y = "Strict accuracy",
    color = "Analysis"
  ) +
  theme_bw()

p

ggsave("LD_and_balanced_strict_accuracy_comparison.pdf", p, width = 10, height = 5)
ggsave("LD_and_balanced_strict_accuracy_comparison.png", p, width = 10, height = 5, dpi = 300)


ld_raw <- bind_rows(lapply(ld_files[file.exists(ld_files)], read_tsv))

facet_levels <- c("0.125", "0.25", "0.5", "0.75", "0.875")
facet_labels <- c(
  "0.125" = "BC2 (Europe)",
  "0.25"  = "BC1 (Europe)",
  "0.5"   = "F1",
  "0.75"  = "BC1 (America)",
  "0.875" = "BC2 (America)"
)

df_sum_ld <- ld_raw %>%
  filter(metric %in% facet_levels) %>%
  mutate(
    class_label = factor(metric, levels = facet_levels),
    N = as.integer(n_per_class)
  ) %>%
  group_by(class_label, N) %>%
  summarise(
    med = median(strict_acc, na.rm = TRUE),
    q25 = quantile(strict_acc, 0.25, na.rm = TRUE),
    q75 = quantile(strict_acc, 0.75, na.rm = TRUE),
    .groups = "drop"
  )


ld_raw <- bind_rows(lapply(ld_files[file.exists(ld_files)], read_tsv))

facet_levels <- c("0.125", "0.25", "0.5", "0.75", "0.875")
facet_labels <- c(
  "0.125" = "BC2 (Europe)",
  "0.25"  = "BC1 (Europe)",
  "0.5"   = "F1",
  "0.75"  = "BC1 (America)",
  "0.875" = "BC2 (America)"
)

df_sum_ld <- ld_raw %>%
  filter(metric %in% facet_levels) %>%
  mutate(
    class_label = factor(metric, levels = facet_levels),
    N = as.integer(n_per_class)
  ) %>%
  group_by(class_label, N) %>%
  summarise(
    med = median(strict_acc, na.rm = TRUE),
    q25 = quantile(strict_acc, 0.25, na.rm = TRUE),
    q75 = quantile(strict_acc, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

p_strict_ld <- ggplot(df_sum_ld, aes(x = N, y = med, group = 1)) +
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
  labs(
    x = "N per class (sqrt scale)",
    y = "Strict accuracy (median; IQR)"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
    strip.text = element_text(color = "black"),
    panel.grid.minor = element_blank()
  )

print(p_strict_ld)

ggsave("LD_pruned_strict_accuracy_mainstyle.pdf", p_strict_ld, width = 10, height = 3)
ggsave("LD_pruned_strict_accuracy_mainstyle.png", p_strict_ld, width = 10, height = 3, dpi = 300)




library(readr)
library(dplyr)
library(ggplot2)

full <- read_tsv("snapclust_resampling_summary_table.tsv") %>%
  mutate(analysis = "Full reference")

balanced <- read_tsv("snapclust_resampling_balanced_parental_reference_summary_table.tsv") %>%
  mutate(analysis = "Balanced parental reference")

df <- bind_rows(full, balanced)

df_plot <- df %>%
  filter(metric != "overall")

p <- ggplot(df_plot, aes(x = n_per_class, y = strict_mean, color = analysis)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ metric, nrow = 2) +
  scale_x_continuous(breaks = sort(unique(df_plot$n_per_class))) +
  labs(
    x = "N per class",
    y = "Strict accuracy",
    color = "Analysis"
  ) +
  theme_bw()


p




library(readr)
library(dplyr)
library(ggplot2)

balanced <- read_tsv("snapclust_resampling_balanced_parental_reference_summary_table.tsv") %>%
  mutate(analysis = "Balanced parental reference")

full <- read_tsv("snapclust_resampling_summary_table.tsv") %>%
  mutate(analysis = "Full reference")

df <- bind_rows(full, balanced)

facet_levels <- c("0.125", "0.25", "0.5", "0.75", "0.875")
facet_labels <- c(
  "0.125" = "BC2 (Europe)",
  "0.25"  = "BC1 (Europe)",
  "0.5"   = "F1",
  "0.75"  = "BC1 (America)",
  "0.875" = "BC2 (America)"
)

df_plot <- df %>%
  filter(metric %in% facet_levels) %>%
  mutate(
    class_label = factor(metric, levels = facet_levels),
    N = as.integer(n_per_class)
  )

df_sum <- df_plot %>%
  group_by(class_label, N, analysis) %>%
  summarise(
    med = median(strict_mean, na.rm = TRUE),
    q25 = quantile(strict_mean, 0.25, na.rm = TRUE),
    q75 = quantile(strict_mean, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

p <- ggplot(df_sum, aes(x = N, y = med, color = analysis, fill = analysis, group = analysis)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  facet_wrap(~ class_label, ncol = 5, labeller = as_labeller(facet_labels)) +
  scale_x_continuous(
    trans  = "sqrt",
    breaks = c(20, 50, 100, 200, 500),
    labels = c("20", "50", "100", "200", "500")
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "N per class (sqrt scale)",
    y = "Strict accuracy (median; IQR)",
    color = "Analysis",
    fill = "Analysis"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
    strip.text = element_text(color = "black"),
    panel.grid.minor = element_blank()
  )

p

ggsave("balanced_vs_full_strict_accuracy.pdf", p, width = 10, height = 3)
ggsave("balanced_vs_full_strict_accuracy.png", p, width = 10, height = 3, dpi = 300)



library(readr)
library(dplyr)
library(ggplot2)

balanced_raw <- read_tsv("snapclust_resampling_balanced_parental_reference_all_runs.tsv")

facet_levels <- c("0.125", "0.25", "0.5", "0.75", "0.875")
facet_labels <- c(
  "0.125" = "BC2 (Europe)",
  "0.25"  = "BC1 (Europe)",
  "0.5"   = "F1",
  "0.75"  = "BC1 (America)",
  "0.875" = "BC2 (America)"
)

df_sum_balanced <- balanced_raw %>%
  filter(metric %in% facet_levels) %>%
  mutate(
    class_label = factor(metric, levels = facet_levels),
    N = as.integer(n_per_class)
  ) %>%
  group_by(class_label, N) %>%
  summarise(
    med = median(strict_acc, na.rm = TRUE),
    q25 = quantile(strict_acc, 0.25, na.rm = TRUE),
    q75 = quantile(strict_acc, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

p_balanced <- ggplot(df_sum_balanced, aes(x = N, y = med, group = 1)) +
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
  labs(
    x = "N per class (sqrt scale)",
    y = "Strict accuracy (median; IQR)"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
    strip.text = element_text(color = "black"),
    panel.grid.minor = element_blank()
  )

p_balanced

ggsave("balanced_parental_reference_strict_accuracy_mainstyle.pdf", p_balanced, width = 10, height = 3)
ggsave("balanced_parental_reference_strict_accuracy_mainstyle.png", p_balanced, width = 10, height = 3, dpi = 300)