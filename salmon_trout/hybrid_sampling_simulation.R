#!/usr/bin/env Rscript

## fixed_diagnostic_snp_hybrid_sampling_simulation.R
##
## Minimal simulation for fixed diagnostic SNPs.
##
## Genotypes are encoded as Trout-allele dosage:
##   0 = Salmon/Salmon
##   1 = Salmon/Trout
##   2 = Trout/Trout
##
## Main outputs:
##   hybrid_class_definitions.csv
##   hybrid_sampling_replicate_results.csv
##   strict_classification_accuracy_by_sampling_depth.png
##
## There is no scenario-level summary table.
## hybrid_sampling_replicate_results.csv contains one row per simulated sampled set
## (= one parameter combination x one replicate).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(readr)
  library(ggplot2)
})

## ---------------------------------------------------------------------------
## User settings
## ---------------------------------------------------------------------------

out_dir <- "fixed_diagnostic_snp_hybrid_sampling_outputs"

n_snps_grid <- c(5, 10, 20, 40)
n_sample_grid <- c(10, 20, 50, 100, 200)
n_replicates <- 10

## Settings for the single figure.
## The figure shows strict classification accuracy as a function of the
## number of individuals resampled from each true hybrid class.
n_per_class_grid <- c(5, 10, 20, 50, 100, 200)
n_snps_for_strict_accuracy_plot <- 20
strict_accuracy_classes <- c(
  "BC2_to_Salmon",
  "BC1_to_Salmon",
  "F1",
  "BC1_to_Trout",
  "BC2_to_Trout"
)

## Probability that hybridisation has occurred in the sampled population.
## If this is 0, the sampled set contains no hybrids.
## If this is 1, hybridisation is always possible at the individual level.
p_pop_hybridisation_grid <- c(1)

## Individual-level probability of being hybrid, conditional on population-level
## hybridisation having occurred. This is the main 0--1 parameter controlling
## hybrid prevalence among sampled fish.
p_individual_hybrid_given_pop_grid <- c(
  0, 0.001, 0.005, 0.01, 0.025, 0.05,
  0.10, 0.25, 0.50, 0.75, 1.00
)

## Among non-hybrid individuals, probability of belonging to the Trout reference
## class. 0.5 = Salmon and Trout reference classes are equally likely among
## non-hybrids. 0 = all non-hybrids are Salmon. 1 = all non-hybrids are Trout.
p_nonhybrid_is_trout <- 0.5

## Hybrid classes. The ancestry columns are written as decimals.
## weight controls the composition among hybrid individuals and is normalised internally.
## F1 and F2 both have Salmon=0.5 and Trout=0.5, but their genotype models differ.
hybrid_scheme <- tribble(
  ~class_id,        ~salmon_ancestry, ~trout_ancestry, ~genotype_model, ~bc_generation, ~weight,
  "F1",             0.5000,           0.5000,          "F1",           NA_integer_,    1,
  "F2",             0.5000,           0.5000,          "F2",           NA_integer_,    1,
  "BC1_to_Salmon",  0.7500,           0.2500,          "BC_to_Salmon", 1L,             1,
  "BC2_to_Salmon",  0.8750,           0.1250,          "BC_to_Salmon", 2L,             1,
  "BC3_to_Salmon",  0.9375,           0.0625,          "BC_to_Salmon", 3L,             1,
  "BC1_to_Trout",   0.2500,           0.7500,          "BC_to_Trout",  1L,             1,
  "BC2_to_Trout",   0.1250,           0.8750,          "BC_to_Trout",  2L,             1,
  "BC3_to_Trout",   0.0625,           0.9375,          "BC_to_Trout",  3L,             1
)

## Optional data imperfections.
missing_rate <- 0.00
genotyping_error_rate <- 0.00

## Hybrid call threshold using q_trout_hat = mean(Trout-allele dosage / 2).
## Individuals between these two reference-parent zones are called hybrid.
hybrid_q_trout_lower <- 0.0625
hybrid_q_trout_upper <- 0.9375

seed <- 1

## ---------------------------------------------------------------------------
## Class scheme
## ---------------------------------------------------------------------------

format_prop <- function(x) {
  sub("\\.?0+$", "", sprintf("%.4f", x))
}

make_ancestry_ratio <- function(salmon_ancestry, trout_ancestry) {
  paste0(
    "Salmon ", format_prop(salmon_ancestry),
    " : Trout ", format_prop(trout_ancestry)
  )
}

make_class_scheme <- function(hybrid_scheme) {
  nonhybrid_scheme <- tribble(
    ~class_id,             ~salmon_ancestry, ~trout_ancestry, ~genotype_model,      ~bc_generation, ~weight,
    "nonhybrid_Salmon",    1.0000,           0.0000,          "nonhybrid_Salmon",  NA_integer_,    NA_real_,
    "nonhybrid_Trout",     0.0000,           1.0000,          "nonhybrid_Trout",   NA_integer_,    NA_real_
  )

  bind_rows(nonhybrid_scheme, hybrid_scheme) |>
    mutate(
      expected_ancestry_ratio = make_ancestry_ratio(salmon_ancestry, trout_ancestry),
      true_is_hybrid = !class_id %in% c("nonhybrid_Salmon", "nonhybrid_Trout"),
      expected_q_trout = trout_ancestry,
      expected_heterozygosity = case_when(
        genotype_model %in% c("nonhybrid_Salmon", "nonhybrid_Trout") ~ 0,
        genotype_model == "F1" ~ 1,
        genotype_model == "F2" ~ 0.5,
        genotype_model %in% c("BC_to_Salmon", "BC_to_Trout") ~ 0.5^bc_generation,
        TRUE ~ NA_real_
      ),
      probability_among_hybrids = if_else(
        true_is_hybrid,
        weight / sum(weight[true_is_hybrid], na.rm = TRUE),
        NA_real_
      )
    ) |>
    select(
      class_id,
      expected_ancestry_ratio,
      salmon_ancestry,
      trout_ancestry,
      genotype_model,
      bc_generation,
      true_is_hybrid,
      expected_q_trout,
      expected_heterozygosity,
      weight,
      probability_among_hybrids
    )
}

## ---------------------------------------------------------------------------
## Genotype simulation
## ---------------------------------------------------------------------------

simulate_true_genotypes <- function(class_id, n_individuals, n_snps, class_scheme) {
  if (n_individuals == 0) {
    return(matrix(numeric(0), nrow = 0, ncol = n_snps))
  }

  cls <- class_scheme |>
    filter(.data$class_id == !!class_id) |>
    slice(1)

  if (nrow(cls) != 1) stop("Unknown class_id: ", class_id)

  genotype_model <- cls$genotype_model
  generation <- cls$bc_generation

  if (genotype_model == "nonhybrid_Salmon") {
    return(matrix(0L, nrow = n_individuals, ncol = n_snps))
  }

  if (genotype_model == "nonhybrid_Trout") {
    return(matrix(2L, nrow = n_individuals, ncol = n_snps))
  }

  if (genotype_model == "F1") {
    return(matrix(1L, nrow = n_individuals, ncol = n_snps))
  }

  if (genotype_model == "F2") {
    return(matrix(
      rbinom(n_individuals * n_snps, size = 2, prob = 0.5),
      nrow = n_individuals,
      ncol = n_snps
    ))
  }

  if (genotype_model == "BC_to_Salmon") {
    return(matrix(
      rbinom(n_individuals * n_snps, size = 1, prob = 0.5^generation),
      nrow = n_individuals,
      ncol = n_snps
    ))
  }

  if (genotype_model == "BC_to_Trout") {
    return(matrix(
      1L + rbinom(n_individuals * n_snps, size = 1, prob = 1 - 0.5^generation),
      nrow = n_individuals,
      ncol = n_snps
    ))
  }

  stop("Unknown genotype_model: ", genotype_model)
}

apply_genotyping_error <- function(G, error_rate) {
  if (error_rate <= 0) return(G)
  if (error_rate >= 0.5) stop("genotyping_error_rate should be < 0.5.")

  dim_G <- dim(G)
  g <- as.integer(as.vector(G))

  observed <- rbinom(length(g), size = g, prob = 1 - error_rate) +
    rbinom(length(g), size = 2L - g, prob = error_rate)

  matrix(observed, nrow = dim_G[1], ncol = dim_G[2])
}

apply_missingness <- function(G, missing_rate) {
  if (missing_rate <= 0) return(G)
  if (missing_rate >= 1) stop("missing_rate should be < 1.")

  G <- matrix(as.numeric(G), nrow = nrow(G), ncol = ncol(G))
  G[matrix(runif(length(G)) < missing_rate, nrow = nrow(G), ncol = ncol(G))] <- NA_real_
  G
}

simulate_genotypes_for_classes <- function(class_ids, n_snps, class_scheme) {
  G <- matrix(NA_real_, nrow = length(class_ids), ncol = n_snps)

  for (cls in unique(class_ids)) {
    idx <- which(class_ids == cls)
    G[idx, ] <- simulate_true_genotypes(cls, length(idx), n_snps, class_scheme) |>
      apply_genotyping_error(genotyping_error_rate) |>
      apply_missingness(missing_rate)
  }

  G
}

sample_class_vector <- function(n_sample,
                                p_pop_hybridisation,
                                p_individual_hybrid_given_pop,
                                class_scheme) {
  pop_hybridisation_occurred <-
    rbinom(1, size = 1, prob = p_pop_hybridisation) == 1

  is_hybrid <- if (pop_hybridisation_occurred) {
    rbinom(n_sample, size = 1, prob = p_individual_hybrid_given_pop) == 1
  } else {
    rep(FALSE, n_sample)
  }

  class_ids <- character(n_sample)

  if (any(is_hybrid)) {
    hybrid_probs <- class_scheme |>
      filter(true_is_hybrid, weight > 0) |>
      transmute(class_id, probability = weight / sum(weight))

    class_ids[is_hybrid] <- sample(
      hybrid_probs$class_id,
      size = sum(is_hybrid),
      replace = TRUE,
      prob = hybrid_probs$probability
    )
  }

  if (any(!is_hybrid)) {
    class_ids[!is_hybrid] <- if_else(
      rbinom(sum(!is_hybrid), size = 1, prob = p_nonhybrid_is_trout) == 1,
      "nonhybrid_Trout",
      "nonhybrid_Salmon"
    )
  }

  list(
    class_ids = class_ids,
    pop_hybridisation_occurred = pop_hybridisation_occurred
  )
}

## ---------------------------------------------------------------------------
## Detection and one-row-per-sampled-set result
## ---------------------------------------------------------------------------

estimate_q_trout <- function(G) {
  rowMeans(G / 2, na.rm = TRUE)
}

call_hybrid_by_q <- function(q_trout_hat) {
  !is.na(q_trout_hat) &
    q_trout_hat > hybrid_q_trout_lower &
    q_trout_hat < hybrid_q_trout_upper
}

make_class_count_string <- function(class_ids) {
  tibble(class_id = class_ids) |>
    count(class_id, name = "n") |>
    arrange(class_id) |>
    mutate(text = paste0(class_id, "=", n)) |>
    pull(text) |>
    paste(collapse = "; ")
}

run_one_replicate <- function(n_snps,
                              n_sample,
                              p_pop_hybridisation,
                              p_individual_hybrid_given_pop,
                              replicate,
                              class_scheme) {
  sampled <- sample_class_vector(
    n_sample = n_sample,
    p_pop_hybridisation = p_pop_hybridisation,
    p_individual_hybrid_given_pop = p_individual_hybrid_given_pop,
    class_scheme = class_scheme
  )

  G <- simulate_genotypes_for_classes(sampled$class_ids, n_snps, class_scheme)
  q_trout_hat <- estimate_q_trout(G)
  hybrid_called_by_q <- call_hybrid_by_q(q_trout_hat)

  individuals <- tibble(
    class_id = sampled$class_ids,
    q_trout_hat = q_trout_hat,
    hybrid_called_by_q = hybrid_called_by_q
  ) |>
    left_join(
      class_scheme |> select(class_id, true_is_hybrid),
      by = "class_id"
    )

  true_n_hybrid <- sum(individuals$true_is_hybrid)
  true_n_nonhybrid <- n_sample - true_n_hybrid
  detected_true_hybrid_n_by_q <- sum(individuals$true_is_hybrid & individuals$hybrid_called_by_q)
  false_positive_n_by_q <- sum(!individuals$true_is_hybrid & individuals$hybrid_called_by_q)
  false_negative_n_by_q <- sum(individuals$true_is_hybrid & !individuals$hybrid_called_by_q)

  tibble(
    replicate = replicate,
    n_snps = n_snps,
    n_sample = n_sample,
    p_pop_hybridisation = p_pop_hybridisation,
    p_individual_hybrid_given_pop = p_individual_hybrid_given_pop,
    p_individual_hybrid_unconditional = p_pop_hybridisation * p_individual_hybrid_given_pop,
    expected_n_hybrid = n_sample * p_pop_hybridisation * p_individual_hybrid_given_pop,
    prob_sample_contains_at_least_one_hybrid_theory =
      p_pop_hybridisation * (1 - (1 - p_individual_hybrid_given_pop)^n_sample),
    pop_hybridisation_occurred = sampled$pop_hybridisation_occurred,
    true_n_hybrid = true_n_hybrid,
    true_n_nonhybrid = true_n_nonhybrid,
    true_any_hybrid = true_n_hybrid > 0,
    detected_hybrid_n_by_q = sum(hybrid_called_by_q),
    detected_true_hybrid_n_by_q = detected_true_hybrid_n_by_q,
    false_positive_n_by_q = false_positive_n_by_q,
    false_negative_n_by_q = false_negative_n_by_q,
    detected_any_true_hybrid_by_q = detected_true_hybrid_n_by_q > 0,
    any_false_positive_by_q = false_positive_n_by_q > 0,
    mean_q_trout_hat = mean(q_trout_hat, na.rm = TRUE),
    min_q_trout_hat = min(q_trout_hat, na.rm = TRUE),
    max_q_trout_hat = max(q_trout_hat, na.rm = TRUE),
    class_counts = make_class_count_string(sampled$class_ids)
  )
}


## ---------------------------------------------------------------------------
## Strict classification accuracy plot
## ---------------------------------------------------------------------------

estimate_heterozygosity <- function(G) {
  rowMeans(G == 1, na.rm = TRUE)
}

classify_by_q_and_heterozygosity <- function(G,
                                             class_scheme,
                                             candidate_classes = class_scheme$class_id) {
  q_trout_hat <- estimate_q_trout(G)
  heterozygosity_hat <- estimate_heterozygosity(G)

  candidate_scheme <- class_scheme |>
    filter(class_id %in% candidate_classes) |>
    select(class_id, expected_q_trout, expected_heterozygosity)

  map2_chr(q_trout_hat, heterozygosity_hat, function(q, h) {
    if (is.na(q) || is.na(h)) return(NA_character_)

    candidate_scheme |>
      mutate(
        distance = sqrt(
          (q - expected_q_trout)^2 +
            (h - expected_heterozygosity)^2
        )
      ) |>
      arrange(distance, class_id) |>
      slice(1) |>
      pull(class_id)
  })
}

run_one_strict_accuracy_replicate <- function(true_class,
                                              n_per_class,
                                              replicate,
                                              class_scheme) {
  G <- simulate_true_genotypes(
    class_id = true_class,
    n_individuals = n_per_class,
    n_snps = n_snps_for_strict_accuracy_plot,
    class_scheme = class_scheme
  ) |>
    apply_genotyping_error(genotyping_error_rate) |>
    apply_missingness(missing_rate)

  predicted_class <- classify_by_q_and_heterozygosity(
    G = G,
    class_scheme = class_scheme
  )

  tibble(
    true_class = true_class,
    n_per_class = n_per_class,
    n_snps = n_snps_for_strict_accuracy_plot,
    replicate = replicate,
    strict_accuracy = mean(predicted_class == true_class, na.rm = TRUE)
  )
}

simulate_strict_classification_accuracy <- function(class_scheme) {
  crossing(
    true_class = strict_accuracy_classes,
    n_per_class = n_per_class_grid,
    replicate = seq_len(n_replicates)
  ) |>
    pmap_dfr(function(true_class, n_per_class, replicate) {
      run_one_strict_accuracy_replicate(
        true_class = true_class,
        n_per_class = n_per_class,
        replicate = replicate,
        class_scheme = class_scheme
      )
    })
}

write_strict_accuracy_plot <- function(class_scheme) {
  accuracy_results <- simulate_strict_classification_accuracy(class_scheme)

  plot_data <- accuracy_results |>
    group_by(true_class, n_per_class, n_snps) |>
    summarise(
      median_accuracy = median(strict_accuracy, na.rm = TRUE),
      q25_accuracy = quantile(strict_accuracy, 0.25, na.rm = TRUE, names = FALSE),
      q75_accuracy = quantile(strict_accuracy, 0.75, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    ) |>
    mutate(
      true_class = factor(true_class, levels = strict_accuracy_classes)
    )

  p <- ggplot(
    plot_data,
    aes(x = n_per_class, y = median_accuracy, group = true_class)
  ) +
    geom_ribbon(
      aes(ymin = q25_accuracy, ymax = q75_accuracy),
      alpha = 0.2
    ) +
    geom_line() +
    geom_point() +
    facet_wrap(~ true_class, nrow = 1) +
    scale_x_sqrt(breaks = n_per_class_grid) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = "Strict classification accuracy across sampling depth",
      subtitle = paste0(
        "Classification based on q_trout_hat and heterozygosity; ",
        n_snps_for_strict_accuracy_plot,
        " fixed diagnostic SNPs"
      ),
      x = "Number of individuals sampled per class (N per class)",
      y = "Strict classification accuracy"
    ) +
    theme_bw()

  ggsave(
    filename = file.path(out_dir, "strict_classification_accuracy_by_sampling_depth.png"),
    plot = p,
    width = 13,
    height = 5,
    dpi = 300
  )
}

## ---------------------------------------------------------------------------
## Run
## ---------------------------------------------------------------------------

set.seed(seed)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

class_scheme <- make_class_scheme(hybrid_scheme)

scenario_grid <- crossing(
  n_snps = n_snps_grid,
  n_sample = n_sample_grid,
  p_pop_hybridisation = p_pop_hybridisation_grid,
  p_individual_hybrid_given_pop = p_individual_hybrid_given_pop_grid,
  replicate = seq_len(n_replicates)
) |>
  mutate(
    simulation_id = row_number(),
    scenario_id = paste0(
      "snps_", n_snps,
      "_sample_", n_sample,
      "_pop_", p_pop_hybridisation,
      "_ind_", p_individual_hybrid_given_pop
    )
  )

message("Running ", nrow(scenario_grid), " simulated sampled sets...")

sample_results <- pmap_dfr(
  scenario_grid,
  function(n_snps,
           n_sample,
           p_pop_hybridisation,
           p_individual_hybrid_given_pop,
           replicate,
           simulation_id,
           scenario_id) {
    run_one_replicate(
      n_snps = n_snps,
      n_sample = n_sample,
      p_pop_hybridisation = p_pop_hybridisation,
      p_individual_hybrid_given_pop = p_individual_hybrid_given_pop,
      replicate = replicate,
      class_scheme = class_scheme
    ) |>
      mutate(
        simulation_id = simulation_id,
        scenario_id = scenario_id,
        .before = replicate
      )
  }
)

write_csv(class_scheme, file.path(out_dir, "hybrid_class_definitions.csv"))
write_csv(sample_results, file.path(out_dir, "hybrid_sampling_replicate_results.csv"))
write_strict_accuracy_plot(class_scheme)

message("Done. Outputs written to: ", normalizePath(out_dir))
