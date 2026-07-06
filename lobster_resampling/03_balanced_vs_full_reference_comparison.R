suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

FULL_TSV <- "snapclust_resampling_all_runs.tsv"
BAL_TSV  <- "snapclust_resampling_balanced_parental_reference_all_runs.tsv"

full <- read_tsv(FULL_TSV, show_col_types = FALSE) |>
  mutate(analysis = "Full reference") |>
  select(-any_of(c("acc_strict", "acc_nearest")))

balanced <- read_tsv(BAL_TSV, show_col_types = FALSE) |>
  mutate(analysis = "Balanced reference") |>
  select(-any_of(c("acc_strict", "acc_nearest")))

class_levels <- c("0", "0.125", "0.25", "0.5", "0.75", "0.875", "1")

class_labels <- c(
  "0"     = "Europe",
  "0.125" = "BC2 (Europe)",
  "0.25"  = "BC1 (Europe)",
  "0.5"   = "F1",
  "0.75"  = "BC1 (America)",
  "0.875" = "BC2 (America)",
  "1"     = "America"
)

common_keys <- inner_join(
  full |> distinct(n_per_class, seed),
  balanced |> distinct(n_per_class, seed),
  by = c("n_per_class", "seed")
)

df <- bind_rows(full, balanced) |>
  semi_join(common_keys, by = c("n_per_class", "seed")) |>
  filter(metric %in% class_levels) |>
  mutate(
    class_label = factor(metric, levels = class_levels),
    N = as.integer(n_per_class),
    analysis = factor(analysis, levels = c("Full reference", "Balanced reference"))
  )

df_summary <- df |>
  group_by(analysis, class_label, N) |>
  summarise(
    strict_med = median(strict_acc, na.rm = TRUE),
    strict_q25 = quantile(strict_acc, 0.25, na.rm = TRUE),
    strict_q75 = quantile(strict_acc, 0.75, na.rm = TRUE),
    post_med = median(median_max_post, na.rm = TRUE),
    post_q25 = quantile(median_max_post, 0.25, na.rm = TRUE),
    post_q75 = quantile(median_max_post, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

base_theme <- theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
    strip.text = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text( hjust = 1),
    legend.position = "bottom"
  )

p_strict <- ggplot(
  df_summary,
  aes(x = N, y = strict_med, group = analysis, linetype = analysis, fill = analysis)
) +
  geom_ribbon(aes(ymin = strict_q25, ymax = strict_q75), alpha = 0.20, color = NA) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.8) +
  facet_wrap(~ class_label, ncol = 7, labeller = as_labeller(class_labels)) +
  scale_x_continuous(
    trans = "sqrt",
    breaks = sort(unique(df_summary$N))
  ) +
  coord_cartesian(ylim = c(0, 1.05)) +
  labs(
    x = NULL,
    y = "Strict accuracy\n(median; IQR)",
    linetype = "Analysis",
    fill = "Analysis",
    tag = "A"
  ) +
  base_theme

p_post <- ggplot(
  df_summary,
  aes(x = N, y = post_med, group = analysis, linetype = analysis, fill = analysis)
) +
  geom_ribbon(aes(ymin = post_q25, ymax = post_q75), alpha = 0.20, color = NA) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.8) +
  facet_wrap(~ class_label, ncol = 7, labeller = as_labeller(class_labels)) +
  scale_x_continuous(
    trans = "sqrt",
    breaks = sort(unique(df_summary$N))
  ) +
  coord_cartesian(ylim = c(0, 1.05)) +
  labs(
    x = "N per class (sqrt scale)",
    y = "Median maximum\nposterior probability\n(median; IQR)",
    linetype = "Analysis",
    fill = "Analysis",
    tag = "B"
  ) +
  base_theme

p_s3 <- p_strict / p_post +
  plot_layout(guides = "collect", heights = c(1, 1)) &
  theme(
    legend.position = "bottom",
    plot.tag = element_text(face = "bold", size = 14)
  )

print(p_s3)

ggsave(
  "Figure_S3_balanced_vs_full_direct_comparison.pdf",
  p_s3,
  width = 15,
  height = 7.5
)

ggsave(
  "Figure_S3_balanced_vs_full_direct_comparison.png",
  p_s3,
  width = 15,
  height = 7.5,
  dpi = 300
)

paired_delta <- full |>
  select(
    n_per_class,
    seed,
    metric,
    strict_acc_full = strict_acc,
    median_max_post_full = median_max_post
  ) |>
  inner_join(
    balanced |>
      select(
        n_per_class,
        seed,
        metric,
        strict_acc_balanced = strict_acc,
        median_max_post_balanced = median_max_post
      ),
    by = c("n_per_class", "seed", "metric")
  ) |>
  filter(metric %in% class_levels) |>
  mutate(
    delta_strict_acc = strict_acc_balanced - strict_acc_full,
    delta_median_max_post = median_max_post_balanced - median_max_post_full
  )

paired_delta_summary <- paired_delta |>
  group_by(n_per_class, metric) |>
  summarise(
    median_delta_strict_acc = median(delta_strict_acc, na.rm = TRUE),
    q25_delta_strict_acc = quantile(delta_strict_acc, 0.25, na.rm = TRUE),
    q75_delta_strict_acc = quantile(delta_strict_acc, 0.75, na.rm = TRUE),
    median_delta_median_max_post = median(delta_median_max_post, na.rm = TRUE),
    q25_delta_median_max_post = quantile(delta_median_max_post, 0.25, na.rm = TRUE),
    q75_delta_median_max_post = quantile(delta_median_max_post, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(paired_delta, "balanced_vs_full_paired_delta_all_runs.tsv")
write_tsv(paired_delta_summary, "balanced_vs_full_paired_delta_summary.tsv")