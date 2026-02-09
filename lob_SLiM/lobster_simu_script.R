#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
})

# =========================================================
# Lobster SLiM outputs -> Figure 2-5 + Table 1
# Files expected in base_dir:
#   all.tsv.popN.tsv
#   all.tsv.sizebin.tsv
#   all.tsv.natal.tsv
# Output:
#   fig2_pop_timeseries.png
#   fig3_patch_boxplots_gen140.png
#   fig4_pyramid_gen200_dispersal.png
#   fig5_natal_composition_gen140.png
#   table1_sex_diff_gen200.csv
# =========================================================

# -----------------------------
# Paths / settings
# -----------------------------
base_dir <- "."
base_tsv <- "all.tsv"

pop_path    <- file.path(base_dir, paste0(base_tsv, ".popN.tsv"))
sizebin_path<- file.path(base_dir, paste0(base_tsv, ".sizebin.tsv"))
natal_path  <- file.path(base_dir, paste0(base_tsv, ".natal.tsv"))

out_dir <- file.path(base_dir, "outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fish_gen <- 51L
reg_gen  <- 101L

GEN_FIG3 <- 140L
GEN_FIG5 <- 140L
GEN_FIG4 <- 200L

# -----------------------------
# Style constants
# -----------------------------
space_sea_palette <- c(
  "#DCB13C",  # Yellow Gold
  "#57BDA2",  # Aqua
  "#2493A2",  # Teal
  "#304A78",  # Medium Blue
  "#2C3259",  # Navy Blue
  "#0A0E28"   # Dark Blue
)

cols_open_nt <- c(
  "Open"    = space_sea_palette[3],  # Teal
  "No Take" = space_sea_palette[1]   # Gold
)

# mpa label map: script naming -> paper naming
mpa_map <- c(
  block_up   = "NT_OP",
  block_down = "OP_NT",
  alt_A      = "NONO",
  alt_B      = "ONON"
)
mpa_levels  <- c("NT_OP", "OP_NT", "NONO", "ONON")
cond_levels <- c("Local", "Dispersal", "Unidirectional")

# -----------------------------
# Helpers
# -----------------------------
stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path, call. = FALSE)
}

save_plot <- function(p, filename, w=14, h=7) {
  ggsave(
    filename = file.path(out_dir, filename),
    plot = p,
    width = w, height = h, dpi = 300
  )
}

normalize_mpa <- function(x) {
  # If already paper-style, keep; otherwise recode.
  y <- recode(as.character(x), !!!mpa_map, .default = as.character(x))
  factor(y, levels = mpa_levels)
}

normalize_cond <- function(x) {
  factor(as.character(x), levels = cond_levels)
}

# =========================================================
# Load data (only needed files)
# =========================================================
stop_if_missing(pop_path)
stop_if_missing(sizebin_path)
stop_if_missing(natal_path)

pop    <- read_tsv(pop_path, show_col_types = FALSE)
sizebin<- read_tsv(sizebin_path, show_col_types = FALSE)
natal  <- read_tsv(natal_path, show_col_types = FALSE, comment = "#")

# =========================================================
# FIGURE 2: Total population size over time (Open vs No Take)
# =========================================================
make_fig2 <- function(pop) {
  pop2 <- pop %>%
    mutate(
      OpenNT = case_when(
        str_detect(area, "^open")    ~ "Open",
        str_detect(area, "^no_take") ~ "No Take",
        TRUE ~ NA_character_
      ),
      mpa  = normalize_mpa(mpa),
      cond = normalize_cond(cond)
    ) %>%
    filter(!is.na(OpenNT))
  
  sumdat <- pop2 %>%
    group_by(cond, mpa, generation, rep, OpenNT) %>%
    summarise(N_rep = sum(N), .groups = "drop") %>%
    group_by(cond, mpa, generation, OpenNT) %>%
    summarise(
      N_mean = mean(N_rep),
      N_sd   = sd(N_rep),
      n_rep  = n(),
      N_se   = N_sd / sqrt(n_rep),
      .groups = "drop"
    ) %>%
    mutate(
      ymin = pmax(0, N_mean - N_se),
      ymax = N_mean + N_se
    )
  
  ggplot(sumdat, aes(x = generation, y = N_mean, color = OpenNT, fill = OpenNT)) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.7, color = NA) +
    geom_line(linewidth = 1.1) +
    geom_vline(xintercept = c(fish_gen, reg_gen), linetype = "dashed",
               color = "black", linewidth = 0.6) +
    facet_grid(cond ~ mpa) +
    scale_color_manual(values = cols_open_nt) +
    scale_fill_manual(values = cols_open_nt) +
    labs(x = "Generation", y = "Total population size", color = NULL, fill = NULL) +
    theme_bw() +
    theme(
      legend.position  = "top",
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text       = element_text(color = "black")
    )
}

p2 <- make_fig2(pop)
save_plot(p2, "fig2_pop_timeseries.png", w=16, h=7)

# =========================================================
# FIGURE 3: Patch-level population size distributions (gen 140)
# =========================================================
make_fig3 <- function(pop, gen = 140L) {
  pop2 <- pop %>%
    mutate(
      pos = as.integer(pos),
      generation = as.integer(generation),
      OpenNT = case_when(
        str_detect(area, "^no_take") ~ "No Take",
        str_detect(area, "^open")    ~ "Open",
        TRUE ~ NA_character_
      ),
      mpa  = normalize_mpa(mpa),
      cond = normalize_cond(cond)
    ) %>%
    filter(!is.na(OpenNT), generation == gen)
  
  boxdat <- pop2 %>%
    group_by(cond, mpa, generation, pos, rep, OpenNT) %>%
    summarise(N = sum(N), .groups = "drop")
  
  ggplot(boxdat, aes(x = factor(pos), y = N, fill = OpenNT)) +
    geom_boxplot() +
    facet_grid(cond ~ mpa) +
    scale_fill_manual(values = cols_open_nt) +
    labs(x = "Patch", y = "Population size", fill = "Area") +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      legend.key       = element_blank()
    )
}

p3 <- make_fig3(pop, gen = GEN_FIG3)
save_plot(p3, "fig3_patch_boxplots_gen140.png", w=16, h=7)

# =========================================================
# FIGURE 4: Population pyramid (gen 200, Dispersal only)
# - Alive individuals from SIZEBIN
# - Open vs No-take
# - <25cm grayscale, >=25 colored by sex
# =========================================================
make_fig4 <- function(sizebin, gen = 200L) {
  sb <- sizebin %>%
    mutate(
      area2 = case_when(
        str_detect(area, "^open")    ~ "Open",
        str_detect(area, "^no_take") ~ "No-take",
        area %in% c("open","no_take") ~ ifelse(area == "open","Open","No-take"),
        TRUE ~ NA_character_
      ),
      mpa  = normalize_mpa(mpa),
      cond = normalize_cond(cond),
      sex  = factor(sex, levels = c("F","M")),
      bin_lo = as.numeric(bin_lo),
      bin_hi = as.numeric(bin_hi),
      bin_mid = (bin_lo + bin_hi) / 2
    ) %>%
    filter(!is.na(area2), generation == gen, cond == "Dispersal")
  
  rep_level <- sb %>%
    group_by(cond, mpa, rep, area2, sex, bin_lo, bin_hi, bin_mid) %>%
    summarise(n_rep = sum(n, na.rm = TRUE), .groups = "drop")
  
  pyr_dat <- rep_level %>%
    group_by(cond, mpa, area2, sex, bin_lo, bin_hi, bin_mid) %>%
    summarise(n_mean = mean(n_rep), .groups = "drop") %>%
    mutate(
      n_plot = ifelse(sex == "M", -n_mean, n_mean),
      size_class = ifelse(bin_mid < 25, "low", "high"),
      fill_group = factor(
        paste0(sex, "_", size_class),
        levels = c("F_low","M_low","F_high","M_high")
      )
    )
  
  ggplot(pyr_dat, aes(x = bin_mid, y = n_plot, fill = fill_group)) +
    geom_col(width = 4.5, color = NA) +
    geom_hline(yintercept = 0, linewidth = 0.9) +
    geom_vline(xintercept = 25, linetype = "dashed", linewidth = 0.8) +
    facet_grid(area2 ~ mpa, scales = "free_y") +
    coord_flip() +
    scale_y_continuous(labels = function(x) abs(x)) +
    scale_fill_manual(
      values = c(
        F_low  = "grey70",
        M_low  = "grey85",
        F_high = "#E64B35",
        M_high = "#4DBBD5"
      ),
      labels = c(
        F_low  = "Female (<25 cm)",
        M_low  = "Male (<25 cm)",
        F_high = "Female (>=25 cm)",
        M_high = "Male (>=25 cm)"
      )
    ) +
    labs(
      x = "Size bin midpoint (cm)",
      y = "Individuals (rep-mean; males left, females right)",
      fill = NULL,
      title = paste0("Population pyramid (alive, sizebin) at generation ", gen),
      subtitle = "Dispersal scenario; <25 cm in grayscale"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey80", linewidth = 0.4),
      strip.background = element_rect(fill = "grey90"),
      legend.position = "top"
    )
}

p4 <- make_fig4(sizebin, gen = GEN_FIG4)
save_plot(p4, "fig4_pyramid_gen200_dispersal.png", w=14, h=8)

# =========================================================
# FIGURE 5: Natal origin composition (gen 140)
# - Stacked proportions within current patch
# - bornPos 0 and 11 pooled as Outside
# =========================================================
make_fig5 <- function(natal, gen = 140L) {
  natal_sum <- natal %>%
    filter(generation == gen) %>%
    group_by(cond, mpa, currPos, bornPos) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    group_by(cond, mpa, currPos) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    mutate(
      cond = normalize_cond(cond),
      mpa  = normalize_mpa(mpa),
      born_grp = case_when(
        bornPos %in% c(0, 11) ~ "Outside",
        TRUE ~ as.character(bornPos)
      ),
      born_grp_f = factor(born_grp, levels = c(as.character(1:10), "Outside"))
    )
  
  # keep "Outside" white; others default ggplot hues (no external palette needed)
  ggplot(natal_sum, aes(x = factor(currPos), y = prop, fill = born_grp_f)) +
    geom_col(width = 0.8, color = "grey30") +
    facet_grid(cond ~ mpa) +
    scale_y_continuous(expand = c(0, 0)) +
    guides(fill = guide_legend(ncol = 1)) +
    labs(x = "Current patch", y = "Proportion", fill = "Natal patch") +
    theme_bw() +
    theme(panel.grid = element_blank())
}

p5 <- make_fig5(natal, gen = GEN_FIG5)
save_plot(p5, "fig5_natal_composition_gen140.png", w=16, h=9)

# =========================================================
# TABLE 1: Sex-specific differences in abundance at gen 200
# - size classes: <25 and 25-35
# - one-sided Wilcoxon signed-rank on (M - F), alternative "less"
# - BH FDR; keep FDR < 0.01
# =========================================================
make_table1 <- function(sizebin, gen = 200L) {
  sb <- sizebin %>%
    mutate(
      area2 = case_when(
        str_detect(area, "^open")    ~ "open",
        str_detect(area, "^no_take") ~ "no_take",
        area %in% c("open","no_take") ~ as.character(area),
        TRUE ~ NA_character_
      ),
      cond = normalize_cond(cond),
      mpa  = normalize_mpa(mpa),
      bin_lo = as.numeric(bin_lo)
    ) %>%
    filter(!is.na(area2), generation == gen) %>%
    mutate(
      size_group = case_when(
        bin_lo < 25 ~ "<25",
        bin_lo >= 25 & bin_lo < 35 ~ "25-35",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(size_group))
  
  dat_sum <- sb %>%
    group_by(cond, mpa, area2, rep, size_group, sex) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    pivot_wider(names_from = sex, values_from = n, values_fill = 0) %>%
    mutate(diff = M - F)
  
  wilcox_res <- dat_sum %>%
    group_by(cond, mpa, area2, size_group) %>%
    summarise(
      n_rep = n(),
      median_diff = median(diff),
      p_value = wilcox.test(diff, alternative = "less", exact = FALSE)$p.value,
      .groups = "drop"
    ) %>%
    mutate(
      direction = case_when(
        median_diff < 0 ~ "M < F",
        median_diff > 0 ~ "M > F",
        TRUE ~ "M = F"
      ),
      p_fdr = p.adjust(p_value, method = "BH")
    ) %>%
    filter(p_fdr < 0.01) %>%
    arrange(size_group, cond, mpa, area2) %>%
    select(cond, mpa, area2, size_group, median_diff, direction, p_fdr)
  
  wilcox_res
}

tab1 <- make_table1(sizebin, gen = GEN_FIG4)
write_csv(tab1, file.path(out_dir, "table1_sex_diff_gen200.csv"))

message("Done. Outputs in: ", normalizePath(out_dir))
