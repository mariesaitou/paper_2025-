# Description: This script generates all figures used in the manuscript 
#   including data wrangling and statistical modeling.



#### 1. Load libraries ####
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(forcats)
library(broom)
library(ComplexUpset)

#### 2. Load and wrangle NanoPlot statistics ####
stats_files <- list.files("nanoplot", pattern = "NanoStats.*\.txt$", full.names = TRUE)

read_stats_file <- function(file) {
  lines <- read_lines(file)
  lines <- lines[!str_detect(lines, "^Metrics dataset")]
  split_lines <- str_split_fixed(lines, "\t+", n = 2)
  tibble(
    metric = str_trim(split_lines[, 1]),
    value = parse_number(str_trim(split_lines[, 2])),
    sample = basename(file) %>%
      str_remove("^NanoStats_?") %>%
      str_remove("_?NanoStats") %>%
      str_remove("\.txt$")
  ) %>%
    filter(metric != "Metrics")
}

stats_long <- map_dfr(stats_files, read_stats_file)

# Assign metadata
metadata <- stats_long %>%
  distinct(sample) %>%
  mutate(
    year = str_extract(sample, "^\d{4}"),
    treatment = case_when(
      str_detect(sample, "Needle") ~ "Needle",
      str_detect(sample, "Vortex") ~ "Vortex",
      str_detect(sample, "Freeze") ~ "Freeze",
      str_detect(sample, "Ctrl|Control") ~ "Control",
      TRUE ~ "Unknown"
    )
  )

# Join metadata and compute log2 enrichment vs control
plot_data <- stats_long %>%
  mutate(filter_status = if_else(str_detect(sample, "filtered"), "Filtered", "Raw"),
         sample_for_join = str_remove(sample, "_filtered$")) %>%
  left_join(metadata, by = c("sample_for_join" = "sample")) %>%
  filter(!is.na(year), !is.na(treatment)) %>%
  mutate(year_filter = paste(year, filter_status, sep = "-"),
         treatment_label = case_when(
           treatment == "Freeze" & year == "2024" ~ "Freeze–Thaw",
           treatment == "Freeze" & year == "2025" ~ "Freeze–Heat",
           TRUE ~ treatment
         ))

ctrl_means <- plot_data %>%
  filter(treatment_label == "Control") %>%
  group_by(year, filter_status, metric) %>%
  summarise(ctrl_mean = mean(value, na.rm = TRUE), .groups = "drop")

plot_enrich <- plot_data %>%
  left_join(ctrl_means, by = c("year", "filter_status", "metric")) %>%
  mutate(log2_enrichment = log2(value / ctrl_mean))

#### 3. Figure 2: Log2 enrichment boxplots ####
metric_labels <- c(
  "mean_qual" = "Mean Quality",
  "number_of_bases" = "Base Numbers",
  "n50" = "N50",
  "number_of_reads" = "Read Count"
)

metric_settings <- list(
  "mean_qual" = list(lim = c(-0.152, 0.137), breaks = c(-0.152, 0, 0.137), labels = c("0.9×", "1×", "1.1×")),
  "number_of_bases" = list(lim = c(-2.1, 3), breaks = c(-1, 0, 1, 1.584), labels = c("0.5×", "1×", "2×", "3×")),
  "n50" = list(lim = c(-1, 1), breaks = c(-1, 0, 1), labels = c("0.5×", "1×", "2×")),
  "number_of_reads" = list(lim = c(-2.1, 4), breaks = c(-1, 0, 1, 1.584, 2), labels = c("0.5×", "1×", "2×", "3×", "4×"))
)

plot_data <- plot_enrich %>%
  filter(metric %in% names(metric_labels)) %>%
  mutate(
    metric_label = factor(metric_labels[metric], levels = metric_labels),
    year_filter_label = str_replace(year_filter, "-", " "),
    treatment_label = factor(treatment_label)
  )

make_sub_plot <- function(metric_id, year_label, show_y = FALSE, show_x = FALSE) {
  data_sub <- plot_data %>%
    filter(metric == metric_id, year_filter_label == year_label)
  setting <- metric_settings[[metric_id]]
  ggplot(data_sub, aes(x = treatment_label, y = log2_enrichment, fill = treatment_label)) +
    geom_boxplot(outlier.shape = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    scale_y_continuous(limits = setting$lim, breaks = setting$breaks, labels = setting$labels) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      x = NULL,
      y = if (show_y) metric_labels[[metric_id]] else NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = if (show_x) element_text(size = 8) else element_blank(),
      axis.title.y = if (show_y) element_text(face = "bold") else element_blank(),
      legend.position = "none"
    )
}

make_column <- function(year_label, show_y = FALSE) {
  plots <- lapply(names(metric_labels), function(m) {
    make_sub_plot(m, year_label, show_y = show_y, show_x = (m == "number_of_reads"))
  })
  wrap_plots(plots, ncol = 1)
}

col1 <- make_column("2024 Raw", show_y = TRUE)
col2 <- make_column("2024 Filtered")
col3 <- make_column("2025 Raw")
col4 <- make_column("2025 Filtered")

header_titles <- c("2024 Raw", "2024 Filtered", "2025 Raw", "2025 Filtered")
headers <- lapply(header_titles, function(title) {
  ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = title, size = 4, fontface = "bold")
})
header_row <- wrap_plots(headers, nrow = 1)

figure2_plot <- header_row / (col1 | col2 | col3 | col4) + 
  plot_layout(heights = c(0.08, 1), guides = "collect") & 
  theme(legend.position = "right")

figure2_plot


#### 4. Figure 3: Forest plot from linear models ####
metrics_for_lm <- c("mean_qual", "n50", "number_of_reads")

lm_results_log2 <- map_dfr(metrics_for_lm, function(m) {
  plot_data %>%
    filter(metric == m, !is.na(value), value > 0) %>%
    mutate(log_value = log2(value)) %>%
    lm(log_value ~ treatment_label + year + filter_status, data = .) %>%
    tidy(conf.int = TRUE) %>%
    mutate(metric = m)
})

forest_data <- lm_results_log2 %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term_clean = case_when(
      term == "treatment_labelNeedle"        ~ "Needle",
      term == "treatment_labelVortex"        ~ "Vortex",
      term == "treatment_labelFreeze–Thaw"   ~ "Freeze–Thaw",
      term == "treatment_labelFreeze–Heat"   ~ "Freeze–Heat",
      term == "filter_statusFiltered"        ~ "Filtered",
      term == "year2025"                     ~ "2025",
      TRUE ~ term
    ),
    metric_clean = case_when(
      metric == "mean_qual"        ~ "Mean Quality",
      metric == "n50"              ~ "Read N50",
      metric == "number_of_reads"  ~ "Read Count",
      TRUE ~ metric
    ),
    dot_size = abs(estimate)
  )

figure3_plot <- ggplot(forest_data, aes(x = estimate, y = term_clean)) +
  geom_point(aes(color = metric_clean, size = dot_size)) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = metric_clean), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ metric_clean, scales = "free_x") +
  scale_size_continuous(range = c(1.5, 6), guide = "none") +
  labs(
    x = "Effect Size (log2 Fold Change)",
    y = NULL,
    title = "Estimated Effects of Fragmentation, Filtering, and Year"
  ) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11),
    legend.position = "none"
  )

figure3_plot


#### 5. Figure 4: SV detection rate by SV length bin and type ####

# Load SV info column from preprocessed VCF data
vcf2024_raw <- read_delim("2024_var_info_8.txt", delim = ";", col_names = FALSE)
vcf2025_raw <- read_delim("2025_var_info_8.txt", delim = ";", col_names = FALSE)

parse_vcf_data <- function(df) {
  df %>%
    mutate(across(everything(), ~ sub(".*?=", "", .))) %>%
    setNames(sub("=.*", "", df[1, ])) %>%
    slice(-1)
}

vcf2024 <- parse_vcf_data(vcf2024_raw)
vcf2025 <- parse_vcf_data(vcf2025_raw)

techs <- c("Ctrl", "Freeze", "Needle", "Vortex", "Gold_standard")

expand_supp_vec <- function(df) {
  df %>%
    mutate(SUPP_VEC = str_pad(SUPP_VEC, width = 5, side = "left", pad = "0")) %>%
    mutate(across(SUPP_VEC, ~ strsplit(., ""))) %>%
    unnest_wider(SUPP_VEC, names_sep = "_") %>%
    rename_with(~ techs, starts_with("SUPP_VEC_")) %>%
    mutate(across(all_of(techs), as.integer))
}

vcf2024 <- vcf2024 %>% expand_supp_vec()
vcf2025 <- vcf2025 %>% expand_supp_vec()

vcf2024 <- vcf2024 %>% mutate(SVLEN = ifelse(SVTYPE == "TRA", 1, abs(as.numeric(SVLEN))))
vcf2025 <- vcf2025 %>% mutate(SVLEN = ifelse(SVTYPE == "TRA", 1, abs(as.numeric(SVLEN))))

detected_2024 <- vcf2024 %>%
  filter(Gold_standard == 1) %>%
  pivot_longer(cols = c(Ctrl, Freeze, Needle, Vortex),
               names_to = "Technology", values_to = "Detected") %>%
  mutate(Year = "2024")

detected_2025 <- vcf2025 %>%
  filter(Gold_standard == 1) %>%
  pivot_longer(cols = c(Ctrl, Freeze, Needle, Vortex),
               names_to = "Technology", values_to = "Detected") %>%
  mutate(Year = "2025")

df_all <- bind_rows(detected_2024, detected_2025)

bin_breaks <- c(0, 100, 1000, 10000, 100000, Inf)
bin_labels <- c("1-100", "101-1k", "1k-10k", "10k-100k", ">100k")

df_all <- df_all %>%
  mutate(SV_bin = cut(SVLEN, breaks = bin_breaks, labels = bin_labels, include.lowest = TRUE))

detection_rate <- df_all %>%
  group_by(Year, Technology, SV_bin, SVTYPE) %>%
  summarise(
    total = n(),
    detected = sum(Detected),
    rate = detected / total,
    .groups = "drop"
  )

shape_palette <- c("Ctrl" = 21, "Freeze" = 22, "Needle" = 23, "Vortex" = 24)
dark_colors <- c("Ctrl" = "black", "Freeze" = "brown", "Needle" = "#377EB8", "Vortex" = "#4DAF4A")

figure4_plot <- ggplot(detection_rate, aes(x = SV_bin, y = rate, shape = Technology, color = Technology)) +
  geom_point(position = position_dodge(width = 0.5),
             size = 3.5, stroke = 1.2, fill = NA) +
  scale_shape_manual(values = shape_palette) +
  scale_color_manual(values = dark_colors) +
  facet_grid(Year ~ SVTYPE) +
  labs(
    title = "Detection Rate by SV Length Bin",
    x = "SV Length Bin",
    y = "Detection Rate"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    legend.position = "right",
    strip.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(color = "black", face = "bold")
  )

figure4_plot


#### 5. Figure 4: Structural Variant Detection Overlaps (UpSet plot) ####

# Convert presence/absence matrix from SUPP_VEC
sets <- c("Ctrl", "Freeze", "Needle", "Vortex", "Gold_standard")

vcf2024[sets] <- vcf2024[sets] == 1
vcf2025[sets] <- vcf2025[sets] == 1

vcf2024$SVTYPE <- as.factor(vcf2024$SVTYPE)
vcf2025$SVTYPE <- as.factor(vcf2025$SVTYPE)

# Compute intersection string
make_intersection_column <- function(df, sets) {
  df %>%
    mutate(intersection = apply(select(., all_of(sets)), 1, function(row) {
      paste0(sets[as.logical(as.integer(row))], collapse = "&")
    }))
}

vcf2024 <- make_intersection_column(vcf2024, sets)
vcf2025 <- make_intersection_column(vcf2025, sets)

# Extract top 10 most frequent intersections
top10_names_2024 <- vcf2024 %>%
  count(intersection, sort = TRUE) %>%
  slice_head(n = 10) %>%
  pull(intersection)

vcf2024_top10 <- vcf2024 %>% filter(intersection %in% top10_names_2024)

# Define color palette
svtype_colors <- c(
  DEL = "#E41A1C", DUP = "#FF7F00", INS = "#4DAF4A",
  INV = "#377EB8", TRA = "#984EA3"
)

# UpSet plot for 2024 (Top10 intersections)
figure4_plot <- upset(
  vcf2024_top10,
  intersect = sets,
  base_annotations = list(
    'Intersection size' = intersection_size(aes(fill = SVTYPE))
  )
) +
  ggtitle("Figure 4: Structural Variant Detection Overlap (2024)") +
  scale_fill_manual(values = svtype_colors)

figure4_plot


#### 6. Figure 5: Detection Rate by SV Length Bin and Type ####

# Convert SVLEN (TRA = 1)
vcf2024 <- vcf2024 %>% mutate(SVLEN = ifelse(SVTYPE == "TRA", 1, abs(as.numeric(SVLEN))))
vcf2025 <- vcf2025 %>% mutate(SVLEN = ifelse(SVTYPE == "TRA", 1, abs(as.numeric(SVLEN))))

# Create detection dataset
get_detected_df <- function(vcf, year) {
  vcf %>%
    filter(Gold_standard == 1) %>%
    pivot_longer(cols = c(Ctrl, Freeze, Needle, Vortex),
                 names_to = "Technology", values_to = "Detected") %>%
    mutate(Year = year)
}

detected_2024 <- get_detected_df(vcf2024, "2024")
detected_2025 <- get_detected_df(vcf2025, "2025")

df_all <- bind_rows(detected_2024, detected_2025)

# SV length bin
bin_breaks <- c(0, 100, 1000, 10000, 100000, Inf)
bin_labels <- c("1-100", "101-1k", "1k-10k", "10k-100k", ">100k")

df_all <- df_all %>%
  mutate(SV_bin = cut(SVLEN, breaks = bin_breaks, labels = bin_labels, include.lowest = TRUE))

# Compute detection rates
detection_rate <- df_all %>%
  group_by(Year, Technology, SV_bin, SVTYPE) %>%
  summarise(
    total = n(),
    detected = sum(Detected),
    rate = detected / total,
    .groups = "drop"
  )

# Plot
shape_palette <- c("Ctrl" = 21, "Freeze" = 22, "Needle" = 23, "Vortex" = 24)
dark_colors <- c("Ctrl" = "black", "Freeze" = "brown", "Needle" = "#377EB8", "Vortex" = "#4DAF4A")

figure5_plot <- ggplot(detection_rate, aes(x = SV_bin, y = rate, shape = Technology, color = Technology)) +
  geom_point(position = position_dodge(width = 0.5),
             size = 3.5, stroke = 1.2, fill = NA) +
  scale_shape_manual(values = shape_palette) +
  scale_color_manual(values = dark_colors) +
  facet_grid(Year ~ SVTYPE) +
  labs(
    title = "Figure 5: Detection Rate by SV Length Bin and SV Type",
    x = "SV Length Bin",
    y = "Detection Rate"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    legend.position = "right",
    strip.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(color = "black", face = "bold")
  )

figure5_plot





library(broom)
library(purrr)
library(dplyr)

metrics_for_lm <- c("mean_qual", "n50", "number_of_reads")

# Linear models with interaction: treatment_label * year
lm_results_interaction <- map_dfr(metrics_for_lm, function(m) {
  plot_data %>%
    filter(metric == m, !is.na(value), value > 0) %>%
    mutate(log_value = log2(value)) %>%
    lm(log_value ~ treatment_label * year + filter_status, data = .) %>%
    tidy(conf.int = TRUE) %>%
    mutate(metric = m)
})
