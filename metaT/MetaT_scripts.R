# R Script to Reproduce Figures and Tables in the MetaT Salmon Paper
# --------------------------------------------------------------
# Figure 2: Taxonomic Composition and Alpha Diversity
# --------------------------------------------------------------

# Required packages
library(tidyverse)
library(vegan)
library(ggplot2)

# Load metadata and taxonomic data
samples_info <- read_csv("samples.csv")
filtered_data <- read_csv("brecken_cleaned.csv")  # Bracken output, family-level

# Normalize taxon counts to make them comparable across samples
scaling_factor <- 200000
normalized_data <- filtered_data %>%
  mutate(across(-taxon_name, ~ round(.x / sum(.x, na.rm = TRUE) * scaling_factor, 0)))

# Compute alpha diversity indices (Shannon and Simpson)
alpha_diversity <- normalized_data %>%
  column_to_rownames("taxon_name") %>%
  t() %>%
  as.data.frame() %>%
  mutate(
    sample = rownames(.),
    Shannon = diversity(., index = "shannon"),
    Simpson = diversity(., index = "simpson")
  ) %>%
  select(sample, Shannon, Simpson) %>%
  inner_join(samples_info, by = c("sample" = "ENA_ID"))

# Reshape for plotting
alpha_diversity_long <- alpha_diversity %>%
  pivot_longer(cols = c(Shannon, Simpson), names_to = "Index", values_to = "Diversity") %>%
  pivot_longer(cols = c(Sex, Status, Tank), names_to = "Variable", values_to = "Group")

# Plot alpha diversity (Figure 2B)
ggplot(alpha_diversity_long, aes(x = Group, y = Diversity, fill = Group)) +
  geom_boxplot() +
  facet_grid(Index ~ Variable, scales = "free", switch = "both") +
  theme_minimal() +
  labs(title = "Alpha Diversity by Group", y = "Diversity Index", x = "") +
  scale_fill_brewer(palette = "Set2")

# Compute relative abundances for family-level barplot (Figure 2A)
relative_abundance <- normalized_data %>%
  pivot_longer(cols = -taxon_name, names_to = "Sample", values_to = "Abundance") %>%
  left_join(samples_info, by = c("Sample" = "ENA_ID")) %>%
  group_by(taxon_name, Tank, Sex, Status) %>%
  summarise(Mean_Abundance = mean(Abundance), .groups = "drop")

# Keep top 10 families based on total abundance
top_families <- normalized_data %>%
  mutate(Total = rowSums(select(., -taxon_name))) %>%
  arrange(desc(Total)) %>%
  slice(1:10) %>%
  pull(taxon_name)

relative_abundance <- relative_abundance %>%
  mutate(taxon_name = ifelse(taxon_name %in% top_families, taxon_name, "Other")) %>%
  group_by(taxon_name, Tank, Sex, Status) %>%
  summarise(Mean_Abundance = sum(Mean_Abundance), .groups = "drop")

# Plot relative taxonomic composition (Figure 2A)
ggplot(relative_abundance, aes(x = interaction(Tank, Sex, Status), y = Mean_Abundance, fill = taxon_name)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  labs(title = "Relative Taxonomic Composition by Group", y = "Relative Abundance (%)", x = "Tank-Sex-Status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")



# --------------------------------------------------------------
# Figure 3: PCA of Microbial Community Structure (Jaccard)
# --------------------------------------------------------------

library(vegan)

# Prepare matrix for Jaccard distance
binary_data <- normalized_data %>%
  mutate(across(-taxon_name, ~ ifelse(. > 0, 1, 0))) %>%
  column_to_rownames("taxon_name") %>%
  t()

# Compute Jaccard distance and perform PCA
jaccard_dist <- vegdist(binary_data, method = "jaccard")
pca_result <- cmdscale(jaccard_dist, k = 2, eig = TRUE)

# Compute variance explained
eig_values <- pca_result$eig
variance_explained <- eig_values / sum(abs(eig_values)) * 100

# Create PCA dataframe
pca_df <- as.data.frame(pca_result$points)
colnames(pca_df) <- c("PC1", "PC2")
pca_df$Sample <- rownames(pca_df)
pca_df <- pca_df %>% left_join(samples_info, by = c("Sample" = "ENA_ID"))

# Plot PCA (Figure 3A)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Status)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, linetype = "dashed") +
  labs(title = "PCA of Microbial Community (Jaccard)",
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)")) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "top")

# --------------------------------------------------------------
# Figure 4: Presence-Absence Histogram and Survival-Associated Taxa
# --------------------------------------------------------------

library(moments)

# Convert to presence-absence and pivot for analysis
master_long <- filtered_data %>%
  pivot_longer(-taxon_name, names_to = "ENA_ID", values_to = "count")

merged_data <- samples_info %>%
  left_join(master_long, by = "ENA_ID") %>%
  mutate(present = ifelse(count > 0, 1, 0))

# Calculate presence ratios per group
status_bacteria_ratio <- merged_data %>%
  group_by(Status, taxon_name) %>%
  summarise(presence_ratio = mean(present), .groups = "drop")

# Pivot for histogram
status_wide <- status_bacteria_ratio %>%
  pivot_wider(names_from = Status, values_from = presence_ratio) %>%
  mutate(difference = `Not Survived` - Survived)

# Plot histogram of presence-absence differences
ggplot(status_wide, aes(x = difference, y = ..count.. + 1, fill = ..x..)) +
  geom_histogram(binwidth = 0.05, color = "black") +
  scale_fill_gradient2(name = "Difference", low = "#0005b4", mid = "white", high = "#b61c1c") +
  theme_minimal() +
  labs(title = "Difference in Bacterial Presence Ratio",
       x = "Presence Ratio Difference (Not Survived - Survived)",
       y = "Count") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")

# Skewness test
skew_value <- skewness(status_wide$difference)
print(paste("Skewness:", skew_value))
t.test(status_wide$difference, mu = 0)


# --------------------------------------------------------------
# Table 2: Top Abundant Families and Genera
# --------------------------------------------------------------

# Calculate mean reads and relative abundance per family
top_families <- normalized_data %>%
  mutate(mean_reads = rowMeans(select(., -taxon_name), na.rm = TRUE)) %>%
  mutate(percent_reads = mean_reads / sum(mean_reads) * 100) %>%
  arrange(desc(mean_reads)) %>%
  slice_head(n = 20)


# Load genus-level Bracken output
genus_data <- read_csv("bracken_genus.csv")

genus_means <- genus_data %>%
  mutate(mean_reads = rowMeans(select(., starts_with("Plate")), na.rm = TRUE)) %>%
  mutate(percent_reads = mean_reads / sum(mean_reads) * 100) %>%
  arrange(desc(mean_reads)) %>%
  slice_head(n = 20)


# --------------------------------------------------------------
# Supplementary Figure S1: Poly(A) Tail Retention
# --------------------------------------------------------------

data <- read.csv("polyA_tail_summary.csv")

data_long <- data %>%
  pivot_longer(cols = starts_with("Percentage_A"),
               names_to = "PolyA_Length",
               values_to = "Ratio") %>%
  mutate(PolyA_Length = as.numeric(gsub("Percentage_A", "", PolyA_Length)),
         Ratio = Ratio / 100)

theoretical <- data.frame(
  PolyA_Length = seq(4, 24, by = 4),
  Probability = (1/4) ^ seq(4, 24, by = 4)
)

ggplot(data_long, aes(x = factor(PolyA_Length), y = Ratio)) +
  geom_boxplot(fill = "gray") +
  geom_line(data = theoretical, aes(x = factor(PolyA_Length), y = Probability, group = 1),
            color = "blue", linetype = "dashed", size = 1) +
  geom_point(data = theoretical, aes(x = factor(PolyA_Length), y = Probability),
             color = "blue", size = 1.5) +
  scale_y_continuous(trans = "sqrt", breaks = c(0, 0.001, 0.01, 0.1, 1), limits = c(0, 1)) +
  theme_minimal() +
  labs(title = "Poly(A) Tail Length Distribution", y = "Proportion of Reads")


# --------------------------------------------------------------
# Supplementary Figure S2: Correlation between MATAM and Bracken
# --------------------------------------------------------------

matam <- read_tsv("contingency_table.tsv")

matam <- matam %>%
  mutate(Taxon = str_extract(Taxonomy, "[^;]+$"),
         SampleID = str_remove(SampleID, "_R[12]_nonmt")) %>%
  group_by(Taxon, SampleID) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  filter(Taxon != "unclassified")

# Assume bracken_wide1 already loaded and processed
bracken_sub <- bracken_wide1 %>%
  filter(taxon_name_genus %in% unique(matam$Taxon)) %>%
  pivot_longer(cols = starts_with("Plate"), names_to = "SampleID", values_to = "Bracken_Abundance") %>%
  rename(Taxon = taxon_name_genus)

merged_data <- matam %>%
  inner_join(bracken_sub, by = c("Taxon", "SampleID")) %>%
  mutate(
    log_MATAM_Abundance = log10(Abundance + 1),
    log_Bracken_Abundance = log10(Bracken_Abundance + 1)
  )

# Calculate correlation
library(broom)
spearman_results <- merged_data %>%
  group_by(Taxon) %>%
  summarise(
    rho = cor(log_MATAM_Abundance, log_Bracken_Abundance, method = "spearman"),
    p_value = cor.test(log_MATAM_Abundance, log_Bracken_Abundance, method = "spearman")$p.value
  )

# Plot
merged_data <- merged_data %>% left_join(spearman_results, by = "Taxon")

ggplot(merged_data, aes(x = log_MATAM_Abundance, y = log_Bracken_Abundance)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", color = "blue") +
  facet_wrap(~ Taxon, scales = "free") +
  theme_minimal() +
  labs(title = "Correlation between MATAM and Bracken Estimates",
       x = "log10(MATAM Abundance + 1)",
       y = "log10(Bracken Abundance + 1)")
