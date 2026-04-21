# =============================================================================
# Ortholog retention analysis (13 species) + Histogram only
# Outputs:
#   1) DEG_with_ortholog_presence.tsv
#   2) ortholog_retention_per_gene.tsv
#   3) ortholog_stats.txt
#   4) Figure5_histogram.png
# =============================================================================


suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(biomaRt)
  library(ggplot2)
})

strip_ver <- function(x) sub("\\.\\d+$", "", x)

# ---- Inputs ------------------------------------------------------------------
deg_file     <- "DEG.csv"  # must include: ensembl_gene_id + category
out_presence <- "DEG_with_ortholog_presence.tsv"

# ---- Fixed 13-species panel --------------------------------------------------
# short codes must match BioMart homolog attribute prefixes
panel <- tibble::tribble(
  ~label,               ~short,
  "Mouse",              "mmusculus",
  "Rat",                "rnorvegicus",
  "Rabbit",             "ocuniculus",
  "Pig",                "sscrofa",
  "Cattle",             "btaurus",
  "Goat",               "chircus",          
  "Sheep",              "oaries",
  "Koala",              "pcinereus",
  "Platypus",           "oanatinus",
  "Chicken",            "ggallus",
  "Xenopus_tropicalis", "xtropicalis",
  "Zebrafish",          "drerio",
  "Medaka",             "olatipes"
)

# ---- Load DEG table -----------------------------------------------------------
deg_tbl <- readr::read_csv(deg_file, show_col_types = FALSE) %>%
  mutate(ensembl_gene_id = strip_ver(ensembl_gene_id))

genes <- unique(deg_tbl$ensembl_gene_id)
stopifnot(length(genes) > 0)

# ---- Connect to Ensembl (human mart) -----------------------------------------
ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# ---- Check available homolog attributes --------------------------------------
attr_names  <- biomaRt::listAttributes(ensembl)$name
avail_short <- sub("_homolog_ensembl_gene$", "", grep("_homolog_ensembl_gene$", attr_names, value = TRUE))

panel <- panel %>% mutate(usable = short %in% avail_short)

if (any(!panel$usable)) {
  message("Dropping unavailable species (not found in BioMart attributes): ",
          paste(panel$label[!panel$usable], collapse = ", "))
}
panel_use <- panel %>% filter(usable)
stopifnot(nrow(panel_use) > 0)

# ---- Presence per species (0/1) ----------------------------------------------
get_presence_one_species <- function(short, label, genes, mart) {
  attr_gene <- sprintf("%s_homolog_ensembl_gene", short)
  
  tb <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", attr_gene),
    filters    = "ensembl_gene_id",
    values     = genes,
    mart       = mart
  ) %>% as.data.table()
  
  setnames(tb, c("ensembl_gene_id", "ortholog_gene"))
  tb[, ensembl_gene_id := strip_ver(ensembl_gene_id)]
  tb[, presence := as.integer(!is.na(ortholog_gene) & ortholog_gene != "")]
  out <- unique(tb[, .(ensembl_gene_id, presence)], by = "ensembl_gene_id")
  
  miss <- setdiff(genes, out$ensembl_gene_id)
  if (length(miss) > 0) {
    out <- rbind(out, data.table(ensembl_gene_id = miss, presence = 0L))
  }
  setnames(out, "presence", label)
  out
}

presence_list <- lapply(seq_len(nrow(panel_use)), function(i) {
  get_presence_one_species(panel_use$short[i], panel_use$label[i], genes, ensembl)
})

presence_wide <- Reduce(function(x, y) merge(x, y, by = "ensembl_gene_id", all = TRUE), presence_list)
for (j in 2:ncol(presence_wide)) presence_wide[[j]][is.na(presence_wide[[j]])] <- 0L

# ---- Join and write presence table -------------------------------------------
final_tbl <- deg_tbl %>%
  left_join(as.data.frame(presence_wide), by = "ensembl_gene_id")

readr::write_tsv(final_tbl, out_presence)

# ---- Per-gene retention ratio -------------------------------------------------
species_cols <- names(presence_wide)[-1]

gene_summary <- final_tbl %>%
  mutate(
    ortholog_ratio = rowSums(across(all_of(species_cols)), na.rm = TRUE) / length(species_cols)
  ) %>%
  select(ensembl_gene_id, category, ortholog_ratio)

readr::write_tsv(gene_summary, "ortholog_retention_per_gene.tsv")

# ---- Stats (saved) ------------------------------------------------------------
kw <- kruskal.test(ortholog_ratio ~ category, data = gene_summary)
pw <- pairwise.wilcox.test(gene_summary$ortholog_ratio, gene_summary$category,
                           p.adjust.method = "bonferroni")

sink("ortholog_stats.txt")
cat("Kruskal–Wallis test\n"); print(kw)
cat("\nPairwise Wilcoxon (Bonferroni)\n"); print(pw)
sink()

# ---- Final figure: histogram only --------------------------------------------
gene_summary$category <- factor(gene_summary$category)

p <- ggplot(gene_summary, aes(x = ortholog_ratio)) +
  geom_histogram(bins = 10, color = "grey40") +
  facet_wrap(~ category, ncol = 1, scales = "free_y") +
  labs(
    x = "Ortholog retention ratio per gene",
    y = "Number of genes",
    title = "Ortholog retention ratio by DEG category"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave("Figure5_histogram.png", p, width = 6.5, height = 7, dpi = 300)

message("Done: 13-species presence + retention ratio + stats + histogram.")






###

# conservation rato
ratio_list <- df2 %>%
  group_by(category) %>%
  summarise(ratios = list(retention_ratio), .groups = "drop")

female_ratios <- ratio_list$ratios[[which(ratio_list$category == "MammaryFemaleOnly")]]
mammary_ratios <- ratio_list$ratios[[which(ratio_list$category == "MammaryGlandOnly")]]
gland_mammary_ratios <- ratio_list$ratios[[which(ratio_list$category == "GlandularPlusMammaryGland")]]

# KS test
ks_female_vs_mammary <- ks.test(female_ratios, mammary_ratios)
ks_female_vs_gland_mammary <- ks.test(female_ratios, gland_mammary_ratios)
ks_mammary_vs_gland_mammary <- ks.test(mammary_ratios, gland_mammary_ratios)

# random
set.seed(1)
n_iter <- 1000

all_retention <- ortho_genes %>%
  mutate(retention_ratio = rowMeans(across(all_of(species_cols)), na.rm = TRUE)) %>%
  pull(retention_ratio)

obs_summary <- df2 %>%
  group_by(category) %>%
  summarise(
    n = n(),
    observed_median = median(retention_ratio, na.rm = TRUE),
    observed_mean = mean(retention_ratio, na.rm = TRUE),
    .groups = "drop"
  )

random_results <- lapply(seq_len(nrow(obs_summary)), function(i) {
  n_i <- obs_summary$n[i]
  cat_i <- obs_summary$category[i]
  
  rand_medians <- replicate(n_iter, {
    median(sample(all_retention, size = n_i, replace = FALSE), na.rm = TRUE)
  })
  
  rand_means <- replicate(n_iter, {
    mean(sample(all_retention, size = n_i, replace = FALSE), na.rm = TRUE)
  })
  
  data.frame(
    category = cat_i,
    iter = seq_len(n_iter),
    random_median = rand_medians,
    random_mean = rand_means
  )
}) %>%
  bind_rows()

empirical_test <- obs_summary %>%
  left_join(
    random_results %>%
      group_by(category) %>%
      summarise(
        empirical_p_median_high = mean(random_median >= first(obs_summary$observed_median[obs_summary$category == category])),
        empirical_p_mean_high   = mean(random_mean   >= first(obs_summary$observed_mean[obs_summary$category == category])),
        random_median_mean = mean(random_median),
        random_mean_mean = mean(random_mean),
        .groups = "drop"
      ),
    by = "category"
  )




random_results$category <- factor(
  random_results$category,
  levels = c(
    "GlandularPlusMammaryGland",
    "MammaryGlandOnly",
    "MammaryFemaleOnly"
  )
)

obs_summary$category <- factor(
  obs_summary$category,
  levels = c(
    "GlandularPlusMammaryGland",
    "MammaryGlandOnly",
    "MammaryFemaleOnly"
  )
)


p_random <- ggplot(random_results, aes(x = category, y = random_median, fill = category)) +
  geom_boxplot(color = "black") +
  geom_point(
    data = obs_summary,
    aes(x = category, y = observed_median),
    shape = 18, size = 3.5, color = "black"
  ) +
  scale_fill_manual(
    values = c(
      GlandularPlusMammaryGland = "#F1605DFF",
      MammaryGlandOnly = "#FEAF77FF",
      MammaryFemaleOnly = "#FCFDBFFF"
    )
  ) +
  scale_x_discrete(
    labels = c(
      GlandularPlusMammaryGland =
        "Epithelial \nsecretory–enriched \n(n = 711)",
      MammaryGlandOnly =
        "Breast-enriched\n(n = 189)",
      MammaryFemaleOnly =
        "Female \nbreast-biased \n(n = 617)"
    )
  ) +
  theme_bw() +
  ylab("Random median retention ratio") +
  theme(
    axis.text.x = element_text(size = 10, lineheight = 0.9),
    panel.grid = element_blank()
  )

