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
