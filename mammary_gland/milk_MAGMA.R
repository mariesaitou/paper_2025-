############################################################
# MAGMA (GWAS Atlas) enrichment analysis for DEG categories
# Outputs:
#   1) magma_best_trait_per_gene.tsv  : best trait per gene (FDR < q_threshold)
#   2) magma_enrichment_all.tsv       : enrichment results for all traits
#   3) magma_enrichment_significant.tsv: subset (FDR < 0.05, positive enrichment)
############################################################


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
})

# ---- Inputs ------------------------------------------------------------------
# Must include: ensembl_gene_id, category
deg_category_file <- "DEG_evo.csv"  # you already write this in your current pipeline

# GWAS Atlas metadata + MAGMA gene-level p-values
info_file  <- "gwasATLAS_v20191115.txt.gz"
magma_file <- "gwasATLAS_v20191115_magma_P.txt.gz"

# ---- Parameters --------------------------------------------------------------
q_threshold <- 1e-4     # same as your script
bg_total    <- 1000L    # random background gene count
set.seed(20251031)

# ---- Helpers -----------------------------------------------------------------
strip_version <- function(x) sub("\\.\\d+$", "", x)

# ---- 1) Load DEG categories --------------------------------------------------
deg_cat <- readr::read_csv(deg_category_file, show_col_types = FALSE) %>%
  mutate(ensembl_gene_id = strip_version(ensembl_gene_id)) %>%
  select(ensembl_gene_id, category) %>%
  distinct()

stopifnot(nrow(deg_cat) > 0)

# ---- 2) Load GWAS Atlas info + MAGMA results --------------------------------
info  <- readr::read_tsv(info_file, show_col_types = FALSE) %>%
  mutate(id = as.integer(id))

magma <- readr::read_tsv(magma_file, show_col_types = FALSE) %>%
  as_tibble()

# MAGMA gene IDs: strip Ensembl version
magma <- magma %>% mutate(GENE = strip_version(GENE))

# ---- 3) Long format + FDR per trait ------------------------------------------
magma_long <- magma %>%
  pivot_longer(-GENE, names_to = "id", values_to = "P") %>%
  mutate(id = as.integer(id)) %>%
  inner_join(select(info, id, SubchapterLevel), by = "id") %>%
  mutate(FDR = p.adjust(P, method = "fdr"))

# ---- 4) Best trait per gene (FDR < q_threshold) ------------------------------
best_trait_per_gene <- magma_long %>%
  filter(FDR < q_threshold) %>%
  group_by(GENE) %>%
  slice_min(FDR, with_ties = FALSE) %>%
  ungroup() %>%
  select(GENE, SubchapterLevel, FDR)

readr::write_tsv(best_trait_per_gene, "magma_best_trait_per_gene.tsv")

# ---- 5) Significant genes in each category ----------------------------------
sig_in_category <- deg_cat %>%
  inner_join(best_trait_per_gene, by = c("ensembl_gene_id" = "GENE")) %>%
  select(category, ensembl_gene_id, SubchapterLevel)

# If nothing passes, still write empty and stop gracefully
if (nrow(sig_in_category) == 0) {
  readr::write_tsv(sig_in_category, "magma_enrichment_significant.tsv")
  stop("No DEG genes matched best_trait_per_gene at the chosen q_threshold.")
}

# ---- 6) Background: random genes from MAGMA universe --------------------------
magma_universe <- unique(magma_long$GENE)
stopifnot(length(magma_universe) >= bg_total)

rand_genes <- sample(magma_universe, bg_total)

rand_best <- magma_long %>%
  filter(GENE %in% rand_genes, FDR < q_threshold) %>%
  group_by(GENE) %>%
  slice_min(FDR, with_ties = FALSE) %>%
  ungroup() %>%
  select(GENE, SubchapterLevel)

bg_counts <- rand_best %>% count(SubchapterLevel, name = "bg_k")
bg_total  <- bg_total

# ---- 7) Fisher enrichment per category x trait -------------------------------
cat_totals <- sig_in_category %>% count(category, name = "M")               # total significant genes in category
cat_trait  <- sig_in_category %>% count(category, SubchapterLevel, name = "k")

# Expand to include all traits observed in background (keeps table rectangular)
trait_universe <- unique(bg_counts$SubchapterLevel)

enrich <- cat_trait %>%
  right_join(
    expand.grid(category = cat_totals$category,
                SubchapterLevel = trait_universe,
                KEEP.OUT.ATTRS = FALSE,
                stringsAsFactors = FALSE),
    by = c("category", "SubchapterLevel")
  ) %>%
  left_join(cat_totals, by = "category") %>%
  left_join(bg_counts,  by = "SubchapterLevel") %>%
  mutate(
    k    = replace_na(k, 0L),
    bg_k = replace_na(bg_k, 0L),
    
    # 2x2 table:
    # a = k (sig in category & trait)
    # b = M-k (sig in category & not trait)
    # c = bg_k (bg & trait)
    # d = bg_total - bg_k (bg & not trait)
    a = as.integer(k),
    b = pmax(as.integer(M) - a, 0L),
    c = as.integer(bg_k),
    d = pmax(bg_total - c, 0L)
  ) %>%
  rowwise() %>%
  mutate(
    valid_fisher = (a + b) > 0 & (c + d) > 0,
    pval = if (valid_fisher) fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value else 1,
    # add pseudocounts to stabilize log2 enrichment
    log2enrich = log2(((a + 0.5) / (as.numeric(M) + 1)) /
                        ((c + 0.5) / (bg_total + 1)))
  ) %>%
  ungroup() %>%
  mutate(
    FDR = p.adjust(pval, method = "fdr"),
    positive = log2enrich > 0
  )

readr::write_tsv(enrich, "magma_enrichment_all.tsv")

# ---- 8) Save significant positive enrichments --------------------------------
enrich_sig <- enrich %>%
  filter(FDR < 0.05, positive, log2enrich > 1)

readr::write_tsv(enrich_sig, "magma_enrichment_significant.tsv")

message("Done: MAGMA enrichment exported.")
