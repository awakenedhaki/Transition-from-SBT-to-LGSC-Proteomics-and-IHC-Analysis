# Loading Dependencies =========================================================
# . Path Handling
library(here)

# . Read/Write Data
library(readr)

# . Gene Set Enrichment Analysis Package
library(fgsea)
library(msigdbr)
library(org.Hs.eg.db)

# . Data Wrangling
# . . Must load after msigdbr & org.Hs.eg.db
library(glue)
library(dplyr)
library(purrr)
library(tibble)
library(janitor)

# Constants ====================================================================
# . Groups of Interest
GROUP_1 <- "LGSC"
GROUP_2 <- "SBT"

# . Input Filename
INPUT <- glue("{GROUP_2}_v_{GROUP_1}.csv")

# . File Paths
DIFFERENTIAL_ABUNDANCE <- here("data", "outputs", "03", INPUT)

# . Output Filename
OUTPUT <- glue("{GROUP_2}_v_{GROUP_1}-gsea.csv")

# . Gene Set for GSEA
# . . C5: Ontology gene sets
# . . C6: Oncogenic signature gene
CATEGORY <- "C5"
# . . NULL if no subcategory
SUBCATEGORY <- "BP"

# . Gene Sets for GSEA
GENE_SETS <- msigdbr(species = "Homo sapiens",
                     category = "C5",
                     subcategory = "BP") |>
  select(gs_name, entrez_gene) |>
  mutate(entrez_gene = as.character(entrez_gene)) |>
  unstack(entrez_gene ~ gs_name)

# Loading Data =================================================================
differential_abundance <- read_csv(DIFFERENTIAL_ABUNDANCE) |>
  select(uniprot, log_fc, deqms_adj_pval)

entrez_ids <- mapIds(org.Hs.eg.db, 
                     keys = differential_abundance$uniprot, 
                     keytype = "UNIPROT", 
                     column = "ENTREZID", 
                     multiVals = "first") |>
  enframe(name = "uniprot", value = "entrezid")

# Ranking Proteins =============================================================
ranked_proteins <- differential_abundance |>
  # Calculate ranking score
  mutate(score = log_fc * -log10(deqms_adj_pval)) |>
  arrange(score) |>
  # Annotate proteins by entrez ID
  left_join(entrez_ids, by = join_by(uniprot)) |>
  # Pull a named vector
  pull(var = log_fc, name = entrezid)

# GSEA =========================================================================
gsea <- 
  # Run Fast Gene Set Enrichment Analysis
  fgsea(pathways = GENE_SETS, stats = ranked_proteins, minSize = 10, maxSize = 500) |>
  # Convert to tibble
  as_tibble() |>
  # Standardized column naming conventions
  clean_names() |>
  # Collapse leading edge genes
  mutate(leading_edge = map_chr(leading_edge, ~paste(.x, collapse = ", ")))

# Saving Data ==================================================================
write_csv(gsea, here("data", "outputs", "04", OUTPUT))
