# ==============================================================================
# Differential Abundance Analysis
# . Refactored from Gian Negri's code
# ==============================================================================

# Loading Dependencies =========================================================
# . Path Handling
library(here)

# . Read/Write Data
library(readr)

# Data Wrangling
library(dplyr)
library(tibble)
library(janitor)
library(magrittr)

# . String Manipulation
library(glue)
library(stringr)

# . Differential Abundance Analysis
library(limma)
library(DEqMS)
library(matrixStats)

# Constants ====================================================================
# . Files Paths
METADATA <- here("data", "metadata.csv")
LOG2_COUNTS_MATRIX <- here("data", "log2_counts_matrix.csv")

# . Factor Order
HISTOTYPES <- c("HGSC", "SBT", "LGSC")

# . Feature of Interest
COLUMN <- "histotype"

# . Groups of Interest
GROUP_1 <- "SBT"
GROUP_2 <- "LGSC"

# . Output Filename
FILENAME <- glue("{GROUP_2}_v_{GROUP_1}.csv")

# Helper Functions =============================================================
select_sample_ids <- function(tbl, col, pattern) {
  tbl %>%
    filter(!!sym(col) == pattern) %>%
    pull(set)
}

# Loading Data =================================================================
metadata <- read_csv(METADATA) %>%
  mutate(histotype = recode(histotype, `mSBT` = "SBT")) %>%
  filter(histotype %in% HISTOTYPES, !exclude) %>%
  select(set, histotype, stage_full, anatomic_site , cellularity, batch)
  
proteome <- read_csv(LOG2_COUNTS_MATRIX)

# Differential Abundance Analysis ==============================================
# . Select desired samples
group_1_sample_ids <- select_sample_ids(metadata, COLUMN, GROUP_1)
group_2_sample_ids <- select_sample_ids(metadata, COLUMN, GROUP_2)

length(group_1_sample_ids)
length(group_2_sample_ids)

# . log2(peptide abundance)
counts_matrix <- proteome %>%
  select(uniprot, group_1_sample_ids, group_2_sample_ids) %>%
  column_to_rownames("uniprot") %>%
  filter(apply(X = .,
               MARGIN = 1,
               FUN = \(row) {
                 length(na.omit(row)) >= ncol(.) * 0.33
               }))

# . Explanatory Variable Levels
labels <- ifelse(colnames(counts_matrix) %in% group_1_sample_ids,
                 "group_1", 
                 "group_2")

# . Design Matrix
design_matrix <- model.matrix(~0 + labels, data = counts_matrix)
colnames(design_matrix) <- gsub("labels", "", colnames(design_matrix))

# . Contrast Matrix
contrast_matrix <- makeContrasts(group_1 - group_2, levels = design_matrix)

# . Fitting linear model
fitted_counts_matrix <- counts_matrix %>%
  lmFit(design = design_matrix) %>%
  contrasts.fit(contrasts = contrast_matrix) %>%
  eBayes()

fitted_counts_matrix$count <- proteome %>%
  filter(uniprot %in% rownames(fitted_counts_matrix)) %>%
  pull(number_of_peptides)

# . Final Output
differentially_abundant_proteins <- spectraCounteBayes(fitted_counts_matrix) %>%
  outputResult(coef_col = 1)  %>%
  as_tibble() %>%
  clean_names() %>%
  rename(uniprot = gene) %>%
  left_join(proteome %>%
              select(uniprot, symbol, description),
            by = "uniprot") %>%
  relocate(c(uniprot, symbol, description), .before = "log_fc") %>%
  rename_at(vars(starts_with("sca")), .funs = ~str_replace(.x, "sca", "deqms")) %>%
  mutate(nlog10_adj_pval = -log10(deqms_adj_pval))

# Saving Results ===============================================================
write_csv(differentially_abundant_proteins, here("data", "outputs", "03", FILENAME))
