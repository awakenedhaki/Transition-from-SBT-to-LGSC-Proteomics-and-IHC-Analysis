# Loading Dependencies =========================================================
# . Path Handling
library(here)

# . Read/Write Data
library(readr)

# Data Wrangling
library(tidyr)
library(dplyr)
library(tibble)

# Clustering
library(ConsensusClusterPlus)

# Constants ====================================================================
# . Files Paths
METADATA <- here("data", "metadata.csv")
LOG2_COUNTS_MATRIX <- here("data", "log2_counts_matrix.csv")

# . Factor Order
HISTOTYPES <- c("HGSC", "SBT", "LGSC")

# . Threshold
MAD_THRESHOLD <- 0.25

# Loading Data =================================================================
metadata <- read_csv(METADATA) |>
  mutate(histotype = recode(histotype, `mSBT` = "SBT")) |>
  filter(histotype %in% HISTOTYPES, !exclude) |>
  select(set, histotype, stage_full, anatomic_site , cellularity, batch)
  
proteome <- read_csv(LOG2_COUNTS_MATRIX) |>
  select(uniprot, (metadata |>
                     pull(set)))

# Data Preprocessing ===========================================================
tranposed_proteome <- proteome |>
  # Transpose
  # . Protein-Sample Matrix => Sample-Protein Matrix
  pivot_longer(cols = starts_with("set"), names_to = "set", values_to = "log2_count") |>
  pivot_wider(names_from = "uniprot", values_from = "log2_count") |>
  select(-set)

variable_proteins <- tranposed_proteome |>
  summarize_all(\(protein) {
    mad(protein, na.rm = TRUE)
  }) |>
  # Filter by MAD threshold
  pivot_longer(cols = everything(), names_to = "uniprot", values_to = "MAD") |>
  slice_max(order_by = MAD, prop = MAD_THRESHOLD) |>
  # Pull UNIPROT identifier vector
  pull(uniprot)

filtered_proteome <- proteome |>
  filter(uniprot %in% variable_proteins) |>
  # Remove rows with NA
  rowwise() |>
  drop_na() |>
  ungroup()

# Consensus Clustering =========================================================
consensus_cluster <- filtered_proteome |>
  column_to_rownames("uniprot") |>
  as.matrix() |>
  ConsensusClusterPlus(maxK = 10, 
                       clusterAlg = "hc",
                       innerLinkage = "ward.D2",
                       finalLinkage = "ward.D2",
                       distance = "euclidean",
                       seed = 123)

N_CLUSTERS <- 3
sample_clusters <- consensus_cluster[[N_CLUSTERS]]$consensusClass |>
  enframe(name = "set", value = "cluster")

# Save Data ====================================================================
write_csv(sample_clusters, here("data", "outputs", "02", "sample_clusters.csv"))
write_csv(filtered_proteome, here("data", "outputs", "02", "filtered_proteome.csv"))
