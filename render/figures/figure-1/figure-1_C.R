# Loading Dependencies =========================================================
# . Path Handling
library(here)

# . Read/Write Data
library(readr)

# . Data Wrangling
library(dplyr)
library(tidyr)
library(tibble)
library(forcats)
library(stringr)

# . Visualization
library(pheatmap)
library(viridisLite)
library(RColorBrewer)

# Constants ====================================================================
# . File Paths
METADATA <- here("data", "metadata.csv")
SAMPLE_CLUSTERS <- here("data", "outputs", "02", "sample_clusters.csv")
FILTERED_PROTEOME <- here("data", "outputs", "02", "filtered_proteome.csv")

# . Heatmap Tree Depth
ROW_TREE_DEPTH <- 3
COLUMN_TREE_DEPTH <- 3

# . Metadata Factor Levels
HISTOTYPE_LEVELS <- c("HGSC", "SBT", "LGSC")
MUTATION_LEVELS <- c("WT", "NRAS", "KRAS", "BRAF", "Not Assayed", "No Material")
FIGO_STAGE_LEVELS <- c("I", "II", "III", "IV", "Unstaged")

# . Visualization Color Palette
COLOR_PALETTE <- list(
  "Mutation" = setNames(brewer.pal(n = 6, name = "Set3"), MUTATION_LEVELS),
  "FIGO Stage" = setNames(brewer.pal(n = 5, name = "Set2"), FIGO_STAGE_LEVELS),
  "Cluster" = c("1" = "#184E77", "2" = "#52B69A",  "3" = "#D9ED92"),
  "Histotype" = c("HGSC" = "#184E77", "SBT" = "#52B69A", "LGSC" = "#D9ED92")
)

# Helper Functions =============================================================
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

# Loading Data =================================================================
filtered_proteome <- read_csv(FILTERED_PROTEOME)

metadata <- read_csv(METADATA, col_select = c(set, histotype, stage_full, mutation))
sample_clusters <- read_csv(SAMPLE_CLUSTERS, col_types = list(cluster = col_factor()))

# Heatmap Column Annotations ===================================================
column_annotations <- metadata |>
  # Adding cluster labels to metadata
  left_join(sample_clusters, by = join_by(set)) |>
  # Removing NA
  filter(!is.na(cluster)) |>
  # Recoding histotype
  mutate(histotype = recode(histotype, `mSBT` = "SBT")) |>
  # Re-level factors
  mutate(histotype = fct_relevel(histotype, !!!HISTOTYPE_LEVELS),
         mutation = fct_relevel(mutation, !!!MUTATION_LEVELS)) |>
  # FIGO Stage Modification
  mutate(stage_full = str_extract(stage_full, pattern = "^([IV]*)"),
         stage_full = ifelse(stage_full == "", "Unstaged", stage_full),
         stage_full = fct_relevel(stage_full, !!!FIGO_STAGE_LEVELS)) |>
  # Renaming features
  rename(Cluster = cluster, Histotype = histotype,
         `FIGO Stage` = stage_full, Mutation = mutation) |>
  column_to_rownames("set")

# Data Preprocessing ===========================================================
filtered_proteome_matrix <- filtered_proteome |>
  column_to_rownames("uniprot") |>
  as.matrix() |>
  t() |>
  scale(scale = TRUE, center = TRUE) |>
  t()

# Quantiled Color Gradient =====================================================
color_breaks <- quantile_breaks(filtered_proteome_matrix, n = 1000)

# Visualization ================================================================
(heatmap <- filtered_proteome_matrix |>
  pheatmap(clustering_method = "ward.D2",
           # Color Parameters
           color = plasma(n = length(color_breaks) - 1),
           breaks = color_breaks,
           annotation_colors = COLOR_PALETTE,
           # Rows Parameters
           show_rownames = FALSE,
           cluster_rows = TRUE,
           cutree_rows = ROW_TREE_DEPTH,
           clustering_distance_rows = "euclidean",
           # Column Parameters
           show_colnames = FALSE,
           cluster_cols = TRUE,
           cutree_cols = COLUMN_TREE_DEPTH,
           clustering_distance_cols = "euclidean",
           annotation_col = column_annotations))

# Save Plot ====================================================================
ggsave(here("figures", "main", "figure-1", "figure-1_C.svg"), plot = heatmap, width = 10, height = 10, dpi = 300)
       