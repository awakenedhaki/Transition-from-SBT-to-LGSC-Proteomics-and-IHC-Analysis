# Loading Dependencies =========================================================
# . Path Handling
library(here)

# . Clean k-means & PCA outputs
library(broom)

# . Read/Write Data
library(readr)

# . Data Wrangling
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(forcats)

# . String Manipulation
library(stringr)

# . Visualization
library(ggplot2)

# Constants ====================================================================
# . Files Paths
METADATA <- here("data", "metadata.csv")
LOG2_COUNTS_MATRIX <- here("data", "log2_counts_matrix.csv")

# . Factor Order
HISTOTYPES <- c("HGSC", "SBT", "LGSC")

# ggplot2 Defaults =============================================================
update_geom_defaults("point", list(size = 4))
theme_set(theme_minimal())

# Loading Data =================================================================
metadata <- read_csv(METADATA) |>
  mutate(histotype = recode(histotype, `mSBT` = "SBT")) |>
  filter(histotype %in% HISTOTYPES, !exclude) |>
  select(set, histotype, stage_full, anatomic_site , cellularity, batch)


read_csv(METADATA) %>%
  mutate(histotype = recode(histotype, `mSBT` = "SBT")) |>
  filter(!(histotype %in% HISTOTYPES)) %>%
  View()
  
 
proteome <- read_csv(LOG2_COUNTS_MATRIX) |>
  select(uniprot, symbol, (metadata |>
                             pull(set)))

# PCA ==========================================================================
pca <- proteome |>
  # Removing HGNC Symbols (No sufficiently unique)
  select(-symbol) |>
  # Labeling rows by UniProt IDs
  column_to_rownames("uniprot") |>
  # Removing rows with missing values
  rowwise() |>
  drop_na() |>
  ungroup() |>
  # Tranpsose
  # . Protein-Sample Matrix => Sample-Protein Matrix
  t() |>
  # Perform PCA
  prcomp(center = TRUE)

components <- pca |>
  # Extracting sample coordinates under the principal components
  tidy(matrix = "samples") |>
  # Rename "row" to "set"
  rename(set = row) |>
  # Annotating samples with metadata
  left_join(metadata, by = join_by(set == set)) |>
  # Selecting the first three PCs
  filter(PC <= 3) |>
  # Renaming PCs
  mutate(PC = str_c("PC", PC)) |>
  # Converting each PC into a column
  pivot_wider(names_from = PC, values_from = value)

eigenvalues <- pca |>
  tidy(matrix = "eigenvalues")

# k-means Clustering ===========================================================
set.seed(123457)
kclusts <- 
  # Create a tibble with k values
  tibble(k = 1:9) |>
  # Perform k-means clustering for each k
  mutate(
    # Perform k-means clustering for each k value
    kclust = map(k, \(x) {
      components |>
        select(set, contains("PC")) |>
        column_to_rownames("set") |>
        kmeans(x, iter.max = 5000)
    }),
    # Cluster assignments
    tidied = map(kclust, tidy),
    # Diagnostic values
    glanced = map(kclust, glance),
    # Metadata annotated
    augmented = map(kclust, augment, components)
  )

clusters <- kclusts |>
  unnest(cols = c(tidied))

clusterings <- kclusts |>
  unnest(cols = c(glanced))

assignments <- kclusts |>
  unnest(cols = c(augmented))

# Diagnostic Visualization =====================================================
# . Elbow Plot
clusterings |>
  ggplot(aes(x = k, y = tot.withinss)) +
    geom_line() +
    geom_point() +
    scale_x_continuous(labels = 1:9, breaks = 1:9)

# . Clustering Assessment
assignments |>
  ggplot(aes(x = PC1, y = PC2)) +
    geom_point(aes(color = .cluster), alpha = 0.8) +
    geom_point(data = clusters, size = 10, shape = "x") +
    facet_wrap(~k) +
    theme(panel.border = element_rect(color = "black", fill = "transparent"))

# k = 3 ========================================================================
cluster_labeled <- kclusts |>
  filter(k == 3) |>
  unnest(col = augmented) |>
  mutate(recoded_cluster = case_match(.cluster,
                                      "1" ~ "3",
                                      "2" ~ "1",
                                      "3" ~ "2")) |>
  mutate(histotype = fct_relevel(histotype, HISTOTYPES))

# . . Proportion of Histotypes per Cluster
cluster_labeled |>
  count(recoded_cluster, histotype) |>
  complete(recoded_cluster, histotype, fill = list(n = 0)) |>
  add_count(recoded_cluster, wt = n) |>
  mutate(prop = (n / nn) * 100) 

# Write Data ===================================================================
cluster_labeled |>
  select(set = row, histotype, stage_full, anatomic_site, cellularity, batch, 
         cluster = recoded_cluster, PC1, PC2) |>
  write_csv(here("data", "outputs", "01", "pca_k_3.csv"))
