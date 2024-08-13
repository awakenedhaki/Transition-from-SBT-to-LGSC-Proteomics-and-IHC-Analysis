# Loading Dependencies =========================================================
# . Path Handling
library(here)

# . Read/Write Data
library(readr)

# . Data Wrangling
library(dplyr)

# . Visualization
library(ggplot2)
library(ggforce)

# Constants ====================================================================
# . File Paths
PCA <- here("data", "outputs", "01", "pca_k_3.csv")

# . Visualizations
SQUARE <- 15
CIRCLE <- 16
TRIANGLE <- 17
ROMBUS <- 18

HISTOTYPES <- list("HGSC" = c("color" = "#184E77", "shape" = TRIANGLE),
                   "SBT"  = c("color" = "#52B69A", "shape" = CIRCLE),
                   "LGSC" = c("color" = "#D9ED92", "shape" = SQUARE))

CLUSTERS <- list("1" = c("color" = "#184E77", "shape" = TRIANGLE),
                 "2" = c("color" = "#52B69A", "shape" = CIRCLE),
                 "3" = c("color" = "#D9ED92", "shape" = SQUARE))

# ggplot2 Defaults =============================================================
update_geom_defaults("point", list(size = 4))
theme_set(theme_bw())

# Helper Functions =============================================================
extract_named_aesthetics <- function(list, aesthetic) {
  attributes <- list |>
    lapply(\(x) {
      x[[aesthetic]]
    }) |>
    unlist()
  
  if (aesthetic == "shape") {
    attributes <- parse_number(attributes)
  }
  
  return(attributes)
}

# Loading Data =================================================================
clustered_pca <- read_csv(PCA,
                          col_select = c(histotype, cluster, PC1, PC2)) |>
  mutate(histotype = fct_relevel(histotype, !!!names(HISTOTYPES)),
         cluster = fct_relevel(as.character(cluster), !!!names(CLUSTERS)))

# Figure 1B ====================================================================
clustered_pca |>
  ggplot(aes(x = PC1, y = PC2)) +
    # Demarcate clusters
    geom_mark_ellipse(aes(fill = cluster), show.legend = FALSE) +
    # Mark samples
    geom_point(aes(shape = cluster, color = histotype)) +
    # Cluster Aesthetics
    scale_shape_manual(values = extract_named_aesthetics(CLUSTERS, "shape")) +
    scale_fill_manual(values = extract_named_aesthetics(CLUSTERS, "color")) +
    # Histotype Aesthetics
    scale_color_manual(values = extract_named_aesthetics(HISTOTYPES, "color")) +
    # Labels
    labs(color = "Histotype", shape = "k-Means Cluster",
         x = "PC1 (15.5%)", y = "PC2 (8.43%)") +
    # Miscellaneous
    theme(panel.border = element_rect(color = "black", fill = "transparent"), 
          axis.title = element_text(size = 12), 
          axis.text = element_text(size = 11),
          legend.position = "bottom")

# Save Plot ====================================================================
ggsave(here("figures", "main", "figure-1", "figure-1_B.svg"), width = 6, height = 6, dpi = 300)
