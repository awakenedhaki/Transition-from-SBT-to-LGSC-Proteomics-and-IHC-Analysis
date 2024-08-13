# Loading Dependencies =========================================================
# . Path Handling
library(here)

# . Read/Write Data
library(readr)

# . Data Wrangling
library(dplyr)
library(tidyr)
library(purrr)
library(forcats)

# . String Manipulation
library(stringr)

# . Visualization
library(ggplot2)

# Constants ====================================================================
# . File Paths
NORMALIZED_CELL_COUNTS <- here("data", "outputs", "00-revisions", "22-002_nuclear_normalized_cell_counts.csv")

# . Histotype Colors
HISTOTYPE_COLORS <- c("HGSC" = "#184E77", "SBT" = "#52B69A", "LGSC" = "#D9ED92")

# ggplot2 Defaults =============================================================
theme_set(theme_bw())

# Loading Data =================================================================
normalized_cell_counts <- read_csv(NORMALIZED_CELL_COUNTS)

# Visualization ================================================================
# . All Histotypes and Cell Types
normalized_cell_counts |>
  pivot_longer(cols = contains("normalized"), 
               names_to = "compartment", 
               values_to = "normalized_counts") |>
  mutate(pseudo_log10_normalized_counts = map_dbl(normalized_counts, \(x) log10(x + 1)),
         histotype = fct_relevel(histotype, !!!names(HISTOTYPE_COLORS))) |>
  ggplot(aes(x = histotype, y = pseudo_log10_normalized_counts, fill = histotype)) +
    geom_boxplot(show.legend = FALSE) +
    facet_grid(compartment ~ cell_type) +
    labs(x = "Histotype", y = "log10(Cell Count / Total Cell Number)") +
    scale_fill_manual(values = HISTOTYPE_COLORS)

# Save Plot ====================================================================
ggsave(here("figures", "revisions", "figure-S7", "figure-S7_A.svg"), 
       plot = last_plot(), width = 8, height = 5, dpi = 300)
