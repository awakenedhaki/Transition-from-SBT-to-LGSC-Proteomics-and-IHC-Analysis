# Loading Dependencies =========================================================
# . Path Handling
library(here)

# . Read/Write Data
library(readr)

# . Data Wrangling
library(dplyr)
library(tidyr)
library(forcats)

# . Visualization
library(ggplot2)

# ggplot2 Defaults =============================================================
theme_set(theme_bw())

# Constants ====================================================================
# . File Paths
PROTEOME <- here("data", "log2_counts_matrix.csv")
METADATA <- here("data", "metadata.csv")

# . Histotype Colors
HISTOTYPE_COLORS <- c("HGSC" = "#184E77", "SBT" = "#52B69A", "LGSC" = "#D9ED92")

# Loading Data =================================================================
# . Loading Metadata
metadata <- read_csv(METADATA) |>
  select(set, histotype) |>
  filter(histotype %in% names(HISTOTYPE_COLORS))

# . Loading Proteome Matrix
ITGA5 <- read_csv(PROTEOME) |>
  select(symbol, contains("set")) |>
  pivot_longer(cols = contains("set"), 
               names_to = "set", 
               values_to = "abundance") |>
  inner_join(metadata, by = join_by(set == set)) |>
  filter(symbol == "ITGA5")

# Visualization ================================================================
ITGA5 |>
  mutate(histotype = fct_relevel(histotype, !!!names(HISTOTYPE_COLORS))) |>
  ggplot(aes(x = histotype, y = abundance, fill = histotype)) +
    geom_boxplot(show.legend = FALSE) +
    scale_fill_manual(values = HISTOTYPE_COLORS)

# Save Plot ====================================================================
ggsave(here("figures", "revisions", "figure-S5", "figure-S5_C.svg"), 
       plot = last_plot(), width = 5, height = 5, dpi = 1200)
