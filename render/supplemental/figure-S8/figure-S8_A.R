# Loading Dependencies =========================================================
# . Path Handling
library(here)

# . Read/Write Data
library(readr)

# . Data Wrangling
library(dplyr)
library(forcats)

# . Visualization
library(ggplot2)

# ggplot2 Defaults =============================================================
theme_set(theme_bw)

# Constants ====================================================================
# . File Paths
CELL_TYPE_DENSITIES <- here("data", "outputs", "06", "cell_type_densities.csv")

# . Histotypes
HISTOTYPE_COLORS <- c("HGSC" = "#184E77", "SBT" = "#52B69A", "LGSC" = "#D9ED92")

# . Cell Types
CELL_TYPES <- c("B cell", "Plasma cell", 
                "CD8- T cell", 
                "CD8+ T cell", 
                "Macrophage")

# Load Data ====================================================================
cell_type_densities <- read_csv(CELL_TYPE_DENSITIES)

# Visualization ================================================================
cell_type_densities |>
  mutate(histotype = ifelse(histotype == "mSBT", "SBT", histotype)) |>
  mutate(histotype = fct_relevel(histotype, !!!names(HISTOTYPE_COLORS)),
         cell_type = fct_relevel(cell_type, !!!CELL_TYPES)) |>
  ggplot(aes(x = histotype, y = pseudo_log10_density, fill = histotype)) + 
    geom_boxplot(show.legend = FALSE) +
    facet_grid(compartment ~ cell_type) +
    labs(x = "Histotype", y = expression(log[10]("Cell Count / Stroma Area"))) +
    scale_fill_manual(values = HISTOTYPE_COLORS) +
    expand_limits(y = c(0, 3.75))

# Save Plot ====================================================================
ggsave(here("figures", "supplemental", "figure-S8", "figure-S8A.svg"), 
       plot = last_plot(), width = 8, height = 5, dpi = 1200)
