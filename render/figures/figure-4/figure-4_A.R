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

# Constants ====================================================================
# . File Paths
PROTEOMIC_CELL_DENSITIES <- here("data", "outputs", "06", "cell_type_densities.csv")

# . Filters
COMPARTMENT <- "stroma"
CELL_TYPES <- c("CD8- T cell", "Macrophage")

# . Histotype Colors
HISTOTYPE_COLORS <- c("SBT" = "#52B69A", "LGSC" = "#D9ED92")

# Loading Data =================================================================
cell_type_densities <- read_csv(PROTEOMIC_CELL_DENSITIES)

# Visualization ================================================================
cell_type_densities |>
  filter(histotype %in% names(HISTOTYPE_COLORS),
         compartment == COMPARTMENT,
         cell_type %in% CELL_TYPES) |>
  mutate(histotype = fct_relevel(histotype, !!!names(HISTOTYPE_COLORS))) |>
  ggplot(aes(x = histotype, y = pseudo_log10_density, fill = histotype)) + 
    geom_boxplot(show.legend = FALSE) +
    facet_wrap(~cell_type) +
    labs(x = "Histotype", y = "stroma") +
    scale_fill_manual(values = HISTOTYPE_COLORS) +
    scale_y_continuous()
    expand_limits(y = c(0, 3.75))

# Save Plot ====================================================================
ggsave(here("figures", "main", "figure-4", "figure-4_A.svg"), 
       plot = last_plot(), width = 8, height = 5, dpi = 1200)
