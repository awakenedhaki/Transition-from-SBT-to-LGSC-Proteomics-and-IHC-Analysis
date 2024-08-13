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
theme_set(theme_bw())

# Constants ====================================================================
# . File Paths
STANFORD_FOXP3_DENSITIES <- here("data", "outputs", "08", "stanford_foxp3_densities.csv")

# . Histotype Colors
HISTOTYPE_COLORS <- c("SBT" = "#52B69A", "LGSC" = "#D9ED92")

# Loading Data =================================================================
foxp3_density <- read_csv(STANFORD_FOXP3_DENSITIES)

# Visualization ================================================================
foxp3_density |>
  mutate(across(contains("density"), \(x) log10(x + 1), .names = "pseudo_log10_{.col}")) |>
  select(histotype, stroma = pseudo_log10_mean_stroma_density) |>
  mutate(histotype = fct_relevel(histotype, !!!names(HISTOTYPE_COLORS)),
         marker = "FOXP3") |>
  ggplot(aes(x = histotype, y = stroma, fill = histotype)) +
    geom_boxplot(show.legend = FALSE) +
    facet_wrap(~marker) +
    scale_fill_manual(values = HISTOTYPE_COLORS) +
    scale_y_continuous(breaks = seq(0, 2.5, 1)) +
    expand_limits(y = c(0, 2.5))

# Save Plot ====================================================================
ggsave(here("figures", "main", "figure-4", "figure-4_F.svg"), 
       plot = last_plot(), width = 5, height = 5, dpi = 1200)
