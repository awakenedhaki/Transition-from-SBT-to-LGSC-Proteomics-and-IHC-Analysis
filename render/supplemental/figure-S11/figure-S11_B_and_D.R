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
PROTEOMIC_FOXP3_DENSITIES <- here("data", "outputs", "08", "22-002_foxp3_densities.csv")
STANFORD_FOXP3_DENSITIES <- here("data", "outputs", "08", "stanford_foxp3_densities.csv")

# . Histotype Colors
HISTOTYPE_COLORS <- c("HGSC" = "#184E77", SBT = "#52B69A", "LGSC" = "#D9ED92")

# Loading Data =================================================================
# . Proteomic Data
proteomic_foxp3_density <- read_csv(PROTEOMIC_FOXP3_DENSITIES)

# . Stanford Data
stanford_foxp3_density <- read_csv(STANFORD_FOXP3_DENSITIES)

# Visualization ================================================================
# . Proteomic Data
proteomic_foxp3_density |>
  select(histotype, tumour = mean_tumour_density) |>
  mutate(histotype = fct_relevel(histotype, !!!names(HISTOTYPE_COLORS)),
         tumour = log10(tumour + 1),
         marker = "FOXP3") |>
  ggplot(aes(x = histotype, y = tumour, fill = histotype)) +
    geom_boxplot(show.legend = FALSE) +
    facet_wrap(~marker) +
    scale_fill_manual(values = HISTOTYPE_COLORS) +
    expand_limits(y = c(0, 2.8))

# . Stanford Data
stanford_foxp3_density |>
  select(histotype, tumour = mean_tumour_density) |>
  mutate(histotype = fct_relevel(histotype, !!!names(HISTOTYPE_COLORS)),
         tumour = log10(tumour + 1),
         marker = "FOXP3") |>
  ggplot(aes(x = histotype, y = tumour, fill = histotype)) +
    geom_boxplot(show.legend = FALSE) +
    facet_wrap(~marker) +
    scale_fill_manual(values = HISTOTYPE_COLORS) +
    expand_limits(y = c(0, 2.8))

# Save Plot ====================================================================
ggsave(here("figures", "supplemental", "figure-S5", "figure-S5_D.svg"), 
       plot = last_plot(), width = 5, height = 5, dpi = 1200)
