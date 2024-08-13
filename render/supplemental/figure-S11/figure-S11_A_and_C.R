# Loading Dependencies =========================================================
# . Path Handling
library(here)

# . Read/Write Data
library(readr)

# . Data Cleaning/Wrangling
library(dplyr)
library(forcats)

# . Visualization
library(ggplot2)

# ggplot2 Defaults =============================================================
theme_set(theme_bw())

# Constants ====================================================================
# . Proteomic Cohort
PROTEOMIC_MACROPHAGE_DENSITIES <- here("data", "outputs", "07", "22-002_macrophage_densities.csv")

# . Stanford Cohort
STANFORD_MACROPHAGE_DENSITIES <- here("data", "outputs", "07", "stanford_macrophage_densities.csv")

# . Histotypes
HISTOTYPE_COLORS <- c("SBT" = "#52B69A", "LGSC" = "#D9ED92")

# Loading Data =================================================================
# . Proteomic Cohort
proteomic_macrophage_densities <- read_csv(PROTEOMIC_MACROPHAGE_DENSITIES) |>
  filter(histotype %in% names(HISTOTYPE_COLORS)) |>
  mutate(histotype = fct_relevel(histotype, !!!names(HISTOTYPE_COLORS)),
         marker = fct_relevel(marker, c("CD68", "CD163"))) |>
  select(histotype, marker, tumour)

# . Stanford Cohort
stanford_macrophage_densities <- read_csv(STANFORD_MACROPHAGE_DENSITIES) |>
  mutate(histotype = case_when(
    histotype == "Serous borderline tumor" ~ "SBT",
    histotype == "Low grade serous carcinoma" ~ "LGSC",
    .default = histotype
  )) |>
  filter(histotype %in% names(HISTOTYPE_COLORS)) |>
  mutate(histotype = fct_relevel(histotype, !!!names(HISTOTYPE_COLORS)),
         marker = fct_relevel(marker, c("CD68", "CD163"))) |>
  select(histotype, marker, tumour)

# Visualization ================================================================
# . Protoemic Cohort
proteomic_macrophage_densities |>
  ggplot(aes(x = histotype, y = tumour, fill = histotype)) +
    geom_boxplot(show.legend = FALSE) +
    facet_wrap(~marker) +
    scale_fill_manual(values = HISTOTYPE_COLORS) 

# . Stanford Cohort
stanford_macrophage_densities |>
  ggplot(aes(x = histotype, y = tumour, fill = histotype)) +
    geom_boxplot(show.legend = FALSE) +
    facet_wrap(~marker) +
    scale_fill_manual(values = HISTOTYPE_COLORS) 
  
# Save Plot ====================================================================
ggsave(here("figures", "supplemental", "figure-S11", "figure-S11_A.svg"), 
       plot = last_plot(), width = 8, height = 5, dpi = 1200)
