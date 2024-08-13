# Loading Dependencies =========================================================
# . Path Handling
library(here)

# . Read/Write Data
library(readr)

# . Data Cleaning/Wrangling
library(dplyr)
library(forcats)
library(janitor)

# . Visualization
library(ggplot2)

# ggplot2 Defaults =============================================================
update_geom_defaults("point", list(size = 3))
theme_set(theme_bw())

# Constants ====================================================================
# . Proteomic Cohort
PROTEOMIC_FAP_SCORES <- here("data", "outputs", "05", "proteomic_fap_h_scores.csv")
PROTEOMIC_MACROPHAGE_DENSITIES <- here("data", "outputs", "07", "22-002_macrophage_densities.csv")

# . Stanford Cohort
STANFORD_FAP_SCORES <- here("data", "outputs", "05", "stanford_fap_h_scores.csv")
STANFORD_MACROPHAGE_DENSITIES <- here("data", "outputs", "07", "stanford_macrophage_densities.csv")

# . Histotypes
HISTOTYPE_COLORS <- c("SBT" = "#52B69A", "LGSC" = "#D9ED92")

# Helper Functions =============================================================
linear_regression <- function(tbl, col) {
  column_name <- as_label(ensym(col))
  
  tbl |>
    mutate({{ col }} := log10({{ col }} + 1)) |>
    ggplot(aes(x = mean_fap_score, y = {{ col }})) +
      geom_point(aes(color = histotype)) +
      geom_smooth(method = "lm") +
      facet_wrap(~marker) +
      labs(x = "Mean FAP Score", y = bquote(log[10](.(column_name) + 1))) +
      scale_color_manual(values = HISTOTYPE_COLORS) +
      theme(legend.position = "bottom")
}

# Loading Data =================================================================
# . Proteomic Cohort
proteomic_fap_scores <- read_csv(PROTEOMIC_FAP_SCORES) |>
  rename(mean_fap_score = mean_score)

proteomic_macrophage_densities <- read_csv(PROTEOMIC_MACROPHAGE_DENSITIES)

# . Stanford Cohort
stanford_fap_scores <- read_csv(STANFORD_FAP_SCORES) |>
  rename(mean_fap_score = mean_score)

stanford_macrophage_densities <- read_csv(STANFORD_MACROPHAGE_DENSITIES) |>
  mutate(histotype = case_when(
    histotype == "Low grade serous carcinoma" ~ "LGSC",
    histotype == "Serous borderline tumor" ~ "SBT",
    .default = histotype
  ))

# Data Merging =================================================================
# . Proteomic Cohort
proteomic_scores <- full_join(proteomic_fap_scores, proteomic_macrophage_densities, 
                              by = join_by(core_id, histotype)) |>
  relocate(core_id, histotype, voa, marker) |>
  filter(
    histotype %in% names(HISTOTYPE_COLORS),
    !is.na(mean_fap_score),
    !is.na(stroma)
  ) |>
  mutate(cohort = "Proteomic")

# . Stanford
stanford_scores <- full_join(stanford_fap_scores, stanford_macrophage_densities, 
                             by = join_by(tma_number, case_code, histotype)) |>
  relocate(case_code, histotype, marker) |>
  filter(
    histotype %in% names(HISTOTYPE_COLORS),
    !is.na(mean_fap_score),
    !is.na(stroma)
  ) |>
  mutate(cohort = "Stanford")

# Visualization ================================================================
# . Proteomic Cohort
proteomic_scores |>
  mutate(marker = fct_relevel(marker, "CD68", "CD163")) |>
  linear_regression(col = stroma) +
  expand_limits(y = c(0.3, 0.7))

# . Stanford Cohort
stanford_scores |>
  mutate(marker = fct_relevel(marker, "CD68", "CD163")) |>
  linear_regression(col = stroma) +
  expand_limits(y = c(0.3, 0.7))

# Save Plot ====================================================================
ggsave(here("figures", "main", "figure-5", "figure-5_A.svg"), 
       plot = last_plot(), width = 8, height = 5, dpi = 1200)
