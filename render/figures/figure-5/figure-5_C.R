# Loading Dependencies =========================================================
# . Path Handling
library(here)

# . Read/Write Data
library(readr)

# . Data Cleaning/Wrangling
library(dplyr)
library(janitor)

# . Visualization
library(ggplot2)

# ggplot2 Defaults =============================================================
update_geom_defaults("point", list(size = 3))
theme_set(theme_bw())

# Constants ====================================================================
# . Proteomic Cohort
PROTEOMIC_FAP_SCORES <- here("data", "outputs", "05", "proteomic_fap_h_scores.csv")
PROTEOMIC_FOXP3_DENSITIES <- here("data", "outputs", "08", "22-002_foxp3_densities.csv")

# . Stanford Cohort
STANFORD_FAP_SCORES <- here("data", "outputs", "05", "stanford_fap_h_scores.csv")
STANFORD_FOXP3_DENSITIES <- here("data", "outputs", "08", "stanford_foxp3_densities.csv")

# . Histotypes
HISTOTYPE_COLORS <- c("SBT" = "#52B69A", "LGSC" = "#D9ED92")

# Helper Functions =============================================================
linear_regression <- function(tbl, col) {
  column_name <- as_label((ensym(col)))
  
  tbl |>
    mutate({{ col }} := log10({{ col }} + 1)) |>
    ggplot(aes(x = mean_fap_score, y = {{ col }})) +
      geom_point(aes(color = histotype)) +
      geom_smooth(method = "lm") +
      labs(x = "Mean FAP Score", y = bquote(log[10](.(column_name) + 1))) +
      scale_color_manual(values = HISTOTYPE_COLORS) +
      theme(legend.position = "bottom")
}

# Loading Data =================================================================
# . Proteomic Cohort
proteomic_fap_scores <- read_csv(PROTEOMIC_FAP_SCORES)
proteomic_foxp3_densities <- read_csv(PROTEOMIC_FOXP3_DENSITIES)

# . Stanford Cohort
stanford_fap_scores <- read_csv(STANFORD_FAP_SCORES)
stanford_foxp3_densities <- read_csv(STANFORD_FOXP3_DENSITIES) |>
  mutate(histotype = case_when(
    histotype == "Low grade serous carcinoma" ~ "LGSC",
    histotype == "Serous borderline tumor" ~ "SBT",
    .default = histotype
  ))

# Data Merging =================================================================
# . Proteomic Cohort
proteomic_scores <- full_join(proteomic_fap_scores, proteomic_foxp3_densities, 
                              by = join_by(core_id, histotype)) |>
  select(core_id, histotype, 
         mean_fap_score = mean_score, 
         stroma = mean_stroma_density) |>
  filter(
    histotype %in% names(HISTOTYPE_COLORS),
    !is.na(mean_fap_score),
    !is.na(stroma)
  ) |>
  mutate(cohort = "Proteomic")

# . Stanford
stanford_scores <- full_join(stanford_fap_scores, stanford_foxp3_densities, 
                             by = join_by(tma_number, case_code, histotype)) |>
  select(case_code, histotype, mean_fap_score = mean_score, stroma = mean_stroma_density) |>
  filter(
    histotype %in% names(HISTOTYPE_COLORS),
    !is.na(mean_fap_score),
    !is.na(stroma)
  ) |>
  mutate(cohort = "Stanford")

# . Combined Scores
combined_scores <- bind_rows(proteomic_scores, stanford_scores)

# Visualization ================================================================
# . Proteomic Cohort
proteomic_scores |>
  linear_regression(col = stroma)

# . Stanford Cohort
stanford_scores |>
  linear_regression(col = stroma)

# . Combined Scores
combined_scores |>
  linear_regression(col = stroma) +
  facet_wrap(~ cohort)

# Save Plot ====================================================================
ggsave(here("figures", "main", "figure-5", "figure-5_C_and_D.svg"), 
       plot = last_plot(), width = 8, height = 5, dpi = 1200)
