# Loading Dependencies =========================================================
# . Path Handling
library(here)

# . Read/Write Data
library(readr)

# . Data Cleaning/Wrangling
library(dplyr)
library(tidyr)
library(purrr)
library(forcats)
library(janitor)

# . Clean Linear Model Output
library(broom)

# Constants ====================================================================
# . Proteomic Cohort
PROTEOMIC_FAP_SCORES <- here("data", "outputs", "05", "proteomic_fap_h_scores.csv")
PROTEOMIC_MACROPHAGE_DENSITIES <- here("data", "outputs", "07", "22-002_macrophage_densities.csv")

# . Stanford Cohort
STANFORD_FAP_SCORES <- here("data", "outputs", "05", "stanford_fap_h_scores.csv")
STANFORD_MACROPHAGE_DENSITIES <- here("data", "outputs", "07", "stanford_macrophage_densities.csv")

# . Histotypes
HISTOTYPES <- c("SBT", "LGSC")

# Helper Functions =============================================================
linear_regression <- function(tbl) {
  tbl |> 
    select(histotype, marker, mean_fap_score, stroma, tumour) |>
    pivot_longer(cols = c(stroma, tumour), names_to = "compartment", values_to = "density") |>
    nest(data = c(histotype, mean_fap_score, density)) |>
    mutate(linear_model = map(data, \(tbl) glance(lm(log10(density + 1) ~ mean_fap_score, data = tbl)))) |>
    unnest(cols = linear_model)
}

correlation_test <- function(tbl) {
  tbl |>
    select(histotype, marker, mean_fap_score, stroma, tumour) |>
    pivot_longer(cols = c(stroma, tumour), names_to = "compartment", values_to = "density") |>
    nest(data = c(histotype, mean_fap_score, density)) |>
    mutate(linear_model = map(data, \(tbl) tidy(cor.test(log10(tbl$density + 1), tbl$mean_fap_score, method = "pearson")))) |>
    unnest(cols = linear_model)
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
    histotype %in% HISTOTYPES,
    !is.na(mean_fap_score),
    !is.na(stroma)
  ) |>
  mutate(cohort = "Proteomic")

# . Stanford
stanford_scores <- full_join(stanford_fap_scores, stanford_macrophage_densities, 
                             by = join_by(tma_number, case_code, histotype)) |>
  relocate(case_code, histotype, marker) |>
  filter(
    histotype %in% HISTOTYPES,
    !is.na(mean_fap_score),
    !is.na(stroma)
  ) |>
  mutate(cohort = "Stanford")

# Linear Model =================================================================
# . Proteomic Cohort
proteomic_lm_results <- proteomic_scores |>
  mutate(marker = fct_relevel(marker, "CD68", "CD163")) |>
  linear_regression() |>
  select(-data)

# . Stanford Cohort
stanford_lm_results <- stanford_scores |>
  mutate(marker = fct_relevel(marker, "CD68", "CD163")) |>
  linear_regression() |>
  select(-data)

# Correlation Test =============================================================
proteomic_cor_results <- proteomic_scores |>
  mutate(marker = fct_relevel(marker, "CD68", "CD163")) |>
  correlation_test() |>
  select(-data)

stanford_cor_results <- stanford_scores |>
  mutate(marker = fct_relevel(marker, "CD68", "CD163")) |>
  correlation_test() |>
  select(-data)

 # Save Results =================================================================
# . Mean FAP Score and FOXP3 Densities
# . . Proteomic Cohort
write_csv(proteomic_scores, file = here("data", "outputs", "09", "22-002_mean_fap_scores_and_macrophage_densities.csv"))

# . . Stanford Cohort
write_csv(stanford_scores, file = here("data", "outputs", "09", "stanford_mean_fap_scores_and_macrophage_densities.csv"))

# . Stroma Results - Linear Model
# . . Proteomic Cohort
proteomic_lm_results |>
  write_csv(file = here("data", "outputs", "09", "22-002_linear_model_results.csv"))

# . . Stanford Cohort
stanford_lm_results |>
  write_csv(file = here("data", "outputs", "09", "stanford_linear_model_results.csv"))

# . Stroma Results - Correlation Test
# . . Proteomic Cohort
proteomic_cor_results |>
  write_csv(file = here("data", "outputs", "09", "22-002_cor_test_results.csv"))

# . . Stanford Cohort
stanford_cor_results |>
  write_csv(file = here("data", "outputs", "09", "stanford_cor_test_results.csv"))
