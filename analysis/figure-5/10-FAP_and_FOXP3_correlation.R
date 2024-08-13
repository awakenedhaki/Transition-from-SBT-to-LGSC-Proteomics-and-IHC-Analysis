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
PROTEOMIC_FOXP3_DENSITIES <- here("data", "outputs", "08", "22-002_foxp3_densities.csv")

# . Stanford Cohort
STANFORD_FAP_SCORES <- here("data", "outputs", "05", "stanford_fap_h_scores.csv")
STANFORD_FOXP3_DENSITIES <- here("data", "outputs", "08", "stanford_foxp3_densities.csv")

# . Histotypes
HISTOTYPES <- c("SBT", "LGSC")

# Helper Functions =============================================================
linear_regression <- function(tbl) {
  tbl |> 
    select(histotype, mean_fap_score, stroma, tumour) |>
    pivot_longer(cols = c(stroma, tumour), names_to = "compartment", values_to = "density") |>
    nest(data = c(histotype, mean_fap_score, density)) |>
    mutate(linear_model = map(data, \(tbl) glance(lm(log10(density + 1) ~ mean_fap_score, data = tbl)))) |>
    unnest(cols = linear_model)
}

correlation_test <- function(tbl) {
  tbl |>
    select(histotype, mean_fap_score, stroma, tumour) |>
    pivot_longer(cols = c(stroma, tumour), names_to = "compartment", values_to = "density") |>
    nest(data = c(histotype, mean_fap_score, density)) |>
    mutate(linear_model = map(data, \(tbl) tidy(cor.test(log10(tbl$density + 1), tbl$mean_fap_score, method = "pearson")))) |>
    unnest(cols = linear_model)
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
         stroma = mean_stroma_density,
         tumour = mean_tumour_density) |>
  filter(
    histotype %in% HISTOTYPES,
    !is.na(mean_fap_score),
    !(is.na(stroma) & is.na(tumour))
  ) |>
  mutate(cohort = "Proteomic")

# . Stanford
stanford_scores <- full_join(stanford_fap_scores, stanford_foxp3_densities, 
                             by = join_by(tma_number, case_code, histotype)) |>
  select(case_code, histotype, 
         mean_fap_score = mean_score, 
         stroma = mean_stroma_density, 
         tumour = mean_tumour_density) |>
  filter(
    histotype %in% HISTOTYPES,
    !is.na(mean_fap_score),
    !(is.na(stroma) & is.na(tumour))
  ) |>
  mutate(cohort = "Stanford")

# Linear Model =================================================================
# . Proteomic Cohort
proteomic_lm_results <- proteomic_scores |>
  select(histotype, mean_fap_score, stroma, tumour) |>
  linear_regression() |>
  select(-data)

# . Stanford Cohort
stanford_lm_results <- stanford_scores |>
  select(histotype, mean_fap_score, stroma, tumour) |>
  linear_regression() |>
  select(-data)

# Correlation Test =============================================================
# . Proteomic Cohort
proteomic_cor_results <- proteomic_scores |>
  select(histotype, mean_fap_score, stroma, tumour) |>
  correlation_test() |>
  select(-data)

# . Stanford Cohort
stanford_cor_results <- stanford_scores |>
  select(histotype, mean_fap_score, stroma, tumour) |>
  correlation_test() |>
  select(-data)

 # Save Results =================================================================
# . Mean FAP Score and FOXP3 Densities
# . . Proteomic Cohort
write_csv(proteomic_scores, file = here("data", "outputs", "10", "22-002_mean_fap_scores_and_foxp3_densities.csv"))

# . . Stanford Cohort
write_csv(stanford_scores, file = here("data", "outputs", "10", "stanford_mean_fap_scores_and_foxp3_densities.csv"))

# . Stroma Results - Linear Model
# . . Proteomic Cohort
proteomic_lm_results |>
  write_csv(file = here("data", "outputs", "10", "22-002_linear_model_results.csv"))

# . . Stanford Cohort
stanford_lm_results |>
  write_csv(file = here("data", "outputs", "10", "stanford_linear_model_results.csv"))

# . Stroma Results - Correlation Test
# . . Proteomic Cohort
proteomic_cor_results |>
  write_csv(file = here("data", "outputs", "10", "22-002_cor_test_results.csv"))

# . . Stanford Cohort
stanford_cor_results |>
  write_csv(file = here("data", "outputs", "10", "stanford_cor_test_results.csv"))