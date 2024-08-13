# Loading Dependencies =========================================================
# . File Handling
library(here)

# . Read/Write Data
library(readr)
library(readxl)

# . Data Wrangling/Cleaning
library(dplyr)
library(tidyr)
library(purrr)
library(janitor)

# . Clean Pairwise Wilcoxon Rank Sum Test
library(broom)

# . String Manipulation
library(glue)
library(stringr)

# Constants ====================================================================
# . Proteomic TMA Map
PROTEOMIC_COHORT <- here("data", "IHC", "22-002.xlsx")
RANGE <- "P1:AB131"

# . . Macrophage Marker Counts
CD68 <- "CD68 (Almira)"
CD163 <- "CD163 (Almira)"

# Relevant Histotypes
HISTOTYPES <- c("SBT", "LGSC")

# Helper Functions =============================================================
pairwise_wilcox_test <- function(tbl, values, group) {
  column_name <- glue("test_{values}")
  
  tbl |>
    mutate({{ column_name }} := map(data, \(tbl) {
      data <- tbl |>
        select({{ group }}, {{ values }}) |>
        unnest(cols = {{ values }})
        
      pairwise.wilcox.test(data[[{{ values }}]], 
                           data[[{{ group }}]], 
                           p.adjust.method = "none",
                           exact = TRUE) |>
        tidy() |>
        clean_names()
    })) |>
    unnest(cols = {{ column_name }}) |>
    nest({{ column_name}} := c(group1, group2, p_value))
}

load_macrophage_density <- function(path, sheet, range) {
  read_xlsx(path, sheet = sheet, range = range) |>
    clean_names() |>
    select(core_id, voa, histotype, contains("density"), notes) |>
    mutate(marker = str_extract(sheet, pattern = "CD\\d+"),
           voa = parse_number(voa))
}

# Loading Data =================================================================
# . Proteomic Cohort
# . . CD68+ Macrophage Counts
cd68 <- load_macrophage_density(PROTEOMIC_COHORT, sheet = CD68, range = RANGE)
  
# . . CD163+ Macrophage density
cd163 <- load_macrophage_density(PROTEOMIC_COHORT, sheet = CD163,range = RANGE)

# Preprocessing Data ===========================================================
macrophages <- bind_rows(cd68, cd163) |>
  # Rename columns
  rename(stroma = stroma_density, tumour = tumour_density) |>
  # Recode histotype + parse VOA
  mutate(histotype = recode(histotype, mSBT = "SBT")) |>
  # Filter out missing notes
  # . If not present, exclude core
  filter(is.na(notes) | !str_detect(notes, "exclude")) |>
  # Filter for relevant histotypes
  filter(histotype %in% HISTOTYPES) |>
  # Mean density per core
  reframe(mean_stroma = mean(stroma, na.rm = TRUE), 
          mean_tumour = mean(tumour, na.rm = TRUE),
          .by = c(core_id, voa, histotype, marker)) |>
  # Pseudo log10 transformation
  mutate(across(contains("mean"), \(x) log10(x + 1), .names = "pseudo_log10_{.col}")) |>
  # Select relevant columns
  select(core_id, voa, histotype, marker, 
         stroma = pseudo_log10_mean_stroma, tumour = pseudo_log10_mean_tumour)

# Wilcoxon Rank Sum Test =======================================================
results <- macrophages |>
  select(marker, histotype, stroma, tumour) |>
  nest(data = c(histotype, stroma, tumour)) |>
  pairwise_wilcox_test(values = "stroma", group = "histotype") |>
  pairwise_wilcox_test(values = "tumour", group = "histotype")

# Save Results =================================================================
# . Macrophage Marker Densities
write_csv(macrophages, file = here("data", "outputs", "07", "22-002_macrophage_densities.csv"))

# . Stroma Results - Pairwise Wilcoxon Rank Sum Test
results |>
  unnest(test_stroma) |>
  select(marker, group1, group2, p_value) |>
  write_csv(file = here("data", "outputs", "07", "22-002_stroma_results.csv"))

# . Tumour Results - Pairwise Wilcoxon Rank Sum Test
results |>
  unnest(test_tumour) |>
  select(marker, group1, group2, p_value) |>
  write_csv(file = here("data", "outputs", "07", "22-002_tumour_results.csv"))
