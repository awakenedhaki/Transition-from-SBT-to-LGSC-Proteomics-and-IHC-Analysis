# Loading Dependencies =========================================================
# . File Handling
library(here)

# . Read/Write Data
library(readr)
library(readxl)

# Clean Fisher Exact Test Outputs
library(broom)

# . Data Wrangling
library(dplyr)
library(tidyr)
library(tibble)
library(forcats)
library(janitor)

# . String Manipulation
library(stringr)

# Constants ====================================================================
# . Proteomic TMA Map
PROTEOMIC_COHORT <- here("data", "IHC", "22-002.xlsx")
PROTOEMIC_SCORE_SHEET <- "FAP (Lars)"
PROTEOMIC_RANGE <- "P1:W131"
PROTEOMIC_HISTOTYPES <- c("HGSC", "SBT", "LGSC")

# . Stanford TMA Map
STANFORD <- here("data", "IHC", "stanford.xlsx")
STANFORD_METADATA <- "Metadata"
STANFORD_SCORE_SHEETS <- c("TMA1 FAP (Lars)", 
                           "TMA2 FAP (Lars)", 
                           "LGSC & Mixed FAP (Lars)")
STANFORD_RANGE <- "O1:R90"
STANFORD_HISTOTYPES <- c("SBT", "LGSC")

#. Whole Sections
WHOLE_SECTIONS <- here("data", "IHC", "whole_tissue_sections.xlsx")

# Helper Functions =============================================================
summarize_fap_score <- function(tbl, .by, .factor_levels) {
  tbl |>
    # Summarize FAP scores
    reframe(mean_score = mean(score, na.rm = TRUE), 
            .by = {{ .by }}) |>
    # Categorize FAP scores
    mutate(h_score = case_when(
      mean_score < 10 ~ "Negative",
      mean_score >= 10 ~ "Positive"
    )) |>
    # Filter for histotypes of interest
    filter(histotype %in% {{ .factor_levels }}) |>
    # Level histotypes
    mutate(histotype = fct_relevel(histotype, !!!{{ .factor_levels }}))
}

fisher_exact_test <- function(tbl) {
  tbl |>
    count(histotype, h_score) |>
    select(histotype, h_score, n) |>
    # Pivot table
    pivot_wider(names_from = histotype, values_from = n) |>
    # Fisher Exact Test
    column_to_rownames("h_score") |> 
    as.matrix() |>
    fisher.test() |>
    tidy()
}

# Loading Metadata =============================================================
# . Stanford Metadata
stanford_metadata <- read_excel(STANFORD, sheet = STANFORD_METADATA) |>
  clean_names() |>
  select(tma_number, 
         row = row_number, 
         column = column_number, 
         case_code, 
         histotype = diagnosis) |>
  mutate_at(.vars = c("row", "column"), .funs = parse_number) |>
  mutate(histotype = case_match(histotype,
    "Low grade serous carcinoma" ~ "LGSC",
    "Serous borderline tumor" ~ "SBT",
    .default = histotype
  )) |>
  mutate(tma_number = case_match(tma_number,
    "LGS/Mixed" ~ "LGSC/Mixed",
    "TMA 1" ~ "TMA1",
    "TMA 2" ~ "TMA2",
    .default = tma_number
  ))

# Loading Data =================================================================
# . Proteomic Cohort
proteomic_cohort <- read_excel(PROTEOMIC_COHORT, 
                               sheet = PROTOEMIC_SCORE_SHEET, 
                               range = PROTEOMIC_RANGE) |>
  clean_names() |>
  mutate_at(.vars = c("voa", "score"), .funs = parse_number) |>
  select(-notes)
  
# . Stanford Cohort
stanford_cohort <- tibble()
for (sheet in STANFORD_SCORE_SHEETS) {
  # Loading Stanford TMA scoring sheet
  score_sheet <- read_excel(STANFORD, sheet = sheet, range = STANFORD_RANGE) |>
    clean_names() |>
    mutate(tma_number = (sheet |>
                           str_remove(pattern = " FAP \\(Lars\\)") |>
                           str_replace_all(pattern = " ", replacement = "") |>
                           str_replace(pattern = "&", replacement = "/"))) |>
    rename(row = row_number, column = column_number) |>
    mutate_at(.vars = c("row", "column", "score"), .funs = parse_number)
  
  # Append to Stanford Cohort tibble
  stanford_cohort <- bind_rows(stanford_cohort, score_sheet)
}

# . Proteomic Whole Sections
proteomic_whole_sections <- read_excel(WHOLE_SECTIONS) |>
    clean_names() |>
    filter(tma == "22-002") |>
    select(-c(tma, notes)) |>
    pivot_longer(cols = contains("region"),
                 names_pattern = "region_(\\d)",
                 names_to = "region",
                 values_to = "score")

# Annotating Data ==============================================================
# . Stanford Cohort
annotated_stanford_cohort <- stanford_cohort |>
  full_join(stanford_metadata, by = join_by(tma_number, row, column)) |>
  select(tma_number, case_code, histotype, score)

# Preprocessing Data ===========================================================
# . Proteomic Cohort
proteomic_fap_h_scores <- proteomic_cohort |>
  summarize_fap_score(.by = c(core_id, histotype),
                      .factor_levels = PROTEOMIC_HISTOTYPES) |>
  # Remove 1 LGSC case
  # . Both cores were adipose tissue
  filter(!is.na(h_score))

# . Stanford Cohort
stanford_fap_h_scores <- annotated_stanford_cohort |>
  summarize_fap_score(.by = c(tma_number, case_code, histotype), 
                      .factor_levels = STANFORD_HISTOTYPES)

# . Proteomic Whole Sections
whole_sections_fap_h_scores <- proteomic_whole_sections |>
  summarize_fap_score(.by = c(slide_id, histotype), 
                      .factor_levels = c("SBT", "LGSC"))

# Fisher Exact Test ============================================================
# . Proteomic Cohort
proteomic_result <- proteomic_fap_h_scores |>
  fisher_exact_test()

# . Stanford Cohort
stanford_result <- stanford_fap_h_scores |>
  fisher_exact_test()

# . Proteomic Whole Sections
whole_sections_result <- whole_sections_fap_h_scores |>
  fisher_exact_test()

# Save Preprocessed Data =======================================================
# . Proteomic Cohort
write_csv(proteomic_fap_h_scores, 
          here("data", "outputs", "05", "proteomic_fap_h_scores.csv"))

# . Stanford Cohort
write_csv(stanford_fap_h_scores, 
          here("data", "outputs", "05", "stanford_fap_h_scores.csv"))

# . Proteomic Whole Sections
write_csv(whole_sections_fap_h_scores, 
          here("data", "outputs", "05", "whole_sections_fap_h_scores.csv"))
