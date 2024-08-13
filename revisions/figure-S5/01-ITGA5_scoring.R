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
# PROTOEMIC_SCORE_SHEET <- "FAP (Lars)"
PROTEOMIC_RANGE <- "P1:W131"
PROTEOMIC_HISTOTYPES <- c("SBT", "LGSC")
PROTEOMIC_SCORE_SHEET <- "ITGA5 (Felix)"
PROTEOMIC_SCORE_RANGE <- "P1:V131"

# . Stanford TMA Map
STANFORD <- here("data", "IHC", "stanford.xlsx")
STANFORD_METADATA <- "Metadata"
STANFORD_SCORE_SHEETS <- c("TMA1 ITGA5 (Felix)", 
                           "TMA2 ITGA5 (Felix)", 
                           "LGSC & Mixed ITGA5 (Felix)")
STANFORD_RANGE <- "O1:R90"
STANFORD_HISTOTYPES <- c("SBT", "LGSC")

# Helper Functions =============================================================
summarize_itga5_score <- function(tbl, .by, .factor_levels) {
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
    complete(histotype, h_score, fill = list(n = 0)) |>
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
  select(tma = tma_number, 
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
  mutate(tma = case_match(tma,
                          "LGS/Mixed" ~ "LGSC/Mixed",
                          "TMA 1" ~ "TMA1",
                          "TMA 2" ~ "TMA2",
                          .default = tma
  ))

# Loading Data =================================================================
# . Proteomic Cohort
proteomic_cohort <- read_excel(PROTEOMIC_COHORT, 
                               sheet = PROTEOMIC_SCORE_SHEET, 
                               range = PROTEOMIC_SCORE_RANGE, 
                               col_types = c("numeric",      # Sector
                                             "numeric",      # Row
                                             "numeric",      # Column
                                             "numeric",      # Core ID
                                             "text",         # VOA
                                             "text",         # Histotype
                                             "numeric")) |>  # Score
  clean_names()

# . Stanford Cohort
stanford_cohort <- tibble()
for (sheet in STANFORD_SCORE_SHEETS) {
  # Loading Stanford TMA scoring sheet
  score_sheet <- read_excel(STANFORD, sheet = sheet, range = STANFORD_RANGE) |>
    clean_names() |>
    mutate(tma = (sheet |>
                    str_remove(pattern = " ITGA5 \\(Felix\\)") |>
                    str_replace_all(pattern = " ", replacement = "") |>
                    str_replace(pattern = "&", replacement = "/"))) |>
    rename(row = row_number, column = column_number) |>
    mutate_at(.vars = c("row", "column", "score"), .funs = parse_number)
  
  # Append to Stanford Cohort tibble
  stanford_cohort <- bind_rows(stanford_cohort, score_sheet)
}

# Annotating Data ==============================================================
# . Stanford Cohort
annotated_stanford_cohort <- stanford_cohort |>
  full_join(stanford_metadata, by = join_by(tma, row, column)) |>
  select(tma, case_code, histotype, score)

# Preprocessing Data ===========================================================
# . Proteomic Cohort
proteomic_itga5_h_scores <- proteomic_cohort |>
  summarize_itga5_score(.by = c(core_id, histotype), 
                        .factor_levels = PROTEOMIC_HISTOTYPES)
# . Stanford Cohort
stanford_itga5_h_scores <- annotated_stanford_cohort |>
  summarize_itga5_score(.by = c(case_code, histotype), 
                        .factor_levels = STANFORD_HISTOTYPES)

# Fisher Exact Test ============================================================
# . Proteomic Cohort
proteomic_result <- proteomic_itga5_h_scores |>
  drop_na() |>
  fisher_exact_test()

# . Stanford Cohort
stanford_result <- stanford_itga5_h_scores |>
  fisher_exact_test()

# Save Preprocessed Data =======================================================
# . Stanford Cohort
write_csv(proteomic_itga5_h_scores, 
          here("data", "outputs", "00-revisions", "proteomic_itga5_h_scores.csv"))

# . Stanford Cohort
write_csv(stanford_itga5_h_scores, 
          here("data", "outputs", "00-revisions", "stanford_itga5_h_scores.csv"))
