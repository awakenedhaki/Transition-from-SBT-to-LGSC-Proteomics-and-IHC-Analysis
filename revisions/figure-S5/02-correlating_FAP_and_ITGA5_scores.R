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

# . Visualization
library(scales)
library(ggplot2)

# ggplot2 Defaults =============================================================
update_geom_defaults("point", list(size = 3))
theme_set(theme_bw())

# Constants ====================================================================
# . Proteomic TMA Map
PROTEOMIC_COHORT <- here("data", "IHC", "22-002.xlsx")
PROTEOMIC_RANGE <- "P1:W131"
PROTEOMIC_HISTOTYPES <- c("SBT", "LGSC")

# . . Scoring Sheets
PROTEOMIC_SCORE_RANGE <- "P1:V131"
PROTOEMIC_FAP_SCORE_SHEET <- "FAP (Lars)"
PROTEOMIC_ITGA5_SCORE_SHEET <- "ITGA5 (Felix)"

# . Stanford TMA Map
STANFORD <- here("data", "IHC", "stanford.xlsx")
STANFORD_METADATA <- "Metadata"
STANFORD_HISTOTYPES <- c("SBT", "LGSC")

# . . Scoring Sheets

STANFORD_RANGE <- "O1:R90"
STANFORD_FAP_SCORE_SHEETS <- c("TMA1 FAP (Lars)", 
                           "TMA2 FAP (Lars)", 
                           "LGSC & Mixed FAP (Lars)")
STANFORD_ITGA5_SCORE_SHEETS <- c("TMA1 ITGA5 (Felix)", 
                           "TMA2 ITGA5 (Felix)", 
                           "LGSC & Mixed ITGA5 (Felix)")

# . Aesthetics
HISTOTYPE_COLORS <- c("SBT" = "#52B69A", "LGSC" = "#D9ED92")

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
# . . FAP
fap_proteomic_cohort <- read_excel(PROTEOMIC_COHORT, 
                                   sheet = PROTOEMIC_FAP_SCORE_SHEET, 
                                   range = PROTEOMIC_RANGE) |>
  clean_names() |>
  mutate_at(.vars = c("voa", "score"), .funs = parse_number) |>
  select(-notes) |>
  rename(fap = score)

# . . ITGA5
itga5_proteomic_cohort <- read_excel(PROTEOMIC_COHORT, 
                                     sheet = PROTEOMIC_ITGA5_SCORE_SHEET, 
                                     range = PROTEOMIC_SCORE_RANGE, 
                                     col_types = c("numeric",      # Sector
                                                   "numeric",      # Row
                                                   "numeric",      # Column
                                                   "numeric",      # Core ID
                                                   "text",         # VOA
                                                   "text",         # Histotype
                                                   "numeric")) |>  # Score
  clean_names() |>
  mutate(voa = parse_number(voa)) |>
  rename(itga5 = score)

# . . Merged Proteomic Scores
proteomic_cohort <- full_join(x = fap_proteomic_cohort,
                              y = itga5_proteomic_cohort,
                              by = join_by(sector, row, column, core_id, voa, histotype))

# . Stanford Cohort
# . . FAP
fap_stanford_cohort <- tibble()
for (sheet in STANFORD_FAP_SCORE_SHEETS) {
  # Loading Stanford TMA scoring sheet
  score_sheet <- read_excel(STANFORD, sheet = sheet, range = STANFORD_RANGE) |>
    clean_names() |>
    mutate(tma = (sheet |>
                    str_remove(pattern = " FAP \\(Lars\\)") |>
                    str_replace_all(pattern = " ", replacement = "") |>
                    str_replace(pattern = "&", replacement = "/"))) |>
    rename(row = row_number, column = column_number) |>
    mutate_at(.vars = c("row", "column", "score"), .funs = parse_number) |>
    select(row, column, tma, fap = score)
  
  # Append to Stanford Cohort tibble
  fap_stanford_cohort <- bind_rows(fap_stanford_cohort, score_sheet)
}

# . . ITGA5
itga5_stanford_cohort <- tibble()
for (sheet in STANFORD_ITGA5_SCORE_SHEETS) {
  # Loading Stanford TMA scoring sheet
  score_sheet <- read_excel(STANFORD, sheet = sheet, range = STANFORD_RANGE) |>
    clean_names() |>
    mutate(tma = (sheet |>
                    str_remove(pattern = " ITGA5 \\(Felix\\)") |>
                    str_replace_all(pattern = " ", replacement = "") |>
                    str_replace(pattern = "&", replacement = "/"))) |>
    rename(row = row_number, column = column_number) |>
    mutate_at(.vars = c("row", "column", "score"), .funs = parse_number) |>
    select(row, column, tma, itga5 = score)
  
  # Append to Stanford Cohort tibble
  itga5_stanford_cohort <- bind_rows(itga5_stanford_cohort, score_sheet)
}

# . . Merged Stanford
stanford_cohort <- full_join(x = fap_stanford_cohort,
                             y = itga5_stanford_cohort,
                             by = join_by(tma, row, column))

# Annotating Data ==============================================================
# . Stanford Cohort
annotated_stanford_cohort <- stanford_cohort |>
  full_join(stanford_metadata, by = join_by(tma, row, column)) |>
  select(tma, case_code, histotype, fap, itga5)

# Combined Scores ==============================================================
combined_scores <- bind_rows(
  (proteomic_cohort |>
     reframe(mean_fap_score = mean(fap, na.rm = TRUE),
             mean_itga5_score = mean(itga5, na.rm = TRUE),
             .by = c(core_id, voa, histotype)) |>
     filter(histotype %in% PROTEOMIC_HISTOTYPES) |>
     select(histotype, mean_fap_score, mean_itga5_score) |>
     mutate(cohort = "Proteomic")),
  (annotated_stanford_cohort |>
     reframe(mean_fap_score = mean(fap, na.rm = TRUE),
             mean_itga5_score = mean(itga5, na.rm = TRUE),
             .by = c(tma, case_code, histotype)) |>
     filter(histotype %in% STANFORD_HISTOTYPES) |>
     select(histotype, mean_fap_score, mean_itga5_score) |>
     mutate(cohort = "Stanford"))
)

# Visualization ================================================================
# . Proteomic Cohort
proteomic_cohort |>
  reframe(mean_fap_score = mean(fap, na.rm = TRUE),
          mean_itga5_score = mean(itga5, na.rm = TRUE),
          .by = c(core_id, voa, histotype)) |>
  filter(histotype %in% PROTEOMIC_HISTOTYPES) |>
  ggplot(aes(x = mean_fap_score, y = mean_itga5_score)) +
    geom_point(aes(color = histotype)) +
    geom_smooth(method = "lm") +
    labs(x = "Mean FAP Score", y = "Mean ITGA5 (CD49e) Score", color = "Histotype") +
    scale_color_manual(values = HISTOTYPE_COLORS)

# . Stanford Cohort
annotated_stanford_cohort |>
  reframe(mean_fap_score = mean(fap, na.rm = TRUE),
          mean_itga5_score = mean(itga5, na.rm = TRUE),
          .by = c(tma, case_code, histotype)) |>
  filter(histotype %in% STANFORD_HISTOTYPES) |>
  ggplot(aes(x = mean_fap_score, y = mean_itga5_score)) +
    geom_point(aes(color = histotype)) +
    geom_smooth(method = "lm") +
    labs(x = "Mean FAP Score", y = "Mean ITGA5 (CD49e) Score", color = "Histotype") +
    scale_color_manual(values = HISTOTYPE_COLORS)

# . Combined Scores
combined_scores |>
  ggplot(aes(x = mean_fap_score, y = mean_itga5_score)) +
    geom_point(aes(color = histotype)) +
    geom_smooth(method = "lm") +
    facet_wrap(~cohort) +
    labs(x = "Mean FAP Score", y = "Mean ITGA5 (CD49e) Score", color = "Histotype") +
    scale_color_manual(values = HISTOTYPE_COLORS) +
    theme(legend.position = "bottom")
    
# Save Plot ====================================================================
ggsave(here("figures", "revisions", "figure-S5", "figure-S5_D.svg"), 
       plot = last_plot(), width = 8, height = 5, dpi = 1200)
