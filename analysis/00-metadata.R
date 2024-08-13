# Loading Dependencies =========================================================
# . Path Handling
library(here)

# . Read/Write Data
library(readr)
library(readxl)

# . Data Wrangling
library(dplyr)
library(tidyr)
library(janitor)

# . String Manipulation
library(stringr)

# Constants ====================================================================
# . File Paths
SAMPLES <- here("data", "metadata", "samples.csv")
CURATION <- here("data", "metadata", "curation.xlsx")
MUTATIONS <- here("data", "metadata", "genotyping", "mutations.csv")
HISTOTYPES <- here("data", "metadata", "validated_histotypes.xlsx")

# . Excel Sheet Names
CURATION_LGSC_SHEET <- "LGSC_final"
CURATION_SBOT_SHEET <- "SBOT_final"
CURATION_HGSC_SHEET <- "HGSC_final"

# . Histotype Abbreviations
HISTOTYPE_ABBREVIATIONS <- c("low-grade serous" = "LGSC",
                             "high-grade serous" = "HGSC",
                             "micropapillary SBOT" = "mSBT",
                             "serous borderline tumour" = "SBT")

# Helper Functions =============================================================
extract_voa_and_block <-  function(tbl, column) {
  tbl |>
    extract({{ column }}, 
            into = c("voa", "block"), 
            regex = "VOA(\\d+)(?:#)?(\\w)", 
            convert = TRUE)
}

# Load Data ====================================================================
# . Sample metadata
samples <- read_csv(SAMPLES) |>
  clean_names() |>
  select(set, sample, acc_num, 
         # Status
         stage_full, grade, status,
         # Treatment
         chemo, rt, 
         # Dates
         age_dx, age_surg, date_surg, 
         # Proteomics
         batch, total_abundance, exclude)

# . Validated histotypes
histotypes <- read_xlsx(HISTOTYPES) |>
  clean_names() |>
  select(-c(x1, exclude)) |>
  rename(voa = sample) |>
  mutate(voa = parse_number(voa))

# . Mutations
mutations <- read_csv(MUTATIONS) |>
  clean_names() |>
  mutate(mutation = case_when(
    mutation == "N/A" ~ "No Material",
    is.na(mutation) ~ "Not Assayed",
    .default = mutation
  ))

# . Sample Collection (ie. macrodissection)
lgsc_curation <- read_xlsx(CURATION, sheet = CURATION_LGSC_SHEET) |>
  clean_names() |>
  select(-c(tumor_type, tumor_bank_patient_id, comments, additional)) |>
  extract_voa_and_block(column = id) |>
  mutate(cellularity = as.character(cellularity))

sbot_curation <- read_xlsx(CURATION, sheet = CURATION_SBOT_SHEET) |>
  clean_names() |>
  select(-c(tumor_type, tumor_bank_patient_id, comment, additional)) |>
  extract_voa_and_block(column = id)

hgsc_curation <- read_xlsx(CURATION, sheet = CURATION_HGSC_SHEET) |>
  clean_names() |>
  select(-c(tumor_type)) |>
  extract_voa_and_block(column = id) |>
  rename(number_of_scrolls = p_number) |>
  mutate(cellularity = as.character(cellularity),
         amount = as.character(amount))

curation <- bind_rows(lgsc_curation, sbot_curation, hgsc_curation) |>
  mutate(sc_md = str_to_lower(sc_md),
         amount = parse_number(amount)) |>
  rename(amount_cm2 = amount)

# Processing Data ==============================================================
processed_metadata <- samples |>
  # Splitting VOA and Block
  extract_voa_and_block(column = sample) |>
  # Adding validated histotypes
  left_join(histotypes, by = join_by(set, voa)) |>
  mutate(histotype = recode(histotype, !!!HISTOTYPE_ABBREVIATIONS)) |>
  # Adding mutations
  left_join(mutations, by = join_by(set, voa, acc_num, block)) |>
  # Add curation details
  left_join(curation, by = join_by(voa, block))

# Save Metadata ================================================================
write_csv(processed_metadata, here("data", "metadata.csv"))
