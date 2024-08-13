# Loading Dependencies =========================================================
# . Path Handling
library(here)

# . Read/Write Data
library(readr)
library(readxl)

# . Data Wrangling
library(dplyr)
library(purrr)
library(tidyr)
library(forcats)
library(janitor)

# . Visualization
library(ggplot2)

# ggplot2 Defaults =============================================================
theme_set(theme_bw())

# Constants ====================================================================
# . Proteomic Metadata
PROTEOMIC_METADATA <- here("data", "IHC", "22-002.xlsx")
SHEET <- "Core IDs"

# . Proteomic Compartment Areas
IMMUNE_PROFILE <- here("data", "mIHC", "immune_profile_multiplex.xlsx")

# . Proteomic Histotypes
PROTEOMIC_HISTOTYPE_COLORS <- c("HGSC" = "#184E77", SBT = "#52B69A", "LGSC" = "#D9ED92")

# . Stanford Metadata
STANFORD <- here("data", "IHC", "stanford.xlsx")
METADATA <- "Metadata"

# . Stanford Compartment Areas
TMA1_AREAS <- here("data", "IHC", "raw", "stanford", "tumour_areas", "TMA1-PanCK-trial1_data.tsv")
TMA2_AREAS <- here("data", "IHC", "raw", "stanford", "tumour_areas", "TMA2-PanCK-trial1_data.tsv")
LGSC_AREAS <- here("data", "IHC", "raw", "stanford", "tumour_areas", "Mixed-TMA-PanCK-trial1_data.tsv")

AREA <- c("TMA 1" = TMA1_AREAS, "TMA 2" = TMA2_AREAS, "LGS/Mixed" = LGSC_AREAS)

# . Stanford Histotypes
STANFORD_HISTOTYPE_COLORS <- c(SBT = "#52B69A", "LGSC" = "#D9ED92")

# Helper Functions =============================================================
clean_tma_areas <- function(tbl, tma) {
  tma_metadata <- stanford_metadata |>
    filter(tma_number == tma)
  
  tbl |>
    filter(class == "Region*") |>
    select(tma_core,
           stroma_area = classifier_stroma_area_um_2,
           tumour_area = classifier_tumor_area_um_2) |>
    mutate_if(.predicate = is.double, .funs = \(x) x / 1e6) |>
    separate_wider_delim(cols = tma_core, delim = "-", names = c("row", "column")) |>
    mutate(column = parse_number(column),
           row = map_dbl(row, \(x) which(LETTERS == x) - 1)) |>
    filter(row != 0) |>
    left_join(tma_metadata, by = join_by(row, column))
}

# Loading Data =================================================================
# . Proteomic Metadata
proteomic_metadata <- read_excel(PROTEOMIC_METADATA, sheet = SHEET) |>
  clean_names() |> 
  mutate(row = case_when(
    sector == 2 ~ row + 1,
    .default = row
  ))

# . Proteomic Compartmnet Areas
proteomic_areas <- read_xlsx(IMMUNE_PROFILE) |>
  clean_names() |>
  select(core_id, voa, 
         row = tma_row, 
         column = tma_column, 
         stroma_area = stroma_area_mm2,
         tumour_area = tumour_area_mm2) |>
  mutate(column = map_dbl(column, \(x) which(LETTERS == x)),
         sector = case_when(
           row < 8 ~ 1,
           row > 8 ~ 2
         )) |>
  mutate(across(contains("area"), parse_number)) |>
  full_join(proteomic_metadata, by = join_by(voa, core_id, sector, row, column))

# . Stanford Metadata
stanford_metadata <- read_xlsx(STANFORD, sheet = METADATA, range = "A1:I268") |>
  clean_names() |>
  select(tma_number, row = row_number, column = column_number, 
         case_code, histotype = diagnosis, exclude) |>
  mutate_at(.vars = c("row", "column"), .funs = parse_number)

# . Stanford Compartment Areas
stanford_areas <- tibble()
for (tma in names(AREA)) {
  stanford_area <- read_tsv(AREA[tma]) |>
    clean_names(replace = c("µ" = "u")) |>
    clean_tma_areas(tma) |>
    filter(!exclude)
  
  stanford_areas <- bind_rows(stanford_areas, stanford_area)
}

# Visualization ================================================================
# . Proteomic Compartment Areas
proteomic_areas |>
  select(histotype, stroma_area, tumour_area) |>
  filter(histotype %in% names(PROTEOMIC_HISTOTYPE_COLORS)) |>
  pivot_longer(cols = c(stroma_area, tumour_area), names_to = "compartment", values_to = "area", names_pattern = "(\\w+)_area") |>
  mutate(histotype = fct_relevel(histotype, !!!names(PROTEOMIC_HISTOTYPE_COLORS))) |>
  ggplot(aes(x = histotype, y = area, fill = histotype)) +
    geom_boxplot(show.legend = FALSE) +
    facet_wrap(~compartment) +
    scale_fill_manual(values = PROTEOMIC_HISTOTYPE_COLORS)

# . Stanford Compartment Areas
stanford_areas |>
  select(histotype, stroma_area, tumour_area) |>
  mutate(histotype = recode(histotype, `Serous borderline tumor` = "SBT", `Low grade serous carcinoma` = "LGSC")) |>
  pivot_longer(cols = c(stroma_area, tumour_area), names_to = "compartment", values_to = "area", names_pattern = "(\\w+)_area") |>
  ggplot(aes(x = histotype, y = area, fill = histotype)) +
    geom_boxplot(show.legend = FALSE) +
    facet_wrap(~compartment) +
    scale_fill_manual(values = STANFORD_HISTOTYPE_COLORS)

# Save Plot ====================================================================
ggsave(here("figures", "supplemental", "figure-S11", "figure-S11_E.svg"), 
       plot = last_plot(), width = 8, height = 5, dpi = 1200)
