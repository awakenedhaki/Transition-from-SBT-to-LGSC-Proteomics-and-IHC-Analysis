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
# . File Paths
STANFORD <- here("data", "IHC", "stanford.xlsx")

# . Tumour Area Paths
TMA1_AREAS <- here("data", "IHC", "raw", "stanford", "tumour_areas", "TMA1-PanCK-trial1_data.tsv")
TMA2_AREAS <- here("data", "IHC", "raw", "stanford", "tumour_areas", "TMA2-PanCK-trial1_data.tsv")
LGSC_AREAS <- here("data", "IHC", "raw", "stanford", "tumour_areas", "Mixed-TMA-PanCK-trial1_data.tsv")

AREA <- c("TMA 1" = TMA1_AREAS, "TMA 2" = TMA2_AREAS, "LGS/Mixed" = LGSC_AREAS)

# . Sheets
# . . Metadata
METADATA <- "Metadata"

# . . FOXP3
TMA1_FOXP3 <- here("data", "IHC", "raw", "stanford", "tregs", "TMA1.txt")
TMA2_FOXP3 <- here("data", "IHC", "raw", "stanford", "tregs", "TMA2.txt")
LGSC_FOXP3 <- here("data", "IHC", "raw", "stanford", "tregs", "LGS-Mixed.txt")

FOXP3 <- c("TMA 1" = TMA1_FOXP3, 
           "TMA 2" = TMA2_FOXP3, 
           "LGS/Mixed" = LGSC_FOXP3)

# Helper Functions =============================================================
clean_tma_areas <- function(tbl, tma) {
  tma_metadata <- metadata |>
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
    mutate(adj_p_value = p.adjust(p_value, method = "BH")) |>
    nest({{ column_name}} := c(group1, group2, p_value, adj_p_value))
}

# Loading Data =================================================================
# . Metadata
metadata <- read_xlsx(STANFORD, sheet = METADATA, range = "A1:I268") |>
  clean_names() |>
  select(tma_number, row = row_number, column = column_number, 
         case_code, histotype = diagnosis, exclude) |>
  mutate_at(.vars = c("row", "column"), .funs = parse_number)

# . Tumour Areas
areas <- tibble()
for (tma in names(AREA)) {
  area <- read_tsv(AREA[tma]) |>
    clean_names(replace = c("µ" = "u")) |>
    clean_tma_areas(tma)
  
  areas <- bind_rows(areas, area)
}

# . FOXP3
counts <- tibble()
for (tma in names(FOXP3)) {
  count <- read_tsv(FOXP3[tma]) |>
    clean_names(replace = c("µ" = "u")) |>
    select(tma_number = image, classification, core = tma_core) |>
    count(tma_number, classification, core) |>
    separate_wider_delim(cols = core, delim = "-", names = c("row", "column")) |>
    mutate(column = parse_number(column),
           row = map_dbl(row, \(x) which(LETTERS == x)),
           tma_number = case_when(
             tma_number == "TMA1" ~ "TMA 1",
             tma_number == "TMA2" ~ "TMA 2",
             .default = tma_number
           )) 
  
  counts <- bind_rows(counts, count)
}

# Preprocessing Data ===========================================================
foxp3_density <- counts |>
  # Annotate FOXP3 counts with Stanford TMA metadata
  full_join(metadata, by = join_by(tma_number, row, column)) |>
  # Annotate FOXP3 counts with tissue compartment areas
  full_join(areas, by = join_by(tma_number, row, column, histotype, case_code, exclude)) |>
  # Exclude cores
  filter(!exclude, !is.na(stroma_area), !is.na(tumour_area), !is.na(classification)) |>
  select(-exclude) |>
  # Relocating columns
  relocate(tma_number, row, column, case_code, histotype, n) |>
  # Pivoting data to wide format
  pivot_wider(names_from = classification, values_from = n) |>
  # Normalizing FOXP3 counts by tissue compartment area
  rename(stroma = Stroma, tumour = Tumor) |>
  mutate(stroma_density = stroma / stroma_area,
         tumour_density = tumour / tumour_area) |>
  # Selecting relevant columns
  select(tma_number, case_code, histotype, stroma_density, tumour_density) |>
  # Calculating the mean density by tissue compartment
  reframe(mean_stroma_density = mean(stroma_density, na.rm = TRUE),
          mean_tumour_density = mean(tumour_density, na.rm = TRUE),
          .by = c(tma_number, case_code, histotype)) |>
  # Abbreviating histotype names
  mutate(histotype = case_when(
    histotype == "Low grade serous carcinoma" ~ "LGSC",
    histotype == "Serous borderline tumor" ~ "SBT",
    .default = histotype
  )) |>
  # Correct TMA names
  mutate(tma_number = case_match(tma_number,
    "LGS/Mixed" ~ "LGSC/Mixed",
    "TMA 1" ~ "TMA1",
    "TMA 2" ~ "TMA2",
    .default = tma_number
  ))
  
# Wilcoxon Rank Sum Test =======================================================
results <- foxp3_density |>
  mutate(across(contains("density"), \(x) log10(x + 1), .names = "pseudo_log10_{.col}")) |>
  select(histotype, stroma = pseudo_log10_mean_stroma_density, tumour = pseudo_log10_mean_tumour_density) |>
  pivot_longer(cols = c(stroma, tumour), names_to = "compartment", values_to = "density") |>
  nest(data = -compartment) |>
  pairwise_wilcox_test(values = "density", group = "histotype")

# Save Results =================================================================
# . FOXP3 Densities
write_csv(foxp3_density, file = here("data", "outputs", "08", "stanford_foxp3_densities.csv"))

# . Stroma Results - Pairwise Wilcoxon Rank Sum Test
results |>
  filter(compartment == "stroma") |>
  unnest(test_density) |>
  select(group1, group2, p_value) |>
  write_csv(file = here("data", "outputs", "08", "stanford_stroma_results.csv"))

# . Tumour Results - Pairwise Wilcoxon Rank Sum Test
results |>
  filter(compartment == "tumour") |>
  unnest(test_density) |>
  select(group1, group2, p_value) |>
  write_csv(file = here("data", "outputs", "08", "stanford_tumour_results.csv"))

