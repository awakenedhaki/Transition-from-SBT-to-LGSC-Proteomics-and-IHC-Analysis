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
FOXP3 <- here("data", "IHC", "raw", "22-002", "tregs", "FOXP3_counts.tsv")

# . Proteomic TMA Map
PROTEOMIC_METADATA <- here("data", "IHC", "22-002.xlsx")
SHEET <- "Core IDs"

# . Proteomic Core Areas
TUMOUR_AREAS <- here("data", "mIHC", "immune_profile_multiplex.xlsx")

# . Histotypes
HISTOTYPES <- c("HGSC", "SBT", "LGSC")

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
    mutate(adj_p_value = p.adjust(p_value, method = "BH")) |>
    nest({{ column_name}} := c(group1, group2, p_value, adj_p_value))
}

# Loading Metadata =============================================================
metadata <- read_excel(PROTEOMIC_METADATA, sheet = SHEET) |>
  clean_names() |> 
  mutate(row = case_when(
    sector == 2 ~ row + 1,
    .default = row
  ))
  
# Loading Data =================================================================
# . Area
area <- read_excel(TUMOUR_AREAS) |>
  # Clean column names
  clean_names() |>
  # Select relevant column
  select(core_id, row = tma_row, column = tma_column, 
         stroma_area = stroma_area_mm2, tumour_area = tumour_area_mm2) |>
  # Extract row and column information
  mutate(column = map_dbl(column, \(x) which(LETTERS == x)),
         sector = case_when(
           row < 8 ~ 1,
           row > 8 ~ 2
         )) |>
  # Relocate columns
  relocate(core_id, sector, row, column) |>
  # Parse numeric values
  mutate_at(vars(stroma_area, tumour_area), parse_number)

# . FOXP3 Counts
counts <- read_delim(FOXP3, delim = "\t", col_select = c("class", "name")) |>
  # Clean column names
  rename(compartment = class, core = name) |>
  # Label explicit missing values as "Stroma"
  replace_na(list(compartment = "Stroma")) |>
  # Count number of cells per core per compartment
  count(core, compartment) |>
  # Assign a value of 0 to missing compartments in counts table
  complete(core, compartment, fill = list(n = 0)) |>
  # Extract sector, row, column information
  extract(col = core, 
          into = c("sector", "row", "column"), 
          regex = "(\\d)(\\d)(\\d{2})", 
          convert = TRUE) |>
  # Adjust row number for sector 2
  mutate(row = case_when(
    sector == 2 ~ row + 8,
    .default = row
  )) |>
  # Merge with metadata to get histotype information
  full_join(metadata, by = join_by(sector, row, column)) |>
  mutate(histotype = ifelse(histotype == "mSBT", "SBT", histotype)) |>
  # Filter out irrelevant histotypes and missing values
  filter(
    histotype %in% HISTOTYPES,
    !is.na(n)
  ) |>
  # Relocate columns
  relocate(core_id, voa, sector, row, column, histotype)

# Preprocessing Data ===========================================================
foxp3_density <- counts |>
  # Convert compartment to lower case
  mutate(compartment = map_chr(compartment, \(x) str_to_lower(x))) |>
  # Rename compartment
  mutate(compartment = ifelse(compartment == "tumor", "tumour", compartment)) |>
  # Pivot data
  pivot_wider(names_from = compartment, values_from = n) |>
  # Merge with area data
  left_join(area, by = join_by(core_id, sector, row, column)) |>
  # Calculate density
  mutate(stroma_density = stroma / stroma_area,
         tumour_density = tumour / tumour_area) |>
  # Remove irrelevant columns
  select(-c(voa, sector, row, column, stroma, tumour, contains("area"))) |>
  # Calculating the mean density by tissue compartment
  reframe(mean_stroma_density = mean(stroma_density, na.rm = TRUE),
          mean_tumour_density = mean(tumour_density, na.rm = TRUE), 
          .by = c(core_id, histotype))

# Wilcoxon Rank Sum Test =======================================================
results <- foxp3_density |>
  mutate(across(contains("density"), \(x) log10(x + 1), .names = "pseudo_log10_{.col}")) |>
  select(histotype, stroma = pseudo_log10_mean_stroma_density, tumour = pseudo_log10_mean_tumour_density) |>
  pivot_longer(cols = c(stroma, tumour), names_to = "compartment", values_to = "density") |>
  nest(data = -compartment) |>
  pairwise_wilcox_test(values = "density", group = "histotype")

 # Save Results =================================================================
# . FOXP3 Densities
write_csv(foxp3_density, file = here("data", "outputs", "08", "22-002_foxp3_densities.csv"))

# . Stroma Results - Pairwise Wilcoxon Rank Sum Test
results |>
  filter(compartment == "stroma") |>
  unnest(test_density) |>
  select(group1, group2, p_value) |>
  write_csv(file = here("data", "outputs", "08", "22-002_stroma_results.csv"))

# . Tumour Results - Pairwise Wilcoxon Rank Sum Test
results |>
  filter(compartment == "tumour") |>
  unnest(test_density) |>
  select(group1, group2, p_value) |>
  write_csv(file = here("data", "outputs", "08", "22-002_tumour_results.csv"))
