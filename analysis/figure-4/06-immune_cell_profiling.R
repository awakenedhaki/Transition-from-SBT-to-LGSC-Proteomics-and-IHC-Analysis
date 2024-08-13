# Loading Dependencies =========================================================
# . Path Handling
library(here)

# . Read/Write Data
library(readr)
library(readxl)

# . Data Wrangling
library(dplyr)
library(tidyr)
library(purrr)
library(forcats)
library(janitor)

# . Clean Pairwise Wilcoxon Rank Sum Test
library(broom)

# . String Manipulation
library(glue)
library(stringr)

# Constants ====================================================================
# . File Paths
METADATA <- here("data", "IHC", "22-002.xlsx")
IMMUNE_PROFILE <- here("data", "mIHC", "immune_profile_multiplex.xlsx")

# . Cell Type Markers
CELL_TYPES <- c("cd68p"        = "Macrophage",
                "cd3p_cd8p"    = "CD8+ T cell", 
                "cd3p_cd8n"    = "CD8- T cell",
                "cd20p_cd79ap" = "B cell",
                "cd20n_cd79ap" = "Plasma cell")

# . Histotype Abbreviations
HISTOTYPES <- c("HGSC" = "#184E77", "SBT" = "#52B69A", "LGSC" = "#D9ED92")

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

# Loading Data =================================================================
metadata <- read_excel(METADATA, sheet = "Core IDs") |>
  clean_names() |>
  mutate(voa = parse_number(voa),
         column = LETTERS[column])

immune_profile <- read_excel(IMMUNE_PROFILE)

# Preprocessing Data ===========================================================
cell_type_densities <- immune_profile |>
  # Replace "+" and "-" with "p" and "n", respectively
  rename_with(.fn = \(colname) {
    colname |>
      str_replace_all(pattern = "\\+", replacement = "p") |>
      str_replace_all(pattern = "-", replacement = "n")
  }) |>
  # Clean column names
  clean_names() |>
  # Exclude noted cores
  filter(is.na(notes)) |>
  # Extract relevant columns
  select(core_id, voa, 
         row = tma_row, 
         column = tma_column, 
         stroma_area = stroma_area_mm2,
         tumour_area = tumour_area_mm2,
         matches("cells$")) |>
  # Label cores by histotype
  mutate(voa = parse_number(voa)) |>
  left_join(metadata, by = join_by(core_id, voa, column)) |>
  relocate(histotype, .after = column) |>
  mutate(histotype = recode(histotype, mSBT = "SBT"),
         histotype = fct_relevel(histotype, !!!names(HISTOTYPES))) |>
  filter(histotype %in% names(HISTOTYPES)) |>
  # Remove Core 55 (adipose tissue)
  filter(core_id != 55) |>
  # Convert character columns to double
  mutate(across(c("stroma_area", "tumour_area", matches("cells$")), parse_number)) |>
  # Extract compartment and cell type markers
  pivot_longer(cols = matches("cells$"),
               names_to = "population",
               values_to = "count") |>
  extract(col = "population", 
          into = c("compartment", "marker"), 
          regex = "(stroma|tumour)_(.*)_cells$") |>
  # Pivot longer compartment area
  pivot_longer(cols = contains("area"),
               names_to = "compartment_area",
               values_to = "area") |>
  mutate(compartment_area = str_remove(compartment_area, "_area")) |>
  filter(compartment == compartment_area) |>
  select(-compartment_area) |>
  # Label cell type names
  mutate(cell_type = recode(marker, !!!CELL_TYPES)) |>
  # Compartment area-based normalization
  mutate(density = count / area) |>
  reframe(density = mean(density, na.rm = TRUE), 
          .by  = c(core_id, voa, histotype, compartment, cell_type)) |>
  mutate(pseudo_log10_density = log10(density + 1))

# Wilcoxon Rank Sum Test =======================================================
results <- cell_type_densities |>
  # Selecting relevant columns
  select(histotype, compartment, cell_type, pseudo_log10_density) |>
  # Nesting density by cell_type
  nest(.by = cell_type) |>
  # Pivoting wider by compartment
  mutate(data = map(data, \(tbl) {
    tbl |>
      pivot_wider(names_from = compartment, 
                  values_from = pseudo_log10_density,
                  values_fn = list)
  })) |>
  # Run pairwise wilcoxon test
  pairwise_wilcox_test(values = "stroma", group = "histotype") |>
  pairwise_wilcox_test(values = "tumour", group = "histotype")

# Save Results =================================================================
# . Cell Type Densities
write_csv(cell_type_densities, file = here("data", "outputs", "06", "cell_type_densities.csv"))

# . Stroma Results - Pairwise Wilcoxon Rank Sum Test
results |>
  unnest(test_stroma) |>
  select(cell_type, group1, group2, p_value, adj_p_value) |>
  write_csv(file = here("data", "outputs", "06", "stroma_results.csv"))

# . Tumour Results - Pairwise Wilcoxon Rank Sum Test
results |>
  unnest(test_tumour) |>
  select(cell_type, group1, group2, p_value, adj_p_value) |>
  write_csv(file = here("data", "outputs", "06", "tumour_results.csv"))
