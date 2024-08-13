# Loading Dependencies =========================================================
# . Path Handling
library(here)

# . Read/Write Data
library(readr)
library(readxl)

# . Clean Wilcoxon Rank-Sum Test
library(broom)

# . Data Wrangling/Cleaning
library(dplyr)
library(tidyr)
library(purrr)
library(janitor)

# . String Manipulation
library(stringr)

# Constants ====================================================================
# . HALO Raw Data
HALO <- here("data", "mIHC", "raw", "total_summary_results.xlsx")

# . Multiplex Data
MULTIPLEX <- here("data", "mIHC", "immune_profile_multiplex.xlsx")

# . Proteomic Metadata
PROTEOMIC_METADATA <- here("data", "IHC", "22-002.xlsx")

# . Cell Types
CELL_TYPES <- c("cd3p_cd8p"    = "CD8+ T cell", 
                "cd3p_cd8n"    = "CD8- T cell",
                "cd68p"        = "Macrophage",
                "cd20p_cd79ap" = "B cell",
                "cd20n_cd79ap" = "Plasma cell")

# . Histotypes
HISTOTYPES <- c("HGSC", "LGSC", "SBT")

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
# . Loading HALO Metadata
total_cell_counts <- read_xlsx(HALO) |>
  clean_names() |>
  select(row = tma_row, column = tma_column, contains("total_cells")) |>
  mutate(column = map_dbl(column, \(x) which(LETTERS == x))) |>
  rename_with(.cols = contains("cells"), .fn = \(x) str_remove(x, "_cells"))
 
# . Loading Proteomic Metadata
histotypes <- read_xlsx(PROTEOMIC_METADATA, sheet = "Core IDs") |>
  clean_names() |>
  mutate(row = case_when(sector == 2 ~ row + 1, .default = row))

# Loading Data =================================================================
# . Loading Multiplex Data
immune_cell_counts <- read_xlsx(MULTIPLEX) |>
  clean_names(replace = c("\\+" = "p", "-" = "n")) |>
  select(core_id, row = tma_row, column = tma_column, matches("cells$")) |>
  mutate(across(.cols = matches("cells$"), .fns = parse_number))

# Preprocessing Data ===========================================================
normalized_cell_counts <- immune_cell_counts |>
  mutate(column = map_dbl(column, \(x) which(LETTERS == x))) |>
  left_join(total_cell_counts, by = join_by(row == row, column == column)) |>
  pivot_longer(cols = matches("cells$"), 
               names_to = c("compartment", "marker"), 
               names_pattern = "(stroma|tumour)_(.*)_cells", 
               values_to = "count") |>
  mutate(cell_type = recode(marker, !!!CELL_TYPES)) |>
  pivot_wider(names_from = compartment, values_from = count) |>
  mutate(normalized_stroma = stroma / stroma_total,
         normalized_tumour = tumour / tumour_total) |>
  left_join(histotypes, by = join_by(core_id, row, column)) |>
  filter(histotype %in% HISTOTYPES) |>
  select(histotype, cell_type, normalized_stroma, normalized_tumour)

# Wilcoxon Rank Sum Test =======================================================
results <- normalized_cell_counts |>
  mutate(across(.cols = contains("normalized"), .fns = \(x) log10(x + 1))) |>
  nest(data = -cell_type) |> 
  pairwise_wilcox_test(values = "normalized_stroma", group = "histotype") |>
  pairwise_wilcox_test(values = "normalized_tumour", group = "histotype")

# Save Results =================================================================
# . Nuclei Normalized Cell Counts
write_csv(normalized_cell_counts, file = here("data", "outputs", "00-revisions", "22-002_nuclear_normalized_cell_counts.csv"))

# . Stroma Results - Pairwise Wilcoxon Rank Sum Test
results |>
  unnest(test_normalized_stroma) |>
  select(cell_type, group1, group2, p_value, adj_p_value) |>
  write_csv(file = here("data", "outputs", "00-revisions", "22-002_nuclear_normalized_stroma_results.csv"))

# . Tumour Results - Pairwise Wilcoxon Rank Sum Test
results |>
  unnest(test_normalized_tumour) |>
  select(cell_type, group1, group2, p_value, adj_p_value) |>
  write_csv(file = here("data", "outputs", "00-revisions", "22-002_nuclear_normalized_tumour_results.csv"))
