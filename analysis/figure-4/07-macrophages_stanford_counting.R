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

# . . CD68
TMA1_CD68 <- "TMA1 CD68 (Almira)"
TMA2_CD68 <- "TMA2 CD68 (Almira)"
LGSC_CD68 <- "LGSC & Mixed CD68 (Almira)"

CD68 <- c("TMA 1" = TMA1_CD68, 
          "TMA 2" = TMA2_CD68, 
          "LGS/Mixed" = LGSC_CD68)

# . . CD163
TMA1_CD163 <- "TMA1 CD163 (Almira)"
TMA2_CD163 <- "TMA2 CD163 (Rodrigo)"
LGSC_CD163 <- "LGSC & Mixed CD163 (Almira)"

CD163 <- c("TMA 1" = TMA1_CD163, 
           "TMA 2" = TMA2_CD163, 
           "LGS/Mixed" = LGSC_CD163)

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

clean_macrophage_counts <- function(tbl, tma, marker) {
  tbl |>
    rename(row = row_number, column = column_number) |>
    extract(stroma_tumor, 
            into = c("stroma", "tumour"), 
            regex = "(\\d+)/(\\d+)", 
            convert = TRUE) |>
    mutate_at(.vars = c("row", "column"), .funs = parse_number) |>
    mutate(marker = {{ marker }}, tma_number = {{ tma }})
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

# Loading Metadata =============================================================
metadata <- read_xlsx(STANFORD, sheet = METADATA, range = "A1:I268") |>
  clean_names() |>
  select(tma_number, row = row_number, column = column_number, 
         case_code, histotype = diagnosis, exclude) |>
  mutate_at(.vars = c("row", "column"), .funs = parse_number)

# Loading Data =================================================================
# . Tumour Areas
areas <- tibble()
for (tma in names(AREA)) {
  area <- read_tsv(AREA[tma]) |>
    clean_names(replace = c("µ" = "u")) |>
    clean_tma_areas(tma)
  
  areas <- bind_rows(areas, area)
}

# . CD68 Counts
cd68 <- tibble()
for (tma in names(CD68)) {
  counts <- read_xlsx(STANFORD, sheet = CD68[tma], range = "O1:Q90") |>
    clean_names() |>
    clean_macrophage_counts(tma = tma, marker = "CD68")
  
  cd68 <- bind_rows(cd68, counts)
}

# . CD163 Counts
cd163 <- tibble()
for (tma in names(CD163)) {
  counts <- read_xlsx(STANFORD, sheet = CD163[tma], range = "O1:Q90") |>
    clean_names() |>
    clean_macrophage_counts(tma = tma, marker = "CD163")
  
  cd163 <- bind_rows(cd163, counts)
}

# . Merge Data
counts <- bind_rows(cd68, cd163) |>
  left_join(areas, by = join_by(tma_number, row, column), relationship = "many-to-many") |>
  filter(!exclude) |>
  select(-c(row, column))

# Preprocessing Data ===========================================================
macrophages <- counts |>
  # Normalizing cell counts by tissue compartment area
  mutate(stroma_density = stroma / stroma_area,
         tumour_density = tumour / tumour_area)|>
  # Pivoting data to long format
  pivot_longer(cols = c(stroma_density, tumour_density), 
               names_to = "compartment",
               names_pattern = "(\\w+)_density",
               values_to = "density") |>
  # Calculating the mean density of each marker by compartment
  reframe(mean_density = mean(density, na.rm = TRUE), 
          .by = c(tma_number, case_code, histotype, marker, compartment)) |>
  # Pseudo log10 transformation
  mutate(pseudo_log10_mean_density = log10(mean_density + 1)) |>
  # Remove rows with missing densities
  filter(!is.na(mean_density)) |>
  select(-mean_density) |>
  # Pivoting data to wide format
  pivot_wider(names_from = compartment, values_from = pseudo_log10_mean_density) |>
  # Correct TMA names
  mutate(tma_number = case_match(tma_number,
    "LGS/Mixed" ~ "LGSC/Mixed",
    "TMA 1" ~ "TMA1",
    "TMA 2" ~ "TMA2",
    .default = tma_number
  ))

# Wilcoxon Rank Sum Test =======================================================
results <- macrophages |>
  select(histotype, marker, stroma, tumour) |>
  nest(data = c(histotype, stroma, tumour)) |>
  pairwise_wilcox_test(values = "stroma", group = "histotype") |>
  pairwise_wilcox_test(values = "tumour", group = "histotype")

# Save Results =================================================================
# . Macrophage Marker Densities
write_csv(macrophages, file = here("data", "outputs", "07", "stanford_macrophage_densities.csv"))

# . Stroma Results - Pairwise Wilcoxon Rank Sum Test
results |>
  unnest(test_stroma) |>
  select(marker, group1, group2, p_value, adj_p_value) |>
  write_csv(file = here("data", "outputs", "07", "stanford_stroma_results.csv"))

# . Tumour Results - Pairwise Wilcoxon Rank Sum Test
results |>
  unnest(test_tumour) |>
  select(marker, group1, group2, p_value, adj_p_value) |>
  write_csv(file = here("data", "outputs", "07", "stanford_tumour_results.csv"))
