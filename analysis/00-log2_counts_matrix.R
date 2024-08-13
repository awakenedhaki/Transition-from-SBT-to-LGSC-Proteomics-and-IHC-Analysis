# Loading Dependencies =========================================================
# . Path Handling
library(here)

# . Data Cleaning
library(janitor)

# Constants ====================================================================
LOG2_COUNTS_MATRIX <- here("data", "proteome", "protein_log2_filtered_2pep_1unique.xlsx")

# Loading Data =================================================================
proteome <- read_excel(LOG2_COUNTS_MATRIX) %>%
  clean_names()

# Save Data ====================================================================
write_csv(proteome, here("data", "log2_counts_matrix.csv"))
