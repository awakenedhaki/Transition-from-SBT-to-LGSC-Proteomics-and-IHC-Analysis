# Loading Dependencies =========================================================
# . Path Handling
library(here)

# . Read/Write Data
library(readr)

# . Data Wrangling
library(dplyr)
library(forcats)

# . Visualization
library(scales)
library(ggplot2)

# Constants ====================================================================
# . File Paths
PROTEOMIC <- here("data", "outputs", "05", "proteomic_fap_h_scores.csv")
STANFORD <- here("data", "outputs", "05", "stanford_fap_h_scores.csv")
WHOLE_SECTIONS <- here("data", "outputs", "05", "whole_sections_fap_h_scores.csv")

# . Aesthetics
PROTEOMIC_HISTOTYPES <- c("HGSC" = "#184E77", "SBT" = "#52B69A", "LGSC" = "#D9ED92")
STANFORD_HISTOTYPES <- c("SBT" = "#52B69A", "LGSC" = "#D9ED92")

# Helper Functions =============================================================
bar_plot <- function(tbl, palette) {
  tbl |>
    ggplot(aes(x = h_score, y = proportion, fill = histotype)) +
      geom_col(color = "black", position = position_dodge(0.9), show.legend = FALSE) +
      geom_text(aes(label = percent(proportion)), vjust = -0.5) +
      facet_wrap(~histotype) +
      labs(x = "FAP Staining Score", y = "Percent", fill = "Histotype") +
      scale_fill_manual(values = palette) +
      scale_y_continuous(label = label_percent()) +
      expand_limits(y = c(0, 1))
}

# ggplot2 Defaults =============================================================
theme_set(theme_bw())

# Loading Data =================================================================
# . Proteomic Cohort
proteomic <- read_csv(PROTEOMIC) |>
  mutate(histotype = fct_relevel(histotype, !!!names(PROTEOMIC_HISTOTYPES)))

# . Stanford Cohort
stanford <- read_csv(STANFORD) |>
  mutate(histotype = fct_relevel(histotype, !!!names(STANFORD_HISTOTYPES)))

# . Proteomic Whole Sections
whole_sections <- read_csv(WHOLE_SECTIONS) |>
  mutate(histotype = fct_relevel(histotype, !!!names(PROTEOMIC_HISTOTYPES)))

# Preprocessing Data ===========================================================
# . Proteomic Cohort
summarized_proteomic_proportions <- proteomic |>
  count(histotype, h_score) |>
  add_count(histotype, wt = n, name = "total") |>
  mutate(proportion = n / total)

# . Stanford Cohort
summarized_stanford_proportions <- stanford |>
  count(histotype, h_score) |>
  add_count(histotype, wt = n, name = "total") |>
  mutate(proportion = n / total)

# . Proteomic Whole Sections
summarized_whole_sections_proportions <- whole_sections |>
  count(histotype, h_score) |>
  add_count(histotype, wt = n, name = "total") |>
  mutate(proportion = n / total)

# Visualization ================================================================
# . Proteomic Cohort
summarized_proteomic_proportions |>
  bar_plot(palette = PROTEOMIC_HISTOTYPES)
    
# . Stanford Cohort
summarized_stanford_proportions |>
  bar_plot(palette = STANFORD_HISTOTYPES)

# . Proteomic Whole Sections
summarized_whole_sections_proportions |>
  bar_plot(palette = STANFORD_HISTOTYPES)

# Save Plot ==================================================================== 
ggsave(here("figures", "main", "figure-3", "figure-3_C.svg"), 
       plot = last_plot(), width = 5, height = 3, dpi = 1200)
