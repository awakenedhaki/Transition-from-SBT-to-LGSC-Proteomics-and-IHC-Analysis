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
theme_set(theme_bw())

# Constants ====================================================================
WHOLE_SECTIONS <- here("data", "IHC", "whole_tissue_sections.xlsx")

HISTOTYPES <- c("SBT"  = "#52B69A",
                "Invasive" = "#D9ED92",
                "LGSC" = "#D9ED92")

# Helper Functions =============================================================
summarize_fap_score <- function(tbl, .by, .factor_levels) {
  tbl |>
    # Summarize FAP scores
    reframe(mean_score = mean(score, na.rm = TRUE), 
            .by = {{ .by }}) |>
    # Categorize FAP scores
    mutate(h_score = case_when(
      mean_score < 10 ~ "Negative",
      mean_score >= 10 ~ "Positive"
    )) |>
    # Filter for histotypes of interest
    filter(histotype %in% {{ .factor_levels }}) |>
    # Level histotypes
    mutate(histotype = fct_relevel(histotype, !!!{{ .factor_levels }}))
}

fisher_exact_test <- function(tbl) {
  tbl |>
    count(histotype, h_score) |>
    select(histotype, h_score, n) |>
    # Pivot table
    pivot_wider(names_from = histotype, values_from = n) |>
    # Fisher Exact Test
    column_to_rownames("h_score") |> 
    as.matrix() |>
    fisher.test() |>
    tidy()
}

# Loading Data =================================================================
whole_sections <- read_excel(WHOLE_SECTIONS) |>
  clean_names() |>
  filter(is.na(tma)) |>
  select(-c(tma, notes)) |>
  pivot_longer(cols = contains("region"),
               names_pattern = "region_(\\d)",
               names_to = "region",
               values_to = "score") |>
  mutate(histotype = recode(histotype, Invasive = "LGSC"))

# Preprocessing Data ===========================================================
whole_sections_fap_h_scores <- whole_sections |>
  summarize_fap_score(.by = c(slide_id, histotype),
                      .factor_levels = c("SBT", "LGSC"))

# Fisher Exact Test ============================================================
whole_sections_fap_h_scores |>
  fisher_exact_test()

# Visualization ================================================================
whole_sections_fap_h_scores |>
  count(histotype, h_score) |>
  add_count(histotype, name = "total", wt = n) |>
  mutate(proportion = n / total) |>
  ggplot(aes(x = h_score, y = proportion, fill = histotype)) +
    geom_col(color = "black", position = position_dodge(0.9), show.legend = FALSE) +
    geom_text(aes(label = percent(proportion)), vjust = -0.5) +
    facet_wrap(~histotype) +
    scale_fill_manual(values = HISTOTYPES) +
    scale_y_continuous(labels = label_percent(scale = 100)) +
    expand_limits(y = c(0, 1)) +
    labs(x = "Histotype", y = "Percent of Cases")

# Save Plot ====================================================================
ggsave(here("figures", "revisions", "figure-S4", "figure-S4_A.svg"), 
       plot = last_plot(), width = 5, height = 3, dpi = 1200)
