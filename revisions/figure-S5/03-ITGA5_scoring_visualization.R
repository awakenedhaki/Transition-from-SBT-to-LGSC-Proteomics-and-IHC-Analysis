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
PROTEOMIC <- here("data", "outputs", "00-revisions", "proteomic_itga5_h_scores.csv")
STANFORD <- here("data", "outputs", "00-revisions", "stanford_itga5_h_scores.csv")

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
      labs(x = "CD49e (ITGA5) Staining Score", y = "Percent", fill = "Histotype") +
      scale_fill_manual(values = palette) +
      scale_y_continuous(label = label_percent()) +
      expand_limits(y = c(0, 1.05))
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

# Visualization ================================================================
# . Proteomic Cohort
summarized_proteomic_proportions |>
  drop_na() |>
  bar_plot(palette = PROTEOMIC_HISTOTYPES)

# . Stanford Cohort
summarized_stanford_proportions |>
  bar_plot(palette = STANFORD_HISTOTYPES)

# Save Plot ==================================================================== 
ggsave(here("figures", "revisions", "figure-S5", "figure-S5_B.svg"), 
       plot = last_plot(), width = 5, height = 3, dpi = 1200)
