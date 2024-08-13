# Loading Dependencies =========================================================
# . Path Handling
library(here)

# . Read/Write Data
library(readr)

# . Data Wrangling
library(dplyr)
library(forcats)

# . Visualization
library(ggplot2)

# Constants ====================================================================
# . File Paths
STANFORD_MACROPHAGE_DENSITIES <- here("data", "outputs", "07", "stanford_macrophage_densities.csv")

# . Histotype Colors
HISTOTYPE_COLORS <- c("SBT" = "#52B69A", "LGSC" = "#D9ED92")

# Loading Data =================================================================
macrophage_densities <- read_csv(STANFORD_MACROPHAGE_DENSITIES)

# Visualization ================================================================
macrophage_densities |>
  mutate(histotype = case_when(histotype == "Serous borderline tumor" ~ "SBT",
                               histotype == "Low grade serous carcinoma" ~ "LGSC",
                               .default = histotype), 
         histotype = fct_relevel(histotype, !!!names(HISTOTYPE_COLORS)),
         marker = fct_relevel(marker, c("CD68", "CD163"))) |>
  ggplot(aes(x = histotype, y = stroma, fill = histotype)) +
    geom_boxplot(show.legend = FALSE) +
    facet_wrap(~marker) +
    scale_fill_manual(values = HISTOTYPE_COLORS) +
    expand_limits(y = c(1, 3.75))
    
# Save Plot ====================================================================
ggsave(here("figures", "main", "figure-4", "figure-4_E.svg"), 
       plot = last_plot(), width = 8, height = 5, dpi = 1200)
