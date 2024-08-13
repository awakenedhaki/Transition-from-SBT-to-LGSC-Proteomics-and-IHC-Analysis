# Loading Dependencies =========================================================
# . Path Handling
library(here)

# . Read/Write Data
library(readr)

# . Data Wrangling
library(dplyr)

# . String Manipulation
library(glue)
library(stringr)

# . Visualization
library(ggplot2)

# Constants ====================================================================
# . Groups of Interest
GROUP_1 <- "LGSC"
GROUP_2 <- "HGSC"

# . Output Filename
FILENAME <- glue("{GROUP_2}_v_{GROUP_1}-fgsea.csv")

# . File Paths
GSEA <- here("data", "outputs", "04", FILENAME)

# . Threshold Values
ALPHA <- 0.05

# . Aesthetics
HISTOTYPE_COLORS <- c("HGSC" = "#184E77", "SBT" = "#52B69A", "LGSC" = "#D9ED92")

# ggplot2 Defaults =============================================================
theme_set(theme_bw())

# Loading Data =================================================================
gsea <- read_csv(GSEA) |>
  select(pathway, highlighted, padj, nes, size, leading_edge)

# Preprocessing Data ===========================================================
filtered_gsea <- gsea |>
  filter(highlighted) |>
  filter(padj < ALPHA) |>
  mutate(pathway = str_replace_all(str_remove(pathway, "GOBP_"), "_", " ")) |>
  mutate(histotype = case_when(
    nes > 0 ~ GROUP_1,
    nes < 0 ~ GROUP_2
  ))

# Visualization ================================================================
filtered_gsea |>
  mutate(pathway = reorder_within(pathway, by = nes, within = histotype)) |>
  ggplot(aes(x = nes, y = pathway, fill = histotype)) +
    geom_col() +
    labs(x = "Normalized Enrichment Score", 
         y = "Biological Process", 
         fill = "Histotype") +
    scale_y_reordered() +
    scale_fill_manual(values = HISTOTYPE_COLORS)

# Save Plot ====================================================================
ggsave(here("figures", "main", "figure-2", "figure-2_C.svg"), 
       plot = last_plot(), width = 7, height = 5, dpi = 300, scale = 1.5)
       