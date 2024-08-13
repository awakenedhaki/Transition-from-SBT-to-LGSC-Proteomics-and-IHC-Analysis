# Loading Dependencies =========================================================
# . Path Handling
library(here)

# . Read/Write Data
library(readr)

# . Data Wrangling
library(glue)
library(purrr)
library(tidyr)

# . Visualization
library(ggplot2)
library(ggrepel)

# Constants ====================================================================
# . Groups of Interest
GROUP_1 <- "SBT"
GROUP_2 <- "mSBT"

# . Output Filename
FILENAME <- glue("{GROUP_2}_v_{GROUP_1}.csv")

# . File Paths
DIFFERENTIAL_EXPRESSION <- here("data", "outputs", "03", FILENAME)

# . Volcano Plot Aesthetics
DASHED_LINE <- 2

# . Thresholds
ALPHA <- 0.05
TOP_N_PROTEINS <- 10
FC_THRESHOLD <- 1.5
DELTA_FOLD_CHANGE_ZERO <- 0.01

# . Histotype Colors
HISTOTYPE_COLORS <- c("HGSC" = "#184E77", SBT = "#52B69A", "LGSC" = "#D9ED92")

# ggplot2 Defaults =============================================================
update_geom_defaults("point", list(size = 3))
update_geom_defaults("hline", list(linetype = DASHED_LINE))
update_geom_defaults("vline", list(linetype = DASHED_LINE))

theme_set(theme_bw())

# Loading Data =================================================================
differential_expression <- read_csv(DIFFERENTIAL_EXPRESSION)

# Volcano Plot =================================================================
ggplot(mapping = aes(x = log_fc, y = nlog10_adj_pval)) +
  # Mark non-significant proteins
  geom_point(color = "lightgrey", 
             data = differential_expression) +
  # Label thresholds
  geom_vline(xintercept = c(-log2(FC_THRESHOLD), log2(FC_THRESHOLD))) +
  geom_hline(yintercept = -log10(ALPHA)) +
  # Aesthetics
  scale_color_manual(values = HISTOTYPE_COLORS) +
  labs(x = expression(log[2]("Fold-Change")),
       y = "-log10(adj. p-value)",
       color = "Histotype")

# Save Plot ====================================================================
ggsave(here("figures", "supplemental", "figure-S3", "figure-S3_A.svg"), plot = last_plot(), 
       width = 5, height = 4, dpi = 1200, scale = 1.5)
