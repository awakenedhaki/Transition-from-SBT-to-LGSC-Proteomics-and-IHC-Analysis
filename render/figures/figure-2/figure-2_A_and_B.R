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
GROUP_1 <- "LGSC"
GROUP_2 <- "SBT"

# . Output Filename
FILENAME <- glue("{GROUP_2}_v_{GROUP_1}.csv")

# . File Paths
DIFFERENTIAL_EXPRESSION <- here("data", "outputs", "03", FILENAME)

# . Volcano Plot Aesthetics
DASHED_LINE <- 2

HISTOTYPE_COLORS <- c("HGSC" = "#184E77", "SBT" = "#52B69A", "LGSC" = "#D9ED92")

# . Thresholds
ALPHA <- 0.05
TOP_N_PROTEINS <- 10
FC_THRESHOLD <- 1.5
DELTA_FOLD_CHANGE_ZERO <- 0.01

# ggplot2 Defaults =============================================================
update_geom_defaults("point", list(size = 3))
update_geom_defaults("hline", list(linetype = DASHED_LINE))
update_geom_defaults("vline", list(linetype = DASHED_LINE))

theme_set(theme_bw())

# Loading Data =================================================================
differential_expression <- read_csv(DIFFERENTIAL_EXPRESSION)

# Subset Data ==================================================================
non_significant_proteins <- differential_expression |>
  filter(abs(log_fc) <= log2(FC_THRESHOLD)) |>
  select(symbol, log_fc, nlog10_adj_pval)
  
significant_proteins <- differential_expression |>
  filter(deqms_adj_pval < ALPHA, abs(log_fc) > log2(FC_THRESHOLD)) |>
  select(symbol, log_fc, nlog10_adj_pval) |>
  # Label proteins by group
  mutate(sign = ifelse(log_fc > 0, GROUP_1, GROUP_2)) |>
  # Group by sign
  nest(data = -sign) |>
  # Select top N proteins for each group
  mutate(top_n_proteins = map(data, \(x) {
    x |>
      slice_max(order_by = abs(log_fc * nlog10_adj_pval), n = TOP_N_PROTEINS)
  }))

# Volcano Plot =================================================================
ggplot(mapping = aes(x = log_fc, y = nlog10_adj_pval)) +
  # Mark non-significant proteins
  geom_point(color = "lightgrey", 
             data = non_significant_proteins) +
  # Mark significant proteins
  geom_point(aes(color = sign), 
             data = (significant_proteins |>
                       unnest(cols = data)),
             show.legend = FALSE) +
  # Annotate top N proteins for each group
  geom_text_repel(aes(label = symbol), 
                  min.segment.length = 0,
                  seed = 123,
                  data = (significant_proteins |>
                          unnest(cols = top_n_proteins))) +
  # Label thresholds
  geom_vline(xintercept = c(-log2(FC_THRESHOLD), log2(FC_THRESHOLD))) +
  geom_hline(yintercept = -log10(ALPHA)) +
  # Aesthetics
  scale_color_manual(values = HISTOTYPE_COLORS) +
  labs(x = expression(log[2]("Fold-Change")),
       y = expression("-log"[10]("adj. p-value")),
       color = "Histotype")

# Save Plot ====================================================================
ggsave(here("figures", "main", "figure-2", "figure-2_B.svg"), plot = last_plot(), 
       width = 5, height = 4, dpi = 1200, scale = 1.5)
