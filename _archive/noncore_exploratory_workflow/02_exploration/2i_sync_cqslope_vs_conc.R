suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

base_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality"
out_dir <- file.path(base_dir, "outputs")
fig_dir <- file.path(base_dir, "exploratory_plots", "02_exploration", "2i_sync_cqslope_vs_conc")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
try(source(file.path(repo_dir, "00_helpers", "plot_prefs.R")), silent = TRUE)
try(source(file.path(repo_dir, "00_helpers", "plot_theme_set.R")), silent = TRUE)

# Load composite synchrony with per-solute site-pair metrics if available
sync_df <- read_csv(file.path(out_dir, "HJA_composite_synchrony.csv"), show_col_types = FALSE)

required_cols <- c("Stream1", "Stream2", "solute", "conc_sync_pair", "cqslope_sync_pair")
has_required <- all(required_cols %in% names(sync_df))

if (!has_required) {
  # Fallback: compute pairwise from Abbott and Wymore inputs if present
  abbott <- read_csv(file.path(out_dir, "HJA_Abbott_synchrony_windows.csv"), show_col_types = FALSE)
  wymore <- read_csv(file.path(out_dir, "HJA_wymore_crosssite_sync.csv"), show_col_types = FALSE)

  # Expect columns: Stream1, Stream2, solute, pearson_r (Abbott); similarity/sync (Wymore)
  has_abbott <- all(c("Stream1", "Stream2", "solute", "pearson_r") %in% names(abbott))
  has_wymore_similarity <- all(c("Stream1", "Stream2", "solute", "similarity") %in% names(wymore))
  has_wymore_sync <- all(c("Stream1", "Stream2", "solute", "sync") %in% names(wymore))
  if (has_abbott && (has_wymore_similarity || has_wymore_sync)) {
    sync_df <- abbott %>%
      select(Stream1, Stream2, solute, conc_sync_pair = pearson_r) %>%
      inner_join(
        (if (has_wymore_similarity) {
           wymore %>% select(Stream1, Stream2, solute, cqslope_sync_pair = similarity)
         } else {
           wymore %>% select(Stream1, Stream2, solute, cqslope_sync_pair = sync)
         }),
        by = c("Stream1", "Stream2", "solute")
      )
  } else {
    stop("Required synchrony inputs not found: need per-solute site-pair metrics")
  }
}

# Clean and ensure ordering symmetry (optional: enforce Stream1 < Stream2)
sync_df <- sync_df %>%
  filter(!is.na(conc_sync_pair), !is.na(cqslope_sync_pair)) %>%
  mutate(solute = factor(solute, levels = if (exists("ALL_SOLUTES")) ALL_SOLUTES else solute)) %>%
  mutate(pair_id = paste(pmin(Stream1, Stream2), pmax(Stream1, Stream2), sep = "_") )

# Plot: cq slope sync vs conc sync per solute-site-pair
p <- ggplot(sync_df, aes(x = conc_sync_pair, y = cqslope_sync_pair)) +
  geom_point(aes(color = solute), alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "gray40", alpha = 0.15) +
  facet_wrap(~ solute, scales = "free") +
  labs(
    x = get_sync_label("conc_sync_pair"),
    y = get_sync_label("cqslope_sync_pair")
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "sync_cqslope_vs_conc_by_solute_pair.png"), p, width = 14, height = 10, dpi = 300)
