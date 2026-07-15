# =============================================================================
# CLUSTER HEATMAP: Site × Solute Grid colored by Modal Cluster
# =============================================================================
# Shows at a glance which sites and solutes behave similarly
# - Rows = Sites
# - Columns = Solutes
# - Color = Most common cluster
# - Annotation = % of years in that cluster (stability)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

rm(list = ls())

# Source helpers
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
source(file.path(repo_dir, "00_helpers", "plot_prefs.R"))

# Paths
base_dir   <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality"
output_dir <- file.path(base_dir, "outputs")
plot_dir   <- file.path(base_dir, "exploratory_plots", "02_exploration", "2d_cluster_heatmap")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
message("Loading data...")
cluster_wy <- readr::read_csv(
  file.path(output_dir, "ClusterStreams_allSolutes_byWaterYear.csv"),
  show_col_types = FALSE
)

# Load solute classification
solute_class <- readr::read_csv(
  file.path(output_dir, "solute_behavior_classification.csv"),
  show_col_types = FALSE
)

# =============================================================================
# CALCULATE MODAL CLUSTER AND STABILITY
# =============================================================================

# For each site-solute combo:
# - Modal cluster = most common cluster across all years
# - Stability = % of years in modal cluster
modal_cluster <- cluster_wy %>%
  group_by(Stream_Name, chemical, Cluster_climRef) %>%
  summarise(n_years = n(), .groups = "drop") %>%
  group_by(Stream_Name, chemical) %>%
  mutate(
    total_years = sum(n_years),
    pct_years = 100 * n_years / total_years
  ) %>%
  slice_max(n_years, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Stream_Name, chemical,
         modal_cluster = Cluster_climRef,
         stability_pct = pct_years,
         n_years_modal = n_years,
         total_years) %>%
  left_join(solute_class %>% select(chemical, behavior_type), by = "chemical") %>%
  mutate(
    Stream_Name = factor(Stream_Name, levels = site_order),
    chemical = factor(chemical, levels = solute_order),
    modal_cluster = factor(modal_cluster, levels = c("1", "2", "3", "4")),
    behavior_type = factor(behavior_type, levels = c("Geogenic", "Transitional", "Biogenic"))
  )

# =============================================================================
# FIGURE 1: Basic heatmap with modal cluster colors
# =============================================================================

p_heatmap_basic <- ggplot(modal_cluster,
                          aes(x = chemical, y = Stream_Name, fill = modal_cluster)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.0f%%", stability_pct)),
            size = 3, color = "white", fontface = "bold") +
  scale_fill_cluster() +
  labs(
    x = "Solute",
    y = "Site",
    title = "Modal Cluster Assignment by Site and Solute",
    subtitle = "Color = Most common cluster | Numbers = % of years in that cluster (stability)",
    fill = "Modal\nCluster"
  ) +
  theme_hja(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    panel.border = element_blank(),
    axis.ticks = element_blank()
  )

ggsave(
  file.path(plot_dir, "modal_cluster_heatmap.png"),
  p_heatmap_basic, width = 10, height = 7, dpi = 300, bg = "white"
)

message("✓ Saved: modal_cluster_heatmap.png")

# =============================================================================
# FIGURE 2: Heatmap with solute type facets
# =============================================================================

p_heatmap_by_type <- ggplot(modal_cluster,
                            aes(x = chemical, y = Stream_Name, fill = modal_cluster)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.0f", stability_pct)),
            size = 2.5, color = "white", fontface = "bold") +
  facet_grid(. ~ behavior_type, scales = "free_x", space = "free_x") +
  scale_fill_cluster() +
  labs(
    x = "Solute",
    y = "Site",
    title = "Modal Cluster Assignment by Solute Behavior Type",
    subtitle = "Geogenic solutes show stronger association with Cluster 1 | Numbers = stability (%)",
    fill = "Modal\nCluster"
  ) +
  theme_hja(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "right",
    panel.border = element_blank(),
    axis.ticks = element_blank()
  )

ggsave(
  file.path(plot_dir, "modal_cluster_heatmap_by_type.png"),
  p_heatmap_by_type, width = 12, height = 7, dpi = 300, bg = "white"
)

message("✓ Saved: modal_cluster_heatmap_by_type.png")

# =============================================================================
# FIGURE 3: Stability heatmap (continuous scale showing % stability)
# =============================================================================

p_stability <- ggplot(modal_cluster,
                      aes(x = chemical, y = Stream_Name, fill = stability_pct)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = modal_cluster),
            size = 4, color = "white", fontface = "bold") +
  scale_fill_gradient(low = "#fee5d9", high = "#a50f15",
                      name = "Stability\n(% years)",
                      limits = c(0, 100)) +
  labs(
    x = "Solute",
    y = "Site",
    title = "Cluster Stability by Site and Solute",
    subtitle = "Color intensity = stability (% of years in modal cluster) | Numbers = modal cluster"
  ) +
  theme_hja(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    panel.border = element_blank(),
    axis.ticks = element_blank()
  )

ggsave(
  file.path(plot_dir, "cluster_stability_heatmap.png"),
  p_stability, width = 10, height = 7, dpi = 300, bg = "white"
)

message("✓ Saved: cluster_stability_heatmap.png")

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

message("\n=== SUMMARY STATISTICS ===\n")

# Overall stability
message("Overall cluster stability:")
summary_stats <- modal_cluster %>%
  summarise(
    mean_stability = mean(stability_pct),
    median_stability = median(stability_pct),
    min_stability = min(stability_pct),
    max_stability = max(stability_pct)
  )
print(summary_stats)

# Stability by solute
message("\nMean stability by solute:")
by_solute <- modal_cluster %>%
  group_by(chemical, behavior_type) %>%
  summarise(
    mean_stability = mean(stability_pct),
    modal_cluster_mode = names(which.max(table(modal_cluster))),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_stability))
print(as.data.frame(by_solute))

# Stability by site
message("\nMean stability by site:")
by_site <- modal_cluster %>%
  group_by(Stream_Name) %>%
  summarise(
    mean_stability = mean(stability_pct),
    modal_cluster_mode = names(which.max(table(modal_cluster))),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_stability))
print(as.data.frame(by_site))

# Most/least stable combinations
message("\nMost stable site-solute combinations:")
most_stable <- modal_cluster %>%
  arrange(desc(stability_pct)) %>%
  head(10) %>%
  select(Stream_Name, chemical, modal_cluster, stability_pct)
print(as.data.frame(most_stable))

message("\nLeast stable site-solute combinations:")
least_stable <- modal_cluster %>%
  arrange(stability_pct) %>%
  head(10) %>%
  select(Stream_Name, chemical, modal_cluster, stability_pct)
print(as.data.frame(least_stable))

message("\n=== COMPLETE ===")
message("Created 3 heatmap figures")
