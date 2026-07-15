# =============================================================================
# CLUSTER CENTROIDS: Clean view of cluster-defining patterns
# =============================================================================
# Shows ONLY the mean behavior per cluster (not all site-solutes)
# - One line per cluster
# - Ribbons show ±1 SD
# - Faceted by solute type (geogenic vs transitional)
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
plot_dir   <- file.path(base_dir, "exploratory_plots", "02_exploration", "2c_cluster_centroids")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
message("Loading cluster data...")
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
# COMPUTE CLUSTER CENTROIDS
# =============================================================================

# Calculate mean ± SD for each cluster-solute-month combination
# This is the centroid representing the "typical" behavior in that cluster
cluster_centroids <- cluster_wy %>%
  select(chemical, Cluster_climRef, `1`:`12`) %>%
  pivot_longer(`1`:`12`, names_to = "month", values_to = "z_norm") %>%
  mutate(month = as.numeric(month)) %>%
  group_by(chemical, Cluster_climRef, month) %>%
  summarise(
    mean_z = mean(z_norm, na.rm = TRUE),
    sd_z = sd(z_norm, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  left_join(solute_class %>% select(chemical, behavior_type), by = "chemical") %>%
  mutate(
    Cluster = factor(Cluster_climRef, levels = c("1", "2", "3", "4")),
    chemical = factor(chemical, levels = solute_order),
    behavior_type = factor(behavior_type, levels = c("Geogenic", "Transitional", "Biogenic"))
  )

# =============================================================================
# FIGURE 1: All solutes, one panel per cluster
# =============================================================================

p_by_cluster <- ggplot(cluster_centroids,
                       aes(x = month, y = mean_z,
                           color = chemical, fill = chemical)) +
  geom_ribbon(aes(ymin = mean_z - sd_z, ymax = mean_z + sd_z),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.0, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  facet_wrap(~ Cluster, ncol = 2,
             labeller = labeller(Cluster = function(x) paste("Cluster", x))) +
  scale_x_continuous(breaks = 1:12, labels = month_labels) +
  scale_color_solute(name = "Solute") +
  scale_fill_solute(name = "Solute") +
  labs(
    x = "Month",
    y = "Normalized Concentration",
    title = "Cluster Centroids: Mean Seasonal Patterns",
    subtitle = "Lines = cluster mean behavior (averaged across all site-years in that cluster) ± SD"
  ) +
  theme_hja(base_size = 11) +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "right",
    legend.key.height = unit(0.4, "cm")
  )

ggsave(
  file.path(plot_dir, "cluster_centroids_by_cluster.png"),
  p_by_cluster, width = 14, height = 10, dpi = 300, bg = "white"
)

message("✓ Saved: cluster_centroids_by_cluster.png")

# =============================================================================
# FIGURE 2: Faceted by solute, comparing clusters
# =============================================================================

p_by_solute <- ggplot(cluster_centroids,
                      aes(x = month, y = mean_z,
                          color = Cluster, fill = Cluster)) +
  geom_ribbon(aes(ymin = mean_z - sd_z, ymax = mean_z + sd_z),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.1, alpha = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  facet_wrap(~ chemical, scales = "free_y", ncol = 3) +
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = month_labels[c(1, 4, 7, 10)]) +
  scale_color_cluster() +
  scale_fill_cluster() +
  labs(
    x = "Month",
    y = "Normalized Concentration",
    title = "Cluster Centroids by Solute",
    subtitle = "Lines = cluster mean behavior (averaged across all site-years in that cluster) ± SD"
  ) +
  theme_hja(base_size = 11) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  file.path(plot_dir, "cluster_centroids_by_solute.png"),
  p_by_solute, width = 14, height = 12, dpi = 300, bg = "white"
)

message("✓ Saved: cluster_centroids_by_solute.png")

# =============================================================================
# FIGURE 3: Faceted by behavior type (geogenic vs transitional)
# =============================================================================

p_by_type <- ggplot(cluster_centroids,
                    aes(x = month, y = mean_z,
                        color = Cluster, fill = Cluster)) +
  geom_ribbon(aes(ymin = mean_z - sd_z, ymax = mean_z + sd_z),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.0, alpha = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  facet_grid(behavior_type ~ chemical, scales = "free_y") +
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = month_labels[c(1, 4, 7, 10)]) +
  scale_color_cluster() +
  scale_fill_cluster() +
  labs(
    x = "Month",
    y = "Normalized Concentration",
    title = "Cluster Centroids by Solute Behavior Type",
    subtitle = "Rows = Behavior type (geogenic vs transitional) | Lines = cluster mean ± SD"
  ) +
  theme_hja(base_size = 9) +
  theme(
    panel.spacing = unit(0.3, "lines"),
    strip.text.y = element_text(size = 9, face = "bold", angle = 0, hjust = 0),
    strip.text.x = element_text(size = 8, face = "bold", angle = 90, hjust = 0),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    legend.position = "bottom"
  )

ggsave(
  file.path(plot_dir, "cluster_centroids_by_behavior_type.png"),
  p_by_type, width = 18, height = 10, dpi = 300, bg = "white"
)

message("✓ Saved: cluster_centroids_by_behavior_type.png")

# =============================================================================
# FIGURE 4: Geogenic solutes only (simpler view)
# =============================================================================

p_geogenic <- cluster_centroids %>%
  filter(behavior_type == "Geogenic") %>%
  ggplot(aes(x = month, y = mean_z,
             color = Cluster, fill = Cluster)) +
  geom_ribbon(aes(ymin = mean_z - sd_z, ymax = mean_z + sd_z),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2, alpha = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  facet_wrap(~ chemical, scales = "free_y", ncol = 2) +
  scale_x_continuous(breaks = 1:12, labels = month_labels) +
  scale_color_cluster() +
  scale_fill_cluster() +
  labs(
    x = "Month",
    y = "Normalized Concentration",
    title = "Cluster Centroids: Geogenic Solutes",
    subtitle = "Geogenic solutes (Ca, Na, K, DSi) are strongly associated with Cluster 1"
  ) +
  theme_hja(base_size = 12) +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  file.path(plot_dir, "cluster_centroids_geogenic_only.png"),
  p_geogenic, width = 10, height = 8, dpi = 300, bg = "white"
)

message("✓ Saved: cluster_centroids_geogenic_only.png")

message("\n=== COMPLETE ===")
message("Created 4 cluster centroid figures")
