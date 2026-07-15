# =============================================================================
# CLUSTER CHARACTERIZATION: NESTED FACETS
# =============================================================================
# Facet by solute (rows) and cluster (columns)
# One line per site within each panel
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

rm(list = ls())

# Source helpers
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
source(file.path(repo_dir, "00_helpers", "plot_prefs.R"))

# Paths
base_dir   <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality"
output_dir <- file.path(base_dir, "outputs")
plot_dir   <- file.path(base_dir, "exploratory_plots", "02_exploration", "2b_clusters", "annual")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
message("Loading data...")
cluster_wy <- readr::read_csv(
  file.path(output_dir, "ClusterStreams_allSolutes_byWaterYear.csv"),
  show_col_types = FALSE
)

# Setup - use existing definitions from plot_prefs.R
# month_labels, solute_order, site_order, cluster_colors, site_colors, cluster_labels
# are all already defined by sourcing plot_prefs.R

# =============================================================================
# PREPARE DATA
# =============================================================================

cluster_avg <- cluster_wy %>%
  select(Stream_Name, chemical, water_year, Cluster_climRef, `1`:`12`) %>%
  pivot_longer(`1`:`12`, names_to = "month", values_to = "z_norm") %>%
  mutate(month = as.numeric(month)) %>%
  group_by(Stream_Name, chemical, Cluster_climRef, month) %>%
  summarise(
    mean_z = mean(z_norm, na.rm = TRUE),
    se_z = sd(z_norm, na.rm = TRUE) / sqrt(n()),
    n_years = n(),
    .groups = "drop"
  ) %>%
  filter(n_years >= 1) %>%
  mutate(Cluster = Cluster_climRef) %>%
  apply_factor_orders() %>%
  add_cluster_labels(cluster_col = "Cluster")

# =============================================================================
# FIGURE 1: Nested facets (Solute × Cluster), colored by site
# =============================================================================

p_nested <- ggplot(cluster_avg,
                   aes(x = month, y = mean_z,
                       group = Stream_Name,
                       color = Stream_Name)) +
  geom_line(alpha = 0.5, linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70", linewidth = 0.3) +
  facet_grid(chemical ~ Cluster,
             scales = "free_y",
             labeller = labeller(Cluster = function(x) paste("Cluster", x))) +
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = month_labels[c(1, 4, 7, 10)]) +
  scale_color_site(name = "Site") +
  labs(
    x = "Month",
    y = "Normalized Concentration",
    title = "Site Behavior by Solute and Cluster",
    subtitle = "Rows = Solute | Columns = Cluster | Lines = Site averages (mean monthly pattern when site-solute appears in that cluster)"
  ) +
  theme_hja(base_size = 10) +
  theme(
    panel.spacing = unit(0.3, "lines"),
    strip.text.y = element_text(size = 9, face = "bold", angle = 0, hjust = 0),
    strip.text.x = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    legend.position = "right",
    legend.key.height = unit(0.4, "cm")
  )

ggsave(
  file.path(plot_dir, "annual_nested_facets_solute_x_cluster.png"),
  p_nested, width = 12, height = 16, dpi = 300, bg = "white"
)

message("✓ Saved: annual_nested_facets_solute_x_cluster.png")

# =============================================================================
# FIGURE 2: Same structure but with ribbons for uncertainty
# =============================================================================

p_nested_ribbon <- ggplot(cluster_avg,
                          aes(x = month, y = mean_z,
                              group = Stream_Name,
                              color = Stream_Name,
                              fill = Stream_Name)) +
  geom_ribbon(aes(ymin = mean_z - se_z, ymax = mean_z + se_z),
              alpha = 0.1, color = NA) +
  geom_line(alpha = 0.6, linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70", linewidth = 0.3) +
  facet_grid(chemical ~ Cluster,
             scales = "free_y",
             labeller = labeller(Cluster = function(x) paste("Cluster", x))) +
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = month_labels[c(1, 4, 7, 10)]) +
  scale_color_site(name = "Site") +
  scale_fill_site(name = "Site") +
  labs(
    x = "Month",
    y = "Normalized Concentration",
    title = "Site Behavior by Solute and Cluster (with SE ribbons)",
    subtitle = "Rows = Solute | Columns = Cluster | Lines = Site averages (mean monthly pattern when site-solute appears in that cluster) ± SE"
  ) +
  theme_hja(base_size = 10) +
  theme(
    panel.spacing = unit(0.3, "lines"),
    strip.text.y = element_text(size = 9, face = "bold", angle = 0, hjust = 0),
    strip.text.x = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    legend.position = "right",
    legend.key.height = unit(0.4, "cm")
  )

ggsave(
  file.path(plot_dir, "annual_nested_facets_with_se.png"),
  p_nested_ribbon, width = 12, height = 16, dpi = 300, bg = "white"
)

message("✓ Saved: annual_nested_facets_with_se.png")

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

message("\n=== SUMMARY ===")

# Count site appearances per cluster
site_cluster_counts <- cluster_avg %>%
  distinct(Stream_Name, chemical, Cluster) %>%
  count(Stream_Name, Cluster, name = "n_solutes") %>%
  pivot_wider(names_from = Cluster, values_from = n_solutes, values_fill = 0,
              names_prefix = "Cluster_")

message("\nSite appearances per cluster (count of solutes):")
print(site_cluster_counts)

# Count solute appearances per cluster
solute_cluster_counts <- cluster_avg %>%
  distinct(Stream_Name, chemical, Cluster) %>%
  count(chemical, Cluster, name = "n_sites") %>%
  pivot_wider(names_from = Cluster, values_from = n_sites, values_fill = 0,
              names_prefix = "Cluster_")

message("\nSolute appearances per cluster (count of sites):")
print(solute_cluster_counts)

message("\n=== COMPLETE ===")
message("All figures saved to:", plot_dir)
