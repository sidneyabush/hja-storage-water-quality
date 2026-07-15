# =============================================================================
# SMALL MULTIPLES BY SITE: See within-site coherence
# =============================================================================
# Shows all solutes for each site, colored by modal cluster
# Helps answer: Do all solutes at a site behave similarly?
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
plot_dir   <- file.path(base_dir, "exploratory_plots", "02_exploration", "2e_cluster_by_site")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
message("Loading data...")
cluster_wy <- readr::read_csv(
  file.path(output_dir, "ClusterStreams_allSolutes_byWaterYear.csv"),
  show_col_types = FALSE
)

# =============================================================================
# CALCULATE SITE-SOLUTE AVERAGES
# =============================================================================

# Average across all years for each site-solute combination
site_solute_avg <- cluster_wy %>%
  select(Stream_Name, chemical, water_year, Cluster_climRef, `1`:`12`) %>%
  pivot_longer(`1`:`12`, names_to = "month", values_to = "z_norm") %>%
  mutate(month = as.numeric(month)) %>%
  group_by(Stream_Name, chemical, month) %>%
  summarise(
    mean_z = mean(z_norm, na.rm = TRUE),
    .groups = "drop"
  )

# Get modal cluster for coloring
modal_cluster <- cluster_wy %>%
  group_by(Stream_Name, chemical, Cluster_climRef) %>%
  summarise(n_years = n(), .groups = "drop") %>%
  group_by(Stream_Name, chemical) %>%
  slice_max(n_years, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Stream_Name, chemical, modal_cluster = Cluster_climRef)

# Combine
site_solute_data <- site_solute_avg %>%
  left_join(modal_cluster, by = c("Stream_Name", "chemical")) %>%
  mutate(
    Stream_Name = factor(Stream_Name, levels = site_order),
    chemical = factor(chemical, levels = solute_order),
    modal_cluster = factor(modal_cluster, levels = c("1", "2", "3", "4"))
  )

# =============================================================================
# FIGURE 1: One panel per site, all solutes
# =============================================================================

p_by_site <- ggplot(site_solute_data,
                    aes(x = month, y = mean_z,
                        group = chemical,
                        color = modal_cluster)) +
  geom_line(linewidth = 1.0, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  facet_wrap(~ Stream_Name, ncol = 3, scales = "free_y") +
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = month_labels[c(1, 4, 7, 10)]) +
  scale_color_cluster() +
  labs(
    x = "Month",
    y = "Mean Normalized Concentration",
    title = "Solute Behavior by Site",
    subtitle = "Lines = solutes colored by their modal cluster | Shows within-site coherence",
    color = "Modal\nCluster"
  ) +
  theme_hja(base_size = 10) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "right"
  )

ggsave(
  file.path(plot_dir, "solutes_by_site_panels.png"),
  p_by_site, width = 14, height = 12, dpi = 300, bg = "white"
)

message("✓ Saved: solutes_by_site_panels.png")

# =============================================================================
# FIGURE 2: Separate panels for watershed groups
# =============================================================================

# Define watershed groups
watershed_groups <- tibble(
  Stream_Name = c("GSMACK", "GSWS08", "GSWS07", "GSWS06",
                  "GSWS02", "GSLOOK", "GSWS01",
                  "GSWS10", "GSWS09"),
  group = c(rep("Group 1: Lower elevation", 4),
            rep("Group 2: Mid elevation", 3),
            rep("Group 3: Upper elevation", 2))
)

site_data_grouped <- site_solute_data %>%
  left_join(watershed_groups, by = "Stream_Name")

p_by_group <- ggplot(site_data_grouped,
                     aes(x = month, y = mean_z,
                         group = interaction(Stream_Name, chemical),
                         color = modal_cluster,
                         linetype = Stream_Name)) +
  geom_line(linewidth = 0.9, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  facet_wrap(~ group, ncol = 1, scales = "free_y") +
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = month_labels[c(1, 4, 7, 10)]) +
  scale_color_cluster() +
  labs(
    x = "Month",
    y = "Mean Normalized Concentration",
    title = "Solute Behavior by Watershed Group",
    subtitle = "Lines = site-solute combinations | Color = modal cluster | Line type = site",
    color = "Modal\nCluster",
    linetype = "Site"
  ) +
  theme_hja(base_size = 11) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "right"
  )

ggsave(
  file.path(plot_dir, "solutes_by_watershed_group.png"),
  p_by_group, width = 14, height = 10, dpi = 300, bg = "white"
)

message("✓ Saved: solutes_by_watershed_group.png")

# =============================================================================
# FIGURE 3: Grid layout (like original spaghetti)
# =============================================================================

p_grid <- ggplot(site_solute_data,
                 aes(x = month, y = mean_z,
                     group = chemical,
                     color = chemical)) +
  geom_line(linewidth = 1.1, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  facet_wrap(~ Stream_Name, ncol = 3, scales = "free_y") +
  scale_x_continuous(breaks = 1:12, labels = month_labels) +
  scale_color_solute(name = "Solute") +
  labs(
    x = "Month",
    y = "Mean Normalized Concentration",
    title = "Mean Solute Patterns by Site",
    subtitle = "Lines = solutes (averaged across all years)",
    color = "Solute"
  ) +
  theme_hja(base_size = 10) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "right",
    legend.key.height = unit(0.4, "cm")
  )

ggsave(
  file.path(plot_dir, "solutes_by_site_colored_by_solute.png"),
  p_grid, width = 14, height = 12, dpi = 300, bg = "white"
)

message("✓ Saved: solutes_by_site_colored_by_solute.png")

message("\n=== COMPLETE ===")
message("Created 3 site-based figures")
