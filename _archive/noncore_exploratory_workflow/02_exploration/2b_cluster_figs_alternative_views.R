# =============================================================================
# ALTERNATIVE CLUSTER VIEWS
# =============================================================================
# Figure 1: Facet by Site × Solute, color by Cluster
# Figure 2: Facet by Site × Cluster, color by Solute
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
plot_dir   <- file.path(base_dir, "exploratory_plots", "02_exploration", "2b_clusters", "annual")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
message("Loading data...")
cluster_wy <- readr::read_csv(
  file.path(output_dir, "ClusterStreams_allSolutes_byWaterYear.csv"),
  show_col_types = FALSE
)

# Prepare data
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
  mutate(
    Cluster = factor(Cluster_climRef, levels = c("1", "2", "3", "4")),
    chemical = factor(chemical, levels = solute_order),
    Stream_Name = factor(Stream_Name, levels = site_order)
  )

# =============================================================================
# FIGURE 1: Site × Solute facets, colored by Cluster
# =============================================================================

p_site_solute <- ggplot(cluster_avg,
                        aes(x = month, y = mean_z,
                            group = Cluster,
                            color = Cluster)) +
  geom_line(alpha = 0.5, linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70", linewidth = 0.3) +
  facet_grid(Stream_Name ~ chemical,
             scales = "free_y",
             labeller = labeller(chemical = label_value)) +
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = month_labels[c(1, 4, 7, 10)]) +
  scale_color_cluster() +
  labs(
    x = "Month",
    y = "Normalized Concentration",
    title = "Cluster Behavior by Site and Solute",
    subtitle = "Rows = Site | Columns = Solute | Lines = Cluster averages (mean monthly pattern when site-solute appears in that cluster)",
    color = "Cluster"
  ) +
  theme_hja(base_size = 9) +
  theme(
    panel.spacing = unit(0.2, "lines"),
    strip.text.y = element_text(size = 8, face = "bold", angle = 0, hjust = 0),
    strip.text.x = element_text(size = 8, face = "bold", angle = 90, hjust = 0),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    legend.position = "bottom"
  )

ggsave(
  file.path(plot_dir, "annual_facet_site_x_solute.png"),
  p_site_solute, width = 18, height = 14, dpi = 300, bg = "white"
)

message("✓ Saved: annual_facet_site_x_solute.png")

# =============================================================================
# FIGURE 2: Site × Cluster facets, colored by Solute
# =============================================================================

p_site_cluster <- ggplot(cluster_avg,
                         aes(x = month, y = mean_z,
                             group = chemical,
                             color = chemical)) +
  geom_line(alpha = 0.5, linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70", linewidth = 0.3) +
  facet_grid(Stream_Name ~ Cluster,
             scales = "free_y",
             labeller = labeller(Cluster = function(x) paste("Cluster", x))) +
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = month_labels[c(1, 4, 7, 10)]) +
  scale_color_solute(name = "Solute") +
  labs(
    x = "Month",
    y = "Normalized Concentration",
    title = "Solute Behavior by Site and Cluster",
    subtitle = "Rows = Site | Columns = Cluster | Lines = Solute averages (mean monthly pattern when site-solute appears in that cluster)"
  ) +
  theme_hja(base_size = 9) +
  theme(
    panel.spacing = unit(0.2, "lines"),
    strip.text.y = element_text(size = 8, face = "bold", angle = 0, hjust = 0),
    strip.text.x = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    legend.position = "bottom",
    legend.key.height = unit(0.3, "cm")
  )

ggsave(
  file.path(plot_dir, "annual_facet_site_x_cluster.png"),
  p_site_cluster, width = 10, height = 14, dpi = 300, bg = "white"
)

message("✓ Saved: annual_facet_site_x_cluster.png")

# =============================================================================
# BONUS FIGURE 3: Cluster × Solute facets, colored by Site
# =============================================================================

p_cluster_solute <- ggplot(cluster_avg,
                            aes(x = month, y = mean_z,
                                group = Stream_Name,
                                color = Stream_Name)) +
  geom_line(alpha = 0.6, linewidth = 0.8) +
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
  theme_hja(base_size = 9) +
  theme(
    panel.spacing = unit(0.2, "lines"),
    strip.text.y = element_text(size = 8, face = "bold", angle = 0, hjust = 0),
    strip.text.x = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    legend.position = "bottom",
    legend.key.height = unit(0.3, "cm")
  )

ggsave(
  file.path(plot_dir, "annual_facet_cluster_x_solute.png"),
  p_cluster_solute, width = 12, height = 16, dpi = 300, bg = "white"
)

message("✓ Saved: annual_facet_cluster_x_solute.png")

message("\n=== COMPLETE ===")
message("Created 3 alternative views:")
message("  1. Site × Solute (color by Cluster)")
message("  2. Site × Cluster (color by Solute)")
message("  3. Cluster × Solute (color by Site)")
