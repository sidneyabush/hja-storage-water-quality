# =============================================================================
# CLUSTER CHARACTERIZATION: BY SOLUTE
# =============================================================================
# Faceted by solute, showing site-cluster combinations
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

# Setup
month_labels <- c("J","F","M","A","M","J","J","A","S","O","N","D")
solute_order <- c("Ca", "Mg", "Na", "K", "DSi", "Cl", "SO4", "DOC", "NH3", "NO3", "PO4")

if (!exists("cluster_colors")) {
  cluster_colors <- c("1" = "#CFA980", "2" = "#98B89F", "3" = "#5E8AA1", "4" = "#526B8E")
}

# =============================================================================
# PREPARE DATA: Average monthly pattern per site-solute-cluster
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
  mutate(
    Cluster = factor(Cluster_climRef, levels = c("1", "2", "3", "4")),
    chemical = factor(chemical, levels = solute_order),
    site_cluster = interaction(Stream_Name, Cluster, sep = " - C")
  )

# =============================================================================
# FIGURE: Faceted by solute, colored by cluster
# =============================================================================

p_by_solute <- ggplot(cluster_avg,
                      aes(x = month, y = mean_z,
                          group = site_cluster,
                          color = Cluster)) +
  geom_line(alpha = 0.6, linewidth = 1.0) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  facet_wrap(~ chemical, scales = "free_y", ncol = 3) +
  scale_x_continuous(breaks = 1:12, labels = month_labels) +
  scale_color_cluster() +
  labs(
    x = "Month",
    y = "Mean Normalized Concentration",
    title = "Site-Cluster Behavior by Solute",
    subtitle = "Lines = Site-cluster combinations (mean monthly pattern when that site-solute appears in that cluster)"
  ) +
  theme_hja(base_size = 12) +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  file.path(plot_dir, "annual_cluster_by_solute_faceted.png"),
  p_by_solute, width = 14, height = 12, dpi = 300, bg = "white"
)

message("✓ Saved: annual_cluster_by_solute_faceted.png")

# =============================================================================
# FIGURE 2: Faceted by solute, with separate lines for each site (colored by site)
# =============================================================================

# Add site colors
site_order <- c("GSWS09", "GSWS10", "GSWS01", "GSLOOK", "GSWS02", "GSWS06", "GSWS07", "GSWS08", "GSMACK")
if (!exists("site_colors")) {
  # Generate site colors if not defined
  site_colors <- setNames(
    scales::hue_pal()(length(site_order)),
    site_order
  )
}

p_by_solute_site <- ggplot(cluster_avg,
                           aes(x = month, y = mean_z,
                               group = site_cluster,
                               color = Stream_Name,
                               linetype = Cluster)) +
  geom_line(alpha = 0.7, linewidth = 1.1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  facet_wrap(~ chemical, scales = "free_y", ncol = 3) +
  scale_x_continuous(breaks = 1:12, labels = month_labels) +
  scale_color_site(name = "Site") +
  scale_linetype_manual(values = c("1" = "solid", "2" = "dashed",
                                    "3" = "dotted", "4" = "dotdash"),
                        name = "Cluster") +
  labs(
    x = "Month",
    y = "Mean Normalized Concentration",
    title = "Site-Cluster Behavior by Solute (colored by site)",
    subtitle = "Lines = Site-cluster combinations (mean monthly pattern when that site-solute appears in that cluster) | Line type = Cluster"
  ) +
  theme_hja(base_size = 12) +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  file.path(plot_dir, "annual_cluster_by_solute_colored_by_site.png"),
  p_by_solute_site, width = 14, height = 12, dpi = 300, bg = "white"
)

message("✓ Saved: annual_cluster_by_solute_colored_by_site.png")

# =============================================================================
# FIGURE 3: Each solute gets its own larger panel with ribbons for uncertainty
# =============================================================================

# Create one plot per solute
solute_plots <- list()

for (sol in solute_order) {
  sol_data <- cluster_avg %>% filter(chemical == sol)

  if (nrow(sol_data) == 0) next

  p <- ggplot(sol_data, aes(x = month, y = mean_z,
                            group = site_cluster,
                            color = Cluster, fill = Cluster)) +
    geom_ribbon(aes(ymin = mean_z - se_z, ymax = mean_z + se_z),
                alpha = 0.15, color = NA) +
    geom_line(alpha = 0.7, linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
    scale_x_continuous(breaks = 1:12, labels = month_labels) +
    scale_color_cluster() +
    scale_fill_cluster() +
    labs(
      title = sol,
      x = "Month",
      y = "Norm. Conc."
    ) +
    theme_hja(base_size = 11) +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      legend.position = "none"
    )

  solute_plots[[sol]] <- p
}

# Combine
combined_solute <- wrap_plots(solute_plots, ncol = 3) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Site-Cluster Behavior by Solute (with SE ribbons)",
    subtitle = "Lines = Site-cluster combinations (mean monthly pattern when that site-solute appears in that cluster) ± SE",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      legend.position = "bottom"
    )
  )

ggsave(
  file.path(plot_dir, "annual_cluster_by_solute_with_ribbons.png"),
  combined_solute, width = 16, height = 14, dpi = 300, bg = "white"
)

message("✓ Saved: annual_cluster_by_solute_with_ribbons.png")

message("\n=== COMPLETE ===")
message("All figures saved to:", plot_dir)
