# =============================================================================
# CLUSTER CHARACTERIZATION: YEARLY AVERAGE PATTERNS
# =============================================================================
# Shows how site-solutes behave WITHIN each cluster
# (averaged across the years they appear in that cluster)
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
# PREPARE DATA: Average monthly pattern per site-solute per cluster
# =============================================================================

# For each site-solute-cluster combination, average across years
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
  filter(n_years >= 1) %>%  # Keep site-solutes even if only 1 year in cluster
  mutate(
    Cluster = factor(Cluster_climRef, levels = c("1", "2", "3", "4")),
    chemical = factor(chemical, levels = solute_order)
  )

# Count site-solutes per cluster
cluster_counts <- cluster_avg %>%
  distinct(Stream_Name, chemical, Cluster) %>%
  count(Cluster, name = "n_site_solutes")

message("\nSite-solute counts per cluster:")
print(cluster_counts)

# =============================================================================
# FIGURE 1: Cluster characterization with yearly averages
# =============================================================================

cluster_panels <- list()

for (cl in levels(cluster_avg$Cluster)) {
  cl_data <- cluster_avg %>% filter(Cluster == cl)
  n_site_solutes <- cluster_counts %>% filter(Cluster == cl) %>% pull(n_site_solutes)

  # Left panel: Average patterns per site-solute (colored by solute)
  p_left <- ggplot(cl_data, aes(x = month, y = mean_z,
                                group = interaction(Stream_Name, chemical),
                                color = chemical)) +
    geom_line(alpha = 0.7, linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
    scale_x_continuous(breaks = 1:12, labels = month_labels) +
    scale_color_solute(name = "Solute") +
    labs(
      x = "Month",
      y = "Mean Normalized Concentration",
      title = paste0("Cluster ", cl, " (n=", n_site_solutes, " site-solutes)"),
      subtitle = "Average behavior when in this cluster",
      color = "Solute"
    ) +
    theme_hja(base_size = 11) +
    theme(
      legend.position = "none"
    )

  # Right panel: Solute composition
  solute_counts <- cl_data %>%
    distinct(Stream_Name, chemical) %>%
    count(chemical) %>%
    complete(chemical = factor(solute_order, levels = solute_order), fill = list(n = 0)) %>%
    mutate(chemical = factor(chemical, levels = solute_order))

  p_right <- ggplot(solute_counts, aes(x = chemical, y = n, fill = chemical)) +
    geom_col(alpha = 0.8) +
    scale_fill_solute(name = "Solute") +
    labs(x = "Solute", y = "Count", title = "Site-solute count") +
    theme_hja(base_size = 11) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
    )

  cluster_panels[[cl]] <- p_left + p_right + plot_layout(widths = c(2, 1))
}

combined <- wrap_plots(cluster_panels, ncol = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Cluster Characterization: Average Yearly Patterns",
    subtitle = "Lines = Site-solute combinations (mean monthly pattern when that site-solute appears in that cluster, averaged across all years)",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

ggsave(
  file.path(plot_dir, "annual_cluster_avg_patterns_per_site_solute.png"),
  combined, width = 14, height = 18, dpi = 300, bg = "white"
)

message("\n✓ Saved: annual_cluster_avg_patterns_per_site_solute.png")

# =============================================================================
# FIGURE 2: Faceted view showing all site-solutes per cluster
# =============================================================================

p_faceted <- ggplot(cluster_avg, aes(x = month, y = mean_z,
                                     group = interaction(Stream_Name, chemical),
                                     color = chemical)) +
  geom_line(alpha = 0.6, linewidth = 1.0) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  facet_wrap(~ Cluster, ncol = 2, labeller = labeller(Cluster = function(x) paste("Cluster", x))) +
  scale_x_continuous(breaks = 1:12, labels = month_labels) +
  scale_color_solute(name = "Solute") +
  labs(
    x = "Month",
    y = "Mean Normalized Concentration",
    title = "Average Site-Solute Behavior Per Cluster",
    subtitle = "Lines = Site-solute combinations (mean monthly pattern when that site-solute appears in that cluster, averaged across all years)",
    color = "Solute"
  ) +
  theme_hja(base_size = 12) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  file.path(plot_dir, "annual_cluster_faceted_avg_patterns.png"),
  p_faceted, width = 12, height = 10, dpi = 300, bg = "white"
)

message("✓ Saved: annual_cluster_faceted_avg_patterns.png")

message("\n=== COMPLETE ===")
message("All figures saved to:", plot_dir)
