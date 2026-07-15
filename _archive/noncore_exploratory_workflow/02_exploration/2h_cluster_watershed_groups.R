# =============================================================================
# CLUSTER PATTERNS BY WATERSHED GROUP
# =============================================================================
# Shows cluster centroids grouped by watershed position
# Group 1 (lower): GSMACK, GSWS08, GSWS07, GSWS06
# Group 2 (mid): GSWS02, GSLOOK, GSWS01
# Group 3 (upper): GSWS10, GSWS09
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
plot_dir   <- file.path(base_dir, "exploratory_plots", "02_exploration", "2h_watershed_groups")
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
# DEFINE WATERSHED GROUPS
# =============================================================================

watershed_groups <- tibble(
  Stream_Name = c("GSMACK", "GSWS08", "GSWS07", "GSWS06",
                  "GSWS02", "GSLOOK", "GSWS01",
                  "GSWS10", "GSWS09"),
  group = c(rep("Group 1: Lower elevation", 4),
            rep("Group 2: Mid elevation", 3),
            rep("Group 3: Upper elevation", 2)),
  group_order = c(rep(1, 4), rep(2, 3), rep(3, 2))
) %>%
  mutate(group = factor(group, levels = c("Group 1: Lower elevation",
                                          "Group 2: Mid elevation",
                                          "Group 3: Upper elevation")))

# =============================================================================
# COMPUTE GROUP CENTROIDS
# =============================================================================

group_centroids <- cluster_wy %>%
  select(Stream_Name, chemical, Cluster_climRef, `1`:`12`) %>%
  left_join(watershed_groups, by = "Stream_Name") %>%
  pivot_longer(`1`:`12`, names_to = "month", values_to = "z_norm") %>%
  mutate(month = as.numeric(month)) %>%
  group_by(group, chemical, Cluster_climRef, month) %>%
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
# FIGURE 1: Centroids by watershed group (faceted by cluster)
# =============================================================================

p_by_cluster <- ggplot(group_centroids,
                       aes(x = month, y = mean_z,
                           color = group, linetype = group)) +
  geom_line(linewidth = 1.0, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  facet_grid(Cluster ~ chemical, scales = "free_y",
             labeller = labeller(Cluster = function(x) paste("Cluster", x))) +
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = month_labels[c(1, 4, 7, 10)]) +
  scale_color_manual(values = c("#08519c", "#3182bd", "#9ecae1"),
                     name = "Watershed\nGroup") +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"),
                        name = "Watershed\nGroup") +
  labs(
    x = "Month",
    y = "Normalized Concentration",
    title = "Cluster Behavior by Watershed Group",
    subtitle = "Rows = Cluster | Columns = Solute | Lines = Watershed group means"
  ) +
  theme_hja(base_size = 9) +
  theme(
    panel.spacing = unit(0.2, "lines"),
    strip.text.y = element_text(size = 8, face = "bold", angle = 0, hjust = 0),
    strip.text.x = element_text(size = 8, face = "bold", angle = 90, hjust = 0),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    legend.position = "right"
  )

ggsave(
  file.path(plot_dir, "centroids_by_watershed_group_and_cluster.png"),
  p_by_cluster, width = 18, height = 12, dpi = 300, bg = "white"
)

message("✓ Saved: centroids_by_watershed_group_and_cluster.png")

# =============================================================================
# FIGURE 2: One panel per watershed group
# =============================================================================

p_by_group <- ggplot(group_centroids,
                     aes(x = month, y = mean_z,
                         color = Cluster, fill = Cluster)) +
  geom_ribbon(aes(ymin = mean_z - sd_z, ymax = mean_z + sd_z),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.0, alpha = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  facet_grid(group ~ chemical, scales = "free_y") +
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = month_labels[c(1, 4, 7, 10)]) +
  scale_color_cluster() +
  scale_fill_cluster() +
  labs(
    x = "Month",
    y = "Normalized Concentration",
    title = "Cluster Centroids by Watershed Group",
    subtitle = "Rows = Watershed group | Columns = Solute | Lines = Cluster means ± SD",
    color = "Cluster",
    fill = "Cluster"
  ) +
  theme_hja(base_size = 9) +
  theme(
    panel.spacing = unit(0.2, "lines"),
    strip.text.y = element_text(size = 9, face = "bold", angle = 0, hjust = 0),
    strip.text.x = element_text(size = 8, face = "bold", angle = 90, hjust = 0),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    legend.position = "right"
  )

ggsave(
  file.path(plot_dir, "centroids_by_group_panels.png"),
  p_by_group, width = 18, height = 10, dpi = 300, bg = "white"
)

message("✓ Saved: centroids_by_group_panels.png")

# =============================================================================
# FIGURE 3: Geogenic solutes only by watershed group
# =============================================================================

p_geogenic <- group_centroids %>%
  filter(behavior_type == "Geogenic") %>%
  ggplot(aes(x = month, y = mean_z,
             color = Cluster, fill = Cluster)) +
  geom_ribbon(aes(ymin = mean_z - sd_z, ymax = mean_z + sd_z),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.1, alpha = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  facet_grid(group ~ chemical, scales = "free_y") +
  scale_x_continuous(breaks = 1:12, labels = month_labels) +
  scale_color_cluster() +
  scale_fill_cluster() +
  labs(
    x = "Month",
    y = "Normalized Concentration",
    title = "Geogenic Solute Patterns by Watershed Group",
    subtitle = "Rows = Watershed group | Columns = Geogenic solutes | All groups show strong Cluster 1 (diluting) behavior",
    color = "Cluster",
    fill = "Cluster"
  ) +
  theme_hja(base_size = 10) +
  theme(
    strip.text = element_text(size = 9, face = "bold"),
    legend.position = "right"
  )

ggsave(
  file.path(plot_dir, "geogenic_by_watershed_group.png"),
  p_geogenic, width = 12, height = 8, dpi = 300, bg = "white"
)

message("✓ Saved: geogenic_by_watershed_group.png")

# =============================================================================
# SUMMARY: Cluster membership by watershed group
# =============================================================================

message("\n=== CLUSTER MEMBERSHIP BY WATERSHED GROUP ===\n")

cluster_counts <- cluster_wy %>%
  left_join(watershed_groups, by = "Stream_Name") %>%
  group_by(group, Cluster_climRef) %>%
  summarise(n_obs = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(pct = 100 * n_obs / sum(n_obs)) %>%
  arrange(group, Cluster_climRef)

print(as.data.frame(cluster_counts))

message("\n=== COMPLETE ===")
message("Created 3 watershed group figures")
