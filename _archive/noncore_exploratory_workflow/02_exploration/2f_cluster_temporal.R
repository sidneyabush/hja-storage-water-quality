# =============================================================================
# CLUSTER MEMBERSHIP OVER TIME: Temporal stability and switching
# =============================================================================
# Shows how cluster membership changes over time for each site-solute
# Helps answer: Are clusters temporally stable or driven by inter-annual variability?
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
plot_dir   <- file.path(base_dir, "exploratory_plots", "02_exploration", "2f_cluster_temporal")
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
# PREPARE DATA
# =============================================================================

temporal_data <- cluster_wy %>%
  select(Stream_Name, chemical, water_year, Cluster_climRef) %>%
  left_join(solute_class %>% select(chemical, behavior_type), by = "chemical") %>%
  mutate(
    Stream_Name = factor(Stream_Name, levels = site_order),
    chemical = factor(chemical, levels = solute_order),
    Cluster = factor(Cluster_climRef, levels = c("1", "2", "3", "4")),
    behavior_type = factor(behavior_type, levels = c("Geogenic", "Transitional", "Biogenic"))
  )

# =============================================================================
# FIGURE 1: All site-solutes over time (tile plot)
# =============================================================================

p_all_temporal <- temporal_data %>%
  mutate(site_solute = paste(Stream_Name, chemical, sep = "-")) %>%
  ggplot(aes(x = water_year, y = site_solute, fill = Cluster)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_cluster() +
  labs(
    x = "Water Year",
    y = "Site-Solute",
    title = "Cluster Membership Over Time",
    subtitle = "Each row = one site-solute combination | Color = cluster assignment each year",
    fill = "Cluster"
  ) +
  theme_hja(base_size = 8) +
  theme(
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    legend.position = "right",
    panel.border = element_blank(),
    axis.ticks = element_blank()
  )

ggsave(
  file.path(plot_dir, "temporal_all_site_solutes.png"),
  p_all_temporal, width = 14, height = 20, dpi = 300, bg = "white"
)

message("✓ Saved: temporal_all_site_solutes.png")

# =============================================================================
# FIGURE 2: Faceted by behavior type
# =============================================================================

p_by_type <- temporal_data %>%
  mutate(site_solute = paste(Stream_Name, chemical, sep = "-")) %>%
  ggplot(aes(x = water_year, y = site_solute, fill = Cluster)) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_grid(behavior_type ~ ., scales = "free_y", space = "free_y") +
  scale_fill_cluster() +
  labs(
    x = "Water Year",
    y = "Site-Solute",
    title = "Cluster Membership Over Time by Solute Type",
    subtitle = "Geogenic solutes show more temporal stability in Cluster 1",
    fill = "Cluster"
  ) +
  theme_hja(base_size = 9) +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "right",
    panel.border = element_blank(),
    axis.ticks = element_blank()
  )

ggsave(
  file.path(plot_dir, "temporal_by_behavior_type.png"),
  p_by_type, width = 14, height = 16, dpi = 300, bg = "white"
)

message("✓ Saved: temporal_by_behavior_type.png")

# =============================================================================
# FIGURE 3: Selected representative site-solutes (line plots)
# =============================================================================

# Select interesting examples: stable vs variable
representative <- bind_rows(
  # Stable geogenic
  temporal_data %>% filter(Stream_Name == "GSLOOK", chemical == "Na") %>% mutate(type = "Stable (Na at GSLOOK)"),
  temporal_data %>% filter(Stream_Name == "GSWS01", chemical == "Ca") %>% mutate(type = "Stable (Ca at GSWS01)"),
  # Variable transitional
  temporal_data %>% filter(Stream_Name == "GSWS02", chemical == "Cl") %>% mutate(type = "Variable (Cl at GSWS02)"),
  temporal_data %>% filter(Stream_Name == "GSMACK", chemical == "SO4") %>% mutate(type = "Variable (SO4 at GSMACK)"),
  # Biogenic
  temporal_data %>% filter(Stream_Name == "GSLOOK", chemical == "NO3") %>% mutate(type = "Biogenic (NO3 at GSLOOK)"),
  temporal_data %>% filter(Stream_Name == "GSWS06", chemical == "DOC") %>% mutate(type = "Biogenic (DOC at GSWS06)")
)

p_examples <- ggplot(representative,
                     aes(x = water_year, y = as.numeric(Cluster),
                         color = Cluster, group = 1)) +
  geom_line(linewidth = 1.0, alpha = 0.7) +
  geom_point(size = 2.5, alpha = 0.9) +
  facet_wrap(~ type, ncol = 2) +
  scale_y_continuous(breaks = 1:4, limits = c(0.5, 4.5)) +
  scale_color_cluster() +
  labs(
    x = "Water Year",
    y = "Cluster",
    title = "Representative Examples of Temporal Cluster Membership",
    subtitle = "Stable vs variable cluster membership patterns over time",
    color = "Cluster"
  ) +
  theme_hja(base_size = 11) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  file.path(plot_dir, "temporal_representative_examples.png"),
  p_examples, width = 12, height = 8, dpi = 300, bg = "white"
)

message("✓ Saved: temporal_representative_examples.png")

# =============================================================================
# FIGURE 4: Cluster transition frequencies
# =============================================================================

# Calculate transitions between consecutive years
transitions <- temporal_data %>%
  arrange(Stream_Name, chemical, water_year) %>%
  group_by(Stream_Name, chemical) %>%
  mutate(
    prev_cluster = lag(Cluster),
    transition = paste(prev_cluster, "→", Cluster)
  ) %>%
  filter(!is.na(prev_cluster)) %>%
  ungroup()

# Count transitions
transition_counts <- transitions %>%
  count(prev_cluster, Cluster, name = "n_transitions") %>%
  complete(prev_cluster = factor(1:4), Cluster = factor(1:4), fill = list(n_transitions = 0))

p_transitions <- ggplot(transition_counts,
                        aes(x = prev_cluster, y = Cluster, fill = n_transitions)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = n_transitions), size = 5, color = "white", fontface = "bold") +
  scale_fill_gradient(low = "#fee5d9", high = "#a50f15",
                      name = "Number of\nTransitions") +
  labs(
    x = "Previous Year Cluster",
    y = "Current Year Cluster",
    title = "Cluster Transition Matrix",
    subtitle = "Diagonal = stable (same cluster) | Off-diagonal = transitions"
  ) +
  theme_hja(base_size = 12) +
  theme(
    legend.position = "right",
    panel.border = element_blank(),
    axis.ticks = element_blank()
  )

ggsave(
  file.path(plot_dir, "cluster_transition_matrix.png"),
  p_transitions, width = 8, height = 7, dpi = 300, bg = "white"
)

message("✓ Saved: cluster_transition_matrix.png")

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

message("\n=== TEMPORAL STABILITY STATISTICS ===\n")

# Calculate stability (% of years same as previous year)
stability_stats <- transitions %>%
  mutate(stable = prev_cluster == Cluster) %>%
  summarise(
    total_transitions = n(),
    n_stable = sum(stable),
    pct_stable = 100 * mean(stable)
  )

message("Overall temporal stability:")
print(stability_stats)

# Stability by solute type
by_type <- transitions %>%
  mutate(stable = prev_cluster == Cluster) %>%
  group_by(behavior_type) %>%
  summarise(
    n_transitions = n(),
    pct_stable = 100 * mean(stable),
    .groups = "drop"
  )

message("\nStability by solute type:")
print(as.data.frame(by_type))

message("\n=== COMPLETE ===")
message("Created 4 temporal cluster figures")
