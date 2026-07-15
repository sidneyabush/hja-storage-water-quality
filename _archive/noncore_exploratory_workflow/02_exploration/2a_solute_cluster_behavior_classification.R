# =============================================================================
# SOLUTE CLUSTER BEHAVIOR CLASSIFICATION
# =============================================================================
# Quantify whether solutes are geogenic, biogenic, or transitional based on
# their cluster membership patterns across all sites and years
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)  # For Shannon diversity
})

rm(list = ls())

# Source helpers
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
source(file.path(repo_dir, "00_helpers", "plot_prefs.R"))

# Paths
base_dir   <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality"
output_dir <- file.path(base_dir, "outputs")

# Load data
message("Loading data...")
cluster_wy <- readr::read_csv(
  file.path(output_dir, "ClusterStreams_allSolutes_byWaterYear.csv"),
  show_col_types = FALSE
)

# =============================================================================
# CALCULATE METRICS PER SOLUTE
# =============================================================================

solute_behavior <- cluster_wy %>%
  group_by(chemical, Cluster_climRef) %>%
  summarise(n_occurrences = n(), .groups = "drop") %>%
  group_by(chemical) %>%
  mutate(
    total_occurrences = sum(n_occurrences),
    pct_in_cluster = 100 * n_occurrences / total_occurrences
  ) %>%
  ungroup()

# Calculate diversity metrics per solute
solute_metrics <- solute_behavior %>%
  group_by(chemical) %>%
  summarise(
    # Total occurrences
    n_total = sum(n_occurrences),

    # Number of clusters occupied
    n_clusters = n(),

    # Dominant cluster and its percentage
    dominant_cluster = Cluster_climRef[which.max(pct_in_cluster)],
    pct_in_dominant = max(pct_in_cluster),

    # Shannon diversity (higher = more evenly distributed across clusters)
    # H' = -sum(p_i * log(p_i)) where p_i is proportion in cluster i
    shannon_diversity = diversity(n_occurrences, index = "shannon"),

    # Effective number of clusters (exp(H'))
    # Interprets diversity as "how many equally-common clusters would give this diversity?"
    effective_n_clusters = exp(shannon_diversity),

    # Simpson's evenness (1 - dominance)
    # Close to 1 = evenly distributed, close to 0 = dominated by one cluster
    evenness = diversity(n_occurrences, index = "simpson"),

    # Coefficient of variation in cluster percentages
    cv_cluster_pct = sd(pct_in_cluster) / mean(pct_in_cluster),

    .groups = "drop"
  ) %>%
  arrange(desc(shannon_diversity))

# =============================================================================
# CLASSIFY SOLUTES: Geogenic, Biogenic, or Transitional
# =============================================================================

# Classification based on cluster membership patterns:
# - Geogenic: Strongly in Cluster 1 (diluting, weathering-derived)
# - Biogenic: Strongly in Clusters 2/3/4 (non-diluting, biologically controlled)
# - Transitional: Mixed between Cluster 1 (diluting) AND Clusters 2/3/4 (non-diluting)

# Key insight: Treat Clusters 2, 3, and 4 as a unified "non-diluting" category

# Calculate % in Cluster 1 for each solute
cluster1_pct <- solute_behavior %>%
  filter(Cluster_climRef == "1") %>%
  select(chemical, pct_cluster1 = pct_in_cluster)

solute_classification <- solute_metrics %>%
  left_join(cluster1_pct, by = "chemical") %>%
  mutate(
    pct_cluster1 = replace_na(pct_cluster1, 0),  # If not in cluster 1 at all
    pct_in_C234 = 100 - pct_cluster1,  # % in Clusters 2, 3, or 4 combined
    behavior_type = case_when(
      # Manual override based on known biogeochemistry at HJA
      # Cl: Conservative tracer, not biologically controlled (Kincaid et al. 2024)
      chemical == "Cl" ~ "Transitional",

      # General rules based on cluster behavior
      # Geogenic: >70% in Cluster 1 (diluting behavior)
      pct_cluster1 > 70 ~ "Geogenic",

      # Biogenic: >70% in Clusters 2/3/4 combined (non-diluting behavior)
      # Equivalently: <30% in Cluster 1
      pct_cluster1 < 30 ~ "Biogenic",

      # Transitional: 30-70% in Cluster 1 (mixed diluting/non-diluting)
      TRUE ~ "Transitional"
    )
  ) %>%
  mutate(
    chemical = factor(chemical, levels = solute_order),
    behavior_type = factor(behavior_type, levels = c("Geogenic", "Transitional", "Biogenic"))
  ) %>%
  arrange(chemical)

# =============================================================================
# PRINT RESULTS
# =============================================================================

message("\n=== SOLUTE CLUSTER BEHAVIOR CLASSIFICATION ===\n")

message("Full metrics table:")
print(as.data.frame(solute_classification), row.names = FALSE)

message("\n--- Classification Summary ---")
classification_summary <- solute_classification %>%
  count(behavior_type, name = "n_solutes")
print(classification_summary)

message("\n--- Geogenic Solutes ---")
geogenic <- solute_classification %>%
  filter(behavior_type == "Geogenic") %>%
  select(chemical, dominant_cluster, pct_in_dominant, shannon_diversity)
print(as.data.frame(geogenic), row.names = FALSE)

message("\n--- Biogenic Solutes ---")
biogenic <- solute_classification %>%
  filter(behavior_type == "Biogenic") %>%
  select(chemical, dominant_cluster, pct_in_dominant, shannon_diversity)
print(as.data.frame(biogenic), row.names = FALSE)

message("\n--- Transitional Solutes ---")
transitional <- solute_classification %>%
  filter(behavior_type == "Transitional") %>%
  select(chemical, dominant_cluster, pct_in_dominant, shannon_diversity, effective_n_clusters)
print(as.data.frame(transitional), row.names = FALSE)

# =============================================================================
# VISUALIZE CLUSTER DISTRIBUTION PER SOLUTE
# =============================================================================

# Create output directory
plot_dir <- file.path(base_dir, "exploratory_plots", "02_exploration", "2a_solute_classification")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

p_distribution <- solute_behavior %>%
  left_join(solute_classification %>% select(chemical, behavior_type), by = "chemical") %>%
  mutate(chemical = factor(chemical, levels = solute_order)) %>%
  ggplot(aes(x = chemical, y = pct_in_cluster, fill = factor(Cluster_climRef))) +
  geom_col(position = "stack") +
  facet_wrap(~ behavior_type, scales = "free_x", ncol = 1) +
  scale_fill_cluster() +
  labs(
    x = "Solute",
    y = "Percentage of Occurrences (%)",
    title = "Solute Cluster Membership Distribution",
    subtitle = "Stacked bars show % of site-year occurrences in each cluster",
    fill = "Cluster"
  ) +
  theme_hja(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

ggsave(
  file.path(plot_dir, "solute_cluster_distribution.png"),
  p_distribution, width = 10, height = 10, dpi = 300, bg = "white"
)

message("\n✓ Saved: solute_cluster_distribution.png")

# =============================================================================
# DIVERSITY METRICS PLOT
# =============================================================================

p_diversity <- solute_classification %>%
  ggplot(aes(x = pct_in_dominant, y = shannon_diversity,
             color = behavior_type, label = chemical)) +
  geom_vline(xintercept = 60, linetype = "dashed", color = "gray50", alpha = 0.5) +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "gray50", alpha = 0.5) +
  geom_point(size = 4, alpha = 0.8) +
  ggrepel::geom_text_repel(size = 3.5, max.overlaps = 20) +
  scale_color_manual(
    values = c("Geogenic" = "#F9D5A7", "Biogenic" = "#74B49B", "Transitional" = "#A75D5D"),
    name = "Behavior Type"
  ) +
  labs(
    x = "% in Dominant Cluster",
    y = "Shannon Diversity Index",
    title = "Solute Behavior Classification",
    subtitle = "Transitional solutes have high diversity or weak dominance"
  ) +
  theme_hja(base_size = 12) +
  theme(legend.position = "right")

ggsave(
  file.path(plot_dir, "solute_diversity_classification.png"),
  p_diversity, width = 10, height = 7, dpi = 300, bg = "white"
)

message("✓ Saved: solute_diversity_classification.png")

# =============================================================================
# SAVE CLASSIFICATION TO CSV
# =============================================================================

dir.create(file.path(output_dir), showWarnings = FALSE, recursive = TRUE)

write_csv(
  solute_classification,
  file.path(output_dir, "solute_behavior_classification.csv")
)

message("✓ Saved: solute_behavior_classification.csv")

message("\n=== COMPLETE ===")
