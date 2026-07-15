# =============================================================================
# CLUSTER ORDINATION: PCA of Monthly Patterns
# =============================================================================
# Visualize cluster separation in reduced dimensional space
# Shows if clusters are well-separated and what drives separation
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(ggrepel)
})

rm(list = ls())

# Source helpers
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
source(file.path(repo_dir, "00_helpers", "plot_prefs.R"))

# Paths
base_dir   <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality"
output_dir <- file.path(base_dir, "outputs")
plot_dir   <- file.path(base_dir, "exploratory_plots", "02_exploration", "2g_cluster_ordination")
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
# PREPARE DATA FOR ORDINATION
# =============================================================================

# Get monthly z-scores for each observation
monthly_matrix <- cluster_wy %>%
  select(Stream_Name, chemical, water_year, Cluster_climRef, `1`:`12`) %>%
  left_join(solute_class %>% select(chemical, behavior_type), by = "chemical")

# Create matrix of monthly values (rows = observations, cols = months)
month_data <- monthly_matrix %>%
  select(`1`:`12`) %>%
  as.matrix()

rownames(month_data) <- paste(monthly_matrix$Stream_Name,
                               monthly_matrix$chemical,
                               monthly_matrix$water_year, sep = "_")

# Remove rows with any NAs
complete_rows <- complete.cases(month_data)
month_data_clean <- month_data[complete_rows, ]
metadata_clean <- monthly_matrix[complete_rows, ]

message(sprintf("Running PCA on %d complete observations...", nrow(month_data_clean)))

# =============================================================================
# RUN PCA
# =============================================================================

pca_result <- rda(month_data_clean, scale = FALSE)  # Already z-scored

# Extract scores
pca_scores <- scores(pca_result, display = "sites", choices = 1:4) %>%
  as.data.frame() %>%
  mutate(
    Stream_Name = metadata_clean$Stream_Name,
    chemical = metadata_clean$chemical,
    water_year = metadata_clean$water_year,
    Cluster = factor(metadata_clean$Cluster_climRef, levels = c("1", "2", "3", "4")),
    behavior_type = factor(metadata_clean$behavior_type,
                          levels = c("Geogenic", "Transitional", "Biogenic"))
  )

# Variance explained
variance_explained <- summary(eigenvals(pca_result))
pc1_var <- round(variance_explained[2, 1] * 100, 1)
pc2_var <- round(variance_explained[2, 2] * 100, 1)
pc3_var <- round(variance_explained[2, 3] * 100, 1)

# =============================================================================
# FIGURE 1: PCA colored by cluster
# =============================================================================

p_pca_cluster <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 2, alpha = 0.6) +
  stat_ellipse(level = 0.68, linewidth = 1, alpha = 0.8) +  # 68% confidence ellipse
  scale_color_cluster() +
  labs(
    x = paste0("PC1 (", pc1_var, "% variance)"),
    y = paste0("PC2 (", pc2_var, "% variance)"),
    title = "PCA of Monthly Concentration Patterns",
    subtitle = "Points = site-solute-year observations | Ellipses = 68% confidence intervals per cluster",
    color = "Cluster"
  ) +
  theme_hja(base_size = 12) +
  theme(legend.position = "right")

ggsave(
  file.path(plot_dir, "pca_by_cluster.png"),
  p_pca_cluster, width = 10, height = 8, dpi = 300, bg = "white"
)

message("✓ Saved: pca_by_cluster.png")

# =============================================================================
# FIGURE 2: PCA colored by solute behavior type
# =============================================================================

p_pca_type <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = behavior_type, shape = Cluster)) +
  geom_point(size = 2.5, alpha = 0.7) +
  scale_color_manual(
    values = c("Geogenic" = "#F9D5A7", "Biogenic" = "#74B49B", "Transitional" = "#A75D5D"),
    name = "Solute Type"
  ) +
  scale_shape_manual(values = c("1" = 16, "2" = 17, "3" = 15, "4" = 18),
                     name = "Cluster") +
  labs(
    x = paste0("PC1 (", pc1_var, "% variance)"),
    y = paste0("PC2 (", pc2_var, "% variance)"),
    title = "PCA by Solute Behavior Type and Cluster",
    subtitle = "Shows relationship between solute type and cluster assignment"
  ) +
  theme_hja(base_size = 12) +
  theme(legend.position = "right")

ggsave(
  file.path(plot_dir, "pca_by_solute_type.png"),
  p_pca_type, width = 10, height = 8, dpi = 300, bg = "white"
)

message("✓ Saved: pca_by_solute_type.png")

# =============================================================================
# FIGURE 3: PCA faceted by behavior type
# =============================================================================

p_pca_facet <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 2, alpha = 0.7) +
  stat_ellipse(level = 0.68, linewidth = 0.8) +
  facet_wrap(~ behavior_type, ncol = 3) +
  scale_color_cluster() +
  labs(
    x = paste0("PC1 (", pc1_var, "% variance)"),
    y = paste0("PC2 (", pc2_var, "% variance)"),
    title = "PCA by Solute Behavior Type",
    subtitle = "Geogenic solutes cluster more tightly than transitional/biogenic",
    color = "Cluster"
  ) +
  theme_hja(base_size = 11) +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  file.path(plot_dir, "pca_faceted_by_type.png"),
  p_pca_facet, width = 14, height = 5, dpi = 300, bg = "white"
)

message("✓ Saved: pca_faceted_by_type.png")

# =============================================================================
# FIGURE 4: PC1 vs PC3 (alternative view)
# =============================================================================

p_pca_13 <- ggplot(pca_scores, aes(x = PC1, y = PC3, color = Cluster)) +
  geom_point(size = 2, alpha = 0.6) +
  stat_ellipse(level = 0.68, linewidth = 1) +
  scale_color_cluster() +
  labs(
    x = paste0("PC1 (", pc1_var, "% variance)"),
    y = paste0("PC3 (", pc3_var, "% variance)"),
    title = "PCA: PC1 vs PC3",
    subtitle = "Alternative ordination view",
    color = "Cluster"
  ) +
  theme_hja(base_size = 12) +
  theme(legend.position = "right")

ggsave(
  file.path(plot_dir, "pca_pc1_vs_pc3.png"),
  p_pca_13, width = 10, height = 8, dpi = 300, bg = "white"
)

message("✓ Saved: pca_pc1_vs_pc3.png")

# =============================================================================
# FIGURE 5: Loading plot (month contributions)
# =============================================================================

loadings <- scores(pca_result, display = "species", choices = 1:2) %>%
  as.data.frame() %>%
  mutate(month = 1:12,
         month_label = month_labels)

p_loadings <- ggplot(loadings, aes(x = PC1, y = PC2, label = month_label)) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "gray30", linewidth = 0.8) +
  geom_text_repel(size = 4, fontface = "bold", max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    x = paste0("PC1 loadings (", pc1_var, "%)"),
    y = paste0("PC2 loadings (", pc2_var, "%)"),
    title = "PCA Loadings: Month Contributions",
    subtitle = "Arrows show how each month contributes to PC axes"
  ) +
  theme_hja(base_size = 12)

ggsave(
  file.path(plot_dir, "pca_loadings.png"),
  p_loadings, width = 8, height = 8, dpi = 300, bg = "white"
)

message("✓ Saved: pca_loadings.png")

# =============================================================================
# STATISTICS
# =============================================================================

message("\n=== PCA SUMMARY ===\n")
message(sprintf("PC1 explains %.1f%% of variance", pc1_var))
message(sprintf("PC2 explains %.1f%% of variance", pc2_var))
message(sprintf("PC3 explains %.1f%% of variance", pc3_var))
message(sprintf("First 3 PCs explain %.1f%% of total variance",
                sum(variance_explained[2, 1:3]) * 100))

message("\n=== COMPLETE ===")
message("Created 5 ordination figures")
