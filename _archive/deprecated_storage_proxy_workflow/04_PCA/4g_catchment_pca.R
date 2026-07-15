# Restored from archive on 2025-12-04

suppressPackageStartupMessages({
  library(tidyverse)
  library(FactoMineR)
  library(factoextra)
  library(corrplot)
  library(patchwork)
})
try(source(file.path("/Users/sidneybush/Documents/GitHub/hja-water-quality", "00_helpers", "plot_theme_set.R")), silent = TRUE)

repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
try(source(file.path(repo_dir, "00_helpers", "plot_prefs.R")), silent = TRUE)

base_dir   <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality"
output_dir <- file.path(base_dir, "outputs")
plot_dir   <- file.path(base_dir, "exploratory_plots", "04_PCA", "4g_catchment")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

theme_clean <- function(base_size = 12) {
  theme_bw(base_size = base_size) + theme(panel.grid = element_blank(), strip.background = element_blank())
}

catchment <- read_csv(file.path(output_dir, "Catchment_site_characteristics.csv"), show_col_types = FALSE)
site_data <- read_csv(file.path(output_dir, "HJA_exploratory_site.csv"), show_col_types = FALSE)
mega <- read_csv(file.path(output_dir, "HJA_mega_90d_windows_CQ_hydro_static.csv"), show_col_types = FALSE)
clusters_modal <- read_csv(file.path(output_dir, "ClusterStreams_allSolutes_modalClusters.csv"), show_col_types = FALSE)
stability <- read_csv(file.path(output_dir, "ClusterStreams_stability_metrics.csv"), show_col_types = FALSE)

# Placeholder: implement PCA and joins as needed
message("Catchment PCA script restored. Fill in analysis steps if needed.")
# =============================================================================
# STEP 04g: CATCHMENT CHARACTERISTICS PCA
# =============================================================================
# Goal: Use PCA to understand how catchment characteristics relate to:
#   - CQ behavior patterns
#   - Cluster membership
#   - Synchrony patterns
#
# ═══════════════════════════════════════════════════════════════════════════
# CRITICAL DISTINCTION (2025 Update):
# ═══════════════════════════════════════════════════════════════════════════
# This PCA relates catchment properties to stream chemistry. Key findings:
#
# 1. CATCHMENT → CQ BEHAVIOR: STRONG (R²m ≈ 7%)
#    - Geology/topography controls WHAT CQ behavior a site exhibits
#    - DSi: Lava1_per (r=-0.84), K: Lava1_per (r=-0.88)
#    - PO4: Pyro_per (r=0.87), Ca/Mg/Na: Ash_Per (r=0.75-0.87)
#    - Clusters show significant CQ slope differences (F=237.5, p<0.001)
#
# 2. CATCHMENT → OUTLET SYNCHRONY: WEAK (R²m < 1%)
#    - Catchment characteristics DON'T predict how similar sites are to outlet
#    - NULL MODEL WINS for outlet synchrony prediction
#    - MTT (Tier 2 isotope) explains 37% of CQ-slope sync with outlet
#
# KEY INSIGHT: Geology controls WHAT behavior, not HOW SIMILAR to outlet
#
# References:
#   - 03_stats/3l_outlet_synchrony_predictors.R (tiered modeling)
#   - 03_stats/3m_catchment_cq_concentration_predictors.R (catchment→CQ)
# ═══════════════════════════════════════════════════════════════════════════
#
# Variables:
#   - Physical: Area, Elevation, Slope, Aspect
#   - Land use: Harvest %, Age
#   - Geology: Lava1, Lava2, Ash, Pyro percentages
#   - Hydrology: DR, MTT, FYw
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(FactoMineR)
  library(factoextra)
  library(corrplot)
  library(patchwork)
})

rm(list = ls())

# Source helpers
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
try(source(file.path(repo_dir, "00_helpers", "plot_prefs.R")), silent = TRUE)

# Paths
base_dir   <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality"
output_dir <- file.path(base_dir, "outputs")
plot_dir   <- file.path(base_dir, "exploratory_plots", "04_PCA", "4g_catchment")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Theme
theme_clean <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(panel.grid = element_blank(), strip.background = element_blank())
}

message("\n=== LOADING DATA ===\n")

# Load catchment characteristics
catchment <- read_csv(file.path(output_dir, "Catchment_site_characteristics.csv"), show_col_types = FALSE)

# Load CQ and cluster data to join behavior metrics
site_data <- read_csv(file.path(output_dir, "HJA_exploratory_site.csv"), show_col_types = FALSE)
mega <- read_csv(file.path(output_dir, "HJA_mega_90d_windows_CQ_hydro_static.csv"), show_col_types = FALSE)
clusters_modal <- read_csv(file.path(output_dir, "ClusterStreams_allSolutes_modalClusters.csv"), show_col_types = FALSE)
stability <- read_csv(file.path(output_dir, "ClusterStreams_stability_metrics.csv"), show_col_types = FALSE)

message("  Catchment data:", nrow(catchment), "sites\n")

# #############################################################################
# BUILD INTEGRATED DATASET
# #############################################################################
message("\n=== BUILDING INTEGRATED DATASET ===\n")

# Site-level CQ behavior summary
site_behavior <- mega %>%
  group_by(Stream_Name) %>%
  summarise(
    mean_cq = mean(cq_slope.x, na.rm = TRUE),
    pct_positive = mean(cq_slope.x > 0, na.rm = TRUE) * 100,
    pct_sync = mean(cq_sync == "sync", na.rm = TRUE) * 100,
    cq_variability = sd(cq_slope.x, na.rm = TRUE),
    mean_storage = mean(Q_dS_range_mm, na.rm = TRUE),
    mean_RBI = mean(RBI, na.rm = TRUE),
    .groups = "drop"
  )

# Site-level stability
site_stability <- stability %>%
  group_by(Stream_Name) %>%
  summarise(mean_stability = mean(stability, na.rm = TRUE), .groups = "drop")

# Combine all data
pca_data <- catchment %>%
  left_join(site_behavior, by = "Stream_Name") %>%
  left_join(site_stability, by = "Stream_Name")

message("  Combined dataset:", nrow(pca_data), "sites\n")

# #############################################################################
# PCA: CATCHMENT CHARACTERISTICS ONLY
# #############################################################################
message("\n=== PCA: CATCHMENT CHARACTERISTICS ===\n")

# Select catchment variables
catchment_vars <- c("Area_km2", "Elevation_mean_m", "Slope_mean", "Harvest", "Age",
                    "Lava1_per", "Lava2_per", "Ash_Per", "Pyro_per",
                    "DR_Overall", "MTT_overall", "Fyw_overall")

# Filter to available and valid
catchment_vars <- catchment_vars[catchment_vars %in% names(pca_data)]
pca_input <- pca_data %>%
  select(Stream_Name, all_of(catchment_vars)) %>%
  column_to_rownames("Stream_Name") %>%
  drop_na()

if (nrow(pca_input) >= 5 && ncol(pca_input) >= 3) {
  # Run PCA
  pca_result <- PCA(pca_input, scale.unit = TRUE, graph = FALSE)
  
  # Scree plot
  p_scree <- fviz_screeplot(pca_result, addlabels = TRUE, 
                             title = "Catchment Characteristics PCA - Variance Explained")
  ggsave(file.path(plot_dir, "01_catchment_screeplot.png"), p_scree, width = 10, height = 7, dpi = 300)
  
  # Variable contributions
  p_var <- fviz_pca_var(pca_result, col.var = "contrib",
                        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                        repel = TRUE,
                        title = "Catchment Variables - PCA Loadings")
  ggsave(file.path(plot_dir, "02_catchment_variable_plot.png"), p_var, width = 10, height = 9, dpi = 300)
  
  # Individual sites
  p_ind <- fviz_pca_ind(pca_result, col.ind = "cos2",
                        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                        repel = TRUE,
                        title = "Sites in Catchment PCA Space")
  ggsave(file.path(plot_dir, "03_catchment_site_plot.png"), p_ind, width = 10, height = 9, dpi = 300)
  
  # Biplot
  p_biplot <- fviz_pca_biplot(pca_result, 
                              col.var = "#2166AC",
                              col.ind = "cos2",
                              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                              repel = TRUE,
                              title = "Catchment Characteristics - Biplot")
  ggsave(file.path(plot_dir, "04_catchment_biplot.png"), p_biplot, width = 12, height = 10, dpi = 300)
  
  # Extract PC scores for regression
  pc_scores <- as.data.frame(pca_result$ind$coord) %>%
    rownames_to_column("Stream_Name") %>%
    left_join(pca_data, by = "Stream_Name")
  
  message("\n  PC1 explains", round(pca_result$eig[1,2], 1), "% of variance\n")
  message("  PC2 explains", round(pca_result$eig[2,2], 1), "% of variance\n")
  
  # Variable contributions
  message("\n  Top PC1 contributors:\n")
  contrib1 <- pca_result$var$contrib[,1] %>% sort(decreasing = TRUE)
  print(head(contrib1, 5))
  
  # Save contributions
  var_contrib <- as.data.frame(pca_result$var$contrib) %>%
    rownames_to_column("Variable")
  write_csv(var_contrib, file.path(output_dir, "04_PCA/catchment_pca_contributions.csv"))
}

# #############################################################################
# RELATE PCA AXES TO CQ BEHAVIOR
# #############################################################################
message("\n=== PC AXES VS CQ BEHAVIOR ===\n")

if (exists("pc_scores") && "mean_cq" %in% names(pc_scores)) {
  # PC1 vs CQ slope
  if (sum(!is.na(pc_scores$Dim.1) & !is.na(pc_scores$mean_cq)) >= 5) {
    cor1 <- cor.test(pc_scores$Dim.1, pc_scores$mean_cq)
    message("  PC1 vs mean CQ: r =", round(cor1$estimate, 3), ", p =", round(cor1$p.value, 4), "\n")
    
    p_pc1_cq <- ggplot(pc_scores, aes(x = Dim.1, y = mean_cq)) +
      geom_point(size = 4, alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE, color = "#2166AC") +
      geom_text(aes(label = Stream_Name), vjust = -1, size = 3) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(x = "Catchment PC1", y = "Mean CQ Slope",
           title = paste0("Catchment PC1 vs CQ Slope (r = ", round(cor1$estimate, 2), ")")) +
      theme_clean()
    
    ggsave(file.path(plot_dir, "05_pc1_vs_cq.png"), p_pc1_cq, width = 10, height = 8, dpi = 300)
  }
  
  # PC1 vs synchrony
  if (sum(!is.na(pc_scores$Dim.1) & !is.na(pc_scores$pct_sync)) >= 5) {
    cor2 <- cor.test(pc_scores$Dim.1, pc_scores$pct_sync)
    message("  PC1 vs % sync: r =", round(cor2$estimate, 3), ", p =", round(cor2$p.value, 4), "\n")
    
    p_pc1_sync <- ggplot(pc_scores, aes(x = Dim.1, y = pct_sync)) +
      geom_point(size = 4, alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE, color = "#2166AC") +
      geom_text(aes(label = Stream_Name), vjust = -1, size = 3) +
      labs(x = "Catchment PC1", y = "% Synchronous",
           title = paste0("Catchment PC1 vs Synchrony (r = ", round(cor2$estimate, 2), ")")) +
      theme_clean()
    
    ggsave(file.path(plot_dir, "06_pc1_vs_synchrony.png"), p_pc1_sync, width = 10, height = 8, dpi = 300)
  }
  
  # PC1 vs stability
  if (sum(!is.na(pc_scores$Dim.1) & !is.na(pc_scores$mean_stability)) >= 5) {
    cor3 <- cor.test(pc_scores$Dim.1, pc_scores$mean_stability)
    message("  PC1 vs stability: r =", round(cor3$estimate, 3), ", p =", round(cor3$p.value, 4), "\n")
    
    p_pc1_stab <- ggplot(pc_scores, aes(x = Dim.1, y = mean_stability)) +
      geom_point(size = 4, alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE, color = "#2166AC") +
      geom_text(aes(label = Stream_Name), vjust = -1, size = 3) +
      labs(x = "Catchment PC1", y = "Mean Cluster Stability",
           title = paste0("Catchment PC1 vs Stability (r = ", round(cor3$estimate, 2), ")")) +
      theme_clean()
    
    ggsave(file.path(plot_dir, "07_pc1_vs_stability.png"), p_pc1_stab, width = 10, height = 8, dpi = 300)
  }
}

# #############################################################################
# PHYSICAL VS GEOLOGY PCA
# #############################################################################
message("\n=== PHYSICAL VS GEOLOGY COMPONENTS ===\n")

# Run separate PCA for physical and geology
physical_vars <- c("Area_km2", "Elevation_mean_m", "Slope_mean", "Harvest", "Age")
geology_vars <- c("Lava1_per", "Lava2_per", "Ash_Per", "Pyro_per")
hydro_vars <- c("DR_Overall", "MTT_overall", "Fyw_overall")

# Identify which sites have what geology
geology_composition <- pca_data %>%
  select(Stream_Name, all_of(geology_vars[geology_vars %in% names(pca_data)])) %>%
  mutate(
    dominant_geology = case_when(
      Lava1_per > 50 ~ "Lava1 Dominant",
      Lava2_per > 50 ~ "Lava2 Dominant",
      Ash_Per > 30 ~ "Ash Rich",
      TRUE ~ "Mixed"
    )
  )

message("\nGeology composition by site:\n")
print(geology_composition)

# Save
write_csv(geology_composition, file.path(output_dir, "04_PCA/site_geology_composition.csv"))

# Plot sites colored by dominant geology in PCA space
if (exists("pc_scores")) {
  pc_scores <- pc_scores %>%
    left_join(geology_composition %>% select(Stream_Name, dominant_geology), by = "Stream_Name")
  
  p_geol <- ggplot(pc_scores, aes(x = Dim.1, y = Dim.2, color = dominant_geology)) +
    geom_point(size = 5) +
    geom_text(aes(label = Stream_Name), vjust = -1, size = 3, color = "black") +
    scale_color_brewer(palette = "Set1", name = "Dominant\nGeology") +
    labs(x = paste0("PC1 (", round(pca_result$eig[1,2], 1), "%)"),
         y = paste0("PC2 (", round(pca_result$eig[2,2], 1), "%)"),
         title = "Sites in Catchment PCA Space by Geology") +
    theme_clean()
  
  ggsave(file.path(plot_dir, "08_pca_by_geology.png"), p_geol, width = 11, height = 9, dpi = 300)
}

# #############################################################################
# CORRELATION MATRIX
# #############################################################################
message("\n=== CATCHMENT CORRELATION MATRIX ===\n")

all_vars <- c(catchment_vars, "mean_cq", "pct_sync", "mean_stability")
all_vars <- all_vars[all_vars %in% names(pca_data)]

corr_data <- pca_data %>%
  select(all_of(all_vars)) %>%
  drop_na()

if (nrow(corr_data) >= 5) {
  corr_mat <- cor(corr_data, use = "pairwise.complete.obs")
  
  png(file.path(plot_dir, "09_catchment_correlations.png"), width = 12, height = 10, units = "in", res = 300)
  corrplot(corr_mat, method = "color", type = "upper",
           addCoef.col = "black", number.cex = 0.6,
           tl.col = "black", tl.srt = 45, tl.cex = 0.7,
           title = "Catchment-Behavior Correlations", mar = c(0,0,2,0))
  dev.off()
}

# #############################################################################
# SUMMARY
# #############################################################################
message("\n╔════════════════════════════════════════════════════════════════╗\n")
message("║  CATCHMENT PCA ANALYSIS COMPLETE                              ║\n")
message("╚════════════════════════════════════════════════════════════════╝\n\n")

message("Outputs saved to:\n")
message("  Plots: ", plot_dir, "\n")
message("  Data:  ", file.path(output_dir, "04_PCA"), "\n\n")
