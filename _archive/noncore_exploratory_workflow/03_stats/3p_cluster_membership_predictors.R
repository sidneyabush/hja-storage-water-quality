# =============================================================================
# 3p_cluster_membership_predictors.R
# =============================================================================
# PREDICTING CLUSTER MEMBERSHIP FROM CATCHMENT PROPERTIES
#
# RESEARCH QUESTION:
#   What catchment characteristics predict which seasonal concentration pattern
#   cluster a site-solute belongs to?
#
# NOTE: Clusters are based on DTW clustering of monthly concentration z-scores,
#       NOT CQ behavior. See plot_prefs.R for cluster label definitions.
#
# HYPOTHESES:
#   - If geology determines cluster: geology predictors should be significant
#   - If cluster is dynamic: temporal/storage variables matter more
#
# APPROACH:
#   1. Multinomial logistic regression (cluster ~ catchment)
#   2. Random forest for variable importance
#   3. Cluster stability analysis
#
# OUTPUT:
#   - Predictor importance rankings
#   - Confusion matrices
#   - Cluster-catchment visualizations
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(nnet)          # Multinomial logistic
  library(randomForest)  # Variable importance
  library(patchwork)
  library(broom)
})

rm(list = ls())

# Source helpers
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
try(source(file.path(repo_dir, "00_helpers", "plot_prefs.R")), silent = TRUE)

# Paths
base_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality"
out_dir <- file.path(base_dir, "outputs")
data_dir <- file.path(base_dir, "data")
fig_dir <- file.path(base_dir, "exploratory_plots", "03_stats", "3p_cluster_membership")
res_dir <- file.path(out_dir, "03_stats")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Theme
theme_clean <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(panel.grid = element_blank(), strip.background = element_blank())
}

message("\n╔════════════════════════════════════════════════════════════════╗\n")
message("║  CLUSTER MEMBERSHIP PREDICTION                                 ║\n")
message("╚════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# 1. LOAD DATA
# =============================================================================

message("=== 1. LOADING DATA ===\n\n")

# Modal cluster assignments (site × solute)
clusters_modal <- read_csv(file.path(out_dir, "ClusterStreams_allSolutes_modalClusters.csv"),
                           show_col_types = FALSE)

# Cluster stability
clusters_stability <- read_csv(file.path(out_dir, "ClusterStreams_stability_metrics.csv"),
                               show_col_types = FALSE)

# Site characteristics (cleaned - no McGuire etc)
site_chars <- read_csv(file.path(out_dir, "Catchment_site_characteristics.csv"),
                       show_col_types = FALSE)

message("  Modal clusters:", nrow(clusters_modal), "site-solute combinations\n")
message("  Stability metrics:", nrow(clusters_stability), "records\n")
message("  Site characteristics:", nrow(site_chars), "sites\n")

# =============================================================================
# 2. PREPARE DATA FOR MODELING
# =============================================================================

message("\n=== 2. PREPARING DATA ===\n\n")

# Join clusters with site characteristics
model_data <- clusters_modal %>%
  rename(solute = chemical) %>%
  left_join(site_chars, by = "Stream_Name") %>%
  left_join(clusters_stability %>% rename(solute = chemical), by = c("Stream_Name", "solute")) %>%
  mutate(
    Cluster = factor(Cluster_mode, levels = 1:4, 
                     labels = c("Baseflow Enriched", "Chemostatic", 
                               "Spring/Early Summer Enriched", "Winter Flushing")),
    solute = factor(solute)
  ) %>%
  filter(!is.na(Cluster))

message("  Model data rows:", nrow(model_data), "\n")
message("  Cluster distribution:\n")
print(table(model_data$Cluster))

# Define predictors
geo_preds <- c("Lava1_per", "Lava2_per", "Ash_Per", "Pyro_per")
topo_preds <- c("Area_km2", "Elevation_mean_m", "Slope_mean")
all_preds <- c(geo_preds, topo_preds)

# Check which predictors are available
available_preds <- intersect(all_preds, names(model_data))
message("\n  Available predictors:", paste(available_preds, collapse = ", "), "\n")

# =============================================================================
# 3. MULTINOMIAL LOGISTIC REGRESSION
# =============================================================================

message("\n=== 3. MULTINOMIAL LOGISTIC REGRESSION ===\n\n")

# Full model with geology + topography
model_formula <- as.formula(paste("Cluster ~", paste(available_preds, collapse = " + ")))

# Use complete cases for consistent comparison
model_data_complete <- model_data %>%
  select(Cluster, all_of(available_preds)) %>%
  filter(complete.cases(.)) %>%
  mutate(Cluster = droplevels(Cluster))

# Fit multinomial model
model_full <- nnet::multinom(model_formula, data = model_data_complete, trace = FALSE)

message("Full model (geology + topography):\n")
print(summary(model_full))

# Calculate pseudo R-squared (McFadden's)
null_model <- nnet::multinom(Cluster ~ 1, data = model_data_complete, trace = FALSE)
pseudo_r2 <- 1 - (logLik(model_full) / logLik(null_model))
message("\nPseudo R² (McFadden):", round(as.numeric(pseudo_r2), 3), "\n")

# Likelihood ratio test
lr_test <- tryCatch({
  anova(null_model, model_full, test = "Chisq")
}, error = function(e) {
  message("\nNote: Likelihood ratio test failed -", e$message, "\n")
  NULL
})
if (!is.null(lr_test)) {
  message("\nLikelihood ratio test:\n")
  print(lr_test)
}

# Coefficients table
coef_table <- tidy(model_full, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(term != "(Intercept)") %>%
  arrange(y.level, desc(abs(estimate - 1)))

message("\nOdds Ratios (exponentiated coefficients):\n")
print(coef_table %>% select(y.level, term, estimate, p.value) %>% 
        mutate(estimate = round(estimate, 3), p.value = round(p.value, 3)), n = 30)

# =============================================================================
# 4. RANDOM FOREST VARIABLE IMPORTANCE
# =============================================================================

message("\n=== 4. RANDOM FOREST VARIABLE IMPORTANCE ===\n\n")

# Prepare data for RF (complete cases only)
rf_data <- model_data %>%
  select(Cluster, all_of(available_preds)) %>%
  filter(complete.cases(.))

message("  RF training data:", nrow(rf_data), "observations\n")

# Fit random forest
set.seed(42)
rf_model <- randomForest(Cluster ~ ., data = rf_data, 
                         importance = TRUE, ntree = 500)

message("\nRandom Forest Results:\n")
print(rf_model)

# Variable importance
var_imp <- importance(rf_model) %>%
  as.data.frame() %>%
  rownames_to_column("Predictor") %>%
  arrange(desc(MeanDecreaseGini))

message("\nVariable Importance (Mean Decrease Gini):\n")
print(var_imp %>% select(Predictor, MeanDecreaseGini) %>% 
        mutate(MeanDecreaseGini = round(MeanDecreaseGini, 2)))

# Confusion matrix
message("\nConfusion Matrix:\n")
print(rf_model$confusion)

# OOB error rate
oob_error <- rf_model$err.rate[nrow(rf_model$err.rate), "OOB"]
message("\nOOB Error Rate:", round(oob_error * 100, 1), "%\n")
message("Classification Accuracy:", round((1 - oob_error) * 100, 1), "%\n")

# =============================================================================
# 5. GEOLOGY-ONLY MODEL
# =============================================================================

message("\n=== 5. GEOLOGY-ONLY MODEL ===\n\n")

geo_available <- intersect(geo_preds, available_preds)

if (length(geo_available) >= 2) {
  geo_formula <- as.formula(paste("Cluster ~", paste(geo_available, collapse = " + ")))
  model_geo <- nnet::multinom(geo_formula, data = model_data, trace = FALSE)
  
  pseudo_r2_geo <- 1 - (logLik(model_geo) / logLik(null_model))
  message("Geology-only Pseudo R²:", round(as.numeric(pseudo_r2_geo), 3), "\n")
  
  # Compare models
  message("\nModel Comparison:\n")
  message("  Full model R²:     ", round(as.numeric(pseudo_r2), 3), "\n")
  message("  Geology-only R²:   ", round(as.numeric(pseudo_r2_geo), 3), "\n")
  message("  Geology explains:  ", 
      round(as.numeric(pseudo_r2_geo) / as.numeric(pseudo_r2) * 100, 1), "% of full model\n")
}

# =============================================================================
# 6. CLUSTER STABILITY PREDICTORS
# =============================================================================

message("\n=== 6. CLUSTER STABILITY ANALYSIS ===\n\n")

# What predicts cluster stability?
stability_data <- model_data %>%
  filter(!is.na(stability))

if (nrow(stability_data) > 10) {
  stab_formula <- as.formula(paste("stability ~", paste(available_preds, collapse = " + ")))
  stab_model <- lm(stab_formula, data = stability_data)
  
  message("Stability Model Summary:\n")
  print(summary(stab_model))
  
  # Correlation with stability
  stab_cors <- map_dfr(available_preds, function(pred) {
    test <- cor.test(stability_data[[pred]], stability_data$stability, 
                     use = "complete.obs")
    tibble(Predictor = pred, r = test$estimate, p = test$p.value)
  }) %>%
    arrange(p)
  
  message("\nStability Correlations:\n")
  print(stab_cors %>% mutate(r = round(r, 3), p = round(p, 3)))
}

# =============================================================================
# 7. VISUALIZATIONS
# =============================================================================

message("\n=== 7. CREATING VISUALIZATIONS ===\n\n")

# Variable importance plot
p1 <- var_imp %>%
  mutate(Predictor = fct_reorder(Predictor, MeanDecreaseGini)) %>%
  ggplot(aes(x = MeanDecreaseGini, y = Predictor)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  labs(
    title = "Random Forest: Variable Importance for Cluster Prediction",
    subtitle = paste0("OOB Accuracy: ", round((1 - oob_error) * 100, 1), "%"),
    x = "Mean Decrease in Gini Impurity",
    y = NULL
  ) +
  theme_clean()

ggsave(file.path(fig_dir, "01_rf_variable_importance.png"), p1, 
       width = 9, height = 6, dpi = 300)

# Cluster by geology biplot
if (all(c("Lava1_per", "Ash_Per") %in% names(model_data))) {
  p2 <- ggplot(model_data, aes(x = Lava1_per, y = Ash_Per, color = Cluster)) +
    geom_point(size = 3, alpha = 0.7) +
    stat_ellipse(level = 0.67, linetype = "dashed") +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = "Cluster Membership by Geology",
      subtitle = "67% concentration ellipses",
      x = "Lava1 (%)",
      y = "Ash (%)"
    ) +
    theme_clean() +
    theme(legend.position = "right")
  
  ggsave(file.path(fig_dir, "02_cluster_by_geology.png"), p2, 
         width = 10, height = 8, dpi = 300)
}

# Cluster distribution by site
p3 <- model_data %>%
  count(Stream_Name, Cluster) %>%
  ggplot(aes(x = Stream_Name, y = n, fill = Cluster)) +
  geom_col(position = "fill") +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Cluster Composition by Site",
    subtitle = "Proportion of solute-windows in each cluster",
    x = NULL,
    y = "Proportion"
  ) +
  theme_clean() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(fig_dir, "03_cluster_composition_by_site.png"), p3, 
       width = 10, height = 6, dpi = 300)

# Stability by cluster
if ("stability" %in% names(model_data)) {
  p4 <- ggplot(model_data %>% filter(!is.na(stability)), 
               aes(x = Cluster, y = stability, fill = Cluster)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
    scale_fill_brewer(palette = "Set1") +
    labs(
      title = "Cluster Stability by Cluster Type",
      subtitle = "Do some CQ behaviors persist more than others?",
      x = NULL,
      y = "Stability (proportion of years in modal cluster)"
    ) +
    theme_clean() +
    theme(legend.position = "none")
  
  ggsave(file.path(fig_dir, "04_stability_by_cluster.png"), p4, 
         width = 8, height = 6, dpi = 300)
}

# =============================================================================
# 8. SUMMARY
# =============================================================================

message("\n╔════════════════════════════════════════════════════════════════╗\n")
message("║  CLUSTER MEMBERSHIP ANALYSIS COMPLETE                         ║\n")
message("╚════════════════════════════════════════════════════════════════╝\n\n")

message("KEY FINDINGS:\n\n")

message("1. Random Forest Classification:\n")
message("   - OOB Accuracy:", round((1 - oob_error) * 100, 1), "%\n")
message("   - Top predictor:", var_imp$Variable[1], "\n")
message("   - Geology importance:", 
    round(sum(var_imp$MeanDecreaseGini[var_imp$Variable %in% geo_preds]) / 
            sum(var_imp$MeanDecreaseGini) * 100, 1), "%\n\n")

message("2. Multinomial Logistic Regression:\n")
message("   - Full model Pseudo R²:", round(as.numeric(pseudo_r2), 3), "\n")
if (exists("pseudo_r2_geo")) {
  message("   - Geology-only Pseudo R²:", round(as.numeric(pseudo_r2_geo), 3), "\n")
}

message("\n3. Interpretation:\n")
if (as.numeric(pseudo_r2) > 0.1) {
  message("   - Catchment properties DO predict cluster membership\n")
  message("   - Geology explains meaningful variance in CQ behavior type\n")
} else {
  message("   - Catchment properties explain LITTLE about cluster membership\n")
  message("   - CQ behavior type may be more dynamic/temporal\n")
}

message("\nOutputs saved to:\n")
message("  Figures:", fig_dir, "\n\n")

# Save results
results <- list(
  var_importance = var_imp,
  pseudo_r2_full = as.numeric(pseudo_r2),
  oob_accuracy = 1 - oob_error,
  coef_table = coef_table
)

saveRDS(results, file.path(res_dir, "cluster_membership_results.rds"))
write_csv(var_imp, file.path(res_dir, "cluster_variable_importance.csv"))
