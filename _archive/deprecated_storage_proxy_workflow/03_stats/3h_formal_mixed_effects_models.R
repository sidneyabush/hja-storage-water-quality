# =============================================================================
# STEP 03h: FORMAL MIXED-EFFECTS MODELS
# =============================================================================
# Purpose: Test storage-synchrony relationships using proper statistical models
#          that account for repeated measures and hierarchical data structure
#
# Key insight from 3b: Overall R² ~ 0.10, but within-cluster R² ~ 0.70
# This script extends that finding with:
#   1. Mixed-effects models accounting for repeated measures
#   2. Model comparison (AIC/BIC)
#   3. Fixed effects tables with β, SE, t, p
#   4. Random effects variance components
#   5. Marginal and conditional R²
#
# Research Questions Addressed:
#   RQ1: Does dynamic storage explain shared chemical behavior?
#   RQ2: Does storage variability weaken synchrony?
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(lmerTest)  # For p-values in lmer
  library(MuMIn)     # For R.squaredGLMM
  library(broom.mixed)  # For tidy model outputs
  library(patchwork)
})

rm(list = ls())

# Source helpers
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
try(source(file.path(repo_dir, "00_helpers", "plot_prefs.R")), silent = TRUE)
try(source(file.path(repo_dir, "00_helpers", "workflow_config.R")), silent = TRUE)

# Paths
paths <- get_project_paths()
out_dir <- paths$out_dir
fig_dir <- file.path(paths$fig_root, "03_stats", "3h_mixed_effects")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "03_stats"), showWarnings = FALSE, recursive = TRUE)

message("\n╔════════════════════════════════════════════════════════════════╗\n")
message("║  FORMAL MIXED-EFFECTS MODELS: Storage → Synchrony             ║\n")
message("╚════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# LOAD DATA
# =============================================================================
message("=== LOADING DATA ===\n")

# Annual data with synchrony metrics
annual <- read_csv(file.path(out_dir, "HJA_master_annual.csv"), show_col_types = FALSE)
site_means <- read_csv(file.path(out_dir, "HJA_master_site_means.csv"), show_col_types = FALSE)
sync_annual <- read_csv(file.path(out_dir, "HJA_composite_synchrony_annual.csv"), show_col_types = FALSE)
clusters_modal <- read_csv(file.path(out_dir, "ClusterStreams_allSolutes_modalClusters.csv"), show_col_types = FALSE)

# Join data
annual_full <- annual %>%
  left_join(sync_annual, by = c("Stream_Name", "solute", "water_year")) %>%
  left_join(
    clusters_modal %>% select(Stream_Name, chemical, Cluster_mode),
    by = c("Stream_Name", "solute" = "chemical")
  ) %>%
  mutate(
    Stream_Name = as.factor(Stream_Name),
    solute = as.factor(solute),
    solute_type = factor(categorize_solute(solute), levels = c("Geogenic", "Biogenic", "Nutrient")),
    Cluster = factor(Cluster_mode, levels = cluster_levels)
  )

# Check storage metric
storage_metric <- intersect(PRIMARY_STORAGE_METRIC, names(annual_full))
if (length(storage_metric) == 0) {
  storage_metric <- intersect(c("Q_dS_range_mm", "WB_dS_range_mm"), names(annual_full))
  if (length(storage_metric) == 0) stop("No storage metric found!")
}
storage_metric <- storage_metric[[1]]
message("  Using storage metric:", storage_metric, "\n")

# =============================================================================
# COMPUTE STORAGE DIVERGENCE (Inter-site SD)
# =============================================================================
message("\n=== COMPUTING STORAGE DIVERGENCE ===\n")

# Annual storage divergence (SD across sites within each year-solute)
annual_divergence <- annual_full %>%
  group_by(water_year, solute) %>%
  summarise(
    storage_sd = sd(.data[[storage_metric]], na.rm = TRUE),
    storage_mean = mean(.data[[storage_metric]], na.rm = TRUE),
    n_sites = sum(!is.na(.data[[storage_metric]])),
    .groups = "drop"
  ) %>%
  filter(n_sites >= 3, is.finite(storage_sd))  # Need at least 3 sites

# Join back to annual data with synchrony
model_data <- annual_full %>%
  left_join(annual_divergence, by = c("water_year", "solute")) %>%
  filter(!is.na(storage_sd), !is.na(Cluster)) %>%
  # Center predictors for easier interpretation
  mutate(
    storage_sd_c = as.numeric(scale(storage_sd, center = TRUE, scale = FALSE)),
    storage_mean_c = as.numeric(scale(storage_mean, center = TRUE, scale = FALSE))
  )

message("  Model dataset:", nrow(model_data), "observations\n")
message("  Sites:", n_distinct(model_data$Stream_Name), "\n")
message("  Solutes:", n_distinct(model_data$solute), "\n")
message("  Years:", n_distinct(model_data$water_year), "\n")
message("  Clusters:", paste(levels(droplevels(model_data$Cluster)), collapse = ", "), "\n")

# =============================================================================
# MODEL 1: CONCENTRATION SYNCHRONY ~ STORAGE DIVERGENCE
# =============================================================================
message("\n=== MODEL 1: CONCENTRATION SYNCHRONY ===\n")

if ("conc_sync_allpairs" %in% names(model_data)) {

  # Prepare data (drop NAs)
  df_conc <- model_data %>%
    filter(!is.na(conc_sync_allpairs), !is.na(storage_sd_c)) %>%
    droplevels()

  message("  N =", nrow(df_conc), "observations\n")

  # Model 1a: Null model (random effects only)
  m1a_null <- lmer(
    conc_sync_allpairs ~ 1 + (1|water_year) + (1|solute) + (1|Stream_Name),
    data = df_conc,
    REML = TRUE
  )

  # Model 1b: Storage divergence only
  m1b_storage <- lmer(
    conc_sync_allpairs ~ storage_sd_c + (1|water_year) + (1|solute) + (1|Stream_Name),
    data = df_conc,
    REML = TRUE
  )

  # Model 1c: Storage divergence + Cluster
  m1c_cluster <- lmer(
    conc_sync_allpairs ~ storage_sd_c + Cluster + (1|water_year) + (1|solute) + (1|Stream_Name),
    data = df_conc,
    REML = TRUE
  )

  # Model 1d: Storage divergence × Cluster interaction
  m1d_interaction <- lmer(
    conc_sync_allpairs ~ storage_sd_c * Cluster + (1|water_year) + (1|solute) + (1|Stream_Name),
    data = df_conc,
    REML = TRUE
  )

  # Model comparison
  anova_conc <- anova(m1a_null, m1b_storage, m1c_cluster, m1d_interaction)
  message("\nModel comparison (Concentration Synchrony):\n")
  print(anova_conc)

  # Best model (lowest AIC)
  best_conc <- list(m1a_null, m1b_storage, m1c_cluster, m1d_interaction)[[which.min(anova_conc$AIC)]]
  message("\nBest model: Model", which.min(anova_conc$AIC), "\n")

  # Fixed effects table
  fixed_conc <- tidy(best_conc, effects = "fixed") %>%
    mutate(across(where(is.numeric), ~round(., 4)))
  message("\nFixed effects:\n")
  print(fixed_conc)

  # Random effects
  random_conc <- tidy(best_conc, effects = "ran_pars") %>%
    mutate(across(where(is.numeric), ~round(., 4)))
  message("\nRandom effects:\n")
  print(random_conc)

  # R-squared
  r2_conc <- r.squaredGLMM(best_conc)
  message("\nR-squared:\n")
  message("  Marginal R² (fixed effects only):", round(r2_conc[1, "R2m"], 3), "\n")
  message("  Conditional R² (fixed + random):", round(r2_conc[1, "R2c"], 3), "\n")

  # Save results
  write_csv(fixed_conc, file.path(out_dir, "03_stats/model1_conc_sync_fixed_effects.csv"))
  write_csv(random_conc, file.path(out_dir, "03_stats/model1_conc_sync_random_effects.csv"))
  write_csv(as_tibble(anova_conc, rownames = "model"),
            file.path(out_dir, "03_stats/model1_conc_sync_comparison.csv"))

  # Diagnostic plots
  df_conc$fitted <- fitted(best_conc)
  df_conc$resid <- residuals(best_conc)

  p_diag1 <- ggplot(df_conc, aes(x = fitted, y = resid)) +
    geom_point(alpha = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(method = "loess", se = FALSE, color = "blue") +
    labs(x = "Fitted values", y = "Residuals",
         title = "Model 1: Concentration Synchrony Diagnostics",
         subtitle = "Residuals vs Fitted") +
    theme_hja()

  p_diag2 <- ggplot(df_conc, aes(sample = resid)) +
    stat_qq() +
    stat_qq_line(color = "red") +
    labs(title = "Q-Q Plot", x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_hja()

  p_diag <- p_diag1 + p_diag2
  ggsave(file.path(fig_dir, "model1_conc_sync_diagnostics.png"), p_diag, width = 14, height = 7, dpi = 300)

  message("\n✓ Model 1 complete: Concentration synchrony\n")
}

# =============================================================================
# MODEL 2: CQ-SLOPE SYNCHRONY ~ STORAGE DIVERGENCE
# =============================================================================
message("\n=== MODEL 2: CQ-SLOPE SYNCHRONY ===\n")

if ("cqslope_sync_allpairs" %in% names(model_data)) {

  # Prepare data
  df_cqslope <- model_data %>%
    filter(!is.na(cqslope_sync_allpairs), !is.na(storage_sd_c)) %>%
    droplevels()

  message("  N =", nrow(df_cqslope), "observations\n")

  # Model 2a: Null model
  m2a_null <- lmer(
    cqslope_sync_allpairs ~ 1 + (1|water_year) + (1|solute) + (1|Stream_Name),
    data = df_cqslope,
    REML = TRUE
  )

  # Model 2b: Storage divergence only
  m2b_storage <- lmer(
    cqslope_sync_allpairs ~ storage_sd_c + (1|water_year) + (1|solute) + (1|Stream_Name),
    data = df_cqslope,
    REML = TRUE
  )

  # Model 2c: Storage divergence + Cluster
  m2c_cluster <- lmer(
    cqslope_sync_allpairs ~ storage_sd_c + Cluster + (1|water_year) + (1|solute) + (1|Stream_Name),
    data = df_cqslope,
    REML = TRUE
  )

  # Model 2d: Storage divergence × Cluster interaction
  m2d_interaction <- lmer(
    cqslope_sync_allpairs ~ storage_sd_c * Cluster + (1|water_year) + (1|solute) + (1|Stream_Name),
    data = df_cqslope,
    REML = TRUE
  )

  # Model comparison
  anova_cqslope <- anova(m2a_null, m2b_storage, m2c_cluster, m2d_interaction)
  message("\nModel comparison (CQ-Slope Synchrony):\n")
  print(anova_cqslope)

  # Best model
  best_cqslope <- list(m2a_null, m2b_storage, m2c_cluster, m2d_interaction)[[which.min(anova_cqslope$AIC)]]
  message("\nBest model: Model", which.min(anova_cqslope$AIC), "\n")

  # Fixed effects table
  fixed_cqslope <- tidy(best_cqslope, effects = "fixed") %>%
    mutate(across(where(is.numeric), ~round(., 4)))
  message("\nFixed effects:\n")
  print(fixed_cqslope)

  # Random effects
  random_cqslope <- tidy(best_cqslope, effects = "ran_pars") %>%
    mutate(across(where(is.numeric), ~round(., 4)))
  message("\nRandom effects:\n")
  print(random_cqslope)

  # R-squared
  r2_cqslope <- r.squaredGLMM(best_cqslope)
  message("\nR-squared:\n")
  message("  Marginal R² (fixed effects only):", round(r2_cqslope[1, "R2m"], 3), "\n")
  message("  Conditional R² (fixed + random):", round(r2_cqslope[1, "R2c"], 3), "\n")

  # Save results
  write_csv(fixed_cqslope, file.path(out_dir, "03_stats/model2_cqslope_sync_fixed_effects.csv"))
  write_csv(random_cqslope, file.path(out_dir, "03_stats/model2_cqslope_sync_random_effects.csv"))
  write_csv(as_tibble(anova_cqslope, rownames = "model"),
            file.path(out_dir, "03_stats/model2_cqslope_sync_comparison.csv"))

  # Diagnostics
  df_cqslope$fitted <- fitted(best_cqslope)
  df_cqslope$resid <- residuals(best_cqslope)

  p_diag1 <- ggplot(df_cqslope, aes(x = fitted, y = resid)) +
    geom_point(alpha = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(method = "loess", se = FALSE, color = "blue") +
    labs(x = "Fitted values", y = "Residuals",
         title = "Model 2: CQ-Slope Synchrony Diagnostics",
         subtitle = "Residuals vs Fitted") +
    theme_hja()

  p_diag2 <- ggplot(df_cqslope, aes(sample = resid)) +
    stat_qq() +
    stat_qq_line(color = "red") +
    labs(title = "Q-Q Plot", x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_hja()

  p_diag <- p_diag1 + p_diag2
  ggsave(file.path(fig_dir, "model2_cqslope_sync_diagnostics.png"), p_diag, width = 14, height = 7, dpi = 300)

  message("\n✓ Model 2 complete: CQ-slope synchrony\n")
}

# =============================================================================
# MODEL 3: STORAGE DIVERGENCE × SOLUTE TYPE
# =============================================================================
message("\n=== MODEL 3: SOLUTE TYPE EFFECTS ===\n")

if ("conc_sync_allpairs" %in% names(model_data) && "solute_type" %in% names(model_data)) {

  df_soltype <- model_data %>%
    filter(!is.na(conc_sync_allpairs), !is.na(storage_sd_c), !is.na(solute_type)) %>%
    droplevels()

  message("  N =", nrow(df_soltype), "observations\n")

  # Model 3: Storage × Solute type
  m3_soltype <- lmer(
    conc_sync_allpairs ~ storage_sd_c * solute_type + (1|water_year) + (1|solute) + (1|Stream_Name),
    data = df_soltype,
    REML = TRUE
  )

  # Fixed effects
  fixed_soltype <- tidy(m3_soltype, effects = "fixed") %>%
    mutate(across(where(is.numeric), ~round(., 4)))
  message("\nFixed effects (Storage × Solute Type):\n")
  print(fixed_soltype)

  # R-squared
  r2_soltype <- r.squaredGLMM(m3_soltype)
  message("\nR-squared:\n")
  message("  Marginal R² (fixed effects only):", round(r2_soltype[1, "R2m"], 3), "\n")
  message("  Conditional R² (fixed + random):", round(r2_soltype[1, "R2c"], 3), "\n")

  # Save results
  write_csv(fixed_soltype, file.path(out_dir, "03_stats/model3_solute_type_fixed_effects.csv"))

  message("\n✓ Model 3 complete: Solute type effects\n")
}

# =============================================================================
# SUMMARY TABLE: ALL MODELS
# =============================================================================
message("\n=== CREATING SUMMARY TABLE ===\n")

# Compile model summaries
model_summary <- tibble(
  Model = c("M1: Conc Sync (best)", "M2: CQ-Slope Sync (best)", "M3: Solute Type"),
  Response = c("conc_sync_allpairs", "cqslope_sync_allpairs", "conc_sync_allpairs"),
  N = c(
    if(exists("df_conc")) nrow(df_conc) else NA,
    if(exists("df_cqslope")) nrow(df_cqslope) else NA,
    if(exists("df_soltype")) nrow(df_soltype) else NA
  ),
  AIC = c(
    if(exists("best_conc")) AIC(best_conc) else NA,
    if(exists("best_cqslope")) AIC(best_cqslope) else NA,
    if(exists("m3_soltype")) AIC(m3_soltype) else NA
  ),
  R2_marginal = c(
    if(exists("r2_conc")) r2_conc[1, "R2m"] else NA,
    if(exists("r2_cqslope")) r2_cqslope[1, "R2m"] else NA,
    if(exists("r2_soltype")) r2_soltype[1, "R2m"] else NA
  ),
  R2_conditional = c(
    if(exists("r2_conc")) r2_conc[1, "R2c"] else NA,
    if(exists("r2_cqslope")) r2_cqslope[1, "R2c"] else NA,
    if(exists("r2_soltype")) r2_soltype[1, "R2c"] else NA
  )
) %>%
  mutate(across(where(is.numeric), ~round(., 3)))

print(model_summary)
write_csv(model_summary, file.path(out_dir, "03_stats/model_summary_all.csv"))

# =============================================================================
# INTERPRETATION GUIDE
# =============================================================================
message("\n╔════════════════════════════════════════════════════════════════╗\n")
message("║  INTERPRETATION GUIDE                                         ║\n")
message("╚════════════════════════════════════════════════════════════════╝\n\n")

message("Marginal R² (R2m): Variance explained by fixed effects only\n")
message("  → How much does storage divergence (and cluster) explain?\n\n")

message("Conditional R² (R2c): Variance explained by fixed + random effects\n")
message("  → How much total variance is explained by the model?\n\n")

message("Difference (R2c - R2m): Variance explained by random effects\n")
message("  → How important are year-to-year, solute-specific, and site-specific variations?\n\n")

message("Model comparison (AIC):\n")
message("  → Lower AIC = better model fit\n")
message("  → ΔAIC > 2 indicates meaningful improvement\n\n")

message("Fixed effects table:\n")
message("  β (estimate): Effect size (slope)\n")
message("  SE (std.error): Standard error of the estimate\n")
message("  t-value: Test statistic\n")
message("  p-value: Significance (p < 0.05 = significant)\n\n")

message("╔════════════════════════════════════════════════════════════════╗\n")
message("║  ANALYSIS COMPLETE                                            ║\n")
message("╚════════════════════════════════════════════════════════════════╝\n\n")

message("Outputs saved to:\n")
message("  Tables: ", file.path(out_dir, "03_stats/"), "\n")
message("  Figures:", fig_dir, "\n\n")

message("Key files:\n")
message("  - model_summary_all.csv (all model comparison)\n")
message("  - model1_conc_sync_fixed_effects.csv\n")
message("  - model2_cqslope_sync_fixed_effects.csv\n")
message("  - model3_solute_type_fixed_effects.csv\n\n")
