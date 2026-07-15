# =============================================================================
# STEP 03l: OUTLET SYNCHRONY PREDICTION MODELS
# =============================================================================
# Goal: Identify catchment characteristics that predict synchrony with outlet
#
# RESPONSE VARIABLES:
#   - conc_sync_outlet: Abbott concentration synchrony with GSLOOK
#   - cqslope_sync_outlet: Wymore CQ-slope synchrony with GSLOOK
#
# PREDICTOR TIERS:
#   Tier 1 (all 8 sites): RBI, Q_dS_range_mm, Elevation, Slope, Pyro_per, Lava1_per
#                         + wet_length_days, wet_start_doy (season timing)
#   Tier 2 (isotope subset ~5-6 sites): + MTT_final, FYw_final, DR_Overall
#
# METHODS:
#   1. Univariate screening (correlations)
#   2. Variance partitioning by predictor group
#   3. A priori hypothesis models (not blind dredging)
#   4. AIC comparison and model averaging
#   5. Leave-one-site-out cross-validation
#   6. Cluster-CQ slope relationship testing
#
# A PRIORI HYPOTHESES:
#   H1: "Flashiness dominates" → sync ~ RBI + Elevation
#   H2: "Storage matters" → sync ~ Q_dS + RBI
#   H3: "Geology controls" → sync ~ Pyro_per + Lava1_per
#   H4: "Clusters capture it" → sync ~ pct_cluster3 + stability
#   H5: "Transit time explains" → sync ~ MTT_final + FYw_final (Tier 2 only)
#   H6: "Combined hydro" → sync ~ RBI + Q_dS + Elevation
#   H7: "Season timing" → sync ~ wet_length_days + wet_start_doy (NEW Dec 2025)
#   H8: "Season + Flashiness" → sync ~ wet_length + wet_start + RBI (NEW Dec 2025)
#
# NOTE ON CQ SLOPES: Including mean CQ slope as predictor would be circular
#       since Wymore sync is derived from CQ slopes. We test cluster-CQ
#       relationship separately but do not include CQ slope in sync models.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(MuMIn)
  library(vegan)      # For variance partitioning
  library(broom)
  library(patchwork)
})

rm(list = ls())

# Source helpers
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
source(file.path(repo_dir, "00_helpers", "workflow_config.R"))
try(source(file.path(repo_dir, "00_helpers", "plot_prefs.R")), silent = TRUE)

# Paths
paths <- get_project_paths()
out_dir <- paths$out_dir
fig_dir <- file.path(paths$fig_root, "03_stats", "3l_outlet_sync")
res_dir <- file.path(out_dir, "03_stats")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

message("\n╔════════════════════════════════════════════════════════════════╗\n")
message("║  OUTLET SYNCHRONY PREDICTION MODELS                            ║\n")
message("╚════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# 1. LOAD DATA
# =============================================================================
message("=== 1. LOADING DATA ===\n\n")

# Data directory (raw data files)
data_dir <- file.path(dirname(out_dir), "data")

# Outlet synchrony metrics (created in 1f)
outlet_sync_annual <- read_csv(file.path(out_dir, "HJA_outlet_synchrony_annual.csv"), 
                                show_col_types = FALSE)
outlet_sync_site <- read_csv(file.path(out_dir, "HJA_outlet_synchrony_site_level.csv"), 
                              show_col_types = FALSE)

# Static site characteristics (one row per site) - from catchment characteristics in data folder
static_chars <- read_csv(file.path(data_dir, "Catchment_Charc.csv"), 
                         show_col_types = FALSE) %>%
  rename(Stream_Name = Site) %>%
  select(Stream_Name, any_of(c("Area_km2", "Elevation_mean_m", "Slope_mean", 
                                "Lava1_per", "Lava2_per", "Ash_Per", "Pyro_per")))

# Isotope data (separate file in data folder)
isotope_data <- tryCatch(
  read_csv(file.path(data_dir, "MTT_FYW.csv"), 
           show_col_types = FALSE) %>%
    rename(Stream_Name = site) %>%
    select(Stream_Name, MTTM, FYWM),
  error = function(e) {
    message("  Note: MTT_FYW.csv not found, skipping isotope data\n")
    NULL
  }
)

damping_data <- tryCatch(
  read_csv(file.path(data_dir, "DampingRatios_2025-07-07.csv"), 
           show_col_types = FALSE) %>%
    rename(Stream_Name = site) %>%
    select(Stream_Name, DR_Overall),
  error = function(e) {
    message("  Note: DampingRatios file not found, skipping damping data\n")
    NULL
  }
)

# Hydrologic metrics (need to average across solutes since they vary slightly)
hydro_metrics <- read_csv(file.path(out_dir, "HJA_master_site_means.csv"), 
                          show_col_types = FALSE) %>%
  group_by(Stream_Name) %>%
  summarise(
    RBI = mean(RBI, na.rm = TRUE),
    Q_dS_range_mm = mean(Q_dS_range_mm, na.rm = TRUE),
    RCS_p = mean(RCS_p, na.rm = TRUE),
    RCS_k = mean(RCS_k, na.rm = TRUE),
    FDC_slope_5_95 = mean(FDC_slope_5_95, na.rm = TRUE),
    .groups = "drop"
  )

# Combine site characteristics
site_chars <- static_chars %>%
  left_join(hydro_metrics, by = "Stream_Name")

if (!is.null(isotope_data)) {
  site_chars <- site_chars %>%
    left_join(isotope_data, by = "Stream_Name")
}

if (!is.null(damping_data)) {
  site_chars <- site_chars %>%
    left_join(damping_data, by = "Stream_Name")
}

# Cluster data
clusters_modal <- read_csv(file.path(out_dir, "ClusterStreams_allSolutes_modalClusters.csv"), 
                           show_col_types = FALSE)
clusters_stability <- read_csv(file.path(out_dir, "ClusterStreams_stability_metrics.csv"), 
                               show_col_types = FALSE)

# Season timing data (annual-varying predictors)
season_bounds <- read_csv(file.path(out_dir, "season_boundaries.csv"), 
                          show_col_types = FALSE) %>%
  mutate(
    wet_length_days = as.numeric(wet_end_date - wet_start_date),
    wet_start_doy = lubridate::yday(wet_start_date)
  ) %>%
  select(water_year, wet_length_days, wet_start_doy)

message("  Season timing data loaded:", nrow(season_bounds), "water years\n")

# Annual CQ data for cluster proportions
annual_data <- read_csv(file.path(out_dir, "HJA_master_annual.csv"), 
                        show_col_types = FALSE)

message("  Outlet sync annual rows:", nrow(outlet_sync_annual), "\n")
message("  Outlet sync site-level rows:", nrow(outlet_sync_site), "\n")
message("  Site characteristics rows:", nrow(site_chars), "\n")
message("  Unique sites:", n_distinct(outlet_sync_site$Stream_Name), "\n\n")

# =============================================================================
# 2. BUILD ANALYSIS DATASET
# =============================================================================
message("=== 2. BUILDING ANALYSIS DATASET ===\n\n")

# Calculate cluster proportions per site-year
cluster_props <- annual_data %>%
  filter(!is.na(Cluster_wy)) %>%
  group_by(Stream_Name, water_year) %>%
  summarise(
    n_solutes = n(),
    pct_cluster1 = mean(Cluster_wy == 1, na.rm = TRUE),
    pct_cluster2 = mean(Cluster_wy == 2, na.rm = TRUE),
    pct_cluster3 = mean(Cluster_wy == 3, na.rm = TRUE),
    pct_cluster4 = mean(Cluster_wy == 4, na.rm = TRUE),
    .groups = "drop"
  )

# Calculate cluster stability per site
cluster_stability_site <- clusters_stability %>%
  group_by(Stream_Name) %>%
  summarise(
    cluster_stability_mean = mean(stability, na.rm = TRUE),
    cluster_stability_sd = sd(stability, na.rm = TRUE),
    .groups = "drop"
  )

# Build site-annual dataset for modeling
model_data_annual <- outlet_sync_annual %>%
  left_join(site_chars, by = "Stream_Name") %>%
  left_join(cluster_props, by = c("Stream_Name", "water_year")) %>%
  left_join(cluster_stability_site, by = "Stream_Name") %>%
  left_join(season_bounds, by = "water_year") %>%  # Add season timing predictors
  filter(Stream_Name != OUTLET_SITE)  # Exclude outlet from predictors

message("  Model data (site-annual):", nrow(model_data_annual), "rows\n")
message("  Sites:", n_distinct(model_data_annual$Stream_Name), "\n")
message("  Solutes:", n_distinct(model_data_annual$solute), "\n")
message("  Years:", n_distinct(model_data_annual$water_year), "\n\n")

# Build site-level dataset
model_data_site <- outlet_sync_site %>%
  left_join(site_chars, by = "Stream_Name") %>%
  left_join(cluster_stability_site, by = "Stream_Name") %>%
  filter(Stream_Name != OUTLET_SITE)

message("  Model data (site-level):", nrow(model_data_site), "rows\n\n")

# =============================================================================
# 3. UNIVARIATE SCREENING
# =============================================================================
message("=== 3. UNIVARIATE SCREENING ===\n\n")

# Function to compute correlations with outlet sync
compute_univariate_cors <- function(df, response_var, predictors) {
  results <- map_dfr(predictors, function(pred) {
    if (!pred %in% names(df)) return(NULL)
    x <- df[[pred]]
    y <- df[[response_var]]
    idx <- is.finite(x) & is.finite(y)
    if (sum(idx) < 5) return(NULL)
    
    test <- cor.test(x[idx], y[idx])
    tibble(
      predictor = pred,
      response = response_var,
      r = test$estimate,
      p = test$p.value,
      n = sum(idx),
      sig = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", p < 0.1 ~ ".", TRUE ~ "")
    )
  })
  results %>% arrange(p)
}

# Tier 1 predictors (all sites) - now includes season timing
tier1_preds <- c("RBI", "Q_dS_range_mm", "Elevation_mean_m", "Slope_mean", 
                 "Area_km2", "Pyro_per", "Lava1_per", "Lava2_per", "Ash_Per",
                 "pct_cluster1", "pct_cluster3", "cluster_stability_mean",
                 "wet_length_days", "wet_start_doy")  # NEW: Season timing predictors

# Tier 2 predictors (isotope subset) - using correct column names
tier2_preds <- c(tier1_preds, "MTTM", "FYWM", "DR_Overall")

# Aggregate to site-year level for more power
site_year_means <- model_data_annual %>%
  group_by(Stream_Name, water_year) %>%
  summarise(
    conc_sync_outlet = mean(conc_sync_outlet, na.rm = TRUE),
    cqslope_sync_outlet = mean(cqslope_sync_outlet, na.rm = TRUE),
    across(any_of(tier2_preds), ~ mean(., na.rm = TRUE)),
    .groups = "drop"
  )

message("--- Univariate correlations (site-year level) ---\n\n")

cors_conc <- compute_univariate_cors(site_year_means, "conc_sync_outlet", tier1_preds)
cors_cqslope <- compute_univariate_cors(site_year_means, "cqslope_sync_outlet", tier1_preds)

message("Concentration sync with outlet:\n")
print(cors_conc %>% mutate(r = round(r, 3), p = round(p, 4)), n = 20)

message("\nCQ-slope sync with outlet:\n")
print(cors_cqslope %>% mutate(r = round(r, 3), p = round(p, 4)), n = 20)

# Save screening results
univar_results <- bind_rows(cors_conc, cors_cqslope)
write_csv(univar_results, file.path(res_dir, "outlet_sync_univariate_screening.csv"))

# =============================================================================
# 4. A PRIORI HYPOTHESIS MODELS
# =============================================================================
message("\n=== 4. A PRIORI HYPOTHESIS MODELS ===\n\n")

# Use site-year data with random effects for site and year
# Scale predictors for coefficient comparison
model_df <- site_year_means %>%
  mutate(across(where(is.numeric) & !matches("sync"), ~ scale(.) %>% as.vector()))

# Define hypothesis models
fit_hypothesis_models <- function(response_var, df) {
  
  # Null model
  m0 <- tryCatch(
    lmer(as.formula(paste(response_var, "~ 1 + (1|Stream_Name)")), data = df, REML = FALSE),
    error = function(e) NULL
  )
  
  # H1: Flashiness dominates
  m1 <- tryCatch(
    lmer(as.formula(paste(response_var, "~ RBI + Elevation_mean_m + (1|Stream_Name)")), 
         data = df, REML = FALSE),
    error = function(e) NULL
  )
  
  # H2: Storage matters
  m2 <- tryCatch(
    lmer(as.formula(paste(response_var, "~ Q_dS_range_mm + RBI + (1|Stream_Name)")), 
         data = df, REML = FALSE),
    error = function(e) NULL
  )
  
  # H3: Geology controls
  m3 <- tryCatch(
    lmer(as.formula(paste(response_var, "~ Pyro_per + Lava1_per + (1|Stream_Name)")), 
         data = df, REML = FALSE),
    error = function(e) NULL
  )
  
  # H4: Clusters capture it
  m4 <- tryCatch(
    lmer(as.formula(paste(response_var, "~ pct_cluster3 + cluster_stability_mean + (1|Stream_Name)")), 
         data = df, REML = FALSE),
    error = function(e) NULL
  )
  
  # H6: Combined hydro
  m6 <- tryCatch(
    lmer(as.formula(paste(response_var, "~ RBI + Q_dS_range_mm + Elevation_mean_m + (1|Stream_Name)")), 
         data = df, REML = FALSE),
    error = function(e) NULL
  )
  
  # H7: Season timing (NEW) - does interannual variation in season timing predict sync?
  m7 <- tryCatch(
    lmer(as.formula(paste(response_var, "~ wet_length_days + wet_start_doy + (1|Stream_Name)")), 
         data = df, REML = FALSE),
    error = function(e) NULL
  )
  
  # H8: Season + Flashiness combined (NEW)
  m8 <- tryCatch(
    lmer(as.formula(paste(response_var, "~ wet_length_days + wet_start_doy + RBI + (1|Stream_Name)")), 
         data = df, REML = FALSE),
    error = function(e) NULL
  )
  
  models <- list(
    "M0: Null" = m0,
    "H1: Flashiness" = m1,
    "H2: Storage" = m2,
    "H3: Geology" = m3,
    "H4: Clusters" = m4,
    "H6: Combined" = m6,
    "H7: Season timing" = m7,
    "H8: Season+Flash" = m8
  )
  
  # Filter out NULL models
  models <- models[!sapply(models, is.null)]
  
  return(models)
}

# Fit models for concentration sync
message("--- Fitting models for concentration sync ---\n")
models_conc <- fit_hypothesis_models("conc_sync_outlet", model_df)

# Fit models for CQ-slope sync
message("--- Fitting models for CQ-slope sync ---\n")
models_cqslope <- fit_hypothesis_models("cqslope_sync_outlet", model_df)

# =============================================================================
# 5. AIC COMPARISON
# =============================================================================
message("\n=== 5. MODEL COMPARISON (AIC) ===\n\n")

compare_models <- function(models, response_name) {
  if (length(models) == 0) return(NULL)
  
  comparison <- map_dfr(names(models), function(nm) {
    m <- models[[nm]]
    if (is.null(m)) return(NULL)
    
    r2 <- tryCatch(MuMIn::r.squaredGLMM(m), error = function(e) c(R2m = NA, R2c = NA))
    
    tibble(
      model = nm,
      response = response_name,
      AIC = AIC(m),
      BIC = BIC(m),
      logLik = as.numeric(logLik(m)),
      R2_marginal = r2[1],
      R2_conditional = r2[2],
      n_fixef = length(fixef(m)) - 1
    )
  }) %>%
    mutate(
      delta_AIC = AIC - min(AIC, na.rm = TRUE),
      weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC), na.rm = TRUE)
    ) %>%
    arrange(delta_AIC)
  
  return(comparison)
}

aic_conc <- compare_models(models_conc, "conc_sync_outlet")
aic_cqslope <- compare_models(models_cqslope, "cqslope_sync_outlet")

message("Concentration sync model comparison:\n")
print(aic_conc %>% mutate(across(where(is.numeric), ~ round(., 3))), n = 10)

message("\nCQ-slope sync model comparison:\n")
print(aic_cqslope %>% mutate(across(where(is.numeric), ~ round(., 3))), n = 10)

# Save AIC comparison
aic_results <- bind_rows(aic_conc, aic_cqslope)
write_csv(aic_results, file.path(res_dir, "outlet_sync_model_comparison.csv"))

# =============================================================================
# 6. BEST MODEL COEFFICIENTS
# =============================================================================
message("\n=== 6. BEST MODEL COEFFICIENTS ===\n\n")

extract_coefficients <- function(models, aic_table) {
  best_model_name <- aic_table$model[1]
  best_model <- models[[best_model_name]]
  
  if (is.null(best_model)) return(NULL)
  
  coef_df <- broom.mixed::tidy(best_model, effects = "fixed")
  
  # Handle column name variations (p.value vs std.error only)
  if ("p.value" %in% names(coef_df)) {
    coef_df <- coef_df %>%
      mutate(
        model = best_model_name,
        sig = case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**", 
                        p.value < 0.05 ~ "*", p.value < 0.1 ~ ".", TRUE ~ "")
      )
  } else {
    # Calculate approximate p-value from z = estimate/std.error
    coef_df <- coef_df %>%
      mutate(
        model = best_model_name,
        z_value = estimate / std.error,
        p_approx = 2 * pnorm(-abs(z_value)),
        sig = case_when(p_approx < 0.001 ~ "***", p_approx < 0.01 ~ "**", 
                        p_approx < 0.05 ~ "*", p_approx < 0.1 ~ ".", TRUE ~ "")
      )
  }
  
  return(coef_df)
}

coef_conc <- extract_coefficients(models_conc, aic_conc)
coef_cqslope <- extract_coefficients(models_cqslope, aic_cqslope)

if (!is.null(coef_conc)) {
  message("Best model for concentration sync:", aic_conc$model[1], "\n")
  print(coef_conc %>% mutate(across(where(is.numeric), ~ round(., 4))))
}

if (!is.null(coef_cqslope)) {
  message("\nBest model for CQ-slope sync:", aic_cqslope$model[1], "\n")
  print(coef_cqslope %>% mutate(across(where(is.numeric), ~ round(., 4))))
}

# =============================================================================
# 7. LEAVE-ONE-SITE-OUT CROSS-VALIDATION
# =============================================================================
message("\n=== 7. LEAVE-ONE-SITE-OUT CROSS-VALIDATION ===\n\n")

loocv_site <- function(models, df, response_var) {
  best_model <- models[[1]]  # Use first (best AIC) model
  if (is.null(best_model)) return(NULL)
  
  sites <- unique(df$Stream_Name)
  
  cv_results <- map_dfr(sites, function(held_out) {
    train_df <- df %>% filter(Stream_Name != held_out)
    test_df <- df %>% filter(Stream_Name == held_out)
    
    if (nrow(test_df) == 0) return(NULL)
    
    # Refit model on training data
    formula_str <- as.character(formula(best_model))
    refit <- tryCatch(
      lmer(formula(best_model), data = train_df, REML = FALSE),
      error = function(e) NULL
    )
    
    if (is.null(refit)) return(NULL)
    
    # Predict on held-out site
    pred <- predict(refit, newdata = test_df, allow.new.levels = TRUE)
    obs <- test_df[[response_var]]
    
    tibble(
      held_out_site = held_out,
      n_test = length(obs),
      mean_obs = mean(obs, na.rm = TRUE),
      mean_pred = mean(pred, na.rm = TRUE),
      rmse = sqrt(mean((obs - pred)^2, na.rm = TRUE)),
      mae = mean(abs(obs - pred), na.rm = TRUE),
      cor = cor(obs, pred, use = "complete.obs")
    )
  })
  
  return(cv_results)
}

cv_conc <- loocv_site(models_conc, model_df, "conc_sync_outlet")
cv_cqslope <- loocv_site(models_cqslope, model_df, "cqslope_sync_outlet")

if (!is.null(cv_conc)) {
  message("LOOCV for concentration sync:\n")
  print(cv_conc %>% mutate(across(where(is.numeric), ~ round(., 3))))
  message("\nOverall RMSE:", round(mean(cv_conc$rmse, na.rm = TRUE), 4), "\n")
  message("Overall MAE:", round(mean(cv_conc$mae, na.rm = TRUE), 4), "\n")
}

if (!is.null(cv_cqslope)) {
  message("\nLOOCV for CQ-slope sync:\n")
  print(cv_cqslope %>% mutate(across(where(is.numeric), ~ round(., 3))))
  message("\nOverall RMSE:", round(mean(cv_cqslope$rmse, na.rm = TRUE), 4), "\n")
  message("Overall MAE:", round(mean(cv_cqslope$mae, na.rm = TRUE), 4), "\n")
}

# Save CV results
cv_results <- bind_rows(
  cv_conc %>% mutate(response = "conc_sync_outlet"),
  cv_cqslope %>% mutate(response = "cqslope_sync_outlet")
)
write_csv(cv_results, file.path(res_dir, "outlet_sync_loocv_results.csv"))

# =============================================================================
# 8. TIER 2 MODELS (ISOTOPE SUBSET)
# =============================================================================
message("\n=== 8. TIER 2 MODELS (ISOTOPE SUBSET) ===\n\n")

# Filter to sites with isotope data (using correct column names: MTTM, FYWM, DR_Overall)
model_df_tier2 <- model_df %>%
  filter(!is.na(MTTM) & !is.na(FYWM) & !is.na(DR_Overall))

n_sites_tier2 <- n_distinct(model_df_tier2$Stream_Name)
message("Sites with isotope data:", n_sites_tier2, "\n")

if (n_sites_tier2 >= 4) {
  # H5: Transit time explains
  m5_conc <- tryCatch(
    lmer(conc_sync_outlet ~ MTTM + FYWM + (1|Stream_Name), 
         data = model_df_tier2, REML = FALSE),
    error = function(e) NULL
  )
  
  m5_cqslope <- tryCatch(
    lmer(cqslope_sync_outlet ~ MTTM + FYWM + (1|Stream_Name), 
         data = model_df_tier2, REML = FALSE),
    error = function(e) NULL
  )
  
  if (!is.null(m5_conc)) {
    message("\nH5 (Transit Time) for concentration sync:\n")
    r2_tier2 <- MuMIn::r.squaredGLMM(m5_conc)
    message("  R² marginal:", round(r2_tier2[1], 3), "\n")
    message("  R² conditional:", round(r2_tier2[2], 3), "\n")
    print(broom.mixed::tidy(m5_conc, effects = "fixed") %>% 
            mutate(across(where(is.numeric), ~ round(., 4))))
  }
  
  if (!is.null(m5_cqslope)) {
    message("\nH5 (Transit Time) for CQ-slope sync:\n")
    r2_tier2 <- MuMIn::r.squaredGLMM(m5_cqslope)
    message("  R² marginal:", round(r2_tier2[1], 3), "\n")
    message("  R² conditional:", round(r2_tier2[2], 3), "\n")
    print(broom.mixed::tidy(m5_cqslope, effects = "fixed") %>% 
            mutate(across(where(is.numeric), ~ round(., 4))))
  }
} else {
  message("Insufficient sites with isotope data for Tier 2 analysis\n")
}

# =============================================================================
# 9. CLUSTER-CQ SLOPE RELATIONSHIP
# =============================================================================
message("\n=== 9. CLUSTER-CQ SLOPE RELATIONSHIP ===\n\n")

# Test: Do clusters differ in mean CQ slope?
cq_cluster <- annual_data %>%
  filter(!is.na(Cluster_wy), !is.na(cq_slope)) %>%
  select(Stream_Name, solute, water_year, Cluster_wy, cq_slope)

message("Testing: Do the 4 modal clusters differ in mean CQ slope?\n\n")

cluster_cq_summary <- cq_cluster %>%
  group_by(Cluster_wy) %>%
  summarise(
    n = n(),
    mean_cq_slope = mean(cq_slope, na.rm = TRUE),
    sd_cq_slope = sd(cq_slope, na.rm = TRUE),
    pct_positive = mean(cq_slope > 0, na.rm = TRUE) * 100,
    .groups = "drop"
  )

print(cluster_cq_summary %>% mutate(across(where(is.numeric), ~ round(., 3))))

# ANOVA test
anova_result <- aov(cq_slope ~ factor(Cluster_wy), data = cq_cluster)
message("\nANOVA: CQ slope ~ Cluster\n")
print(summary(anova_result))

# Post-hoc Tukey
if (summary(anova_result)[[1]][["Pr(>F)"]][1] < 0.05) {
  message("\nTukey HSD post-hoc:\n")
  tukey_result <- TukeyHSD(anova_result)
  print(tukey_result)
}

# Save cluster-CQ results
write_csv(cluster_cq_summary, file.path(res_dir, "cluster_cq_slope_summary.csv"))

# =============================================================================
# 10. VARIANCE PARTITIONING
# =============================================================================
message("\n=== 10. VARIANCE PARTITIONING ===\n\n")

# Prepare site-level means for variance partitioning
site_means <- model_df %>%
  group_by(Stream_Name) %>%
  summarise(
    conc_sync_outlet = mean(conc_sync_outlet, na.rm = TRUE),
    cqslope_sync_outlet = mean(cqslope_sync_outlet, na.rm = TRUE),
    across(any_of(tier1_preds), ~ mean(., na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  filter(complete.cases(.))

if (nrow(site_means) >= 5) {
  message("Variance partitioning with", nrow(site_means), "sites\n\n")
  
  # Define predictor groups
  hydro_vars <- c("RBI", "Q_dS_range_mm")
  topo_vars <- c("Elevation_mean_m", "Slope_mean")
  geo_vars <- c("Pyro_per", "Lava1_per")
  cluster_vars <- c("pct_cluster3", "cluster_stability_mean")
  
  # Check which variables are available
  hydro_avail <- intersect(hydro_vars, names(site_means))
  topo_avail <- intersect(topo_vars, names(site_means))
  geo_avail <- intersect(geo_vars, names(site_means))
  cluster_avail <- intersect(cluster_vars, names(site_means))
  
  if (length(hydro_avail) > 0 && length(topo_avail) > 0 && length(geo_avail) > 0) {
    # Variance partitioning for concentration sync
    Y <- site_means$conc_sync_outlet
    X_hydro <- site_means[, hydro_avail, drop = FALSE]
    X_topo <- site_means[, topo_avail, drop = FALSE]
    X_geo <- site_means[, geo_avail, drop = FALSE]
    
    vp <- tryCatch(
      varpart(Y, X_hydro, X_topo, X_geo),
      error = function(e) NULL
    )
    
    if (!is.null(vp)) {
      message("Variance partitioning for concentration sync:\n")
      print(vp)
      
      # Save variance partitioning plot
      png(file.path(fig_dir, "variance_partitioning_conc_sync.png"), width = 800, height = 600)
      plot(vp, digits = 2, Xnames = c("Hydrology", "Topography", "Geology"))
      dev.off()
    }
  }
} else {
  message("Insufficient complete cases for variance partitioning\n")
}

# =============================================================================
# 11. SUMMARY INTERPRETATION
# =============================================================================
message("\n╔════════════════════════════════════════════════════════════════╗\n")
message("║                    SUMMARY INTERPRETATION                      ║\n")
message("╚════════════════════════════════════════════════════════════════╝\n\n")

message("OUTLET SYNCHRONY PREDICTION RESULTS:\n\n")

if (!is.null(aic_conc) && nrow(aic_conc) > 0) {
  best_conc <- aic_conc$model[1]
  best_r2_conc <- aic_conc$R2_marginal[1]
  message("1. CONCENTRATION SYNC WITH OUTLET:\n")
  message("   Best model:", best_conc, "\n")
  message("   R² (marginal):", round(best_r2_conc, 3), "\n")
  message("   AIC weight:", round(aic_conc$weight[1], 3), "\n\n")
}

if (!is.null(aic_cqslope) && nrow(aic_cqslope) > 0) {
  best_cqslope <- aic_cqslope$model[1]
  best_r2_cqslope <- aic_cqslope$R2_marginal[1]
  message("2. CQ-SLOPE SYNC WITH OUTLET:\n")
  message("   Best model:", best_cqslope, "\n")
  message("   R² (marginal):", round(best_r2_cqslope, 3), "\n")
  message("   AIC weight:", round(aic_cqslope$weight[1], 3), "\n\n")
}

message("3. CLUSTER-CQ RELATIONSHIP:\n")
message("   Clusters", ifelse(summary(anova_result)[[1]][["Pr(>F)"]][1] < 0.05, "DO", "do NOT"),
    "significantly differ in mean CQ slope\n\n")

message("4. CROSS-VALIDATION:\n")
if (!is.null(cv_conc)) {
  message("   Conc sync LOOCV RMSE:", round(mean(cv_conc$rmse, na.rm = TRUE), 4), "\n")
}
if (!is.null(cv_cqslope)) {
  message("   CQ-slope sync LOOCV RMSE:", round(mean(cv_cqslope$rmse, na.rm = TRUE), 4), "\n")
}

message("\n=== ANALYSIS COMPLETE ===\n")
message("Results saved to:", res_dir, "\n")
message("Figures saved to:", fig_dir, "\n")
