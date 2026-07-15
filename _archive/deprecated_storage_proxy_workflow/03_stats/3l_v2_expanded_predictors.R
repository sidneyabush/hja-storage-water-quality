# =============================================================================
# STEP 03l_v2: EXPANDED OUTLET SYNCHRONY PREDICTION MODELS
# =============================================================================
# Goal: Test comprehensive predictor set for outlet synchrony prediction
#
# EXPANDED PREDICTOR CATEGORIES:
#
#   STORAGE/HYDROLOGIC (7 predictors):
#     - Q_dS_range_mm     : Within-window storage variability (PRIMARY)
#     - DS_sum_mean       : Annual seasonal drawdown (mm) - DIFFERENT from Q_dS!
#     - Q5norm_mean       : Low-flow conditions (5th percentile normalized)
#     - RBI_mean          : Flashiness index
#     - fdc_slope_mean    : Flow duration curve slope
#     - mean_bf_mean      : Mean baseflow fraction
#     - recession_curve_slope_mean : Recession behavior
#
#   LITHOLOGY (4 predictors):
#     - Lava1_per, Lava2_per : % Basaltic lava
#     - Ash_Per             : % Volcanic ash
#     - Pyro_per            : % Pyroclastic deposits
#
#   TOPOGRAPHY (4 predictors):
#     - Slope_mean         : Mean basin slope
#     - Elevation_mean_m   : Mean elevation  
#     - Area_km2           : Drainage area
#     - Aspec_Mean_deg     : Mean aspect
#
#   LAND USE/DISTURBANCE (5 predictors):
#     - Harvest            : % Harvested
#     - Age                : Stand age
#     - Landslide_Young    : % Young landslide terrain
#     - Landslide_Mod      : % Moderate landslide
#     - Landslide_Old      : % Old landslide
#
#   WATER AGE (3 predictors - Tier 2):
#     - DR_Overall         : Damping ratio (attenuation)
#     - MTTM               : Mean transit time
#     - FYWM               : Young water fraction
#
# CIRCULARITY NOTE:
#   DS_sum_mean (annual drawdown) is NOT circular with Q_dS_range_mm because:
#   - Q_dS captures storage variability WITHIN each window
#   - DS_sum captures seasonal depletion over ENTIRE dry season
#   These are complementary measures of different storage processes.
#
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(MuMIn)
  library(corrplot)
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
data_dir <- file.path(dirname(out_dir), "data")
fig_dir <- file.path(paths$fig_root, "03_stats", "3l_v2_expanded")
res_dir <- file.path(out_dir, "03_stats")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

message("\n╔════════════════════════════════════════════════════════════════╗\n")
message("║  EXPANDED OUTLET SYNCHRONY PREDICTION MODELS (v2)              ║\n")
message("╚════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# 1. LOAD ALL DATA SOURCES
# =============================================================================
message("=== 1. LOADING DATA ===\n\n")

# Core site average data with storage metrics and catchment characteristics
site_avgs <- read_csv(file.path(out_dir, "HJA_master_site_means.csv"),
                      show_col_types = FALSE) %>%
  rename(Stream_Name = site)

message("  Site averages:", nrow(site_avgs), "sites with", ncol(site_avgs), "variables\n")

# Outlet synchrony metrics
outlet_sync_annual <- read_csv(file.path(out_dir, "HJA_outlet_synchrony_annual.csv"), 
                                show_col_types = FALSE)
outlet_sync_site <- read_csv(file.path(out_dir, "HJA_outlet_synchrony_site_level.csv"), 
                              show_col_types = FALSE)

message("  Outlet sync (annual):", nrow(outlet_sync_annual), "rows\n")
message("  Outlet sync (site-level):", nrow(outlet_sync_site), "rows\n")

# Isotope data (MTT, Fyw)
isotope_data <- tryCatch(
  read_csv(file.path(data_dir, "MTT_FYW.csv"), show_col_types = FALSE) %>%
    rename(Stream_Name = site) %>%
    select(Stream_Name, MTTM, FYWM),
  error = function(e) { message("  Note: MTT_FYW.csv not found\n"); NULL }
)

# Cluster data
clusters_modal <- read_csv(file.path(out_dir, "ClusterStreams_allSolutes_modalClusters.csv"), 
                           show_col_types = FALSE)
clusters_stability <- read_csv(file.path(out_dir, "ClusterStreams_stability_metrics.csv"), 
                               show_col_types = FALSE)

# =============================================================================
# 2. BUILD COMPREHENSIVE PREDICTOR DATASET
# =============================================================================
message("\n=== 2. BUILDING PREDICTOR DATASET ===\n\n")

# Define predictor categories for variance partitioning
STORAGE_PREDS <- c("Q_dS_range_mm", "DS_sum_mean", "Q5norm_mean", "RBI_mean", 
                   "fdc_slope_mean", "mean_bf_mean", "recession_curve_slope_mean")

LITHOLOGY_PREDS <- c("Lava1_per", "Lava2_per", "Ash_Per", "Pyro_per")

TOPOGRAPHY_PREDS <- c("Slope_mean", "Elevation_mean_m", "Area_km2", "Aspec_Mean_deg")

LANDUSE_PREDS <- c("Harvest", "Age", "Landslide_Young", "Landslide_Mod", 
                   "Landslide_Old", "Landslide_Total")

WATER_AGE_PREDS <- c("DR_Overall", "MTTM", "FYWM")

ALL_PREDS <- c(STORAGE_PREDS, LITHOLOGY_PREDS, TOPOGRAPHY_PREDS, 
               LANDUSE_PREDS, WATER_AGE_PREDS)

# Get window-level Q_dS from hydro data and average by site
hydro_metrics <- tryCatch({
  read_csv(file.path(out_dir, "HJA_master_site_means.csv"), show_col_types = FALSE) %>%
    group_by(Stream_Name) %>%
    summarise(
      Q_dS_range_mm = mean(Q_dS_range_mm, na.rm = TRUE),
      .groups = "drop"
    )
}, error = function(e) {
  tibble(Stream_Name = character(), Q_dS_range_mm = numeric())
})

# Cluster stability per site
cluster_stability_site <- clusters_stability %>%
  group_by(Stream_Name) %>%
  summarise(
    cluster_stability_mean = mean(stability, na.rm = TRUE),
    .groups = "drop"
  )

# Combine all site-level predictors
site_preds <- site_avgs %>%
  left_join(hydro_metrics, by = "Stream_Name") %>%
  left_join(cluster_stability_site, by = "Stream_Name")

if (!is.null(isotope_data)) {
  site_preds <- site_preds %>%
    left_join(isotope_data, by = "Stream_Name")
}

# Check which predictors are available
available_preds <- intersect(ALL_PREDS, names(site_preds))
message("  Available predictors:", length(available_preds), "/", length(ALL_PREDS), "\n")
message("  Missing:", setdiff(ALL_PREDS, available_preds), "\n\n")

# Build analysis dataset
model_data <- outlet_sync_site %>%
  filter(Stream_Name != OUTLET_SITE) %>%
  left_join(site_preds, by = "Stream_Name")

message("  Sites for modeling:", n_distinct(model_data$Stream_Name), "\n")
message("  Solutes:", n_distinct(model_data$solute), "\n\n")

# =============================================================================
# 3. PREDICTOR CORRELATION STRUCTURE
# =============================================================================
message("=== 3. PREDICTOR CORRELATIONS ===\n\n")

# Create correlation matrix for available predictors
pred_matrix <- model_data %>%
  select(Stream_Name, any_of(available_preds)) %>%
  distinct() %>%
  select(-Stream_Name) %>%
  select(where(~ sum(!is.na(.)) >= 5))

cor_matrix <- cor(pred_matrix, use = "pairwise.complete.obs")

# Save correlation plot
png(file.path(fig_dir, "predictor_correlations.png"), width = 12, height = 10, units = "in", res = 150)
corrplot(cor_matrix, method = "color", type = "lower", 
         tl.col = "black", tl.cex = 0.8,
         addCoef.col = "black", number.cex = 0.6,
         title = "Predictor Correlation Structure",
         mar = c(0, 0, 2, 0))
dev.off()

message("  High correlations (|r| > 0.7):\n")
high_cors <- which(abs(cor_matrix) > 0.7 & abs(cor_matrix) < 1, arr.ind = TRUE)
if (nrow(high_cors) > 0) {
  for (i in seq_len(nrow(high_cors))) {
    r <- cor_matrix[high_cors[i, 1], high_cors[i, 2]]
    message(sprintf("    %s ~ %s: r = %.2f\n", 
                rownames(cor_matrix)[high_cors[i, 1]],
                colnames(cor_matrix)[high_cors[i, 2]], r))
  }
}
message("\n")

# =============================================================================
# 4. UNIVARIATE SCREENING (EXPANDED)
# =============================================================================
message("=== 4. UNIVARIATE SCREENING ===\n\n")

# Aggregate to site-level means
site_means <- model_data %>%
  group_by(Stream_Name) %>%
  summarise(
    conc_sync = mean(conc_sync_outlet_mean, na.rm = TRUE),
    cqslope_sync = mean(cqslope_sync_outlet_mean, na.rm = TRUE),
    across(any_of(available_preds), ~ mean(., na.rm = TRUE)),
    .groups = "drop"
  )

# Compute correlations
univar_results <- map_dfr(available_preds, function(pred) {
  if (!pred %in% names(site_means)) return(NULL)
  
  map_dfr(c("conc_sync", "cqslope_sync"), function(resp) {
    x <- site_means[[pred]]
    y <- site_means[[resp]]
    idx <- is.finite(x) & is.finite(y)
    if (sum(idx) < 4) return(NULL)
    
    test <- cor.test(x[idx], y[idx])
    tibble(
      predictor = pred,
      response = resp,
      r = test$estimate,
      p = test$p.value,
      n = sum(idx),
      sig = case_when(p < 0.01 ~ "**", p < 0.05 ~ "*", p < 0.1 ~ ".", TRUE ~ "")
    )
  })
}) %>%
  arrange(response, p)

# Categorize predictors
univar_results <- univar_results %>%
  mutate(
    category = case_when(
      predictor %in% STORAGE_PREDS ~ "Storage/Hydro",
      predictor %in% LITHOLOGY_PREDS ~ "Lithology",
      predictor %in% TOPOGRAPHY_PREDS ~ "Topography",
      predictor %in% LANDUSE_PREDS ~ "Land Use",
      predictor %in% WATER_AGE_PREDS ~ "Water Age",
      TRUE ~ "Other"
    )
  )

message("--- TOP PREDICTORS FOR CONCENTRATION SYNC ---\n")
univar_results %>%
  filter(response == "conc_sync") %>%
  select(category, predictor, r, p, sig) %>%
  mutate(r = round(r, 3), p = round(p, 3)) %>%
  print(n = 15)

message("\n--- TOP PREDICTORS FOR CQ-SLOPE SYNC ---\n")
univar_results %>%
  filter(response == "cqslope_sync") %>%
  select(category, predictor, r, p, sig) %>%
  mutate(r = round(r, 3), p = round(p, 3)) %>%
  print(n = 15)

# Save results
write_csv(univar_results, file.path(res_dir, "expanded_univariate_screening.csv"))

# =============================================================================
# 5. CATEGORY-BASED COMPARISON
# =============================================================================
message("\n=== 5. CATEGORY-BASED ANALYSIS ===\n\n")

# Best predictor in each category
best_by_category <- univar_results %>%
  group_by(response, category) %>%
  slice_min(p, n = 1) %>%
  ungroup() %>%
  arrange(response, p)

message("Best predictor by category:\n")
print(best_by_category %>% 
        select(response, category, predictor, r, p, sig) %>%
        mutate(r = round(r, 3), p = round(p, 3)), n = 20)

# =============================================================================
# 6. EXPANDED A PRIORI HYPOTHESES
# =============================================================================
message("\n=== 6. A PRIORI HYPOTHESIS MODELS ===\n\n")

# Annual-level data for mixed models
model_annual <- outlet_sync_annual %>%
  filter(Stream_Name != OUTLET_SITE) %>%
  left_join(site_preds, by = "Stream_Name") %>%
  group_by(Stream_Name, water_year) %>%
  summarise(
    conc_sync = mean(conc_sync_outlet, na.rm = TRUE),
    cqslope_sync = mean(cqslope_sync_outlet, na.rm = TRUE),
    across(any_of(available_preds), ~ mean(., na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric) & !matches("sync|year"), ~ scale(.) %>% as.vector()))

# Define hypothesis models
hypotheses <- list(
  "H0_null" = "1",
  "H1_flashiness" = "RBI_mean + Slope_mean",
  "H2_storage_window" = "Q_dS_range_mm + RBI_mean",
  "H3_storage_seasonal" = "DS_sum_mean + Q5norm_mean",  # NEW: seasonal drawdown
  "H4_lithology" = "Pyro_per + Lava1_per + Ash_Per",
  "H5_topography" = "Slope_mean + Elevation_mean_m + Area_km2",
  "H6_landuse" = "Harvest + Landslide_Total",
  "H7_combined_hydro" = "RBI_mean + DS_sum_mean + Slope_mean",
  "H8_water_age" = "DR_Overall"  # Tier 2
)

# Fit models for each response
fit_hypothesis <- function(resp_var, df, hyp_list) {
  results <- map_dfr(names(hyp_list), function(hyp_name) {
    formula_str <- paste0(resp_var, " ~ ", hyp_list[[hyp_name]], " + (1|Stream_Name)")
    
    model <- tryCatch(
      lmer(as.formula(formula_str), data = df, REML = FALSE),
      error = function(e) NULL
    )
    
    if (is.null(model)) return(NULL)
    
    r2 <- tryCatch(r.squaredGLMM(model), error = function(e) c(R2m = NA, R2c = NA))
    
    tibble(
      hypothesis = hyp_name,
      response = resp_var,
      formula = hyp_list[[hyp_name]],
      AIC = AIC(model),
      BIC = BIC(model),
      R2m = r2[1],
      R2c = r2[2],
      n = nrow(model.frame(model))
    )
  })
  
  # Calculate delta AIC and weights
  if (nrow(results) > 0) {
    results <- results %>%
      mutate(
        dAIC = AIC - min(AIC, na.rm = TRUE),
        weight = exp(-0.5 * dAIC) / sum(exp(-0.5 * dAIC), na.rm = TRUE)
      ) %>%
      arrange(AIC)
  }
  
  results
}

# Fit for concentration sync
message("--- CONCENTRATION SYNC MODELS ---\n")
hyp_conc <- fit_hypothesis("conc_sync", model_annual, hypotheses)
print(hyp_conc %>% 
        select(hypothesis, formula, R2m, R2c, dAIC, weight) %>%
        mutate(across(where(is.numeric), ~ round(., 3))), n = 15)

message("\n--- CQ-SLOPE SYNC MODELS ---\n")
hyp_cqslope <- fit_hypothesis("cqslope_sync", model_annual, hypotheses)
print(hyp_cqslope %>% 
        select(hypothesis, formula, R2m, R2c, dAIC, weight) %>%
        mutate(across(where(is.numeric), ~ round(., 3))), n = 15)

# Save hypothesis comparison
bind_rows(hyp_conc, hyp_cqslope) %>%
  write_csv(file.path(res_dir, "expanded_hypothesis_comparison.csv"))

# =============================================================================
# 7. KEY COMPARISONS: STORAGE METRICS
# =============================================================================
message("\n=== 7. STORAGE METRIC COMPARISON ===\n\n")

message("Comparing Q_dS (window-level) vs DS_sum (seasonal drawdown):\n\n")

# Get correlations with each other and with sync
storage_comp <- site_means %>%
  select(Stream_Name, conc_sync, cqslope_sync, 
         any_of(c("Q_dS_range_mm", "DS_sum_mean", "Q5norm_mean", "RBI_mean")))

message("Inter-correlation of storage metrics:\n")
storage_cors <- cor(storage_comp %>% select(-Stream_Name), use = "pairwise.complete.obs")
print(round(storage_cors, 3))

message("\nKey insight: Q_dS_range_mm vs DS_sum_mean correlation:\n")
if (all(c("Q_dS_range_mm", "DS_sum_mean") %in% names(storage_comp))) {
  test <- cor.test(storage_comp$Q_dS_range_mm, storage_comp$DS_sum_mean)
  message(sprintf("  r = %.3f, p = %.4f\n", test$estimate, test$p.value))
  message("  If low correlation → they capture different storage aspects\n")
  message("  If high correlation → may be redundant\n")
}

# =============================================================================
# 8. SUMMARY AND RECOMMENDATIONS
# =============================================================================
message("\n╔════════════════════════════════════════════════════════════════╗\n")
message("║  SUMMARY: EXPANDED PREDICTOR ANALYSIS                          ║\n")
message("╚════════════════════════════════════════════════════════════════╝\n\n")

message("BEST PREDICTORS BY CATEGORY:\n")
print(best_by_category %>% 
        filter(p < 0.1) %>%
        select(response, category, predictor, r, sig), n = 20)

message("\nMODEL COMPARISON SUMMARY:\n")
message("  Best for conc_sync:", hyp_conc$hypothesis[1], 
    sprintf("(R²m = %.3f, weight = %.3f)\n", hyp_conc$R2m[1], hyp_conc$weight[1]))
message("  Best for cqslope_sync:", hyp_cqslope$hypothesis[1],
    sprintf("(R²m = %.3f, weight = %.3f)\n", hyp_cqslope$R2m[1], hyp_cqslope$weight[1]))

message("\nOutputs saved to:\n")
message("  Figures:", fig_dir, "\n")
message("  Results:", res_dir, "\n")
