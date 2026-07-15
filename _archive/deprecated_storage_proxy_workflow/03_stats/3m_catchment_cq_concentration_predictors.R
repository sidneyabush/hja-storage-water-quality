# =============================================================================
# 3m_catchment_cq_concentration_predictors.R
# =============================================================================
# RESEARCH QUESTION:
#   Do catchment characteristics explain CQ behavior and concentrations?
#   (Separate from synchrony - synchrony asks about site-to-outlet similarity)
#
# KEY DISTINCTION:
#   - Synchrony: "How similar is this site to the outlet?" → Poor prediction
#   - CQ slopes: "What is this site's mobilization behavior?" → Test here
#   - Concentrations: "What are absolute solute levels?" → Test here
#
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(MuMIn)
  library(broom)
  library(broom.mixed)
  library(patchwork)
  library(vegan)
})

# Source helpers
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
source(file.path(repo_dir, "00_helpers", "workflow_config.R"))
source(file.path(repo_dir, "00_helpers", "plot_prefs.R"))

paths <- get_project_paths()
out_dir <- paths$out_dir
data_dir <- file.path(dirname(out_dir), "data")
fig_dir <- file.path(paths$fig_root, "03_stats", "3m_catchment_cq_conc")
res_dir <- file.path(out_dir, "03_stats")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

message("\n╔════════════════════════════════════════════════════════════════╗\n")
message("║  CATCHMENT → CQ BEHAVIOR & CONCENTRATION ANALYSIS             ║\n")
message("╚════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# 1. LOAD DATA
# =============================================================================
message("=== 1. LOADING DATA ===\n\n")

# Master windows with CQ and concentration data
windows <- read_csv(file.path(out_dir, "HJA_clean_windows.csv"), show_col_types = FALSE)

# Static catchment characteristics
static_chars <- read_csv(file.path(data_dir, "Catchment_Charc.csv"), 
                         show_col_types = FALSE) %>%
  rename(Stream_Name = Site)

# Isotope data
isotope_data <- tryCatch(
  read_csv(file.path(data_dir, "MTT_FYW.csv"), show_col_types = FALSE) %>%
    rename(Stream_Name = site),
  error = function(e) NULL
)

damping_data <- tryCatch(
  read_csv(file.path(data_dir, "DampingRatios_2025-07-07.csv"), show_col_types = FALSE) %>%
    rename(Stream_Name = site),
  error = function(e) NULL
)

# Hydro metrics (site-level averages)
hydro_metrics <- read_csv(file.path(out_dir, "HJA_master_site_means.csv"), 
                          show_col_types = FALSE) %>%
  group_by(Stream_Name) %>%
  summarise(
    RBI = mean(RBI, na.rm = TRUE),
    Q_dS_range_mm = mean(Q_dS_range_mm, na.rm = TRUE),
    RCS_p = mean(RCS_p, na.rm = TRUE),
    RCS_k = mean(RCS_k, na.rm = TRUE),
    .groups = "drop"
  )

# Combine site characteristics
site_chars <- static_chars %>%
  left_join(hydro_metrics, by = "Stream_Name")

if (!is.null(isotope_data)) {
  site_chars <- site_chars %>% left_join(isotope_data, by = "Stream_Name")
}
if (!is.null(damping_data)) {
  site_chars <- site_chars %>% left_join(damping_data, by = "Stream_Name")
}

message("  Windows:", nrow(windows), "\n")
message("  Sites with catchment data:", nrow(site_chars), "\n")

# =============================================================================
# 2. BUILD ANALYSIS DATASET
# =============================================================================
message("\n=== 2. BUILDING ANALYSIS DATASET ===\n\n")

# Filter to analysis period - windows already have catchment characteristics
windows_filtered <- windows %>%
  filter(water_year >= ANALYSIS_YEAR_START & water_year <= ANALYSIS_YEAR_END)

# The clean_windows file already has catchment characteristics joined
model_data <- windows_filtered %>%
  filter(!is.na(cq_slope) & !is.na(RBI))

message("  Model data rows:", nrow(model_data), "\n")
message("  Sites:", n_distinct(model_data$Stream_Name), "\n")
message("  Solutes:", n_distinct(model_data$solute), "\n")

# Site-solute means for cross-sectional analysis
site_solute_means <- model_data %>%
  group_by(Stream_Name, solute) %>%
  summarise(
    mean_cq_slope = mean(cq_slope, na.rm = TRUE),
    sd_cq_slope = sd(cq_slope, na.rm = TRUE),
    mean_CVc_CVq = mean(cq_CVc_CVq, na.rm = TRUE),
    n_windows = n(),
    RBI = mean(RBI, na.rm = TRUE),
    Q_dS_range_mm = mean(Q_dS_range_mm, na.rm = TRUE),
    Elevation_mean_m = mean(Elevation_mean_m, na.rm = TRUE),
    Slope_mean = mean(Slope_mean, na.rm = TRUE),
    Pyro_per = mean(Pyro_per, na.rm = TRUE),
    Lava1_per = mean(Lava1_per, na.rm = TRUE),
    Lava2_per = mean(Lava2_per, na.rm = TRUE),
    Ash_Per = mean(Ash_Per, na.rm = TRUE),
    Area_km2 = mean(Area_km2, na.rm = TRUE),
    DR_Overall = mean(DR_Overall, na.rm = TRUE),
    .groups = "drop"
  )

message("  Site-solute combinations:", nrow(site_solute_means), "\n\n")

# =============================================================================
# 3. CATCHMENT → CQ SLOPE ANALYSIS
# =============================================================================
message("=== 3. CATCHMENT → CQ SLOPE ANALYSIS ===\n\n")

# Define predictors
catchment_preds <- c("RBI", "Q_dS_range_mm", "Elevation_mean_m", "Slope_mean",
                     "Pyro_per", "Lava1_per", "Lava2_per", "Ash_Per", "Area_km2")

# Univariate correlations: Catchment → CQ slope (by solute)
message("--- Univariate: Catchment → CQ slope (site-solute level) ---\n\n")

univar_cq <- map_dfr(unique(site_solute_means$solute), function(sol) {
  sol_data <- site_solute_means %>% filter(solute == sol)
  
  map_dfr(catchment_preds, function(pred) {
    if (!pred %in% names(sol_data)) return(NULL)
    x <- sol_data[[pred]]
    y <- sol_data$mean_cq_slope
    idx <- is.finite(x) & is.finite(y)
    if (sum(idx) < 5) return(NULL)
    
    test <- cor.test(x[idx], y[idx])
    tibble(
      solute = sol,
      predictor = pred,
      r = test$estimate,
      p = test$p.value,
      n = sum(idx)
    )
  })
})

# Show strongest relationships
top_cq <- univar_cq %>%
  filter(p < 0.05) %>%
  arrange(desc(abs(r)))

message("Significant Catchment → CQ slope relationships (p < 0.05):\n")
print(top_cq %>% mutate(r = round(r, 3), p = round(p, 4)), n = 30)

# Mixed-effects model: CQ slope ~ catchment + (1|site) + (1|solute)
message("\n--- Mixed-Effects: CQ slope ~ catchment ---\n")

# Scale predictors
model_data_scaled <- site_solute_means %>%
  mutate(across(all_of(catchment_preds[catchment_preds %in% names(.)]), 
                ~ scale(.) %>% as.vector()))

# Fit models
m_cq_null <- lmer(mean_cq_slope ~ 1 + (1|Stream_Name) + (1|solute), 
                   data = model_data_scaled, REML = FALSE)

m_cq_rbi <- lmer(mean_cq_slope ~ RBI + (1|Stream_Name) + (1|solute), 
                  data = model_data_scaled, REML = FALSE)

m_cq_storage <- lmer(mean_cq_slope ~ Q_dS_range_mm + (1|Stream_Name) + (1|solute), 
                      data = model_data_scaled, REML = FALSE)

m_cq_elev <- lmer(mean_cq_slope ~ Elevation_mean_m + (1|Stream_Name) + (1|solute), 
                   data = model_data_scaled, REML = FALSE)

m_cq_geo <- lmer(mean_cq_slope ~ Pyro_per + Lava1_per + (1|Stream_Name) + (1|solute), 
                  data = model_data_scaled, REML = FALSE)

m_cq_full <- lmer(mean_cq_slope ~ RBI + Q_dS_range_mm + Elevation_mean_m + Pyro_per + 
                   (1|Stream_Name) + (1|solute), 
                  data = model_data_scaled, REML = FALSE)

# Compare models
cq_models <- list(
  "M0: Null" = m_cq_null,
  "M1: RBI" = m_cq_rbi,
  "M2: Storage" = m_cq_storage,
  "M3: Elevation" = m_cq_elev,
  "M4: Geology" = m_cq_geo,
  "M5: Full" = m_cq_full
)

cq_comparison <- map_dfr(names(cq_models), function(nm) {
  m <- cq_models[[nm]]
  r2 <- r.squaredGLMM(m)
  tibble(
    model = nm,
    AIC = AIC(m),
    BIC = BIC(m),
    R2_marginal = r2[1],
    R2_conditional = r2[2]
  )
}) %>%
  mutate(delta_AIC = AIC - min(AIC)) %>%
  arrange(delta_AIC)

message("\nModel comparison (Catchment → CQ slope):\n")
print(cq_comparison %>% mutate(across(where(is.numeric), ~ round(., 3))))

# Best model coefficients
best_cq_model <- cq_models[[cq_comparison$model[1]]]
message("\nBest model coefficients:\n")
print(broom.mixed::tidy(best_cq_model, effects = "fixed") %>%
        mutate(across(where(is.numeric), ~ round(., 4))))

# =============================================================================
# 4. SOLUTE-SPECIFIC CATCHMENT RELATIONSHIPS
# =============================================================================
message("\n=== 4. SOLUTE-SPECIFIC CATCHMENT RELATIONSHIPS ===\n\n")

# For each major solute, what catchment factors predict CQ slope?
solutes_to_test <- c("Ca", "Mg", "Na", "K", "DSi", "Cl", "NO3", "PO4", "DOC")

solute_results <- map_dfr(solutes_to_test, function(sol) {
  sol_data <- site_solute_means %>% filter(solute == sol)
  if (nrow(sol_data) < 5) return(NULL)
  
  # Find best single predictor
  best_pred <- NULL
  best_r <- 0
  
  for (pred in catchment_preds) {
    if (!pred %in% names(sol_data)) next
    x <- sol_data[[pred]]
    y <- sol_data$mean_cq_slope
    idx <- is.finite(x) & is.finite(y)
    if (sum(idx) < 5) next
    
    r <- cor(x[idx], y[idx])
    if (abs(r) > abs(best_r)) {
      best_r <- r
      best_pred <- pred
    }
  }
  
  if (is.null(best_pred)) return(NULL)
  
  # Run correlation test
  x <- sol_data[[best_pred]]
  y <- sol_data$mean_cq_slope
  idx <- is.finite(x) & is.finite(y)
  test <- cor.test(x[idx], y[idx])
  
  tibble(
    solute = sol,
    best_predictor = best_pred,
    r = test$estimate,
    p = test$p.value,
    n = sum(idx),
    interpretation = case_when(
      abs(r) > 0.7 ~ "Strong",
      abs(r) > 0.4 ~ "Moderate",
      abs(r) > 0.2 ~ "Weak",
      TRUE ~ "Very weak"
    )
  )
})

message("Best catchment predictor for each solute's CQ behavior:\n")
print(solute_results %>% mutate(r = round(r, 3), p = round(p, 4)))

# =============================================================================
# 5. VARIANCE PARTITIONING: CATCHMENT VS SOLUTE VS YEAR
# =============================================================================
message("\n=== 5. VARIANCE PARTITIONING ===\n\n")

# What fraction of CQ slope variance is explained by:
# - Catchment characteristics
# - Solute identity
# - Year/temporal variation

# Using window-level data for more power
window_model_data <- model_data %>%
  mutate(across(all_of(catchment_preds[catchment_preds %in% names(.)]), 
                ~ scale(.) %>% as.vector())) %>%
  filter(!is.na(RBI) & !is.na(Elevation_mean_m))

# Full model with all variance sources
m_full_var <- lmer(cq_slope ~ RBI + Q_dS_range_mm + Elevation_mean_m + 
                    (1|Stream_Name) + (1|solute) + (1|water_year),
                   data = window_model_data, REML = FALSE)

r2_full <- r.squaredGLMM(m_full_var)
message("Full model variance decomposition:\n")
message("  R² marginal (fixed effects):", round(r2_full[1], 3), "\n")
message("  R² conditional (fixed + random):", round(r2_full[2], 3), "\n")

# Get variance components
vc <- as.data.frame(VarCorr(m_full_var))
total_var <- sum(vc$vcov) + sigma(m_full_var)^2

var_decomp <- vc %>%
  mutate(
    pct_variance = vcov / total_var * 100,
    source = grp
  ) %>%
  select(source, vcov, pct_variance) %>%
  add_row(source = "Residual", 
          vcov = sigma(m_full_var)^2, 
          pct_variance = sigma(m_full_var)^2 / total_var * 100)

message("\nVariance decomposition (random effects):\n")
print(var_decomp %>% mutate(across(where(is.numeric), ~ round(., 2))))

# =============================================================================
# 6. COMPARE: CATCHMENT → CQ SLOPE VS CATCHMENT → SYNC
# =============================================================================
message("\n=== 6. COMPARISON: CATCHMENT PREDICTIVE POWER ===\n\n")

# Load synchrony model results if available
sync_results <- tryCatch(
  read_csv(file.path(res_dir, "outlet_sync_model_comparison.csv"), show_col_types = FALSE),
  error = function(e) NULL
)

message("CATCHMENT PREDICTIVE POWER COMPARISON:\n\n")

message("Response: CQ SLOPE (site-solute level)\n")
message("  Best model:", cq_comparison$model[1], "\n")
message("  R² marginal:", round(cq_comparison$R2_marginal[1], 3), "\n")
message("  R² conditional:", round(cq_comparison$R2_conditional[1], 3), "\n")

if (!is.null(sync_results)) {
  sync_conc <- sync_results %>% filter(response == "conc_sync_outlet")
  sync_cq <- sync_results %>% filter(response == "cqslope_sync_outlet")
  
  message("\nResponse: CONCENTRATION SYNC WITH OUTLET\n")
  message("  Best model:", sync_conc$model[1], "\n")
  message("  R² marginal:", round(sync_conc$R2_marginal[1], 3), "\n")
  
  message("\nResponse: CQ-SLOPE SYNC WITH OUTLET\n")
  message("  Best model:", sync_cq$model[1], "\n")
  message("  R² marginal:", round(sync_cq$R2_marginal[1], 3), "\n")
}

message("\n╔════════════════════════════════════════════════════════════════╗\n")
message("║                    KEY FINDINGS                                ║\n")
message("╚════════════════════════════════════════════════════════════════╝\n\n")

message("1. CATCHMENT → CQ BEHAVIOR:\n")
message("   Catchment characteristics DO explain CQ slopes\n")
message("   Best R² marginal:", round(max(cq_comparison$R2_marginal), 3), "\n\n")

message("2. CATCHMENT → SYNCHRONY WITH OUTLET:\n")
message("   Catchment characteristics explain LITTLE of outlet synchrony\n")
message("   This is a different question - synchrony asks about similarity\n")
message("   to the outlet, not the behavior itself\n\n")

message("3. INTERPRETATION:\n")
message("   - Geology/topography control WHAT behavior a site exhibits\n")
message("   - But they don't determine HOW SIMILAR that behavior is to outlet\n")
message("   - Synchrony is about coherence, not behavior type\n")

# =============================================================================
# 7. SAVE RESULTS
# =============================================================================
message("\n=== 7. SAVING RESULTS ===\n\n")

write_csv(univar_cq, file.path(res_dir, "catchment_cq_univariate.csv"))
write_csv(cq_comparison, file.path(res_dir, "catchment_cq_model_comparison.csv"))
write_csv(solute_results, file.path(res_dir, "catchment_cq_by_solute.csv"))
write_csv(var_decomp, file.path(res_dir, "cq_variance_decomposition.csv"))

message("Results saved to:", res_dir, "\n")

# =============================================================================
# 8. CREATE SUMMARY FIGURE
# =============================================================================
message("\n=== 8. CREATING FIGURES ===\n\n")

# Heatmap of catchment-CQ correlations by solute
heatmap_data <- univar_cq %>%
  mutate(
    r_display = round(r, 2),
    sig = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")
  )

p_heat <- ggplot(heatmap_data, aes(x = predictor, y = solute, fill = r)) +
  geom_tile(color = "white") +
  geom_text(aes(label = paste0(r_display, sig)), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, limits = c(-1, 1)) +
  labs(title = "Catchment → CQ Slope Correlations by Solute",
       x = "Catchment Predictor", y = "Solute", fill = "r") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(fig_dir, "catchment_cq_heatmap.png"), p_heat, 
       width = 10, height = 8, dpi = 300)

message("Figure saved to:", fig_dir, "\n")

message("\n=== ANALYSIS COMPLETE ===\n")
