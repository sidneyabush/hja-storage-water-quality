# =============================================================================
# STEP 03n: PAIRWISE SYNCHRONY PREDICTION MODELS
# =============================================================================
# Goal: Test whether DIFFERENCES in catchment properties predict pairwise sync
#
# CONTRAST WITH OUTLET-CENTRIC (3l):
#   Outlet: Does Site X sync with GSLOOK? → Null model wins
#   Pairwise: Do similar sites sync with EACH OTHER? → Testing now
#
# HYPOTHESIS:
#   Sites with similar catchment characteristics should be more synchronous
#   → Predict sync(i,j) from |catchment_i - catchment_j| (absolute difference)
#   → Or from mean(catchment_i, catchment_j) (shared environment)
#
# SEASON TIMING ANALYSIS (Dec 2025):
#   Tests whether interannual variation in wet season timing predicts sync
#   → wet_length_days: Duration of wet season (varies ~62-291 days)
#   → wet_start_doy: Day of year wet season begins (varies ~DOY 59-356)
#
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(MuMIn)
  library(broom)
  library(corrplot)
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
fig_dir <- file.path(paths$fig_root, "03_stats", "3n_pairwise")
res_dir <- file.path(out_dir, "03_stats")
update_dir <- file.path(dirname(out_dir), "updates_12042025")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

message("\n╔════════════════════════════════════════════════════════════════╗\n")
message("║  PAIRWISE SYNCHRONY PREDICTION MODELS                          ║\n")
message("╚════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# 1. LOAD DATA
# =============================================================================
message("=== 1. LOADING DATA ===\n\n")

# Site averages with all predictors (now generated in 01_data_prep/1h)
site_avgs <- read_csv(file.path(out_dir, "HJA_master_site_means.csv"),
                      show_col_types = FALSE)

# Pairwise synchrony (Abbott method) - OVERALL aggregated
sync_pairs <- tryCatch(
  read_csv(file.path(out_dir, "HJA_abbott_conc_sitepairs_overall.csv"), show_col_types = FALSE),
  error = function(e) {
    # Try alternate names
    files <- list.files(out_dir, pattern = "sitepair|SitePair|abbott.*pair", ignore.case = TRUE)
    if (length(files) > 0) {
      read_csv(file.path(out_dir, files[1]), show_col_types = FALSE)
    } else {
      NULL
    }
  }
)

# Pairwise synchrony - ANNUAL level (for season timing analysis)
sync_pairs_annual <- tryCatch(
  read_csv(file.path(out_dir, "HJA_pair_sync_metrics.csv"), show_col_types = FALSE) %>%
    filter(time_scale == "annual") %>%
    select(solute, water_year, Stream1, Stream2, Abbott_S, prop_sync_wymore, is_outlet_pair) %>%
    rename(site1 = Stream1, site2 = Stream2),
  error = function(e) NULL
)

# Isotope data
isotope_data <- tryCatch(
  read_csv(file.path(data_dir, "MTT_FYW.csv"), show_col_types = FALSE) %>%
    rename(Stream_Name = site),
  error = function(e) NULL
)

# Season timing data (annual-varying predictors)
season_bounds <- read_csv(file.path(out_dir, "season_boundaries.csv"), 
                          show_col_types = FALSE) %>%
  mutate(
    wet_length_days = as.numeric(wet_end_date - wet_start_date),
    wet_start_doy = lubridate::yday(wet_start_date)
  ) %>%
  select(water_year, wet_length_days, wet_start_doy)

message("  Season timing data loaded:", nrow(season_bounds), "water years\n")

if (is.null(sync_pairs)) {
  message("ERROR: Could not find pairwise synchrony file\n")
  message("Looking for files in:", out_dir, "\n")
  message("Available files:", list.files(out_dir, pattern = "sync|Sync"), "\n")
  quit(status = 1)
}

message("  Pairwise sync (overall) rows:", nrow(sync_pairs), "\n")
if (!is.null(sync_pairs_annual)) {
  message("  Pairwise sync (annual) rows:", nrow(sync_pairs_annual), "\n")
}
message("  Columns:", paste(names(sync_pairs)[1:min(10, ncol(sync_pairs))], collapse = ", "), "\n\n")

# =============================================================================
# 2. IDENTIFY SITE PAIR COLUMNS
# =============================================================================
message("=== 2. PREPARING PAIRWISE DATA ===\n\n")

# Standardize column names
if ("site1" %in% names(sync_pairs) && "site2" %in% names(sync_pairs)) {
  # Already correct
} else if ("Stream_Name1" %in% names(sync_pairs)) {
  sync_pairs <- sync_pairs %>% rename(site1 = Stream_Name1, site2 = Stream_Name2)
} else {
  site_cols <- names(sync_pairs)[grepl("site|stream", names(sync_pairs), ignore.case = TRUE)]
  if (length(site_cols) >= 2) {
    sync_pairs <- sync_pairs %>% rename(site1 = !!site_cols[1], site2 = !!site_cols[2])
  }
}

# Find synchrony column
sync_col <- names(sync_pairs)[grepl("abbott|sync|r_conc", names(sync_pairs), ignore.case = TRUE)][1]
if (!is.na(sync_col) && sync_col != "synchrony") {
  sync_pairs <- sync_pairs %>% rename(synchrony = !!sync_col)
}

message("  Site columns: site1, site2\n")
message("  Sync column: synchrony\n")
message("  Unique site pairs:", nrow(sync_pairs %>% distinct(site1, site2)), "\n\n")

# =============================================================================
# 3. BUILD PAIRWISE PREDICTOR DIFFERENCES
# =============================================================================
message("=== 3. COMPUTING PAIRWISE PREDICTOR DIFFERENCES ===\n\n")

# Combine site averages with isotope data
all_preds <- site_avgs
if (!is.null(isotope_data)) {
  all_preds <- all_preds %>% left_join(isotope_data, by = "Stream_Name")
}

# Key predictors to test (based on outlet analysis)
key_preds <- c(
  # Storage/Hydro
  "RBI_mean", "Q5norm_mean", "DS_sum_mean", "mean_bf_mean", "fdc_slope_mean",
  "recession_curve_slope_mean", "S_annual_mm_mean",
  # Lithology
  "Lava1_per", "Lava2_per", "Ash_Per", "Pyro_per",
  # Topography
  "Slope_mean", "Elevation_mean_m", "Area_km2",
  # Land use
  "Harvest", "Landslide_Total",
  # Water age
  "DR_Overall", "MTTM"
)

# Filter to available predictors
available_preds <- intersect(key_preds, names(all_preds))
message("  Available predictors:", length(available_preds), "\n")

# Create pairwise difference and mean datasets
create_pair_features <- function(df, site1_col, site2_col, pred_cols, site_data) {
  
  # Get predictor values for site1
  site1_data <- site_data %>%
    select(Stream_Name, all_of(pred_cols)) %>%
    rename_with(~ paste0(.x, "_1"), -Stream_Name)
  
  # Get predictor values for site2
  site2_data <- site_data %>%
    select(Stream_Name, all_of(pred_cols)) %>%
    rename_with(~ paste0(.x, "_2"), -Stream_Name)
  
  # Join to pairs
  result <- df %>%
    left_join(site1_data, by = c("site1" = "Stream_Name")) %>%
    left_join(site2_data, by = c("site2" = "Stream_Name"))
  
  # Compute differences and means
  for (pred in pred_cols) {
    col1 <- paste0(pred, "_1")
    col2 <- paste0(pred, "_2")
    
    if (col1 %in% names(result) && col2 %in% names(result)) {
      # Absolute difference (dissimilarity)
      result[[paste0("diff_", pred)]] <- abs(result[[col1]] - result[[col2]])
      # Mean (shared environment)
      result[[paste0("mean_", pred)]] <- (result[[col1]] + result[[col2]]) / 2
    }
  }
  
  result
}

# Apply to sync pairs
pair_data <- create_pair_features(sync_pairs, "site1", "site2", available_preds, all_preds)

message("  Pair data columns:", ncol(pair_data), "\n")
message("  Complete cases:", sum(complete.cases(pair_data %>% select(synchrony, starts_with("diff_")))), "\n\n")

# =============================================================================
# 4. UNIVARIATE SCREENING: DIFFERENCES
# =============================================================================
message("=== 4. UNIVARIATE SCREENING: DO DIFFERENCES PREDICT SYNC? ===\n\n")

# Aggregate to site-pair level (average across solutes/years if needed)
if ("solute" %in% names(pair_data)) {
  pair_means <- pair_data %>%
    group_by(site1, site2) %>%
    summarise(
      synchrony = mean(synchrony, na.rm = TRUE),
      across(starts_with(c("diff_", "mean_")), ~ mean(., na.rm = TRUE)),
      n_obs = n(),
      .groups = "drop"
    )
} else {
  pair_means <- pair_data %>%
    group_by(site1, site2) %>%
    summarise(
      synchrony = mean(synchrony, na.rm = TRUE),
      across(starts_with(c("diff_", "mean_")), ~ mean(., na.rm = TRUE)),
      n_obs = n(),
      .groups = "drop"
    )
}

message("  Unique pairs for analysis:", nrow(pair_means), "\n\n")

# Screen difference predictors
diff_preds <- names(pair_means)[grepl("^diff_", names(pair_means))]
mean_preds <- names(pair_means)[grepl("^mean_", names(pair_means))]

screen_predictors <- function(df, response, predictors, type_label) {
  results <- map_dfr(predictors, function(pred) {
    x <- df[[pred]]
    y <- df[[response]]
    idx <- is.finite(x) & is.finite(y)
    if (sum(idx) < 5) return(NULL)
    
    test <- cor.test(x[idx], y[idx])
    tibble(
      Type = type_label,
      Predictor = pred,
      r = test$estimate,
      p = test$p.value,
      n = sum(idx),
      sig = case_when(p < 0.01 ~ "**", p < 0.05 ~ "*", p < 0.1 ~ ".", TRUE ~ "")
    )
  }) %>%
    arrange(p)
  
  results
}

# Screen both types
diff_results <- screen_predictors(pair_means, "synchrony", diff_preds, "Difference")
mean_results <- screen_predictors(pair_means, "synchrony", mean_preds, "Mean")

all_results <- bind_rows(diff_results, mean_results) %>%
  mutate(
    Variable = gsub("^(diff_|mean_)", "", Predictor),
    Direction = case_when(
      Type == "Difference" & r < 0 ~ "Similar sites sync MORE",
      Type == "Difference" & r > 0 ~ "Similar sites sync LESS",
      Type == "Mean" & r > 0 ~ "Higher values = more sync",
      Type == "Mean" & r < 0 ~ "Higher values = less sync"
    )
  )

message("--- DIFFERENCE PREDICTORS (|site_i - site_j|) ---\n")
message("(Negative r = similar sites sync more)\n\n")
print(diff_results %>% mutate(r = round(r, 3), p = round(p, 3)) %>% head(15))

message("\n--- MEAN PREDICTORS (mean of site_i and site_j) ---\n")
message("(Tests if pairs in similar environments sync more)\n\n")
print(mean_results %>% mutate(r = round(r, 3), p = round(p, 3)) %>% head(15))

# =============================================================================
# 5. KEY COMPARISON: OUTLET VS PAIRWISE
# =============================================================================
message("\n=== 5. COMPARING OUTLET VS PAIRWISE RESULTS ===\n\n")

# Summarize significant results
sig_diff <- diff_results %>% filter(p < 0.1)
sig_mean <- mean_results %>% filter(p < 0.1)

message("SIGNIFICANT DIFFERENCE PREDICTORS (p < 0.1):\n")
if (nrow(sig_diff) > 0) {
  print(sig_diff %>% select(Predictor, r, p, sig) %>% mutate(r = round(r, 3), p = round(p, 3)))
} else {
  message("  None\n")
}

message("\nSIGNIFICANT MEAN PREDICTORS (p < 0.1):\n")
if (nrow(sig_mean) > 0) {
  print(sig_mean %>% select(Predictor, r, p, sig) %>% mutate(r = round(r, 3), p = round(p, 3)))
} else {
  message("  None\n")
}

# =============================================================================
# 6. MIXED EFFECTS MODELS
# =============================================================================
message("\n=== 6. MIXED EFFECTS MODELS ===\n\n")

# If we have enough data, fit mixed models with site random effects
if (nrow(pair_data) > 50 && "solute" %in% names(pair_data)) {
  
  # Prepare data with scaled predictors
  model_df <- pair_data %>%
    filter(!is.na(synchrony)) %>%
    mutate(across(starts_with("diff_"), ~ scale(.) %>% as.vector()))
  
  # Null model
  m0 <- lmer(synchrony ~ 1 + (1|site1) + (1|site2), data = model_df, REML = FALSE)
  
  # Best difference predictor model
  if (nrow(diff_results) > 0) {
    best_diff <- diff_results$Predictor[1]
    formula_str <- paste0("synchrony ~ ", best_diff, " + (1|site1) + (1|site2)")
    m1 <- tryCatch(
      lmer(as.formula(formula_str), data = model_df, REML = FALSE),
      error = function(e) NULL
    )
    
    if (!is.null(m1)) {
      message("Model comparison:\n")
      message("  Null model AIC:", round(AIC(m0), 1), "\n")
      message("  Best diff model AIC:", round(AIC(m1), 1), "(", best_diff, ")\n")
      message("  Delta AIC:", round(AIC(m1) - AIC(m0), 1), "\n")
      
      r2 <- r.squaredGLMM(m1)
      message("  R²m (fixed):", round(r2[1], 3), "\n")
      message("  R²c (total):", round(r2[2], 3), "\n")
    }
  }
}

# =============================================================================
# 7. CREATE COMPARISON FIGURE
# =============================================================================
message("\n=== 7. CREATING FIGURES ===\n\n")

# Figure: Best predictor relationships
if (nrow(diff_results) > 0) {
  
  # Top difference predictor
  top_diff <- diff_results$Predictor[1]
  
  p1 <- ggplot(pair_means, aes_string(x = top_diff, y = "synchrony")) +
    geom_point(alpha = 0.6, size = 3, color = "#2166AC") +
    geom_smooth(method = "lm", se = TRUE, color = "#B2182B") +
    labs(
      x = gsub("diff_", "Δ ", top_diff),
      y = "Pairwise Synchrony",
      title = paste0("Pairwise Sync vs. Catchment Difference"),
      subtitle = sprintf("r = %.3f, p = %.3f", diff_results$r[1], diff_results$p[1])
    ) +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  ggsave(file.path(fig_dir, "pairwise_top_difference.png"), p1, width = 8, height = 6, dpi = 300)
  
  # Copy to update folder
  if (dir.exists(update_dir)) {
    ggsave(file.path(update_dir, "fig7_pairwise_sync_vs_catchment_diff.png"), p1, 
           width = 8, height = 6, dpi = 300)
  }
}

# Summary bar plot
if (nrow(all_results) > 0) {
  plot_data <- all_results %>%
    filter(p < 0.2) %>%
    mutate(
      abs_r = abs(r),
      Predictor = gsub("^(diff_|mean_)", "", Predictor),
      sig_color = case_when(p < 0.05 ~ "p < 0.05", p < 0.1 ~ "p < 0.1", TRUE ~ "p >= 0.1")
    )
  
  if (nrow(plot_data) > 0) {
    p2 <- ggplot(plot_data, aes(x = r, y = reorder(Predictor, abs_r), fill = sig_color)) +
      geom_col() +
      geom_vline(xintercept = 0, linetype = "dashed") +
      facet_wrap(~Type, scales = "free_y") +
      scale_fill_manual(values = c("p < 0.05" = "#B2182B", "p < 0.1" = "#EF8A62", "p >= 0.1" = "gray70")) +
      labs(
        x = "Correlation with Pairwise Synchrony",
        y = "",
        fill = "Significance",
        title = "Pairwise Synchrony: What Predicts Site-to-Site Agreement?",
        subtitle = "Left: Catchment differences | Right: Shared environment (mean)"
      ) +
      theme_bw() +
      theme(panel.grid = element_blank())
    
    ggsave(file.path(fig_dir, "pairwise_predictor_screening.png"), p2, width = 12, height = 8, dpi = 300)
    
    if (dir.exists(update_dir)) {
      ggsave(file.path(update_dir, "fig8_pairwise_predictor_comparison.png"), p2, 
             width = 12, height = 8, dpi = 300)
    }
  }
}

# =============================================================================
# 8. SEASON TIMING ANALYSIS (ANNUAL LEVEL)
# =============================================================================
message("\n=== 8. SEASON TIMING AS PREDICTOR (ANNUAL LEVEL) ===\n\n")

# This section tests whether interannual variation in wet season timing
# predicts pairwise synchrony

if (!is.null(sync_pairs_annual) && nrow(sync_pairs_annual) > 0) {
  
  # Filter to non-outlet pairs
  annual_nonoutlet <- sync_pairs_annual %>%
    filter(!is_outlet_pair) %>%
    left_join(season_bounds, by = "water_year")
  
  message("  Annual non-outlet pairs with season data:", nrow(annual_nonoutlet), "\n")
  
  # Aggregate to pair-year level (average across solutes)
  annual_pair_year <- annual_nonoutlet %>%
    group_by(site1, site2, water_year, wet_length_days, wet_start_doy) %>%
    summarise(
      abbott_sync = mean(Abbott_S, na.rm = TRUE),
      wymore_sync = mean(prop_sync_wymore, na.rm = TRUE),
      n_solutes = n(),
      .groups = "drop"
    )
  
  message("  Pair-years for analysis:", nrow(annual_pair_year), "\n\n")
  
  # Test season timing as predictors
  message("--- Season timing correlations with pairwise sync ---\n\n")
  
  # Abbott sync vs season timing
  if (sum(is.finite(annual_pair_year$abbott_sync)) > 10) {
    cor_length_abbott <- cor.test(annual_pair_year$wet_length_days, annual_pair_year$abbott_sync,
                                   use = "complete.obs")
    cor_start_abbott <- cor.test(annual_pair_year$wet_start_doy, annual_pair_year$abbott_sync,
                                  use = "complete.obs")
    
    message("Abbott sync:\n")
    message(sprintf("  wet_length_days: r = %.3f, p = %.4f %s\n", 
                cor_length_abbott$estimate, cor_length_abbott$p.value,
                ifelse(cor_length_abbott$p.value < 0.05, "*", "")))
    message(sprintf("  wet_start_doy:   r = %.3f, p = %.4f %s\n",
                cor_start_abbott$estimate, cor_start_abbott$p.value,
                ifelse(cor_start_abbott$p.value < 0.05, "*", "")))
  }
  
  # Wymore sync vs season timing
  if (sum(is.finite(annual_pair_year$wymore_sync)) > 10) {
    cor_length_wymore <- cor.test(annual_pair_year$wet_length_days, annual_pair_year$wymore_sync,
                                   use = "complete.obs")
    cor_start_wymore <- cor.test(annual_pair_year$wet_start_doy, annual_pair_year$wymore_sync,
                                  use = "complete.obs")
    
    message("\nWymore sync:\n")
    message(sprintf("  wet_length_days: r = %.3f, p = %.4f %s\n",
                cor_length_wymore$estimate, cor_length_wymore$p.value,
                ifelse(cor_length_wymore$p.value < 0.05, "*", "")))
    message(sprintf("  wet_start_doy:   r = %.3f, p = %.4f %s\n",
                cor_start_wymore$estimate, cor_start_wymore$p.value,
                ifelse(cor_start_wymore$p.value < 0.05, "*", "")))
  }
  
  # Save season timing results
  season_results <- tibble(
    Predictor = c("wet_length_days", "wet_start_doy", "wet_length_days", "wet_start_doy"),
    Response = c("Abbott_sync", "Abbott_sync", "Wymore_sync", "Wymore_sync"),
    r = c(cor_length_abbott$estimate, cor_start_abbott$estimate,
          cor_length_wymore$estimate, cor_start_wymore$estimate),
    p = c(cor_length_abbott$p.value, cor_start_abbott$p.value,
          cor_length_wymore$p.value, cor_start_wymore$p.value),
    sig = ifelse(c(cor_length_abbott$p.value, cor_start_abbott$p.value,
                   cor_length_wymore$p.value, cor_start_wymore$p.value) < 0.05, "*", "")
  )
  
  write_csv(season_results, file.path(res_dir, "pairwise_season_timing_predictors.csv"))
  
  # Create visualization
  p_season <- annual_pair_year %>%
    pivot_longer(cols = c(wet_length_days, wet_start_doy), 
                 names_to = "season_metric", values_to = "season_value") %>%
    pivot_longer(cols = c(abbott_sync, wymore_sync),
                 names_to = "sync_type", values_to = "sync_value") %>%
    mutate(
      season_metric = recode(season_metric, 
                              "wet_length_days" = "Wet Season Length (days)",
                              "wet_start_doy" = "Wet Season Start (DOY)"),
      sync_type = recode(sync_type,
                          "abbott_sync" = "Abbott Sync",
                          "wymore_sync" = "Wymore Sync")
    ) %>%
    ggplot(aes(x = season_value, y = sync_value)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_smooth(method = "lm", se = TRUE, color = "steelblue") +
    facet_grid(sync_type ~ season_metric, scales = "free") +
    labs(
      x = "Season Metric Value",
      y = "Pairwise Synchrony",
      title = "Season Timing as Predictor of Pairwise Synchrony",
      subtitle = "Annual-level analysis: Does interannual variation in wet season timing affect synchrony?"
    ) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "gray90"))
  
  ggsave(file.path(fig_dir, "season_timing_vs_pairwise_sync.png"), p_season, 
         width = 10, height = 8, dpi = 300)
  
  if (dir.exists(update_dir)) {
    ggsave(file.path(update_dir, "fig9_season_timing_pairwise.png"), p_season,
           width = 10, height = 8, dpi = 300)
  }
  
  message("\n  Season timing analysis saved.\n")
  
} else {
  message("  No annual pairwise data available for season timing analysis.\n")
}

# =============================================================================
# 9. SUMMARY
# =============================================================================
message("\n╔════════════════════════════════════════════════════════════════╗\n")
message("║  PAIRWISE VS OUTLET-CENTRIC COMPARISON                         ║\n")
message("╚════════════════════════════════════════════════════════════════╝\n\n")

message("OUTLET-CENTRIC (from 3l):\n")
message("  - Catchment properties do NOT predict sync with outlet\n")
message("  - Null model wins (AIC weight ≈ 1.0)\n")
message("  - Exception: DR_Overall for CQ-slope sync (r = -0.96)\n\n")

message("PAIRWISE (this analysis):\n")
if (nrow(sig_diff) > 0 || nrow(sig_mean) > 0) {
  message("  - Significant predictors found!\n")
  if (nrow(sig_diff) > 0) {
    message("  - Best DIFFERENCE predictor:", sig_diff$Predictor[1], 
        sprintf("(r = %.3f, p = %.3f)\n", sig_diff$r[1], sig_diff$p[1]))
  }
  if (nrow(sig_mean) > 0) {
    message("  - Best MEAN predictor:", sig_mean$Predictor[1],
        sprintf("(r = %.3f, p = %.3f)\n", sig_mean$r[1], sig_mean$p[1]))
  }
} else {
  message("  - No significant predictors (p < 0.1)\n")
  message("  - Similar to outlet-centric: catchment doesn't predict sync\n")
}

message("\nINTERPRETATION:\n")
if (nrow(sig_diff) > 0 && sig_diff$r[1] < 0) {
  message("  Sites with SIMILAR catchment properties sync MORE with each other\n")
  message("  → Catchment controls pairwise sync even though outlet sync is universal\n")
} else if (nrow(sig_diff) == 0 && nrow(sig_mean) == 0) {
  message("  Catchment properties don't predict sync at ANY scale\n")
  message("  → Synchrony emerges from shared climate/hydrology, not static catchment\n")
}

# Save results
write_csv(all_results, file.path(res_dir, "pairwise_sync_predictors.csv"))
if (dir.exists(update_dir)) {
  write_csv(all_results, file.path(update_dir, "09_pairwise_sync_predictors.csv"))
}

message("\nOutputs saved to:\n")
message("  ", fig_dir, "\n")
message("  ", res_dir, "\n")
if (dir.exists(update_dir)) message("  ", update_dir, "\n")
