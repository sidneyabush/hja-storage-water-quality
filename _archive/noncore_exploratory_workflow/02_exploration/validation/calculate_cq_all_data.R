#!/usr/bin/env Rscript
# =============================================================================
# Calculate CQ Analysis Using ALL Data (No Windowing)
# =============================================================================
# Purpose: Calculate site-level CQ slopes using all C-Q pairs (no rolling windows)
#          for comparison with windowed approach (1c_CQ_Rolling_Analysis)
#
# This provides a "true" site-average CQ relationship calculated from the full
# dataset, not averaged from rolling windows.
#
# Output: HJA_CQ_all_data_site_means.csv, HJA_CQ_all_data_seasonal.csv
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

rm(list = ls())

# Paths
base_dir      <- "/Users/sidneybush/Library/CloudStorage/Box-Box"
project_dir   <- file.path(base_dir, "Sidney_Bush", "HJA_Water_Quality")
out_dir       <- file.path(project_dir, "outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# TEMPORAL RANGE (match primary analyses)
# =============================================================================
ANALYSIS_YEAR_START <- 1979
ANALYSIS_YEAR_END   <- 2018

message("[calculate_cq_all_data] CQ All-Data Analysis (no windowing)...")
message("  Temporal range: ", ANALYSIS_YEAR_START, "-", ANALYSIS_YEAR_END)

# =============================================================================
# LOAD DATA
# =============================================================================

# Load seasons
season_df <- readr::read_csv(file.path(out_dir, "HJA_daily_Q_with_seasons.csv"),
                              show_col_types = FALSE) %>%
  mutate(date = as.Date(date), Stream_Name = as.character(Stream_Name)) %>%
  select(Stream_Name, date, water_year, hydrologic_season)

# Load CQ master
cq_master <- readr::read_csv(file.path(out_dir, "HJA_CQ_master.csv"),
                              show_col_types = FALSE) %>%
  rename(date = Date, Qcms = Q_cms, solute = variable) %>%
  mutate(date = as.Date(date), Stream_Name = as.character(Stream_Name)) %>%
  left_join(season_df, by = c("Stream_Name", "date")) %>%
  filter(!is.na(water_year),
         water_year >= ANALYSIS_YEAR_START,
         water_year <= ANALYSIS_YEAR_END,
         !is.na(Qcms), !is.na(value),
         Qcms > 0, value > 0) %>%
  mutate(logQ = log10(Qcms), logC = log10(value))


# =============================================================================
# FIT CQ RELATIONSHIPS - SITE LEVEL (ALL DATA)
# =============================================================================

MIN_OBS <- 10  # Minimum observations to fit CQ relationship

site_cq_all <- cq_master %>%
  group_by(Stream_Name, solute) %>%
  summarise(
    n_obs = n(),
    n_years = n_distinct(water_year),

    # CQ slope (log-log regression)
    cq_slope_all = {
      if (n() >= MIN_OBS &&
          length(unique(logQ)) >= 2 &&
          length(unique(logC)) >= 2 &&
          var(logQ, na.rm = TRUE) > 1e-12 &&
          var(logC, na.rm = TRUE) > 1e-12) {
        coef(lm(logC ~ logQ))[["logQ"]]
      } else {
        NA_real_
      }
    },

    # R-squared
    cq_r2_all = {
      if (n() >= MIN_OBS && !is.na(cq_slope_all)) {
        summary(lm(logC ~ logQ))$r.squared
      } else {
        NA_real_
      }
    },

    # Intercept
    cq_intercept_all = {
      if (n() >= MIN_OBS && !is.na(cq_slope_all)) {
        coef(lm(logC ~ logQ))[["(Intercept)"]]
      } else {
        NA_real_
      }
    },

    # CV(C) and CV(Q)
    CV_C = sd(value, na.rm = TRUE) / mean(value, na.rm = TRUE),
    CV_Q = sd(Qcms, na.rm = TRUE) / mean(Qcms, na.rm = TRUE),

    # CV ratio
    cq_CVc_CVq_all = CV_C / CV_Q,

    # CQ behavior classification (same as windowed approach)
    cq_behavior_all = case_when(
      is.na(cq_slope_all) ~ NA_character_,
      abs(cq_slope_all) < 0.05 ~ "chemostatic",
      cq_slope_all > 0.05 ~ "mobilization",
      cq_slope_all < -0.05 ~ "dilution",
      TRUE ~ "chemostatic"
    ),

    # Mean concentration and discharge
    mean_C = mean(value, na.rm = TRUE),
    mean_Q = mean(Qcms, na.rm = TRUE),

    # Ranges
    range_logC = max(logC, na.rm = TRUE) - min(logC, na.rm = TRUE),
    range_logQ = max(logQ, na.rm = TRUE) - min(logQ, na.rm = TRUE),

    .groups = "drop"
  )



# Summary by behavior
behavior_summary <- site_cq_all %>%
  filter(!is.na(cq_behavior_all)) %>%
  count(cq_behavior_all) %>%
  arrange(desc(n))

for (i in seq_len(nrow(behavior_summary))) {
  message("    ", behavior_summary$cq_behavior_all[i], ": ",
          behavior_summary$n[i], " (",
          round(100 * behavior_summary$n[i] / sum(behavior_summary$n), 1), "%)")
}

# =============================================================================
# FIT CQ RELATIONSHIPS - SEASONAL (ALL DATA WITHIN SEASON)
# =============================================================================
seasonal_cq_all <- cq_master %>%
  filter(!is.na(hydrologic_season)) %>%
  group_by(Stream_Name, solute, hydrologic_season) %>%
  summarise(
    n_obs = n(),

    # CQ slope
    cq_slope_all = {
      if (n() >= MIN_OBS &&
          length(unique(logQ)) >= 2 &&
          length(unique(logC)) >= 2 &&
          var(logQ, na.rm = TRUE) > 1e-12 &&
          var(logC, na.rm = TRUE) > 1e-12) {
        coef(lm(logC ~ logQ))[["logQ"]]
      } else {
        NA_real_
      }
    },

    # R-squared
    cq_r2_all = {
      if (n() >= MIN_OBS && !is.na(cq_slope_all)) {
        summary(lm(logC ~ logQ))$r.squared
      } else {
        NA_real_
      }
    },

    # CV ratio
    CV_C = sd(value, na.rm = TRUE) / mean(value, na.rm = TRUE),
    CV_Q = sd(Qcms, na.rm = TRUE) / mean(Qcms, na.rm = TRUE),
    cq_CVc_CVq_all = CV_C / CV_Q,

    # Behavior
    cq_behavior_all = case_when(
      is.na(cq_slope_all) ~ NA_character_,
      abs(cq_slope_all) < 0.05 ~ "chemostatic",
      cq_slope_all > 0.05 ~ "mobilization",
      cq_slope_all < -0.05 ~ "dilution",
      TRUE ~ "chemostatic"
    ),

    .groups = "drop"
  )

message("  Calculated seasonal CQ slopes for ", nrow(seasonal_cq_all), " combinations")

# =============================================================================
# SAVE OUTPUTS
# =============================================================================

readr::write_csv(site_cq_all, file.path(out_dir, "HJA_CQ_all_data_site_means.csv"))
readr::write_csv(seasonal_cq_all, file.path(out_dir, "HJA_CQ_all_data_seasonal.csv"))

message("  Saved: HJA_CQ_all_data_site_means.csv")
message("  Saved: HJA_CQ_all_data_seasonal.csv")
message("\n[calculate_cq_all_data] Complete")
