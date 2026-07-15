#!/usr/bin/env Rscript
# =============================================================================
# 1h: Combine Master + Clean Tables (FINAL DATA AGGREGATION)
# =============================================================================
# This script combines ALL data sources and creates final analysis files:
#
# DATA SOURCES (clearly labeled):
#   1. WINDOWED chemistry data (wet=75d, dry=150d rolling windows, then aggregated)
#      - CQ metrics: MUST use windowed (matches temporal scale of samples)
#      - Source: 1c_CQ_Rolling_Analysis
#
#   2. STORAGE-PAPER data (added in 1j_import_storage_framework.R)
#      - Finalized storage-paper metrics and Figure 7 framework axes
#      - Source: 05_Storage_Manuscript/final_workflow
#
#   3. STATIC data (catchment characteristics, isotope metrics, etc.)
#      - Topography, geology, land use, MTT, damping ratios
#      - Source: data/Catchment_Charc.csv, data/MTT_FYW.csv, etc.
#
# OUTPUTS:
#   Master files (comprehensive, all columns):
#     - HJA_master_rolling_windows.csv
#     - HJA_master_seasonal_WINDOWED.csv
#     - HJA_master_annual_WINDOWED.csv
#     - HJA_master_site_means_WINDOWED.csv
#
#   Clean files (focused on key predictors for analysis):
#     - HJA_clean_windows.csv
#     - HJA_clean_seasonal_WINDOWED.csv
#     - HJA_clean_annual_WINDOWED.csv
#     - HJA_clean_site_means_WINDOWED.csv
#
# NOTE: All windowed aggregates are clearly labeled as "WINDOWED" in filenames
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(stringr)
})

rm(list = ls())

# =============================================================================
# SETUP
# =============================================================================

base_dir    <- "/Users/sidneybush/Library/CloudStorage/Box-Box"
project_dir <- file.path(base_dir, "Sidney_Bush", "HJA_Water_Quality")
out_dir     <- file.path(project_dir, "outputs")
data_dir    <- file.path(project_dir, "data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

repo_dir <- Sys.getenv(
  "HJA_WQ_REPO_DIR",
  unset = "/Users/sidneybush/Documents/GitHub/hja-water-quality"
)
source(file.path(repo_dir, "00_helpers", "workflow_config.R"))


# Standardize stream IDs
standardize_stream_ids <- function(df) {
  df %>% mutate(across(any_of(c("Stream_Name","Stream1","Stream2","site","Site")),
                       ~ case_when(. %in% c("GSWSMC","GSMACK") ~ "GSMACK", TRUE ~ as.character(.))))
}

# =============================================================================
# LOAD CLUSTER DATA
# =============================================================================
cluster_old_all <- readr::read_csv(file.path(out_dir, "ClusterStreams_allSolutes.csv"),
                                   show_col_types = FALSE) %>%
  standardize_stream_ids()

cluster_wy_all  <- readr::read_csv(file.path(out_dir, "ClusterStreams_allSolutes_byWaterYear.csv"),
                                   show_col_types = FALSE) %>%
  standardize_stream_ids()

cluster_old_solute <- cluster_old_all %>%
  group_by(Stream_Name, chemical) %>%
  summarise(Cluster_oldSolute = dplyr::first(Cluster_climRef), .groups = "drop") %>%
  transmute(Stream_Name, solute = chemical, Cluster_oldSolute)

cluster_wy_solute <- cluster_wy_all %>%
  transmute(Stream_Name, solute = chemical, water_year = water_year,
            Cluster_wy = Cluster_climRef, dist_oldRef = dist_climRef)

cluster_mode_wy <- cluster_wy_solute %>%
  group_by(Stream_Name, solute) %>%
  summarise(Cluster_mode_wy = names(sort(table(Cluster_wy), decreasing = TRUE))[1],
            .groups = "drop")

# =============================================================================
# LOAD WINDOWED CHEMISTRY DATA
# =============================================================================
# CQ rolling window results
cq_windows_all <- readr::read_csv(file.path(out_dir, "CQ_rolling_window_results.csv"),
                                  show_col_types = FALSE) %>%
  standardize_stream_ids()

# Composite synchrony (site-solute level, aggregated across all pairs)
composite_sync <- readr::read_csv(file.path(out_dir, "HJA_composite_synchrony.csv"),
                                  show_col_types = FALSE) %>%
  standardize_stream_ids()

# Build the active window table directly from CQ windows. Storage metrics are
# joined later from the finalized storage-paper outputs in step 1j.
rolling_windows <- cq_windows_all %>%
  filter(comparison_type == "cqslope_CVcCVq") %>%
  transmute(
    Stream_Name,
    window_center = as.Date(window_center),
    water_year = as.integer(water_year),
    hydrologic_season,
    solute = solute1,
    cq_slope_windowed = slope_x,
    cq_behavior = cq_behavior,
    cq_CVc_CVq_windowed = CV,
    cq_sync = sync
  ) %>%
  filter(!is.na(water_year),
         !is.na(solute),
         water_year >= ANALYSIS_YEAR_START,
         water_year <= ANALYSIS_YEAR_END)

# =============================================================================
# DEFINE METRIC COLUMNS BY DATA SOURCE
# =============================================================================

# WINDOWED CQ METRICS (must use windowed to match sample temporal scale)
cq_windowed_cols <- c("cq_slope_windowed","cq_CVc_CVq_windowed")

# Active water-quality window metrics. Storage metrics are appended in step 1j.
windowed_cols <- intersect(cq_windowed_cols, names(rolling_windows))

# =============================================================================
# AGGREGATE WINDOWED DATA
# =============================================================================

# SEASONAL AGGREGATION (windowed → seasonal average)
seasonal_windowed <- rolling_windows %>%
  filter(!is.na(hydrologic_season), !is.na(water_year), !is.na(solute)) %>%
  group_by(Stream_Name, solute, water_year, hydrologic_season) %>%
  summarise(
    across(all_of(windowed_cols), ~ mean(., na.rm = TRUE)),
    # Compute CQ behavior proportions
    n_windows = n(),
    prop_mobilizing = sum(cq_behavior == "mobilizing", na.rm = TRUE) / n(),
    prop_diluting = sum(cq_behavior == "diluting", na.rm = TRUE) / n(),
    prop_chemostatic = sum(cq_behavior == "chemostatic", na.rm = TRUE) / n(),
    # Also keep legacy names for backward compatibility
    prop_enrich = prop_mobilizing,
    prop_dilute = prop_diluting,
    prop_chemostat = prop_chemostatic,
    .groups = "drop"
  ) %>%
  mutate(data_source = "WINDOWED",
         temporal_scale = "seasonal") %>%
  left_join(cluster_wy_solute, by = c("Stream_Name","solute","water_year")) %>%
  left_join(cluster_old_solute, by = c("Stream_Name","solute")) %>%
  left_join(cluster_mode_wy, by = c("Stream_Name","solute"))

# ANNUAL AGGREGATION (windowed → annual average)
annual_windowed <- rolling_windows %>%
  filter(!is.na(water_year), !is.na(solute)) %>%
  group_by(Stream_Name, solute, water_year) %>%
  summarise(
    across(all_of(windowed_cols), ~ mean(., na.rm = TRUE)),
    # Compute CQ behavior proportions
    n_windows = n(),
    prop_mobilizing = sum(cq_behavior == "mobilizing", na.rm = TRUE) / n(),
    prop_diluting = sum(cq_behavior == "diluting", na.rm = TRUE) / n(),
    prop_chemostatic = sum(cq_behavior == "chemostatic", na.rm = TRUE) / n(),
    # Also keep legacy names for backward compatibility
    prop_enrich = prop_mobilizing,
    prop_dilute = prop_diluting,
    prop_chemostat = prop_chemostatic,
    .groups = "drop"
  ) %>%
  mutate(data_source = "WINDOWED",
         temporal_scale = "annual") %>%
  left_join(cluster_wy_solute, by = c("Stream_Name","solute","water_year")) %>%
  left_join(cluster_old_solute, by = c("Stream_Name","solute")) %>%
  left_join(cluster_mode_wy, by = c("Stream_Name","solute"))

# SITE-LEVEL MEANS (windowed → site longterm average)
# Include prop_* columns in the averaging
prop_cols <- c("prop_mobilizing", "prop_diluting", "prop_chemostatic",
               "prop_enrich", "prop_dilute", "prop_chemostat")
site_mean_cols <- c(windowed_cols, prop_cols)
site_mean_cols <- intersect(site_mean_cols, names(annual_windowed))

site_means_windowed <- annual_windowed %>%
  group_by(Stream_Name, solute) %>%
  summarise(across(all_of(site_mean_cols), ~ mean(., na.rm = TRUE)), .groups = "drop") %>%
  mutate(data_source = "WINDOWED",
         temporal_scale = "site_longterm") %>%
  left_join(cluster_old_solute, by = c("Stream_Name","solute")) %>%
  left_join(cluster_mode_wy, by = c("Stream_Name","solute")) %>%
  left_join(composite_sync, by = c("Stream_Name","solute"))

# =============================================================================
# LOAD STATIC CATCHMENT CHARACTERISTICS
# =============================================================================

static_damping <- readr::read_csv(file.path(data_dir, "DampingRatios_2025-07-07.csv"),
                                  show_col_types = FALSE) %>%
  rename(Stream_Name = site) %>% standardize_stream_ids()

static_mtt_fyw <- readr::read_csv(file.path(data_dir, "MTT_FYW.csv"),
                                  show_col_types = FALSE) %>%
  rename(Stream_Name = site) %>% standardize_stream_ids()

static_chars   <- readr::read_csv(file.path(data_dir, "Catchment_Charc.csv"),
                                  show_col_types = FALSE) %>%
  rename(Stream_Name = Site) %>% standardize_stream_ids()

# Load "real" averaged data (stream temp, isotope metrics)
# These are site-level averages calculated directly (not from rolling windows)
static_ave_metrics <- readr::read_csv(file.path(data_dir, "HJA_Ave_StorageMetrics_CatCharacter.csv"),
                                      show_col_types = FALSE) %>%
  rename(Stream_Name = site) %>%
  standardize_stream_ids() %>%
  select(-any_of("...1")) %>%
  # Keep only "real" averaged metrics (stream temp, isotope data)
  select(Stream_Name, matches("JulyM|JST_AT|^Segu$|^McGuire$|^Ortega$|_err$|^DR_Overall$|DR__err"))

# Combine static data
static_site <- static_chars %>%
  mutate(Stream_Name = str_squish(as.character(Stream_Name))) %>%
  left_join(static_damping, by = "Stream_Name") %>%
  left_join(static_mtt_fyw, by = "Stream_Name") %>%
  left_join(static_ave_metrics, by = "Stream_Name")


# =============================================================================
# BUILD MASTER DATASETS
# =============================================================================

# MASTER ROLLING WINDOWS (window-level, includes WINDOWED CQ metrics + static)
mega_window_df <- rolling_windows %>%
  left_join(static_site, by = "Stream_Name") %>%
  left_join(cluster_wy_solute, by = c("Stream_Name","solute","water_year")) %>%
  left_join(cluster_old_solute, by = c("Stream_Name","solute")) %>%
  left_join(cluster_mode_wy, by = c("Stream_Name","solute")) %>%
  filter(!is.na(solute)) %>%
  mutate(data_source = "WINDOWED")

readr::write_csv(mega_window_df, file.path(out_dir, "HJA_master_rolling_windows.csv"))

# MASTER SEASONAL (seasonal aggregates, WINDOWED metrics + static)
seasonal_windowed_with_static <- seasonal_windowed %>%
  left_join(static_site, by = "Stream_Name")

readr::write_csv(seasonal_windowed_with_static,
                 file.path(out_dir, "HJA_master_seasonal_WINDOWED.csv"))
# Also save without WINDOWED suffix for backward compatibility
readr::write_csv(seasonal_windowed_with_static,
                 file.path(out_dir, "HJA_master_seasonal.csv"))

# MASTER ANNUAL (annual aggregates, WINDOWED metrics + static)
annual_windowed_with_static <- annual_windowed %>%
  left_join(static_site, by = "Stream_Name")

readr::write_csv(annual_windowed_with_static,
                 file.path(out_dir, "HJA_master_annual_WINDOWED.csv"))
# Also save without WINDOWED suffix for backward compatibility
readr::write_csv(annual_windowed_with_static,
                 file.path(out_dir, "HJA_master_annual.csv"))

# MASTER SITE MEANS (site-level, WINDOWED metrics + static + synchrony)
master_site_windowed <- static_site %>%
  left_join(site_means_windowed, by = "Stream_Name")

readr::write_csv(master_site_windowed,
                 file.path(out_dir, "HJA_master_site_means_WINDOWED.csv"))
# Also save without WINDOWED suffix for backward compatibility
readr::write_csv(master_site_windowed,
                 file.path(out_dir, "HJA_master_site_means.csv"))

# =============================================================================
# BUILD CLEAN DATASETS (focused on key predictors)
# =============================================================================

# CLEAN WINDOWS
clean_cols_windows <- c("Stream_Name","window_center","water_year","hydrologic_season",
                        "solute","data_source",
                        # CQ metrics (WINDOWED)
                        "cq_slope_windowed","cq_behavior","cq_CVc_CVq_windowed","cq_sync",
                        # Tier 1: Catchment characteristics
                        "Area_km2","Elevation_mean_m","Slope_mean",
                        "Lava1_per","Lava2_per","Ash_Per","Pyro_per",
                        "Age","Harvest","Landslide_Total",
                        # Tier 2: Isotope metrics
                        "DR_Overall","MTT_final","FYw_final",
                        # Clusters
                        "Cluster_oldSolute","Cluster_wy","dist_oldRef","Cluster_mode_wy")

clean_windows <- mega_window_df %>% select(any_of(clean_cols_windows))
readr::write_csv(clean_windows, file.path(out_dir, "HJA_clean_windows.csv"))

# CLEAN TEMPORAL (seasonal/annual)
clean_cols_temporal <- c("Stream_Name","solute","water_year","hydrologic_season",
                         "data_source","temporal_scale",
                         # CQ metrics (WINDOWED)
                         "cq_slope_windowed","cq_CVc_CVq_windowed",
                         # Catchment characteristics
                         "Area_km2","Elevation_mean_m","Slope_mean",
                         "Lava1_per","Lava2_per","Ash_Per","Pyro_per",
                         "Age","Harvest","Landslide_Total",
                         # Tier 2
                         "DR_Overall","MTT_final","FYw_final",
                         # Clusters
                         "Cluster_oldSolute","Cluster_wy","dist_oldRef","Cluster_mode_wy")

clean_seasonal <- seasonal_windowed_with_static %>% select(any_of(clean_cols_temporal))
clean_annual <- annual_windowed_with_static %>% select(any_of(clean_cols_temporal))

readr::write_csv(clean_seasonal, file.path(out_dir, "HJA_clean_seasonal_WINDOWED.csv"))
readr::write_csv(clean_annual, file.path(out_dir, "HJA_clean_annual_WINDOWED.csv"))


# CLEAN SITE MEANS
clean_cols_site <- c("Stream_Name","solute","data_source","temporal_scale",
                     # CQ metrics (WINDOWED)
                     "cq_slope_windowed","cq_CVc_CVq_windowed",
                     # Synchrony metrics (WINDOWED, aggregated)
                     "conc_sync_allpairs","conc_sync_outlet",
                     "cqslope_sync_allpairs","cqslope_sync_outlet",
                     "wymore_crosssite_allpairs","wymore_crosssite_outlet",
                     "wymore_cvcq_consistency",
                     # Tier 1: Catchment characteristics
                     "Area_km2","Elevation_mean_m","Slope_mean",
                     "Lava1_per","Lava2_per","Ash_Per","Pyro_per",
                     "Age","Harvest","Landslide_Total",
                     # Tier 2: Isotope metrics
                     "DR_Overall","MTT_final","FYw_final",
                     # Real averaged metrics (stream temp, etc.)
                     "JulyM_ST_mean","JulyM_AT_mean","JST_AT_mean",
                     # Clusters
                     "Cluster_oldSolute","Cluster_mode_wy")

clean_site_means <- master_site_windowed %>% select(any_of(clean_cols_site))
readr::write_csv(clean_site_means, file.path(out_dir, "HJA_clean_site_means_WINDOWED.csv"))
