# =============================================================================
# 1g: Composite Synchrony
# =============================================================================
# Aggregates Abbott and Wymore synchrony into long-term composite metrics
# Outputs: HJA_composite_synchrony_annual.csv, HJA_composite_synchrony.csv
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

rm(list = ls())

# Output to Box, not repo
repo_dir <- Sys.getenv(
  "HJA_WQ_REPO_DIR",
  unset = "/Users/sidneybush/Documents/GitHub/hja-water-quality"
)
source(file.path(repo_dir, "00_helpers", "workflow_config.R"))
paths <- get_project_paths()
out_dir <- paths$out_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

standardize_stream_ids <- function(df) {
  df %>% mutate(across(any_of(c("Stream_Name","Stream1","Stream2","site","Site")), 
                       ~ case_when(. %in% c("GSWSMC","GSMACK") ~ "GSMACK", TRUE ~ as.character(.))))
}

message("[1g] Composite synchrony...")

# Load Abbott synchrony
sync_raw <- readr::read_csv(file.path(out_dir, "HJA_Abbott_synchrony_windows.csv"), show_col_types = FALSE) %>%
  standardize_stream_ids() %>% 
  mutate(solute = as.character(solute), Stream1 = as.character(Stream1), Stream2 = as.character(Stream2))

outlet <- "GSLOOK"

# Annual synchrony
sync_annual <- sync_raw %>% filter(time_scale == "annual")
sync_annual_long <- sync_annual %>% 
  select(solute, water_year, synchrony_type, Stream1, Stream2, synchrony, pearson_r) %>%
  pivot_longer(cols = c(Stream1, Stream2), names_to = "which_stream", values_to = "Stream_Name")

annual_allpairs <- sync_annual_long %>% 
  group_by(Stream_Name, solute, water_year, synchrony_type) %>% 
  summarise(sync_allpairs = mean(synchrony, na.rm = TRUE), 
            r_allpairs = mean(pearson_r, na.rm = TRUE), .groups = "drop")

annual_allpairs_wide <- annual_allpairs %>% 
  mutate(sync_type = case_when(synchrony_type == "concentration" ~ "conc", 
                               synchrony_type == "cq_slope" ~ "cqslope", 
                               TRUE ~ synchrony_type)) %>%
  select(-synchrony_type) %>% 
  pivot_wider(names_from = sync_type, values_from = c(sync_allpairs, r_allpairs), 
              names_glue = "{sync_type}_{.value}")

sync_annual_outlet <- sync_annual %>% 
  filter(Stream1 == outlet | Stream2 == outlet) %>% 
  mutate(Stream_Name = if_else(Stream1 == outlet, Stream2, Stream1)) %>%
  group_by(Stream_Name, solute, water_year, synchrony_type) %>% 
  summarise(sync_outlet = mean(synchrony, na.rm = TRUE), 
            r_outlet = mean(pearson_r, na.rm = TRUE), .groups = "drop")

annual_outlet_wide <- sync_annual_outlet %>% 
  mutate(sync_type = case_when(synchrony_type == "concentration" ~ "conc", 
                               synchrony_type == "cq_slope" ~ "cqslope", 
                               TRUE ~ synchrony_type)) %>%
  select(-synchrony_type) %>% 
  pivot_wider(names_from = sync_type, values_from = c(sync_outlet, r_outlet), 
              names_glue = "{sync_type}_{.value}")

composite_annual <- annual_allpairs_wide %>% 
  full_join(annual_outlet_wide, by = c("Stream_Name","solute","water_year")) %>% 
  arrange(Stream_Name, solute, water_year)

readr::write_csv(composite_annual, file.path(out_dir, "HJA_composite_synchrony_annual.csv"))

# Long-term synchrony
expected_cols <- c("conc_sync_allpairs","conc_r_allpairs","conc_sync_outlet","conc_r_outlet",
                   "cqslope_sync_allpairs","cqslope_r_allpairs","cqslope_sync_outlet","cqslope_r_outlet")
for (nm in expected_cols) if (!nm %in% names(composite_annual)) composite_annual[[nm]] <- NA_real_

longterm_sync <- composite_annual %>% 
  group_by(Stream_Name, solute) %>% 
  summarise(
    conc_sync_allpairs    = mean(conc_sync_allpairs,    na.rm = TRUE),
    conc_r_allpairs       = mean(conc_r_allpairs,       na.rm = TRUE),
    conc_sync_outlet      = mean(conc_sync_outlet,      na.rm = TRUE),
    conc_r_outlet         = mean(conc_r_outlet,         na.rm = TRUE),
    cqslope_sync_allpairs = mean(cqslope_sync_allpairs, na.rm = TRUE),
    cqslope_r_allpairs    = mean(cqslope_r_allpairs,    na.rm = TRUE),
    cqslope_sync_outlet   = mean(cqslope_sync_outlet,   na.rm = TRUE),
    cqslope_r_outlet      = mean(cqslope_r_outlet,      na.rm = TRUE),
    n_years_sync          = dplyr::n_distinct(water_year[is.finite(conc_sync_allpairs) | is.finite(cqslope_sync_allpairs)]),
    .groups = "drop"
  )

# Load Wymore quadrant data
cq_quad <- readr::read_csv(file.path(out_dir, "CQ_rolling_window_results.csv"), show_col_types = FALSE) %>%
  filter(comparison_type == "cqslope_CVcCVq") %>% 
  transmute(Stream_Name = Stream_Name, solute = solute1, wymore_quad = quadrant, wymore_sync = sync)

wymore_cvcq_summary <- cq_quad %>% 
  filter(!is.na(wymore_sync)) %>% 
  group_by(Stream_Name, solute) %>% 
  summarise(wymore_cvcq_consistency = mean(wymore_sync == "sync", na.rm = TRUE), .groups = "drop")

# Load Wymore cross-site synchrony
wymore_crosssite <- readr::read_csv(file.path(out_dir, "HJA_wymore_crosssite_sync.csv"), show_col_types = FALSE) %>%
  standardize_stream_ids()

# Active annual pair table used by the storage-metric synchrony analysis.
# This replaces the old static HJA_pair_sync_metrics.csv with a regenerated
# product from current Abbott concentration synchrony and Wymore C-Q agreement.
abbott_pair_annual <- sync_raw %>%
  filter(time_scale == "annual", synchrony_type == "concentration") %>%
  transmute(
    solute,
    time_scale = "annual",
    water_year = as.integer(water_year),
    hydrologic_season = "annual",
    Stream1,
    Stream2,
    n_windows_pair_abb = n_windows_pair,
    sync_magnitude = synchrony,
    Abbott_S = synchrony,
    sync_direction = pearson_r,
    is_outlet_pair
  )

wymore_pair_annual <- wymore_crosssite %>%
  filter(!is.na(water_year), !is.na(sync)) %>%
  mutate(
    water_year = as.integer(water_year),
    quadrant = as.character(quadrant)
  ) %>%
  group_by(solute, water_year, Stream1, Stream2, is_outlet_pair) %>%
  summarise(
    n_windows_pair_wym = n(),
    n_I = sum(quadrant == "Q1", na.rm = TRUE),
    n_II = sum(quadrant == "Q2", na.rm = TRUE),
    n_III = sum(quadrant == "Q3", na.rm = TRUE),
    n_IV = sum(quadrant == "Q4", na.rm = TRUE),
    prop_I = mean(quadrant == "Q1", na.rm = TRUE),
    prop_II = mean(quadrant == "Q2", na.rm = TRUE),
    prop_III = mean(quadrant == "Q3", na.rm = TRUE),
    prop_IV = mean(quadrant == "Q4", na.rm = TRUE),
    prop_sync_wymore = mean(sync, na.rm = TRUE),
    prop_async_wymore = mean(!sync, na.rm = TRUE),
    prop_dual_mobilizing = prop_I,
    prop_dual_diluting = prop_III,
    prop_mixed = prop_II + prop_IV,
    behavior_bias = prop_dual_mobilizing - prop_dual_diluting,
    dominant_sync_type = c("dual_mobilizing", "mixed_Q2", "dual_diluting", "mixed_Q4")[
      which.max(c(prop_I, prop_II, prop_III, prop_IV))
    ],
    .groups = "drop"
  ) %>%
  mutate(
    time_scale = "annual",
    hydrologic_season = "annual"
  )

pair_sync_metrics <- full_join(
  abbott_pair_annual,
  wymore_pair_annual,
  by = c("solute", "time_scale", "water_year", "hydrologic_season", "Stream1", "Stream2", "is_outlet_pair")
) %>%
  mutate(comparison_group = "annual_all_windows", r_slope = NA_real_) %>%
  select(
    solute,
    time_scale,
    water_year,
    hydrologic_season,
    Stream1,
    Stream2,
    comparison_group,
    n_windows_pair_wym,
    n_windows_pair_abb,
    sync_magnitude,
    Abbott_S,
    sync_direction,
    r_slope,
    n_I,
    n_II,
    n_III,
    n_IV,
    prop_I,
    prop_II,
    prop_III,
    prop_IV,
    prop_sync_wymore,
    prop_async_wymore,
    prop_dual_mobilizing,
    prop_dual_diluting,
    prop_mixed,
    behavior_bias,
    dominant_sync_type,
    is_outlet_pair
  ) %>%
  arrange(solute, water_year, Stream1, Stream2)

readr::write_csv(pair_sync_metrics, file.path(out_dir, "HJA_pair_sync_metrics.csv"))

wymore_crosssite_long <- wymore_crosssite %>%
  select(solute, window_center, water_year, hydrologic_season, Stream1, Stream2, is_outlet_pair,
         quadrant, sync) %>%
  pivot_longer(cols = c(Stream1, Stream2), names_to = "which_stream", values_to = "Stream_Name")

# Calculate overall sync (Q1 + Q3) and separate Q1 vs Q3
wymore_crosssite_allpairs <- wymore_crosssite_long %>%
  filter(!is.na(sync)) %>%
  group_by(Stream_Name, solute) %>%
  summarise(
    wymore_crosssite_allpairs = mean(sync, na.rm = TRUE),  # Proportion synchronized (Q1 or Q3)
    wymore_Q1_proportion = mean(quadrant == "Q1", na.rm = TRUE),  # Dual mobilizing
    wymore_Q3_proportion = mean(quadrant == "Q3", na.rm = TRUE),  # Dual diluting
    wymore_dominant_quadrant = names(sort(table(quadrant), decreasing = TRUE))[1],  # Most common
    .groups = "drop"
  )

wymore_crosssite_outlet <- wymore_crosssite_long %>%
  filter(!is.na(sync), is_outlet_pair) %>%
  group_by(Stream_Name, solute) %>%
  summarise(
    wymore_crosssite_outlet = mean(sync, na.rm = TRUE),
    wymore_Q1_proportion_outlet = mean(quadrant == "Q1", na.rm = TRUE),
    wymore_Q3_proportion_outlet = mean(quadrant == "Q3", na.rm = TRUE),
    .groups = "drop"
  )

# Combine all UNIVARIATE synchrony metrics (site-solute level)
composite_longterm <- longterm_sync %>%
  left_join(wymore_cvcq_summary, by = c("Stream_Name","solute")) %>%
  left_join(wymore_crosssite_allpairs, by = c("Stream_Name","solute")) %>%
  left_join(wymore_crosssite_outlet, by = c("Stream_Name","solute")) %>%
  arrange(Stream_Name, solute)

readr::write_csv(composite_longterm, file.path(out_dir, "HJA_composite_synchrony.csv"))

# =============================================================================
# MULTIVARIATE SYNCHRONY AGGREGATION (site-pair level, no solute dimension)
# =============================================================================
# Load multivariate Abbott synchrony (pairwise, annual)
sync_mv_raw <- readr::read_csv(file.path(out_dir, "HJA_Abbott_synchrony_multivariate.csv"), show_col_types = FALSE) %>%
  standardize_stream_ids() %>%
  mutate(synchrony_type = as.character(synchrony_type),
         Stream1 = as.character(Stream1),
         Stream2 = as.character(Stream2))

# Aggregate to site-pair level (average across years)
# For each site, calculate average sync with all other sites
sync_mv_annual <- sync_mv_raw %>%
  select(synchrony_type, water_year, Stream1, Stream2, synchrony, pearson_r)

# Pivot longer to get one row per site per year per sync type
sync_mv_annual_long <- sync_mv_annual %>%
  pivot_longer(cols = c(Stream1, Stream2), names_to = "which_stream", values_to = "Stream_Name")

# Calculate average synchrony with all other sites
sync_mv_allpairs <- sync_mv_annual_long %>%
  group_by(Stream_Name, synchrony_type, water_year) %>%
  summarise(
    sync_allpairs = mean(synchrony, na.rm = TRUE),
    r_allpairs = mean(pearson_r, na.rm = TRUE),
    .groups = "drop"
  )

# Pivot wider to get separate columns for each synchrony type
sync_mv_allpairs_wide <- sync_mv_allpairs %>%
  pivot_wider(
    names_from = synchrony_type,
    values_from = c(sync_allpairs, r_allpairs),
    names_glue = "{synchrony_type}_{.value}"
  )

# Calculate outlet-specific multivariate synchrony
sync_mv_outlet <- sync_mv_annual %>%
  filter(Stream1 == outlet | Stream2 == outlet) %>%
  mutate(Stream_Name = if_else(Stream1 == outlet, Stream2, Stream1)) %>%
  group_by(Stream_Name, synchrony_type, water_year) %>%
  summarise(
    sync_outlet = mean(synchrony, na.rm = TRUE),
    r_outlet = mean(pearson_r, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = synchrony_type,
    values_from = c(sync_outlet, r_outlet),
    names_glue = "{synchrony_type}_{.value}"
  )

# Combine annual multivariate synchrony
composite_mv_annual <- sync_mv_allpairs_wide %>%
  full_join(sync_mv_outlet, by = c("Stream_Name", "water_year")) %>%
  arrange(Stream_Name, water_year)

readr::write_csv(composite_mv_annual, file.path(out_dir, "HJA_composite_synchrony_multivariate_annual.csv"))

# Long-term site-level multivariate synchrony (average across years)
composite_mv_longterm <- composite_mv_annual %>%
  group_by(Stream_Name) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

readr::write_csv(composite_mv_longterm, file.path(out_dir, "HJA_composite_synchrony_multivariate.csv"))

# =============================================================================
# MULTIVARIATE WYMORE SYNCHRONY AGGREGATION
# =============================================================================
# Load multivariate Wymore synchrony (pairwise, window-level)
wymore_mv_raw <- readr::read_csv(file.path(out_dir, "HJA_wymore_synchrony_multivariate.csv"), show_col_types = FALSE) %>%
  standardize_stream_ids() %>%
  mutate(synchrony_type = as.character(synchrony_type),
         Stream1 = as.character(Stream1),
         Stream2 = as.character(Stream2))

# Aggregate to annual level (average across windows within year)
wymore_mv_annual <- wymore_mv_raw %>%
  group_by(synchrony_type, water_year, Stream1, Stream2, is_outlet_pair) %>%
  summarise(
    prop_sync = mean(prop_sync, na.rm = TRUE),
    mean_similarity = mean(mean_similarity, na.rm = TRUE),
    prop_Q1 = mean(prop_Q1, na.rm = TRUE),
    prop_Q3 = mean(prop_Q3, na.rm = TRUE),
    n_windows = n(),
    .groups = "drop"
  )

# Pivot longer to get one row per site per year per sync type
wymore_mv_annual_long <- wymore_mv_annual %>%
  pivot_longer(cols = c(Stream1, Stream2), names_to = "which_stream", values_to = "Stream_Name")

# Calculate average synchrony with all other sites (allpairs)
wymore_mv_allpairs <- wymore_mv_annual_long %>%
  group_by(Stream_Name, synchrony_type, water_year) %>%
  summarise(
    prop_sync_allpairs = mean(prop_sync, na.rm = TRUE),
    similarity_allpairs = mean(mean_similarity, na.rm = TRUE),
    prop_Q1_allpairs = mean(prop_Q1, na.rm = TRUE),
    prop_Q3_allpairs = mean(prop_Q3, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # Dominant quadrant for this site-year-synctype
    dominant_quadrant_allpairs = case_when(
      prop_Q1_allpairs > prop_Q3_allpairs ~ "Dual Mobilizing (Q1)",
      prop_Q3_allpairs > prop_Q1_allpairs ~ "Dual Diluting (Q3)",
      TRUE ~ "Mixed"
    )
  )

# Pivot wider
wymore_mv_allpairs_wide <- wymore_mv_allpairs %>%
  pivot_wider(
    names_from = synchrony_type,
    values_from = c(prop_sync_allpairs, similarity_allpairs, prop_Q1_allpairs, prop_Q3_allpairs, dominant_quadrant_allpairs),
    names_glue = "{synchrony_type}_{.value}"
  )

# Calculate outlet-specific multivariate Wymore synchrony
wymore_mv_outlet_annual <- wymore_mv_annual %>%
  filter(Stream1 == outlet | Stream2 == outlet) %>%
  mutate(Stream_Name = if_else(Stream1 == outlet, Stream2, Stream1)) %>%
  group_by(Stream_Name, synchrony_type, water_year) %>%
  summarise(
    prop_sync_outlet = mean(prop_sync, na.rm = TRUE),
    similarity_outlet = mean(mean_similarity, na.rm = TRUE),
    prop_Q1_outlet = mean(prop_Q1, na.rm = TRUE),
    prop_Q3_outlet = mean(prop_Q3, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    dominant_quadrant_outlet = case_when(
      prop_Q1_outlet > prop_Q3_outlet ~ "Dual Mobilizing (Q1)",
      prop_Q3_outlet > prop_Q1_outlet ~ "Dual Diluting (Q3)",
      TRUE ~ "Mixed"
    )
  ) %>%
  pivot_wider(
    names_from = synchrony_type,
    values_from = c(prop_sync_outlet, similarity_outlet, prop_Q1_outlet, prop_Q3_outlet, dominant_quadrant_outlet),
    names_glue = "{synchrony_type}_{.value}"
  )

# Combine annual Wymore multivariate
wymore_mv_composite_annual <- wymore_mv_allpairs_wide %>%
  full_join(wymore_mv_outlet_annual, by = c("Stream_Name", "water_year")) %>%
  arrange(Stream_Name, water_year)

# Add to existing Abbott multivariate annual
composite_mv_annual_combined <- composite_mv_annual %>%
  full_join(wymore_mv_composite_annual, by = c("Stream_Name", "water_year"))

readr::write_csv(composite_mv_annual_combined,
                 file.path(out_dir, "HJA_composite_synchrony_multivariate_annual.csv"))

# Long-term site-level (average across years)
composite_mv_longterm_combined <- composite_mv_annual_combined %>%
  group_by(Stream_Name) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

readr::write_csv(composite_mv_longterm_combined,
                 file.path(out_dir, "HJA_composite_synchrony_multivariate.csv"))
