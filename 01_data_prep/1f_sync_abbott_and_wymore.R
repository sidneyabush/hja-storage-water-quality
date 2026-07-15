# =============================================================================
# 1f: Abbott + Wymore Synchrony Calculations
# =============================================================================
# Computes pairwise synchrony using Abbott et al. and Wymore et al. methods
# 
# OUTPUT METRICS:
#   - Abbott concentration synchrony: "Do sites have high/low concentrations together?"
#   - Wymore CQ-slope synchrony: "Do sites mobilize solutes the same way?"
#
# DEPRECATED (removed Dec 2025):
#   - Abbott CQ-slope synchrony: Anti-correlated with other sync metrics (r=-0.14)
#
# Outputs: HJA_Abbott_synchrony_windows.csv, HJA_wymore_crosssite_sync.csv,
#          HJA_outlet_synchrony_site_level.csv (NEW - for predictive modeling)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(data.table)
})

rm(list = ls())

# Source shared workflow settings and plot preferences
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
source(file.path(repo_dir, "00_helpers", "workflow_config.R"))
source(file.path(repo_dir, "00_helpers", "plot_prefs.R"))

# Output to Box, not repo
base_dir    <- "/Users/sidneybush/Library/CloudStorage/Box-Box"
project_dir <- file.path(base_dir, "Sidney_Bush", "HJA_Water_Quality")
out_dir     <- file.path(project_dir, "outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

WINDOW_DAYS <- 90
HALF_WIN    <- WINDOW_DAYS %/% 2
MIN_PAIRED_WINDOWS <- 2
outlet_sites <- c("GSLOOK", "LO5")

# Helper function for pairwise synchrony (UNIVARIATE - one solute at a time)
compute_pairwise_synchrony <- function(df, value_col, time_scale = c("annual","seasonal")) {
  value_col  <- rlang::ensym(value_col)
  time_scale <- match.arg(time_scale)
  if (time_scale == "annual") {
    grouped <- df %>% group_by(solute, water_year)
  } else {
    grouped <- df %>% filter(!is.na(hydrologic_season)) %>% group_by(solute, water_year, hydrologic_season)
  }
  out <- grouped %>% group_modify(~ {
    chunk <- .x
    sites <- sort(unique(chunk$Stream_Name))
    if (length(sites) < 2) return(tibble())
    pairs <- t(combn(sites, 2)) %>% as_tibble() %>% rename(Stream1 = V1, Stream2 = V2)
    purrr::map_dfr(seq_len(nrow(pairs)), function(i) {
      s1 <- pairs$Stream1[i]; s2 <- pairs$Stream2[i]
      x_df <- chunk %>% filter(Stream_Name == s1) %>% select(window_center, value = !!value_col)
      y_df <- chunk %>% filter(Stream_Name == s2) %>% select(window_center, value = !!value_col)
      xy <- inner_join(x_df, y_df, by = "window_center", suffix = c("_1","_2"))
      if (nrow(xy) < MIN_PAIRED_WINDOWS) return(tibble(Stream1=s1,Stream2=s2,n_windows_pair=nrow(xy),synchrony=NA_real_,pearson_r=NA_real_))
      x_vals <- xy$value_1; y_vals <- xy$value_2
      if (!any(is.finite(x_vals)) || !any(is.finite(y_vals))) return(tibble(Stream1=s1,Stream2=s2,n_windows_pair=nrow(xy),synchrony=NA_real_,pearson_r=NA_real_))
      mx <- mean(x_vals, na.rm = TRUE); my <- mean(y_vals, na.rm = TRUE)
      dev_x <- x_vals - mx; dev_y <- y_vals - my
      prod_abs <- abs(dev_x * dev_y)
      synchrony_val <- sum(prod_abs, na.rm = TRUE) / (length(prod_abs) - 1)
      pearson_r <- suppressWarnings(cor(x_vals, y_vals, use = "complete.obs"))
      tibble(Stream1=s1,Stream2=s2,n_windows_pair=length(prod_abs),synchrony=synchrony_val,pearson_r=pearson_r)
    })
  }) %>% ungroup() %>% mutate(time_scale = time_scale)
  if (!"hydrologic_season" %in% names(out)) out <- out %>% mutate(hydrologic_season = NA_character_)
  out
}

# Helper function for MULTIVARIATE synchrony (all solutes or subset together)
compute_multivariate_synchrony <- function(df, value_col, solute_subset = NULL, time_scale = c("annual","seasonal")) {
  value_col  <- rlang::ensym(value_col)
  time_scale <- match.arg(time_scale)

  # Filter to solute subset if provided
  if (!is.null(solute_subset)) {
    df <- df %>% filter(solute %in% solute_subset)
  }

  if (time_scale == "annual") {
    grouped <- df %>% group_by(water_year)
  } else {
    grouped <- df %>% filter(!is.na(hydrologic_season)) %>% group_by(water_year, hydrologic_season)
  }

  out <- grouped %>% group_modify(~ {
    chunk <- .x
    # Get all site-solute-window combinations
    sites <- sort(unique(chunk$Stream_Name))
    if (length(sites) < 2) return(tibble())

    pairs <- t(combn(sites, 2)) %>% as_tibble() %>% rename(Stream1 = V1, Stream2 = V2)

    purrr::map_dfr(seq_len(nrow(pairs)), function(i) {
      s1 <- pairs$Stream1[i]; s2 <- pairs$Stream2[i]

      # Get data for both sites (all solutes)
      x_df <- chunk %>% filter(Stream_Name == s1) %>%
        select(solute, window_center, value = !!value_col)
      y_df <- chunk %>% filter(Stream_Name == s2) %>%
        select(solute, window_center, value = !!value_col)

      # Join on both solute AND window_center to align multivariate observations
      xy <- inner_join(x_df, y_df, by = c("solute", "window_center"), suffix = c("_1","_2"))

      if (nrow(xy) < MIN_PAIRED_WINDOWS) {
        return(tibble(Stream1=s1, Stream2=s2, n_windows_pair=0, n_solutes=0,
                     synchrony=NA_real_, pearson_r=NA_real_))
      }

      # Group by window_center to calculate multivariate synchrony
      window_sync <- xy %>%
        group_by(window_center) %>%
        summarise(
          n_solutes_win = n(),
          # For each window, calculate sum of abs(dev_x * dev_y) across all solutes
          sync_win = {
            if (n() < 2) {
              NA_real_
            } else {
              mx <- mean(value_1, na.rm = TRUE)
              my <- mean(value_2, na.rm = TRUE)
              dev_x <- value_1 - mx
              dev_y <- value_2 - my
              sum(abs(dev_x * dev_y), na.rm = TRUE)
            }
          },
          .groups = "drop"
        ) %>%
        filter(is.finite(sync_win))

      if (nrow(window_sync) < MIN_PAIRED_WINDOWS) {
        return(tibble(Stream1=s1, Stream2=s2, n_windows_pair=nrow(window_sync),
                     n_solutes=length(unique(xy$solute)),
                     synchrony=NA_real_, pearson_r=NA_real_))
      }

      # Average across windows
      synchrony_val <- mean(window_sync$sync_win, na.rm = TRUE)

      # Calculate multivariate correlation (average of univariate correlations)
      solute_cors <- xy %>%
        group_by(solute) %>%
        summarise(cor_val = suppressWarnings(cor(value_1, value_2, use = "complete.obs")), .groups = "drop")
      pearson_r <- mean(solute_cors$cor_val, na.rm = TRUE)

      tibble(Stream1=s1, Stream2=s2,
             n_windows_pair=nrow(window_sync),
             n_solutes=length(unique(xy$solute)),
             synchrony=synchrony_val,
             pearson_r=pearson_r)
    })
  }) %>% ungroup() %>% mutate(time_scale = time_scale)

  if (!"hydrologic_season" %in% names(out)) out <- out %>% mutate(hydrologic_season = NA_character_)
  out
}

# Load CQ master
cq_master <- readr::read_csv(file.path(out_dir, "HJA_CQ_master.csv"), show_col_types = FALSE) %>%
  rename(date = Date, Qcms = Q_cms, solute = variable) %>%
  mutate(date = as.Date(date), Stream_Name = as.character(Stream_Name), solute = as.character(solute))

# Prepare concentration samples
chem_samples <- cq_master %>% filter(!is.na(value)) %>% select(Stream_Name, solute, date, concentration = value)
conc_stats <- chem_samples %>% group_by(solute) %>% summarise(global_mean = mean(concentration, na.rm = TRUE), global_sd = sd(concentration, na.rm = TRUE), .groups = "drop")
chem_scaled <- chem_samples %>% left_join(conc_stats, by = "solute") %>% mutate(conc_scaled_global = (concentration - global_mean) / global_sd) %>% filter(is.finite(conc_scaled_global))

# Load rolling window CQ results (from 1c)
roll_windows <- readr::read_csv(file.path(out_dir, "windows_wet75_dry150", "roll_pts_keep_windows_wet75_dry150.csv"), show_col_types = FALSE) %>%
  rename(solute = variable, cq_slope = slope) %>%
  mutate(window_center = as.Date(window_center), Stream_Name = as.character(Stream_Name), solute = as.character(solute))

# Load seasons (from 1b) and join
seasons <- readr::read_csv(file.path(out_dir, "HJA_daily_Q_with_seasons.csv"), show_col_types = FALSE) %>%
  mutate(date = as.Date(date), Stream_Name = as.character(Stream_Name)) %>%
  select(Stream_Name, date, water_year, hydrologic_season)

# Join windows with seasons
clean_windows <- roll_windows %>%
  left_join(seasons, by = c("Stream_Name", "window_center" = "date")) %>%
  mutate(hydrologic_season = ifelse(is.na(hydrologic_season), NA_character_, as.character(hydrologic_season)))

# APPLY TEMPORAL FILTER FROM WORKFLOW CONFIG
clean_windows <- clean_windows %>%
  filter(!is.na(water_year), water_year >= ANALYSIS_YEAR_START, water_year <= ANALYSIS_YEAR_END)

windows_core <- clean_windows %>% select(Stream_Name, solute, window_center, water_year, hydrologic_season, cq_slope) %>% filter(!is.na(solute), !is.na(window_center))

# Build window intervals
windows_intervals <- windows_core %>% distinct(Stream_Name, solute, window_center, water_year, hydrologic_season) %>% 
  mutate(window_start = window_center - days(HALF_WIN), window_end = window_center + days(HALF_WIN))

# Join samples to windows using data.table overlap join
samples_dt <- as.data.table(chem_scaled %>% select(Stream_Name, solute, date, conc_scaled_global))
windows_dt <- as.data.table(windows_intervals)
samples_dt[, `:=`(start = date, end = date)]
windows_dt[, `:=`(start = window_start, end = window_end)]
setkey(samples_dt, Stream_Name, solute, start, end)
setkey(windows_dt, Stream_Name, solute, start, end)
conc_win_joined <- foverlaps(samples_dt, windows_dt, by.x = c("Stream_Name","solute","start","end"), by.y = c("Stream_Name","solute","start","end"), nomatch = 0L)

conc_windows <- conc_win_joined %>% as_tibble() %>% 
  transmute(Stream_Name, solute, window_center, water_year, hydrologic_season, conc_scaled_global) %>%
  group_by(Stream_Name, solute, window_center, water_year, hydrologic_season) %>% 
  summarise(conc_scaled_win = mean(conc_scaled_global, na.rm = TRUE), n_samples_win = n(), .groups = "drop")

# Prepare slope windows
slope_windows <- windows_core %>% filter(!is.na(cq_slope)) %>% select(Stream_Name, solute, window_center, water_year, hydrologic_season, cq_slope)
slope_stats <- slope_windows %>% group_by(solute) %>% summarise(slope_mean = mean(cq_slope, na.rm = TRUE), slope_sd = sd(cq_slope, na.rm = TRUE), .groups = "drop")
slope_windows_scaled <- slope_windows %>% left_join(slope_stats, by = "solute") %>% mutate(cq_slope_scaled = (cq_slope - slope_mean) / slope_sd) %>% filter(is.finite(cq_slope_scaled))

# =============================================================================
# UNIVARIATE ABBOTT SYNCHRONY (one solute at a time)
# =============================================================================
# NOTE: Abbott CQ-slope synchrony is DEPRECATED - it was anti-correlated with
#       other sync metrics (r=-0.14 with conc sync) and measured something different
sync_conc_annual <- compute_pairwise_synchrony(conc_windows, conc_scaled_win, time_scale = "annual") %>% mutate(synchrony_type = "concentration")
sync_conc_seasonal <- compute_pairwise_synchrony(conc_windows, conc_scaled_win, time_scale = "seasonal") %>% mutate(synchrony_type = "concentration")

# Combine concentration synchrony results (CQ-slope sync removed)
synchrony_all <- bind_rows(sync_conc_annual, sync_conc_seasonal) %>%
  mutate(time_scale = factor(time_scale, levels = c("annual","seasonal")),
         solute = as.character(solute), Stream1 = as.character(Stream1), Stream2 = as.character(Stream2),
         is_outlet_pair = Stream1 %in% outlet_sites | Stream2 %in% outlet_sites) %>%
  relocate(solute, synchrony_type, time_scale, water_year, hydrologic_season, Stream1, Stream2, is_outlet_pair, n_windows_pair, synchrony, pearson_r) %>%
  arrange(synchrony_type, solute, time_scale, water_year, hydrologic_season, Stream1, Stream2)

readr::write_csv(synchrony_all, file.path(out_dir, "HJA_Abbott_synchrony_windows.csv"))

# =============================================================================
# MULTIVARIATE ABBOTT SYNCHRONY (all solutes or subsets together)
# =============================================================================
# Solute groups defined in plot_prefs.R (GEOGENIC_SOLUTES, BIOGENIC_SOLUTES, NUTRIENT_SOLUTES)

# 1. All solutes together
sync_mv_all_annual <- compute_multivariate_synchrony(
  conc_windows, conc_scaled_win, solute_subset = NULL, time_scale = "annual"
) %>% mutate(synchrony_type = "concentration_mv_all")

# 2. Biogenic solutes only (DOC, NH3, NO3, PO4)
sync_mv_bio_annual <- compute_multivariate_synchrony(
  conc_windows, conc_scaled_win, solute_subset = BIOGENIC_SOLUTES, time_scale = "annual"
) %>% mutate(synchrony_type = "concentration_mv_biogenic")

# 3. Geogenic solutes only (Ca, Mg, Na, K, Cl, SO4)
sync_mv_geo_annual <- compute_multivariate_synchrony(
  conc_windows, conc_scaled_win, solute_subset = GEOGENIC_SOLUTES, time_scale = "annual"
) %>% mutate(synchrony_type = "concentration_mv_geogenic")

# 4. Nutrient solutes only (DSi - weathering-derived but biologically mediated)
sync_mv_nutrient_annual <- compute_multivariate_synchrony(
  conc_windows, conc_scaled_win, solute_subset = NUTRIENT_SOLUTES, time_scale = "annual"
) %>% mutate(synchrony_type = "concentration_mv_nutrient")

# Combine multivariate results
synchrony_mv <- bind_rows(sync_mv_all_annual, sync_mv_bio_annual, sync_mv_geo_annual, sync_mv_nutrient_annual) %>%
  mutate(time_scale = "annual",
         Stream1 = as.character(Stream1),
         Stream2 = as.character(Stream2),
         is_outlet_pair = Stream1 %in% outlet_sites | Stream2 %in% outlet_sites) %>%
  relocate(synchrony_type, time_scale, water_year, Stream1, Stream2, is_outlet_pair,
           n_windows_pair, n_solutes, synchrony, pearson_r) %>%
  arrange(synchrony_type, water_year, Stream1, Stream2)

readr::write_csv(synchrony_mv, file.path(out_dir, "HJA_Abbott_synchrony_multivariate.csv"))

# =============================================================================
# UNIVARIATE WYMORE SYNCHRONY (one solute at a time)
# =============================================================================
wymore_crosssite <- slope_windows %>% select(Stream_Name, solute, window_center, water_year, hydrologic_season, cq_slope) %>% filter(!is.na(cq_slope))

wymore_pairs <- wymore_crosssite %>% group_by(solute, window_center, water_year, hydrologic_season) %>% group_modify(~ {
  chunk <- .x
  sites <- sort(unique(chunk$Stream_Name))
  if (length(sites) < 2) return(tibble())
  pairs <- t(combn(sites, 2)) %>% as_tibble() %>% rename(Stream1 = V1, Stream2 = V2)
  purrr::map_dfr(seq_len(nrow(pairs)), function(i) {
    s1 <- pairs$Stream1[i]; s2 <- pairs$Stream2[i]
    slope1 <- chunk$cq_slope[chunk$Stream_Name == s1]
    slope2 <- chunk$cq_slope[chunk$Stream_Name == s2]
    if (length(slope1) == 0 || length(slope2) == 0) return(tibble())
    slope1 <- slope1[1]; slope2 <- slope2[1]
    quadrant <- case_when(
      slope1 > 0 & slope2 > 0 ~ "Q1", 
      slope1 > 0 & slope2 < 0 ~ "Q2", 
      slope1 < 0 & slope2 < 0 ~ "Q3", 
      slope1 < 0 & slope2 > 0 ~ "Q4", 
      TRUE ~ NA_character_
    )
    # Calculate continuous similarity metrics
    slope_similarity = 1 - abs(slope1 - slope2) / (abs(slope1) + abs(slope2) + 1e-6)  # Ranges 0-1, higher = more similar
    slope_agreement = sign(slope1) == sign(slope2)  # Same direction (both +ve or both -ve)
    tibble(Stream1=s1, Stream2=s2,
           slope1=slope1, slope2=slope2,
           slope_similarity=slope_similarity,
           slope_agreement=slope_agreement,
           quadrant=quadrant,
           sync=quadrant %in% c("Q1","Q3"))
  })
}) %>% ungroup() %>% 
  mutate(is_outlet_pair = Stream1 %in% outlet_sites | Stream2 %in% outlet_sites) %>%
  relocate(solute, window_center, water_year, hydrologic_season, Stream1, Stream2, is_outlet_pair,
           slope1, slope2, slope_similarity, slope_agreement, quadrant, sync) %>%
  arrange(solute, window_center, Stream1, Stream2)

readr::write_csv(wymore_pairs, file.path(out_dir, "HJA_wymore_crosssite_sync.csv"))

# =============================================================================
# MULTIVARIATE WYMORE SYNCHRONY (all solutes or subsets together)
# =============================================================================
# Solute groups defined in plot_prefs.R (BIO_SOLUTES, GEO_SOLUTES)

# Function to calculate multivariate Wymore sync
compute_wymore_multivariate <- function(slope_df, solute_subset = NULL) {
  # Filter to solute subset if provided
  if (!is.null(solute_subset)) {
    slope_df <- slope_df %>% filter(solute %in% solute_subset)
  }

  # Group by time window (not solute) and calculate proportion synchronized
  slope_df %>%
    group_by(window_center, water_year, hydrologic_season) %>%
    group_modify(~ {
      chunk <- .x
      sites <- sort(unique(chunk$Stream_Name))
      if (length(sites) < 2) return(tibble())

      pairs <- t(combn(sites, 2)) %>% as_tibble() %>% rename(Stream1 = V1, Stream2 = V2)

      purrr::map_dfr(seq_len(nrow(pairs)), function(i) {
        s1 <- pairs$Stream1[i]; s2 <- pairs$Stream2[i]

        # Get all solutes for both sites at this window
        slopes1 <- chunk %>% filter(Stream_Name == s1) %>% select(solute, cq_slope)
        slopes2 <- chunk %>% filter(Stream_Name == s2) %>% select(solute, cq_slope)

        # Join to get paired slopes
        slopes_paired <- inner_join(slopes1, slopes2, by = "solute", suffix = c("_1", "_2"))

        if (nrow(slopes_paired) < 2) {
          return(tibble(Stream1=s1, Stream2=s2, n_solutes=nrow(slopes_paired),
                       prop_sync=NA_real_, mean_similarity=NA_real_,
                       prop_Q1=NA_real_, prop_Q3=NA_real_, dominant_quadrant=NA_character_))
        }

        # Calculate for each solute: are they synchronized (Q1 or Q3)?
        slopes_paired <- slopes_paired %>%
          mutate(
            quadrant = case_when(
              cq_slope_1 > 0 & cq_slope_2 > 0 ~ "Q1",
              cq_slope_1 > 0 & cq_slope_2 < 0 ~ "Q2",
              cq_slope_1 < 0 & cq_slope_2 < 0 ~ "Q3",
              cq_slope_1 < 0 & cq_slope_2 > 0 ~ "Q4",
              TRUE ~ NA_character_
            ),
            sync = quadrant %in% c("Q1", "Q3"),
            similarity = 1 - abs(cq_slope_1 - cq_slope_2) / (abs(cq_slope_1) + abs(cq_slope_2) + 1e-6)
          )

        # Calculate proportion synchronized across all solutes
        prop_sync <- mean(slopes_paired$sync, na.rm = TRUE)
        mean_similarity <- mean(slopes_paired$similarity, na.rm = TRUE)

        # Track Q1 vs Q3 separately
        prop_Q1 <- mean(slopes_paired$quadrant == "Q1", na.rm = TRUE)
        prop_Q3 <- mean(slopes_paired$quadrant == "Q3", na.rm = TRUE)

        # Dominant quadrant (Q1 = mobilizing, Q3 = diluting)
        dominant_quad <- ifelse(prop_Q1 > prop_Q3, "Q1",
                               ifelse(prop_Q3 > prop_Q1, "Q3", "Mixed"))

        tibble(Stream1=s1, Stream2=s2, n_solutes=nrow(slopes_paired),
               prop_sync=prop_sync, mean_similarity=mean_similarity,
               prop_Q1=prop_Q1, prop_Q3=prop_Q3, dominant_quadrant=dominant_quad)
      })
    }) %>%
    ungroup()
}

# 1. All solutes together
wymore_mv_all <- compute_wymore_multivariate(wymore_crosssite, solute_subset = NULL) %>%
  mutate(synchrony_type = "wymore_mv_all")

# 2. Biogenic solutes only (DOC, NH3, NO3, PO4)
wymore_mv_bio <- compute_wymore_multivariate(wymore_crosssite, solute_subset = BIOGENIC_SOLUTES) %>%
  mutate(synchrony_type = "wymore_mv_biogenic")

# 3. Geogenic solutes only (Ca, Mg, Na, K, Cl, SO4)
wymore_mv_geo <- compute_wymore_multivariate(wymore_crosssite, solute_subset = GEOGENIC_SOLUTES) %>%
  mutate(synchrony_type = "wymore_mv_geogenic")

# 4. Nutrient solutes only (DSi - weathering-derived but biologically mediated)
wymore_mv_nutrient <- compute_wymore_multivariate(wymore_crosssite, solute_subset = NUTRIENT_SOLUTES) %>%
  mutate(synchrony_type = "wymore_mv_nutrient")

# Combine multivariate Wymore results
wymore_mv <- bind_rows(wymore_mv_all, wymore_mv_bio, wymore_mv_geo, wymore_mv_nutrient) %>%
  mutate(Stream1 = as.character(Stream1),
         Stream2 = as.character(Stream2),
         is_outlet_pair = Stream1 %in% outlet_sites | Stream2 %in% outlet_sites) %>%
  relocate(synchrony_type, window_center, water_year, hydrologic_season,
           Stream1, Stream2, is_outlet_pair, n_solutes, prop_sync, mean_similarity) %>%
  arrange(synchrony_type, window_center, Stream1, Stream2)

readr::write_csv(wymore_mv, file.path(out_dir, "HJA_wymore_synchrony_multivariate.csv"))

# =============================================================================
# OUTLET-SPECIFIC SYNCHRONY METRICS (for predictive modeling)
# =============================================================================
# Create site-level metrics: sync of each site with GSLOOK outlet

# Abbott concentration sync with outlet (site-annual level)
abbott_outlet <- synchrony_all %>%
  filter(is_outlet_pair, time_scale == "annual") %>%
  mutate(
    headwater_site = case_when(
      Stream1 %in% outlet_sites ~ Stream2,
      Stream2 %in% outlet_sites ~ Stream1,
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(headwater_site)) %>%
  group_by(headwater_site, solute, water_year) %>%
  summarise(
    conc_sync_outlet = mean(synchrony, na.rm = TRUE),
    conc_pearson_outlet = mean(pearson_r, na.rm = TRUE),
    n_windows = sum(n_windows_pair, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(Stream_Name = headwater_site)

# Wymore CQ-slope sync with outlet (site-annual level)
wymore_outlet <- wymore_pairs %>%
  filter(is_outlet_pair) %>%
  mutate(
    headwater_site = case_when(
      Stream1 %in% outlet_sites ~ Stream2,
      Stream2 %in% outlet_sites ~ Stream1,
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(headwater_site)) %>%
  group_by(headwater_site, solute, water_year) %>%
  summarise(
    cqslope_sync_outlet = mean(sync, na.rm = TRUE),  # Proportion of windows in sync (windowed)
    cqslope_similarity_outlet = mean(slope_similarity, na.rm = TRUE),  # Mean slope similarity (windowed, 0-1)
    cqslope_agreement_outlet = mean(slope_agreement, na.rm = TRUE),  # Proportion with same sign (windowed)
    n_windows_wymore = n(),
    .groups = "drop"
  ) %>%
  rename(Stream_Name = headwater_site)

# Combine outlet sync metrics
outlet_sync_annual <- abbott_outlet %>%
  full_join(wymore_outlet, by = c("Stream_Name", "solute", "water_year"))

readr::write_csv(outlet_sync_annual, file.path(out_dir, "HJA_outlet_synchrony_annual.csv"))

# Create site-level summary (for modeling with catchment traits)
outlet_sync_site <- outlet_sync_annual %>%
  group_by(Stream_Name, solute) %>%
  summarise(
    conc_sync_outlet_mean = mean(conc_sync_outlet, na.rm = TRUE),
    conc_sync_outlet_sd = sd(conc_sync_outlet, na.rm = TRUE),
    conc_pearson_outlet_mean = mean(conc_pearson_outlet, na.rm = TRUE),
    cqslope_sync_outlet_mean = mean(cqslope_sync_outlet, na.rm = TRUE),
    cqslope_sync_outlet_sd = sd(cqslope_sync_outlet, na.rm = TRUE),
    n_years = n_distinct(water_year),
    .groups = "drop"
  )

readr::write_csv(outlet_sync_site, file.path(out_dir, "HJA_outlet_synchrony_site_level.csv"))
