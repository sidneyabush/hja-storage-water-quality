# =============================================================================
# SEASON DEFINITION: 
# H.J. Andrews - Define seasons by hydrograph shape and discharge dynamics
# =============================================================================
# HYBRID APPROACH:
# - Uses Q extrema (min/max) to define SEARCH WINDOWS
# - Uses slope detection to identify PROCESS TRANSITIONS
#
# Wet season: From wetting start through peak until drying starts
# Dry season: From drying start through recession until next wetting start
#
# Wet start rule (UPDATED):
#   Wet season starts on the date of the FIRST SIGNIFICANT PEAK in Q
#   at or after August 1 of the PREVIOUS calendar year.
# =============================================================================

library(tidyverse)
library(lubridate)
library(zoo)

rm(list = ls())

# =============================================================================
# Paths
# =============================================================================
data_path <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality/data"
raw_dir   <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality/raw_data"

# outputs folder at same level as "data"
output_dir <- file.path(dirname(data_path), "outputs")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

first_existing <- function(paths, label) {
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0) {
    stop(label, " not found. Checked:\n", paste(paths, collapse = "\n"))
  }
  existing[[1]]
}

parse_hja_date <- function(x) {
  x_chr <- as.character(x)
  out <- suppressWarnings(lubridate::ymd(x_chr))
  missing <- is.na(out)
  if (any(missing)) {
    out[missing] <- suppressWarnings(lubridate::mdy(x_chr[missing]))
  }
  as.Date(out)
}

# =============================================================================
# Parameters
# =============================================================================

# Extrema & Slope windows
extreme_window_days <- 14
q_window_days       <- 14
slope_window_days   <- 21

# Peak search
winter_months             <- c(11, 12, 1, 2, 3, 4)
drying_search_start_month <- 5

# Wet-start detection:
# "Significant peak" = Q >= max(peak_factor * min_Q_value, peak_abs_min)
# Peak_factor = 2.0 means first peak that's 2x the summer minimum triggers wet season
# This uses "first significant peak" not "highest peak" - captures earliest wetting
# NOTE: These parameters are NOT validated against ground-truth observations
# See 1b_validate_hydro_seasons.R for sensitivity analysis
peak_factor  <- 2.0   # relative to minimum Q (validated in hydro_seasons_robustness_review.md)
peak_abs_min <- 0.3   # absolute floor for a peak (mm/d) - prevents noise in high-flow systems

# Slope thresholds
# Wetting slope threshold currently UNUSED - wet start defined by peak detection instead
# Could implement slope-based wetting detection in future iteration
wetting_slope_threshold <- 0.003  # LEFT FOR REFERENCE - not currently used
drying_slope_threshold  <- -0.003

# Enforce dry gap between water years
dry_buffer_days <- 90   # days before next wet start to force previous wet to end

# Filter incomplete site-years
min_days_per_year <- 300

# =============================================================================
# Load daily discharge
# =============================================================================
discharge_file <- first_existing(
  file.path(raw_dir, c("HF00402_v14.csv", "HF00402_v15.csv")),
  "HF00402 discharge file"
)

discharge_raw <- read_csv(discharge_file,
                          show_col_types = FALSE) %>%
  mutate(
    date  = parse_hja_date(DATE),
    month = month(date)
  ) %>%
  select(Stream_Name = SITECODE, date, month, Qcms = MEAN_Q) %>%
  filter(!is.na(Qcms)) %>%
  mutate(
    # *** IMPORTANT: normalize names to match CQ + DS workflows ***
    Stream_Name = case_when(
      Stream_Name == "GSLOOK_FULL"             ~ "GSLOOK",
      Stream_Name %in% c("GSWSMC","GSWSMC_FULL") ~ "GSMACK",
      TRUE                                     ~ Stream_Name
    ),
    water_year = if_else(month >= 10, year(date) + 1L, year(date))
  )

# =============================================================================
# Basin-median discharge for season detection
# =============================================================================
global_daily <- discharge_raw %>%
  group_by(date) %>%
  summarise(Qcms = median(Qcms, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    month      = month(date),
    water_year = if_else(month >= 10, year(date) + 1L, year(date))
  )

# =============================================================================
# Core season detection function
# =============================================================================
detect_seasons <- function(df) {
  
  yrs <- sort(unique(df$water_year))
  out <- vector("list", length(yrs)); names(out) <- yrs
  
  for (yr in yrs) {
    
    wy_start   <- as.Date(paste0(yr - 1, "-10-01"))
    wy_end     <- as.Date(paste0(yr,     "-09-30"))
    data_start <- as.Date(paste0(yr - 1, "-07-01"))
    
    sub <- df %>%
      filter(date >= data_start, date <= wy_end) %>%
      arrange(date)
    
    if (nrow(sub) < 100) {
      out[[yr]] <- tibble(
        water_year     = yr,
        min_Q_date     = as.Date(NA),
        wet_start_date = as.Date(NA),
        peak_Q_date    = as.Date(NA),
        wet_end_date   = as.Date(NA),
        min_Q_value    = NA_real_,
        peak_Q_value   = NA_real_
      )
      next
    }
    
    # Rolling stats
    sub <- sub %>%
      mutate(
        Q_roll_min = rollapply(Qcms, extreme_window_days, min,
                               fill = NA, align = "center"),
        Q_roll_max = rollapply(Qcms, extreme_window_days, max,
                               fill = NA, align = "center"),
        Q_smooth   = rollapply(Qcms, q_window_days, mean,
                               fill = NA, align = "right"),
        Q_slope = c(
          rep(NA_real_, slope_window_days),
          (Q_smooth[(1 + slope_window_days):n()] -
             Q_smooth[1:(n() - slope_window_days)]) / slope_window_days
        )
      )
    
    # ========================================================================
    # 1. Summer minimum Q (late-summer baseflow)
    # ========================================================================
    min_search_start <- as.Date(paste0(yr - 1, "-07-01"))
    min_search_end   <- as.Date(paste0(yr - 1, "-11-30"))
    
    min_Q <- sub %>%
      filter(date >= min_search_start,
             date <= min_search_end,
             !is.na(Q_roll_min),
             Qcms == Q_roll_min) %>%
      arrange(Qcms, date) %>%
      slice_head(n = 1)
    
    if (nrow(min_Q) == 0) {
      out[[yr]] <- tibble(
        water_year     = yr,
        min_Q_date     = as.Date(NA),
        wet_start_date = as.Date(NA),
        peak_Q_date    = as.Date(NA),
        wet_end_date   = as.Date(NA),
        min_Q_value    = NA_real_,
        peak_Q_value   = NA_real_
      )
      next
    }
    
    min_Q_date  <- min_Q$date
    min_Q_value <- min_Q$Qcms
    
    # ========================================================================
    # 2. Wet season start: FIRST SIGNIFICANT PEAK as early as Aug 1
    # ========================================================================
    aug1 <- as.Date(paste0(yr - 1, "-08-01"))
    
    wet_window <- sub %>%
      filter(date >= aug1) %>%
      arrange(date)
    
    # default fallback
    wet_start_date <- min_Q_date
    
    if (nrow(wet_window) > 0) {
      
      signif_thresh <- max(peak_factor * min_Q_value, peak_abs_min)
      
      candidate_peaks <- wet_window %>%
        filter(!is.na(Q_roll_max),
               Qcms == Q_roll_max,
               Qcms >= signif_thresh) %>%
        arrange(date)
      
      if (nrow(candidate_peaks) > 0) {
        # Rule: wet season starts at the FIRST significant peak
        wet_start_date <- candidate_peaks$date[1]
      }
    }
    
    # CONSTRAINT: wet_start_date must be >= min_Q_date (summer min occurs first)
    # This prevents illogical scenarios where peak is detected before minimum
    if (!is.na(wet_start_date) && !is.na(min_Q_date) && wet_start_date < min_Q_date) {
      wet_start_date <- min_Q_date
    }
    
    # ========================================================================
    # 3. Peak Q (Nov–Jun)
    # ========================================================================
    peak <- sub %>%
      filter(date >= wy_start,
             date <= as.Date(paste0(yr, "-06-30"))) %>%
      arrange(desc(Qcms), date) %>%
      slice_head(n = 1)
    
    if (nrow(peak) == 0) {
      out[[yr]] <- tibble(
        water_year     = yr,
        min_Q_date     = as.Date(min_Q_date),
        wet_start_date = as.Date(wet_start_date),
        peak_Q_date    = as.Date(NA),
        wet_end_date   = as.Date(NA),
        min_Q_value    = min_Q_value,
        peak_Q_value   = NA_real_
      )
      next
    }
    
    peak_Q_date  <- peak$date
    peak_Q_value <- peak$Qcms
    
    # ========================================================================
    # 4. Wet season end (start of consistent drying after peak)
    # ========================================================================
    drying_month <- ifelse(peak_Q_value > 11, 6, 5)
    
    post_peak <- sub %>%
      filter(date > peak_Q_date,
             month(date) >= drying_month,
             date <= wy_end)
    
    end_threshold_absolute    <- max(min_Q_value * 20, 0.5)
    resume_threshold_absolute <- max(min_Q_value * 30, 0.7)
    
    wet_end_date <- wy_end
    
    if (nrow(post_peak) > 0) {
      
      candidate_ends <- post_peak %>%
        filter(!is.na(Q_slope),
               Q_slope <= drying_slope_threshold,
               Qcms   <  end_threshold_absolute) %>%
        arrange(date)
      
      if (nrow(candidate_ends) > 0) {
        for (j in seq_len(nrow(candidate_ends))) {
          d_end <- candidate_ends$date[j]
          
          lookahead <- sub %>%
            filter(date > d_end,
                   date <= d_end + days(120),
                   date <= wy_end)
          
          if (nrow(lookahead) == 0) {
            wet_end_date <- d_end
            break
          }
          
          if (max(lookahead$Qcms, na.rm = TRUE) < resume_threshold_absolute) {
            wet_end_date <- d_end
            break
          }
        }
      }
    }
    
    # Save for this year
    out[[yr]] <- tibble(
      water_year     = yr,
      min_Q_date     = min_Q_date,
      wet_start_date = wet_start_date,
      peak_Q_date    = peak_Q_date,
      wet_end_date   = wet_end_date,
      min_Q_value    = min_Q_value,
      peak_Q_value   = peak_Q_value
    )
  }
  
  bind_rows(out)
}

# =============================================================================
# Run season detection
# =============================================================================
season_boundaries_raw <- detect_seasons(global_daily)

# VALIDATION: Check chronological consistency before enforcing dry gaps
# Peak should be between wet_start and wet_end
validation_chronology <- season_boundaries_raw %>%
  filter(!is.na(wet_start_date), !is.na(peak_Q_date), !is.na(wet_end_date)) %>%
  mutate(
    peak_before_start = peak_Q_date < wet_start_date,
    peak_after_end = peak_Q_date > wet_end_date
  ) %>%
  filter(peak_before_start | peak_after_end)

if (nrow(validation_chronology) > 0) {
  warning("Found ", nrow(validation_chronology), " years where peak_Q_date is outside wet season bounds")
  message("Check these water years:\n")
  print(validation_chronology %>% select(water_year, wet_start_date, peak_Q_date, wet_end_date))
}

# Chronology clean-up: enforce a dry gap between years
season_boundaries <- season_boundaries_raw %>%
  arrange(water_year) %>%
  mutate(next_wet_start = lead(wet_start_date)) %>%
  mutate(
    wet_end_date = if_else(
      !is.na(next_wet_start) &
        !is.na(wet_end_date) &
        wet_end_date >= next_wet_start,
      next_wet_start - days(dry_buffer_days),
      wet_end_date
    )
  ) %>%
  select(-next_wet_start)

# POST-DRY-GAP VALIDATION: Check that enforcing dry gap didn't create logical errors
validation_after_gap <- season_boundaries %>%
  filter(!is.na(wet_start_date), !is.na(peak_Q_date), !is.na(wet_end_date)) %>%
  mutate(
    peak_after_end = peak_Q_date > wet_end_date,
    wet_end_after_start = wet_end_date >= wet_start_date
  ) %>%
  filter(peak_after_end | !wet_end_after_start)

if (nrow(validation_after_gap) > 0) {
  warning("CRITICAL: Dry gap enforcement created logical errors")
  print(validation_after_gap %>% select(water_year, wet_start_date, peak_Q_date, wet_end_date))
}

# =============================================================================
# Build date → season map (Wet & Dry intervals)
# =============================================================================
season_intervals <- season_boundaries %>%
  arrange(water_year) %>%
  mutate(next_wet_start = lead(wet_start_date))

expand_season_dates <- function(df, start_col, end_col, season_label) {
  if (nrow(df) == 0) {
    return(tibble(date = as.Date(character()), hydrologic_season = character()))
  }
  purrr::map2_dfr(df[[start_col]], df[[end_col]], function(start_date, end_date) {
    if (is.na(start_date) || is.na(end_date) || start_date > end_date) {
      return(tibble(date = as.Date(character()), hydrologic_season = character()))
    }
    tibble(
      date = seq.Date(as.Date(start_date), as.Date(end_date), by = "day"),
      hydrologic_season = season_label
    )
  })
}

wet_dates <- season_intervals %>%
  filter(!is.na(wet_start_date), !is.na(wet_end_date)) %>%
  expand_season_dates("wet_start_date", "wet_end_date", "Wet")

dry_dates <- season_intervals %>%
  filter(!is.na(wet_end_date),
         !is.na(next_wet_start),
         wet_end_date < next_wet_start) %>%
  mutate(
    dry_start_date = wet_end_date + days(1),
    dry_end_date = next_wet_start - days(1)
  ) %>%
  expand_season_dates("dry_start_date", "dry_end_date", "Dry")

season_date_map <- bind_rows(wet_dates, dry_dates) %>%
  distinct(date, .keep_all = TRUE)

# =============================================================================
# Assign seasons to daily discharge (per site)
# =============================================================================
Q_with_seasons_Qextrema <- discharge_raw %>%
  left_join(season_date_map, by = "date")

# Drop incomplete site–years (per site)
site_year_full <- Q_with_seasons_Qextrema %>%
  group_by(Stream_Name, water_year) %>%
  summarise(
    n_days = n_distinct(date),
    .groups = "drop"
  ) %>%
  mutate(full_year = n_days >= min_days_per_year)

Q_with_seasons_Qextrema <- Q_with_seasons_Qextrema %>%
  inner_join(
    site_year_full %>% filter(full_year),
    by = c("Stream_Name", "water_year")
  )

# (Optional) drop GSWSMA / GSWSMF here to match other scripts
Q_with_seasons_Qextrema <- Q_with_seasons_Qextrema %>%
  filter(!Stream_Name %in% c("GSWSMA", "GSWSMF"))

# =============================================================================
# Save outputs
# =============================================================================
write_csv(Q_with_seasons_Qextrema,
          file.path(output_dir, "HJA_daily_Q_with_seasons.csv"))

write_csv(season_boundaries,
          file.path(output_dir, "season_boundaries.csv"))
