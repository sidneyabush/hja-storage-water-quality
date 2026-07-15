# =============================================================================
# HYDRO SEASONS VALIDATION & DIAGNOSTICS
# 1b_validate_hydro_seasons.R
# =============================================================================
# Purpose: Comprehensive validation of my season detection algorithm
# - Sensitivity analysis on key parameters (peak_factor, peak_abs_min)
# - Visual diagnostics (plot detected boundaries on actual discharge)
# - Assertions to catch logical inconsistencies
# - Comparison of first-peak vs alternative definitions
#
# Run this AFTER 1b_define_hydro_seasons.R produces baseline seasons
# =============================================================================

library(tidyverse)
library(lubridate)
library(assertthat)

rm(list = ls())

# Paths
data_path <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality/data"
output_dir <- file.path(dirname(data_path), "outputs")
diag_dir <- file.path(output_dir, "season_diagnostics")
dir.create(diag_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# SECTION 1: LOAD BASELINE DATA
# =============================================================================
message("Loading baseline season definitions...\n")

season_boundaries <- read_csv(
  file.path(output_dir, "season_boundaries.csv"),
  show_col_types = FALSE
)

Q_with_seasons <- read_csv(
  file.path(output_dir, "HJA_daily_Q_with_seasons.csv"),
  show_col_types = FALSE
)

# Load raw discharge for reference
raw_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality/raw_data"

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

discharge_file <- first_existing(
  file.path(raw_dir, c("HF00402_v14.csv", "HF00402_v15.csv")),
  "HF00402 discharge file"
)

discharge_raw <- read_csv(
  discharge_file,
  show_col_types = FALSE
) %>%
  mutate(
    date = parse_hja_date(DATE),
    month = month(date)
  ) %>%
  select(Stream_Name = SITECODE, date, Qcms = MEAN_Q) %>%
  filter(!is.na(Qcms)) %>%
  mutate(
    Stream_Name = case_when(
      Stream_Name == "GSLOOK_FULL" ~ "GSLOOK",
      Stream_Name %in% c("GSWSMC", "GSWSMC_FULL") ~ "GSMACK",
      TRUE ~ Stream_Name
    ),
    water_year = if_else(month(date) >= 10, year(date) + 1L, year(date))
  )

# Basin median for season detection comparison
global_daily <- discharge_raw %>%
  group_by(date) %>%
  summarise(Qcms = median(Qcms, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    month = month(date),
    water_year = if_else(month >= 10, year(date) + 1L, year(date))
  )

message("✓ Loaded baseline data\n")
message("  - Season boundaries:", nrow(season_boundaries), "water years\n")
message("  - Q with seasons:", nrow(Q_with_seasons), "daily records\n")

# =============================================================================
# SECTION 2: ASSERTIONS & DATA QUALITY
# =============================================================================
message("\nRunning assertions on baseline seasons...\n")

# 2.1 Chronological order
sb <- season_boundaries %>% arrange(water_year)
assert_that(
  all(sb$wet_start_date >= sb$min_Q_date, na.rm = TRUE),
  msg = "ERROR: wet_start_date before min_Q_date"
)
message("  ✓ wet_start_date >= min_Q_date\n")

assert_that(
  all(sb$peak_Q_date >= sb$wet_start_date, na.rm = TRUE),
  msg = "ERROR: peak_Q_date before wet_start_date"
)
message("  ✓ peak_Q_date >= wet_start_date\n")

assert_that(
  all(sb$wet_end_date >= sb$peak_Q_date, na.rm = TRUE),
  msg = "ERROR: wet_end_date before peak_Q_date"
)
message("  ✓ wet_end_date >= peak_Q_date\n")

# 2.2 Gap between water years
sb_with_lag <- sb %>%
  arrange(water_year) %>%
  mutate(
    next_wet_start = lead(wet_start_date),
    gap_days = as.numeric(next_wet_start - wet_end_date)
  )

min_gap <- min(sb_with_lag$gap_days, na.rm = TRUE)
message("  ✓ Minimum dry gap:", min_gap, "days\n")

if (min_gap < 30) {
  warning("WARNING: At least one year has <30 day dry gap. Check manually.")
}

# 2.3 No missing values in key columns
assert_that(
  !any(is.na(season_boundaries$water_year)),
  msg = "water_year has NAs"
)
message("  ✓ No missing water_years\n")

na_summary <- season_boundaries %>%
  summarise(
    across(everything(), ~sum(is.na(.x)), .names = "na_{.col}")
  ) %>%
  pivot_longer(everything()) %>%
  filter(value > 0)

if (nrow(na_summary) > 0) {
  message("  ⚠ Missing values detected:\n")
  print(na_summary)
} else {
  message("  ✓ No missing values in season boundaries\n")
}

# =============================================================================
# SECTION 3: SENSITIVITY ANALYSIS
# =============================================================================
message("\n" %+% strrep("=", 70) %+% "\n")
message("SENSITIVITY ANALYSIS: Wet-Start Parameters\n")
message(strrep("=", 70) %+% "\n")

# Re-implement season detection with parameter flexibility
detect_seasons_param <- function(df, peak_factor, peak_abs_min,
                                 extreme_window_days = 14,
                                 q_window_days = 14,
                                 slope_window_days = 21,
                                 drying_slope_threshold = -0.003) {
  
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
        water_year = yr, min_Q_date = as.Date(NA), wet_start_date = as.Date(NA),
        peak_Q_date = as.Date(NA), wet_end_date = as.Date(NA),
        min_Q_value = NA_real_, peak_Q_value = NA_real_
      )
      next
    }
    
    sub <- sub %>%
      mutate(
        Q_roll_min = zoo::rollapply(Qcms, extreme_window_days, min, fill = NA, align = "center"),
        Q_roll_max = zoo::rollapply(Qcms, extreme_window_days, max, fill = NA, align = "center"),
        Q_smooth = zoo::rollapply(Qcms, q_window_days, mean, fill = NA, align = "right"),
        Q_slope = c(
          rep(NA_real_, slope_window_days),
          (Q_smooth[(1 + slope_window_days):n()] - Q_smooth[1:(n() - slope_window_days)]) / slope_window_days
        )
      )
    
    # Min Q
    min_Q <- sub %>%
      filter(date >= as.Date(paste0(yr - 1, "-07-01")),
             date <= as.Date(paste0(yr - 1, "-11-30")),
             !is.na(Q_roll_min), Qcms == Q_roll_min) %>%
      arrange(Qcms, date) %>%
      slice_head(n = 1)
    
    if (nrow(min_Q) == 0) {
      out[[yr]] <- tibble(
        water_year = yr, min_Q_date = as.Date(NA), wet_start_date = as.Date(NA),
        peak_Q_date = as.Date(NA), wet_end_date = as.Date(NA),
        min_Q_value = NA_real_, peak_Q_value = NA_real_
      )
      next
    }
    
    min_Q_date <- min_Q$date
    min_Q_value <- min_Q$Qcms
    
    # Wet start with parameters
    aug1 <- as.Date(paste0(yr - 1, "-08-01"))
    wet_window <- sub %>% filter(date >= aug1) %>% arrange(date)
    wet_start_date <- min_Q_date
    
    if (nrow(wet_window) > 0) {
      signif_thresh <- max(peak_factor * min_Q_value, peak_abs_min)
      candidate_peaks <- wet_window %>%
        filter(!is.na(Q_roll_max), Qcms == Q_roll_max, Qcms >= signif_thresh) %>%
        arrange(date)
      
      if (nrow(candidate_peaks) > 0) {
        wet_start_date <- candidate_peaks$date[1]
      }
    }
    
    # Peak
    peak <- sub %>%
      filter(date >= wy_start, date <= as.Date(paste0(yr, "-06-30"))) %>%
      arrange(desc(Qcms), date) %>%
      slice_head(n = 1)
    
    if (nrow(peak) == 0) {
      out[[yr]] <- tibble(
        water_year = yr, min_Q_date = min_Q_date, wet_start_date = wet_start_date,
        peak_Q_date = as.Date(NA), wet_end_date = as.Date(NA),
        min_Q_value = min_Q_value, peak_Q_value = NA_real_
      )
      next
    }
    
    peak_Q_date <- peak$date
    peak_Q_value <- peak$Qcms
    
    # Wet end (simplified - just find first drying)
    drying_month <- ifelse(peak_Q_value > 11, 6, 5)
    post_peak <- sub %>%
      filter(date > peak_Q_date, month(date) >= drying_month, date <= wy_end)
    
    end_threshold_absolute <- max(min_Q_value * 20, 0.5)
    wet_end_date <- wy_end
    
    if (nrow(post_peak) > 0) {
      candidate_ends <- post_peak %>%
        filter(!is.na(Q_slope), Q_slope <= drying_slope_threshold,
               Qcms < end_threshold_absolute) %>%
        arrange(date)
      
      if (nrow(candidate_ends) > 0) {
        wet_end_date <- candidate_ends$date[1]
      }
    }
    
    out[[yr]] <- tibble(
      water_year = yr, min_Q_date = min_Q_date, wet_start_date = wet_start_date,
      peak_Q_date = peak_Q_date, wet_end_date = wet_end_date,
      min_Q_value = min_Q_value, peak_Q_value = peak_Q_value
    )
  }
  
  bind_rows(out)
}

# Run sensitivity analysis
message("\nTesting parameter combinations...\n")

sensitivity_grid <- expand_grid(
  peak_factor = c(1.5, 1.75, 2.0, 2.25, 2.5),
  peak_abs_min = c(0.2, 0.3, 0.4, 0.5)
)

sensitivity_results <- sensitivity_grid %>%
  mutate(
    seasons = pmap(list(peak_factor, peak_abs_min),
                   ~detect_seasons_param(global_daily, .x, .y)),
    n_detected = map_int(seasons, ~nrow(filter(.x, !is.na(wet_start_date))))
  )

# Summary
message("\nSensitivity Summary:\n")
print(sensitivity_results %>% select(peak_factor, peak_abs_min, n_detected))

# Calculate shifts
baseline_starts <- season_boundaries %>%
  filter(!is.na(wet_start_date)) %>%
  select(water_year, baseline_wet_start = wet_start_date)

shift_analysis <- sensitivity_results %>%
  filter(peak_factor != 2.0 | peak_abs_min != 0.3) %>%
  mutate(
    shifts = map2(seasons, peak_factor, function(s, pf) {
      s %>%
        inner_join(baseline_starts, by = "water_year") %>%
        filter(!is.na(wet_start_date)) %>%
        mutate(
          shift_days = as.numeric(wet_start_date - baseline_wet_start)
        )
    })
  ) %>%
  tidyr::unnest(shifts, names_repair = "minimal")

# Statistics
message("\nWet-start date shifts from baseline (peak_factor=2.0, peak_abs_min=0.3):\n")
shift_summary <- shift_analysis %>%
  group_by(peak_factor, peak_abs_min) %>%
  summarise(
    mean_shift = mean(shift_days, na.rm = TRUE),
    sd_shift = sd(shift_days, na.rm = TRUE),
    min_shift = min(shift_days, na.rm = TRUE),
    max_shift = max(shift_days, na.rm = TRUE),
    .groups = "drop"
  )

print(shift_summary)

# Save sensitivity results
write_csv(shift_summary,
          file.path(diag_dir, "sensitivity_analysis_summary.csv"))

message("\n✓ Sensitivity analysis complete\n")

# =============================================================================
# SECTION 4: VISUAL DIAGNOSTICS
# =============================================================================
message("\n" %+% strrep("=", 70) %+% "\n")
message("GENERATING VISUAL DIAGNOSTICS\n")
message(strrep("=", 70) %+% "\n")

# 4.1 Plot season boundaries for each site
message("\nGenerating per-site boundary plots...\n")

all_sites <- discharge_raw %>%
  distinct(Stream_Name) %>%
  filter(!Stream_Name %in% c("GSWSMA", "GSWSMF")) %>%
  pull(Stream_Name)

for (site in all_sites) {
  
  site_data <- global_daily %>%
    filter(year(date) >= 2000) %>%  # Show recent years only for clarity
    mutate(site = "Basin Median")
  
  site_seasons <- season_boundaries %>%
    filter(water_year >= 2000)
  
  # Create plot
  p <- ggplot(site_data, aes(x = date, y = Qcms)) +
    geom_line(color = "steelblue", linewidth = 0.5) +
    
    # Wet/dry transitions
    geom_vline(aes(xintercept = wet_start_date), data = site_seasons,
               color = "green", linetype = "dashed", linewidth = 0.8, alpha = 0.7) +
    geom_vline(aes(xintercept = wet_end_date), data = site_seasons,
               color = "orange", linetype = "dashed", linewidth = 0.8, alpha = 0.7) +
    
    # Peaks and minima
    geom_point(aes(x = peak_Q_date, y = peak_Q_value), data = site_seasons,
               color = "red", size = 2, shape = 17) +
    geom_point(aes(x = min_Q_date, y = min_Q_value), data = site_seasons,
               color = "blue", size = 2, shape = 16) +
    
    scale_y_log10() +
    labs(
      title = paste("Hydro Season Boundaries:", site),
      subtitle = "Green = wet start, Orange = wet end, Red △ = peak, Blue ● = summer min",
      x = "Date",
      y = "Discharge (mm/d, log scale)"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(
    file.path(diag_dir, paste0("boundaries_", site, ".png")),
    p, width = 14, height = 6, dpi = 100
  )
  
  message("  ✓", site, "\n")
}

# 4.2 Distribution of key dates
message("\nGenerating date distribution plots...\n")

p_dists <- season_boundaries %>%
  filter(!is.na(wet_start_date), !is.na(wet_end_date)) %>%
  pivot_longer(cols = c(min_Q_date, wet_start_date, peak_Q_date, wet_end_date),
               names_to = "event", values_to = "date") %>%
  mutate(
    event = factor(event, levels = c("min_Q_date", "wet_start_date", "peak_Q_date", "wet_end_date"),
                   labels = c("Summer Min", "Wet Start", "Peak", "Wet End")),
    doy = yday(date)
  ) %>%
  ggplot(aes(x = doy, fill = event)) +
  geom_histogram(binwidth = 10, color = "black", alpha = 0.7) +
  facet_wrap(~event, nrow = 2) +
  scale_x_continuous(breaks = c(1, 92, 183, 274, 365),
                     labels = c("Jan 1", "Apr 1", "Jul 1", "Oct 1", "Dec 31")) +
  labs(
    title = "Distribution of Hydrologic Events (All Water Years)",
    x = "Day of Year",
    y = "Count"
  ) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(
  file.path(diag_dir, "event_distributions.png"),
  p_dists, width = 12, height = 8, dpi = 100
)

message("  ✓ Event distributions\n")

# 4.3 Time series of event dates
p_timeline <- season_boundaries %>%
  filter(!is.na(wet_start_date)) %>%
  select(water_year, wet_start_date, peak_Q_date, wet_end_date) %>%
  pivot_longer(cols = -water_year, names_to = "event", values_to = "date") %>%
  mutate(
    event = factor(event, levels = c("wet_start_date", "peak_Q_date", "wet_end_date"),
                   labels = c("Wet Start", "Peak", "Wet End"))
  ) %>%
  ggplot(aes(x = water_year, y = yday(date), color = event, shape = event)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_line(alpha = 0.5) +
  scale_y_continuous(breaks = c(1, 92, 183, 274, 365),
                     labels = c("Jan 1", "Apr 1", "Jul 1", "Oct 1", "Dec 31")) +
  labs(
    title = "Timeline of Key Events by Water Year",
    x = "Water Year",
    y = "Day of Year",
    color = "Event",
    shape = "Event"
  ) +
  theme_bw()

ggsave(
  file.path(diag_dir, "timeline_events.png"),
  p_timeline, width = 12, height = 6, dpi = 100
)

message("  ✓ Timeline of events\n")

# =============================================================================
# SECTION 5: FIRST PEAK VS HIGHEST PEAK COMPARISON
# =============================================================================
message("\n" %+% strrep("=", 70) %+% "\n")
message("ANALYZING: First Peak vs Highest Peak Logic\n")
message(strrep("=", 70) %+% "\n")

# Using baseline parameters (peak_factor=2.0, peak_abs_min=0.3)
current_params <- season_boundaries %>%
  filter(!is.na(min_Q_value)) %>%
  mutate(signif_thresh = pmax(2.0 * min_Q_value, 0.3))

# Get all peaks and compare
peak_comparison <- list()

for (yr in unique(global_daily$water_year)) {
  
  aug1 <- as.Date(paste0(yr - 1, "-08-01"))
  wy_end <- as.Date(paste0(yr, "-09-30"))
  
  year_data <- global_daily %>%
    filter(date >= aug1, date <= wy_end) %>%
    arrange(date) %>%
    mutate(
      Q_roll_max = zoo::rollapply(Qcms, 14, max, fill = NA, align = "center")
    )
  
  # Get baseline min_Q for this year
  min_q_val <- current_params %>%
    filter(water_year == yr) %>%
    pull(min_Q_value)
  
  if (length(min_q_val) == 0 || is.na(min_q_val)) next
  
  signif_thresh <- max(2.0 * min_q_val, 0.3)
  
  # Find all significant peaks
  all_peaks <- year_data %>%
    filter(!is.na(Q_roll_max), Qcms == Q_roll_max, Qcms >= signif_thresh) %>%
    arrange(date) %>%
    mutate(
      peak_type = if_else(row_number() == 1, "First", "Subsequent"),
      water_year = yr
    ) %>%
    select(water_year, date, Qcms, peak_type)
  
  if (nrow(all_peaks) > 0) {
    peak_comparison[[yr]] <- all_peaks
  }
}

peak_comp_df <- bind_rows(peak_comparison)

if (nrow(peak_comp_df) > 0) {
  
  message("\nPeak Detection Analysis:\n")
  message("  Total water years with peaks:", n_distinct(peak_comp_df$water_year), "\n")
  multi_peak_n <- peak_comp_df %>%
    group_by(water_year) %>%
    filter(n() > 1) %>%
    ungroup() %>%
    distinct(water_year) %>%
    nrow()
  message("  Water years with >1 peak:", multi_peak_n, "\n")
  
  # Multi-peak years
  multi_peak_years <- peak_comp_df %>%
    group_by(water_year) %>%
    filter(n() > 1) %>%
    ungroup() %>%
    arrange(water_year, date)
  
  if (nrow(multi_peak_years) > 0) {
    message("\n  Water years with multiple significant peaks:\n")
    print(multi_peak_years %>%
      select(water_year, date, Qcms, peak_type) %>%
      arrange(water_year, date))
    
    # Save for review
    write_csv(multi_peak_years,
              file.path(diag_dir, "multi_peak_years.csv"))
    
    message("\n  ⚠ Consider if 'first peak' misses true seasonal transitions\n")
  }
}

# =============================================================================
# SECTION 6: SUMMARY REPORT
# =============================================================================
message("\n" %+% strrep("=", 70) %+% "\n")
message("VALIDATION SUMMARY\n")
message(strrep("=", 70) %+% "\n")

summary_stats <- season_boundaries %>%
  filter(!is.na(min_Q_value)) %>%
  summarise(
    n_years = n(),
    min_Q_mean = mean(min_Q_value),
    min_Q_sd = sd(min_Q_value),
    peak_Q_mean = mean(peak_Q_value),
    peak_Q_sd = sd(peak_Q_value),
    wet_duration_mean = mean(as.numeric(wet_end_date - wet_start_date), na.rm = TRUE),
    wet_duration_sd = sd(as.numeric(wet_end_date - wet_start_date), na.rm = TRUE)
  )

message("\nBaseline Season Statistics:\n")
print(summary_stats)

message("\n" %+% strrep("=", 70) %+% "\n")
message("DIAGNOSTIC FILES SAVED TO:\n")
message(diag_dir, "\n")
message(strrep("=", 70) %+% "\n\n")

# List all outputs
message("Output files:\n")
diag_files <- list.files(diag_dir, full.names = FALSE)
for (f in diag_files) {
  message("  •", f, "\n")
}

message("\n✓ Validation complete\n")
