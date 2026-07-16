# Fixed hydro-season sensitivity runner — avoids unnest duplicate-name error
library(tidyverse)
library(lubridate)
library(zoo)

data_path <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality/data"
raw_dir   <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality/raw_data"
output_dir <- file.path(dirname(data_path), "outputs")
diag_dir <- file.path(output_dir, "season_diagnostics")
dir.create(diag_dir, showWarnings = FALSE, recursive = TRUE)

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

# Load basin median daily discharge
discharge_file <- first_existing(
  file.path(raw_dir, c("HF00402_v14.csv", "HF00402_v15.csv")),
  "HF00402 discharge file"
)

discharge_raw <- read_csv(discharge_file, show_col_types = FALSE) %>%
  mutate(date = parse_hja_date(DATE), month = month(date)) %>%
  select(Stream_Name = SITECODE, date, Qcms = MEAN_Q) %>%
  filter(!is.na(Qcms)) %>%
  mutate(Stream_Name = case_when(
    Stream_Name == "GSLOOK_FULL" ~ "GSLOOK",
    Stream_Name %in% c("GSWSMC","GSWSMC_FULL") ~ "GSMACK",
    TRUE ~ Stream_Name
  ),
  water_year = if_else(month(date) >= 10, year(date) + 1L, year(date)))

global_daily <- discharge_raw %>%
  group_by(date) %>%
  summarise(Qcms = median(Qcms, na.rm = TRUE), .groups = "drop") %>%
  mutate(month = month(date), water_year = if_else(month >= 10, year(date) + 1L, year(date)))

# Load baseline season boundaries
baseline_file <- file.path(output_dir, "season_boundaries.csv")
if (!file.exists(baseline_file)) stop("Baseline season_boundaries.csv not found in outputs; run 1b_define_hydro_seasons.R first")
baseline <- read_csv(baseline_file, show_col_types = FALSE) %>%
  select(water_year, baseline_wet_start = wet_start_date)

# season detection function (parameterized)
detect_seasons_param <- function(df, peak_factor, peak_abs_min,
                                 extreme_window_days = 14,
                                 q_window_days = 14,
                                 slope_window_days = 21,
                                 drying_slope_threshold = -0.003) {
  yrs <- sort(unique(df$water_year))
  out <- vector("list", length(yrs)); names(out) <- yrs
  for (yr in yrs) {
    wy_start <- as.Date(paste0(yr - 1, "-10-01"))
    wy_end <- as.Date(paste0(yr, "-09-30"))
    data_start <- as.Date(paste0(yr - 1, "-07-01"))
    sub <- df %>% filter(date >= data_start, date <= wy_end) %>% arrange(date)
    if (nrow(sub) < 100) {
      out[[as.character(yr)]] <- tibble(water_year = yr, wet_start_date = as.Date(NA))
      next
    }
    sub <- sub %>% mutate(
      Q_roll_min = zoo::rollapply(Qcms, extreme_window_days, min, fill = NA, align = "center"),
      Q_roll_max = zoo::rollapply(Qcms, extreme_window_days, max, fill = NA, align = "center"),
      Q_smooth = zoo::rollapply(Qcms, q_window_days, mean, fill = NA, align = "right"),
      Q_slope = c(rep(NA_real_, slope_window_days),
                  (Q_smooth[(1 + slope_window_days):n()] - Q_smooth[1:(n() - slope_window_days)]) / slope_window_days)
    )
    # summer minimum
    min_Q <- sub %>% filter(date >= as.Date(paste0(yr - 1, "-07-01")), date <= as.Date(paste0(yr - 1, "-11-30")), !is.na(Q_roll_min), Qcms == Q_roll_min) %>% arrange(Qcms, date) %>% slice_head(n = 1)
    if (nrow(min_Q) == 0) { out[[as.character(yr)]] <- tibble(water_year = yr, wet_start_date = as.Date(NA)); next }
    min_Q_date <- min_Q$date; min_Q_value <- min_Q$Qcms
    # wet start
    aug1 <- as.Date(paste0(yr - 1, "-08-01"))
    wet_window <- sub %>% filter(date >= aug1) %>% arrange(date)
    wet_start_date <- min_Q_date
    if (nrow(wet_window) > 0) {
      signif_thresh <- max(peak_factor * min_Q_value, peak_abs_min)
      candidate_peaks <- wet_window %>% filter(!is.na(Q_roll_max), Qcms == Q_roll_max, Qcms >= signif_thresh) %>% arrange(date)
      if (nrow(candidate_peaks) > 0) wet_start_date <- candidate_peaks$date[1]
    }
    if (!is.na(wet_start_date) && !is.na(min_Q_date) && wet_start_date < min_Q_date) wet_start_date <- min_Q_date
    out[[as.character(yr)]] <- tibble(water_year = yr, wet_start_date = wet_start_date)
  }
  bind_rows(out)
}

# Sensitivity grid
sensitivity_grid <- expand_grid(peak_factor = c(1.5, 1.75, 2.0, 2.25, 2.5), peak_abs_min = c(0.2, 0.3, 0.4, 0.5))

results <- sensitivity_grid %>% mutate(
  seasons = pmap(list(peak_factor, peak_abs_min), ~detect_seasons_param(global_daily, ..1, ..2)),
  n_detected = map_int(seasons, ~sum(!is.na(.x$wet_start_date)))
)

# For each combo compute shift summary relative to baseline
summaries <- map2_dfr(results$seasons, seq_len(nrow(results)), function(seas, i) {
  pf <- results$peak_factor[i]; pa <- results$peak_abs_min[i]
  joined <- seas %>% inner_join(baseline, by = "water_year") %>% filter(!is.na(wet_start_date), !is.na(baseline_wet_start)) %>% mutate(shift_days = as.numeric(wet_start_date - baseline_wet_start))
  tib <- tibble(
    peak_factor = pf,
    peak_abs_min = pa,
    n_detected = nrow(seas %>% filter(!is.na(wet_start_date))),
    mean_shift = mean(joined$shift_days, na.rm = TRUE),
    sd_shift = sd(joined$shift_days, na.rm = TRUE),
    min_shift = ifelse(all(is.na(joined$shift_days)), NA_real_, min(joined$shift_days, na.rm = TRUE)),
    max_shift = ifelse(all(is.na(joined$shift_days)), NA_real_, max(joined$shift_days, na.rm = TRUE))
  )
  tib
})

write_csv(summaries, file.path(diag_dir, "sensitivity_analysis_summary_fixed.csv"))
message("Wrote sensitivity summary to:", file.path(diag_dir, "sensitivity_analysis_summary_fixed.csv"), "\n")

print(summaries)

message("Sensitivity run complete.\n")
