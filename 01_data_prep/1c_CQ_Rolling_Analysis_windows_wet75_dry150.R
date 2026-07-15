# Rolling window analysis - asymmetric windows (wet=75 days, dry=150 days)
suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})
rm(list = ls())

base_dir      <- "/Users/sidneybush/Library/CloudStorage/Box-Box"
project_dir   <- file.path(base_dir, "Sidney_Bush", "HJA_Water_Quality")
root_out_dir  <- file.path(project_dir, "outputs")
out_dir       <- file.path(root_out_dir, "windows_wet75_dry150")
dir.create(root_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load seasons and CQ master
season_df <- readr::read_csv(file.path(project_dir, "outputs", "HJA_daily_Q_with_seasons.csv"), show_col_types = FALSE) %>% mutate(date = as.Date(date), Stream_Name = as.character(Stream_Name)) %>% select(Stream_Name, date, water_year, hydrologic_season)
season_lookup <- season_df %>% distinct(Stream_Name, date, water_year, hydrologic_season)

cq_master_raw <- readr::read_csv(file.path(project_dir, "outputs", "HJA_CQ_master.csv"), show_col_types = FALSE) %>% rename(date = Date, Qcms = Q_cms) %>% mutate(date = as.Date(date), Stream_Name = as.character(Stream_Name))

cq_master <- cq_master_raw %>% left_join(season_df, by = c("Stream_Name", "date")) %>% filter(!is.na(hydrologic_season), !is.na(Qcms), !is.na(value), Qcms > 0, value > 0) %>% mutate(logQ = log10(Qcms), logC = log10(value))

fit_cq <- function(df, min_obs = 4, var_eps = 1e-12) {
  n <- nrow(df)
  if (n < min_obs) return(tibble(slope = NA_real_))
  if (length(unique(df$logQ)) < 2 || length(unique(df$logC)) < 2) return(tibble(slope = NA_real_))
  vq <- var(df$logQ, na.rm = TRUE); vc <- var(df$logC, na.rm = TRUE)
  if (!is.finite(vq) || !is.finite(vc) || vq < var_eps || vc < var_eps) return(tibble(slope = NA_real_))
  m <- lm(logC ~ logQ, data = df)
  tibble(slope = coef(m)[["logQ"]])
}

# Parameters: asymmetric half-window by season
half_wet  <- 75 %/% 2
half_dry  <- 150 %/% 2
MIN_OBS_WIN <- 5

# Conservative guards copied from main script
FIRST_PASS_B_MAX           <- 5.0
FIRST_PASS_CVR_MAX         <- 10.0
FIRST_PASS_LOGQ_SPREAD_MIN <- 0.10
B_ABS_MAX        <- 3.0
CVR_MIN          <- 0.0
CVR_MAX          <- 5.5
MIN_OBS_SUPPORT  <- 8
LOGQ_SPREAD_MIN  <- 0.10

roll_pts_raw <- cq_master %>% arrange(Stream_Name, variable, date) %>% group_by(Stream_Name, variable) %>% tidyr::nest() %>% mutate(res = purrr::map(data, function(df) {
  centers <- sort(unique(df$date))
  purrr::map_dfr(centers, function(d) {
    # determine season for the center date (use season_lookup)
    site_name <- if ("Stream_Name" %in% names(df)) df$Stream_Name[1] else NA_character_
    s <- season_lookup %>% filter(Stream_Name == site_name, date == d) %>% pull(hydrologic_season)
    half <- if (length(s) == 0 || is.na(s) || is.null(s)) 90 else if (s == "Wet") half_wet else half_dry
    wdf <- df %>% filter(date >= d - days(half), date <= d + days(half))
    n_obs <- sum(is.finite(wdf$Qcms) & is.finite(wdf$value))
    if (n_obs < MIN_OBS_WIN) return(tibble())
    rng_logQ <- if (sum(is.finite(wdf$logQ)) >= 2) diff(range(wdf$logQ[wdf$logQ %>% is.finite()])) else NA_real_
    fit <- fit_cq(wdf, min_obs = MIN_OBS_WIN)
    if (!is.finite(fit$slope)) return(tibble())
    mc <- mean(wdf$value, na.rm = TRUE); sc <- sd(wdf$value, na.rm = TRUE)
    mq <- mean(wdf$Qcms,  na.rm = TRUE); sq <- sd(wdf$Qcms,  na.rm = TRUE)
    cv_c <- if (is.finite(mc) && mc > 0) sc/mc else NA_real_
    cv_q <- if (is.finite(mq) && mq > 0) sq/mq else NA_real_
    cvr  <- if (is.finite(cv_c) && is.finite(cv_q) && cv_q > 0) cv_c/cv_q else NA_real_
    if (!is.finite(cvr) || cvr <= 0) return(tibble())
    if (abs(fit$slope) > FIRST_PASS_B_MAX)          return(tibble())
    if (is.finite(cvr) && cvr > FIRST_PASS_CVR_MAX) return(tibble())
    if (!is.finite(rng_logQ) || rng_logQ < FIRST_PASS_LOGQ_SPREAD_MIN) return(tibble())
    tibble(window_center = d, slope = fit$slope, CVc_CVq = cvr, n_obs = n_obs, range_logQ = rng_logQ)
  })
})) %>% select(-data) %>% unnest(res) %>% ungroup()

roll_pts_keep <- roll_pts_raw %>% mutate(weak_support  = (n_obs < MIN_OBS_SUPPORT) | (!is.finite(range_logQ)) | (range_logQ < LOGQ_SPREAD_MIN), extreme_slope = is.finite(slope)   & (abs(slope) > B_ABS_MAX), extreme_cvr   = is.finite(CVc_CVq) & (CVc_CVq < CVR_MIN | CVc_CVq > CVR_MAX), keep = !(weak_support & (extreme_slope | extreme_cvr))) %>% filter(keep) %>% filter(!is.na(variable))

# Log dropped solutes
all_solutes <- cq_master %>% distinct(variable) %>% pull(variable)
present_solutes <- roll_pts_keep %>% distinct(variable) %>% pull(variable)
missing_solutes <- setdiff(all_solutes, present_solutes)
write_csv(tibble(dropped_solute = missing_solutes), file.path(out_dir, "dropped_solutes_windows_wet75_dry150.csv"))

write_csv(roll_pts_keep, file.path(out_dir, "roll_pts_keep_windows_wet75_dry150.csv"))

classify_quad_slope_CVcCVq <- function(slope, CV) {
  case_when(is.na(slope) | is.na(CV) ~ NA_character_, CV >= 1 & slope > 0.1 ~ "1", CV <  1 & slope > 0.1 ~ "2", CV <  1 & slope < -0.1 ~ "3", CV >= 1 & slope < -0.1 ~ "4", TRUE ~ NA_character_)
}
classify_quad_slope_slope <- function(slope_x, slope_y) {
  case_when(is.na(slope_x) | is.na(slope_y) ~ NA_character_, slope_x > 0.1 & slope_y > 0.1 ~ "1", slope_x < -0.1 & slope_y > 0.1 ~ "2", slope_x < -0.1 & slope_y < -0.1 ~ "3", slope_x > 0.1 & slope_y < -0.1 ~ "4", TRUE ~ NA_character_)
}
get_sync <- function(quadrant) { case_when(is.na(quadrant) ~ NA_character_, quadrant %in% c("1","3") ~ "sync", quadrant %in% c("2","4") ~ "async", TRUE ~ NA_character_) }
cq_behavior_from_slope <- function(slope) { case_when(is.na(slope) ~ NA_character_, slope > 0.1 ~ "mobilizing", slope < -0.1 ~ "diluting", TRUE ~ "chemostatic") }

slope_CVcCVq_long <- roll_pts_keep %>% mutate(comparison_type = "cqslope_CVcCVq", solute1 = variable, solute2 = NA_character_, comparison_name = variable, slope_x = slope, slope_y = NA_real_, CV = CVc_CVq, quadrant = classify_quad_slope_CVcCVq(slope, CVc_CVq), sync = get_sync(quadrant), cq_behavior = cq_behavior_from_slope(slope_x)) %>% select(Stream_Name, window_center, comparison_type, solute1, solute2, comparison_name, quadrant, sync, cq_behavior, slope_x, slope_y, CV)

all_solutes_present <- roll_pts_keep %>% distinct(variable) %>% pull(variable)
if (length(all_solutes_present) >= 2) {
  pairs <- t(combn(all_solutes_present, 2)) %>% as_tibble() %>% rename(sol1 = V1, sol2 = V2)
} else {
  pairs <- tibble(sol1 = character(), sol2 = character())
}

if (nrow(pairs) > 0) {
  slope_slope_long <- pairs %>% rowwise() %>% do({
    x_sol <- .$sol1; y_sol <- .$sol2
    x_data <- roll_pts_keep %>% filter(variable == x_sol) %>% select(Stream_Name, window_center, slope_x = slope)
    y_data <- roll_pts_keep %>% filter(variable == y_sol) %>% select(Stream_Name, window_center, slope_y = slope)
    inner_join(x_data, y_data, by = c("Stream_Name","window_center")) %>% mutate(comparison_type = "cqslope_cqslope", solute1 = x_sol, solute2 = y_sol, comparison_name = paste(x_sol, y_sol, sep = "_"), CV = NA_real_, quadrant = classify_quad_slope_slope(slope_x, slope_y), sync = get_sync(quadrant), cq_behavior = NA_character_) %>% select(Stream_Name, window_center, comparison_type, solute1, solute2, comparison_name, quadrant, sync, cq_behavior, slope_x, slope_y, CV)
  }) %>% ungroup()
} else {
  slope_slope_long <- tibble()
}

slope_CVcCVq_long <- slope_CVcCVq_long %>% left_join(season_lookup, by = c("Stream_Name", "window_center" = "date"))
slope_slope_long <- slope_slope_long %>% left_join(season_lookup, by = c("Stream_Name", "window_center" = "date"))

cq_quadrants_long <- bind_rows(slope_CVcCVq_long, slope_slope_long) %>% arrange(Stream_Name, window_center, comparison_type, comparison_name)

variant_quadrants_path <- file.path(out_dir, "CQ_rolling_window_results_windows_wet75_dry150.csv")
general_quadrants_path <- file.path(root_out_dir, "CQ_rolling_window_results.csv")
readr::write_csv(cq_quadrants_long, variant_quadrants_path)
readr::write_csv(cq_quadrants_long, general_quadrants_path)

cq_solute_base <- slope_CVcCVq_long %>% transmute(Stream_Name, window_center, water_year, hydrologic_season, solute = solute1, cq_slope = slope_x, CVc_CVq = CV, cv_quadrant = quadrant, cv_sync = sync, cq_behavior = cq_behavior)

pair_long <- slope_slope_long %>% transmute(Stream_Name, window_center, solute = solute1, partner_solute = solute2, slope_pair_quadrant = quadrant, slope_pair_sync = sync)

if (nrow(pair_long) > 0) {
  pair_wide <- tidyr::pivot_wider(pair_long, names_from = partner_solute, values_from = c(slope_pair_quadrant, slope_pair_sync), names_glue = "{partner_solute}_{.value}")
} else {
  pair_wide <- tibble()
}

if (nrow(pair_wide) > 0) {
  cq_solute_window_summary <- cq_solute_base %>% left_join(pair_wide, by = c("Stream_Name","window_center","solute"))
} else {
  cq_solute_window_summary <- cq_solute_base
}

variant_solute_path <- file.path(out_dir, "CQ_solute_window_summary_windows_wet75_dry150.csv")
general_solute_path <- file.path(root_out_dir, "CQ_solute_window_summary.csv")
readr::write_csv(cq_solute_window_summary, variant_solute_path)
readr::write_csv(cq_solute_window_summary, general_solute_path)

counts <- roll_pts_keep %>% group_by(variable) %>% summarise(n_windows = n(), .groups = "drop")
write_csv(counts, file.path(out_dir, "n_windows_per_solute_windows_wet75_dry150.csv"))

