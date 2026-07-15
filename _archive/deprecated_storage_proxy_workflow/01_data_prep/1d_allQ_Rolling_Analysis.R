# =============================================================================
# HJA Rolling Hydrologic & Storage Metrics (ASYMMETRIC windows)
#
# Uses window_center dates from CQ analysis for alignment
# Window sizes: Wet=75d, Dry=150d (matches 1c_CQ_Rolling_Analysis.R)
#
# Outputs: HJA_rolling_hydro_storage_90d.csv  
#  - RCS: recession p, k, RCS_n
#  - Q_dS_*: Kirchner/Staudinger discharge–storage ranges (mm)
#  - WB_dS_*: Water-balance storage (mm): net / min / max / range
#  - FDC_slope_5_95, Q_FDC01/50/99, Q05
#  - RBI: Richards–Baker Index
#  - water_year, hydrologic_season
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(zoo)
  library(progress)
})

rm(list = ls())

# =============================================================================
# Paths
# =============================================================================
base_dir    <- "/Users/sidneybush/Library/CloudStorage/Box-Box"
project_dir <- file.path(base_dir, "Sidney_Bush", "HJA_Water_Quality")
data_dir    <- file.path(project_dir, "data")
out_dir     <- file.path(project_dir, "outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Parameters - ASYMMETRIC windows (matches 1c_CQ_Rolling_Analysis.R)
# =============================================================================
WINDOW_DAYS_WET  <- 75
WINDOW_DAYS_DRY  <- 150
half_win_wet     <- WINDOW_DAYS_WET %/% 2
half_win_dry     <- WINDOW_DAYS_DRY %/% 2
MIN_REC_POINTS   <- 5
MIN_Q_POINTS_WIN <- 10
FDC_PROB_RANGE   <- c(5, 95)

# =============================================================================
# Load Q + seasons
# =============================================================================
q_with_seasons <- read_csv(file.path(out_dir, "HJA_daily_Q_with_seasons.csv"),
                           show_col_types = FALSE) %>%
  mutate(
    date        = as.Date(date),
    Stream_Name = as.character(Stream_Name)
  ) %>%
  filter(!is.na(hydrologic_season)) %>%
  select(Stream_Name, date, Qcms, water_year, hydrologic_season)

season_lookup <- q_with_seasons %>%
  distinct(Stream_Name, date, water_year, hydrologic_season)

# =============================================================================
# Drainage areas → convert Q to mm/day
# =============================================================================
da_df <- read_csv(file.path(data_dir, "drainage_area.csv"),
                  show_col_types = FALSE) %>%
  transmute(
    Stream_Name = as.character(SITECODE),
    DA_m2       = DA_M2
  )

q_daily <- q_with_seasons %>%
  left_join(da_df, by = "Stream_Name") %>%
  filter(!is.na(DA_m2), is.finite(Qcms)) %>%
  mutate(
    # Qcms = m^3/s; convert to mm/day over basin
    Q_mm_d = Qcms / DA_m2 * 86400 * 1000
  )

# =============================================================================
# Daily P/Q/ET water balance file
# =============================================================================
wb_daily <- read_csv(file.path(data_dir, "daily_water_balance_ET_Hamon-Zhang_coeff_interp.csv"),
                     show_col_types = FALSE) %>%
  mutate(
    date        = parse_date_time(DATE, orders = c("Ymd","Y-m-d","mdy","dmy")) %>% as.Date(),
    month       = month(date),
    water_year  = if_else(month >= 10, year(date) + 1L, year(date)) %>% as.integer(),
    Stream_Name = case_when(
      SITECODE == "GSLOOK_FULL" ~ "GSLOOK",
      TRUE                      ~ as.character(SITECODE)
    )
  ) %>%
  rename(
    P_mm_d    = P_mm_d,
    Q_mm_d_wb = Q_mm_d,   # Q in mm/d from WB bookkeeping
    ET_mm_d   = ET_mm_d
  ) %>%
  filter(!Stream_Name %in% c("COLD","LONGER"))

# =============================================================================
# Merge Q-based and WB-based daily
# =============================================================================
daily_all <- q_daily %>%
  left_join(
    wb_daily %>% select(Stream_Name, date, P_mm_d, Q_mm_d_wb, ET_mm_d),
    by = c("Stream_Name","date")
  ) %>%
  mutate(
    # WB_dS_daily_mm = P - Q - ET from the WB file (mm/day)
    WB_dS_daily_mm = if_else(
      !is.na(P_mm_d) & !is.na(Q_mm_d_wb) & !is.na(ET_mm_d),
      P_mm_d - Q_mm_d_wb - ET_mm_d,
      NA_real_
    )
  ) %>%
  # keep only (Stream_Name, date) pairs that actually exist in the seasons file
  semi_join(q_with_seasons %>% select(Stream_Name, date),
            by = c("Stream_Name","date"))

# Drop watershed 3 becasue we don't have chemistry data for that site:
daily_all <- daily_all %>%
  filter(Stream_Name != "GSWS03")

# =============================================================================
# Helper Methods: FDC, RCS/Q_dS, RBI, Q05, WB_dS
# =============================================================================
compute_fdc <- function(Q) {
  Qpos <- Q[Q > 0 & is.finite(Q)]
  if (length(Qpos) < MIN_Q_POINTS_WIN) {
    return(tibble(exceedance = numeric(0), Q = numeric(0)))
  }
  Qs <- sort(Qpos, decreasing = TRUE)
  n  <- length(Qs)
  tibble(
    exceedance = (seq_len(n) - 0.44) / (n + 0.12) * 100,
    Q          = Qs
  )
}

get_Q_at_prob <- function(fdc, prob) {
  if (nrow(fdc) < 2) return(NA_real_)
  approx(fdc$exceedance, fdc$Q, xout = prob)$y
}

calc_FDC_slope <- function(fdc, pr = c(5, 95)) {
  if (nrow(fdc) < 2) return(NA_real_)
  d <- fdc %>%
    filter(exceedance >= pr[1],
           exceedance <= pr[2],
           Q > 0)
  if (nrow(d) < 5) return(NA_real_)
  coef(lm(log10(Q) ~ exceedance, data = d))[["exceedance"]]
}

calc_RCS_QdS <- function(df) {
  df_r <- df %>%
    arrange(date) %>%
    mutate(
      dt   = as.numeric(difftime(date, lag(date), units = "days")),
      dQ   = (lag(Q_mm_d) - Q_mm_d) / dt,
      rain = if_else(!is.na(P_mm_d), P_mm_d > 0, FALSE)
    ) %>%
    filter(!is.na(dQ), dQ > 0, Q_mm_d > 0, !rain)
  
  if (nrow(df_r) < MIN_REC_POINTS) {
    return(tibble(
      RCS_p            = NA_real_,
      RCS_k            = NA_real_,
      Q_dS_range_mm    = NA_real_,
      Q_dS_high_med_mm = NA_real_,
      Q_dS_med_low_mm  = NA_real_,
      Q_FDC01_mm_d     = NA_real_,
      Q_FDC50_mm_d     = NA_real_,
      Q_FDC99_mm_d     = NA_real_,
      Q_max_mm_d       = max(df$Q_mm_d, na.rm = TRUE),
      Q_min_mm_d       = suppressWarnings(min(df$Q_mm_d[df$Q_mm_d > 0], na.rm = TRUE)),
      RCS_n            = nrow(df_r)
    ))
  }
  
  fit <- lm(log(dQ) ~ log(Q_mm_d), data = df_r)
  p   <- coef(fit)[["log(Q_mm_d)"]]
  k   <- exp(coef(fit)[["(Intercept)"]])
  
  fdc  <- compute_fdc(df$Q_mm_d)
  Q01  <- get_Q_at_prob(fdc, 1)
  Q50  <- get_Q_at_prob(fdc, 50)
  Q99  <- get_Q_at_prob(fdc, 99)
  Qmax <- max(df$Q_mm_d, na.rm = TRUE)
  Qmin <- suppressWarnings(min(df$Q_mm_d[df$Q_mm_d > 0], na.rm = TRUE))
  
  dS <- function(Qu, Ql, k, p) {
    if (any(is.na(c(Qu, Ql, k, p)))) return(NA_real_)
    if (p >= 2 || k <= 0) return(NA_real_)
    (Qu^(2 - p) - Ql^(2 - p)) / (k * (2 - p))
  }
  
  tibble(
    RCS_p            = p,
    RCS_k            = k,
    Q_dS_range_mm    = dS(Qmax, Qmin, k, p),
    Q_dS_high_med_mm = dS(Q01,  Q50, k, p),
    Q_dS_med_low_mm  = dS(Q50,  Q99, k, p),
    Q_FDC01_mm_d     = Q01,
    Q_FDC50_mm_d     = Q50,
    Q_FDC99_mm_d     = Q99,
    Q_max_mm_d       = Qmax,
    Q_min_mm_d       = Qmin,
    RCS_n            = nrow(df_r)
  )
}

calc_RBI <- function(Q) {
  if (length(Q) < 3 || all(is.na(Q))) return(NA_real_)
  d <- tibble(Q = Q) %>%
    mutate(dQ = Q - lag(Q)) %>%
    filter(!is.na(dQ))
  total_Q <- sum(d$Q, na.rm = TRUE)
  if (!is.finite(total_Q) || total_Q <= 0) return(NA_real_)
  sum(abs(d$dQ), na.rm = TRUE) / total_Q
}

calc_Q05 <- function(Q) {
  if (length(Q) < 3 || all(is.na(Q))) return(NA_real_)
  quantile(Q, 0.05, na.rm = TRUE)
}

calc_WB_dS_stats <- function(WB) {
  if (all(is.na(WB))) {
    return(tibble(
      WB_dS_net_mm   = NA_real_,
      WB_dS_min_mm   = NA_real_,
      WB_dS_max_mm   = NA_real_,
      WB_dS_range_mm = NA_real_
    ))
  }
  cum <- cumsum(replace_na(WB, 0))
  cum <- cum - dplyr::first(cum)
  wb_min <- min(cum, na.rm = TRUE)
  wb_max <- max(cum, na.rm = TRUE)
  tibble(
    WB_dS_net_mm   = sum(WB, na.rm = TRUE),
    WB_dS_min_mm   = wb_min,
    WB_dS_max_mm   = wb_max,
    WB_dS_range_mm = wb_max - wb_min
  )
}

streams <- sort(unique(daily_all$Stream_Name))

# =============================================================================
# Load CQ window centers as CANONICAL dates (from 1c_CQ_Rolling_Analysis.R)
# This ensures storage metrics align with chemistry sample windows
# =============================================================================
cq_windows <- read_csv(file.path(out_dir, "CQ_solute_window_summary.csv"),
                       show_col_types = FALSE) %>%
  select(Stream_Name, window_center, hydrologic_season) %>%
  distinct() %>%
  mutate(window_center = as.Date(window_center))

message("Using ", nrow(cq_windows), " canonical CQ window centers for storage calculations")

pb <- progress_bar$new(
  total      = length(streams),
  format     = "  computing rolling windows [:bar] :current/:total streams (:percent) eta=:eta :stream",
  clear      = FALSE,
  width      = 80,
  show_after = 0,    # show immediately
  force      = TRUE  # show even if output not recognized as a tty
)

# Draw the empty bar without printing the R6 object
invisible(pb$tick(0))

rolling_list <- vector("list", length(streams))

for (i in seq_along(streams)) {
  s <- streams[i]
  
  df_site <- daily_all %>%
    filter(Stream_Name == s) %>%
    arrange(date)
  
  # Use ONLY CQ window centers for this stream (canonical dates)
  site_cq_windows <- cq_windows %>%
    filter(Stream_Name == s)
  
  if (nrow(site_cq_windows) == 0) {
    rolling_list[[i]] <- tibble()
    pb$tick(tokens = list(stream = s))
    next
  }
  
  res_site <- purrr::map_dfr(seq_len(nrow(site_cq_windows)), function(j) {
    d <- site_cq_windows$window_center[j]
    season <- site_cq_windows$hydrologic_season[j]
    
    # Use asymmetric half-window based on season
    half_win <- if (!is.na(season) && season == "Wet") half_win_wet else half_win_dry
    
    w <- df_site %>%
      filter(date >= d - days(half_win),
             date <= d + days(half_win))
    
    if (nrow(w) < 30) return(tibble())
    
    tibble(
      Stream_Name      = s,
      window_center    = d,
      window_size_days = half_win * 2,
      RBI              = calc_RBI(w$Q_mm_d),
      FDC_slope_5_95   = calc_FDC_slope(compute_fdc(w$Q_mm_d), FDC_PROB_RANGE),
      Q05_mm_d         = calc_Q05(w$Q_mm_d)
    ) %>%
      bind_cols(
        calc_RCS_QdS(w),
        calc_WB_dS_stats(w$WB_dS_daily_mm)
      )
  })
  
  rolling_list[[i]] <- res_site
  
  # Update progress bar, showing current stream name
  pb$tick(tokens = list(stream = s))
}

rolling_raw <- bind_rows(rolling_list)

# =============================================================================
# Save
# =============================================================================
write_csv(rolling_raw, file.path(out_dir, "HJA_rolling_hydro_storage_90d.csv"))
