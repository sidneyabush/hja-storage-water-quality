# =======================================================================
# HJA C–Q Master and Monthly (with MDL-based ½ DL substitution for zeros)
# =======================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(stringr)
})

rm(list = ls())

# =============================================================================
# PATHS
# =============================================================================
repo_dir <- Sys.getenv(
  "HJA_WQ_REPO_DIR",
  unset = "/Users/sidneybush/Documents/GitHub/hja-water-quality"
)
source(file.path(repo_dir, "00_helpers", "workflow_config.R"))

paths <- get_project_paths()
raw_dir  <- paths$raw_dir
data_dir <- paths$data_dir

# exploratory plots (goes with other exploratory plots, not under data)
base_box_dir <- dirname(data_dir)  # /HJA_Water_Quality
plot_dir <- file.path(base_box_dir, "exploratory_plots", "01_data_prep", "conc_sample_bias")
write_prep_diagnostic_plots <- tolower(Sys.getenv("HJA_WRITE_PREP_DIAGNOSTIC_PLOTS", "false")) %in%
  c("1", "true", "yes", "y")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
if (write_prep_diagnostic_plots) {
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
}

# =============================================================================
# OUTPUT DIRECTORY (same level as "data")
# =============================================================================
output_dir <- file.path(dirname(data_dir), "outputs")
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
# VARIABLE CONFIGURATION (final selected solutes)
# =============================================================================
variable_rename_map <- c(
  "SI"   = "DSi",
  "PO4P" = "PO4",
  "NH3N" = "NH3",
  "NO3N" = "NO3",
  "NA"   = "Na",
  "K"    = "K",
  "CA"   = "Ca",
  "MG"   = "Mg",
  "SO4S" = "SO4",  # same values, CCAL data reported as SO4 (ion)
  "CL"   = "Cl",
  "DOC"  = "DOC"
)
keep_vars <- unname(variable_rename_map)

# =============================================================================
# CCAL detection limits (mg/L)
# Source: CCAL Analytical Detection Levels (accessed 2025-11-12)
# =============================================================================
USE_ML <- FALSE  # FALSE = use MDL, TRUE = use ML

ccal_dl <- tibble::tribble(
  ~variable, ~MDL_mgL, ~ML_mgL,
  "PO4", 0.001, 0.003,
  "NH3", 0.003, 0.009,
  "NO3", 0.001, 0.003,
  "Na",  0.010, 0.030,
  "K",   0.030, 0.100,
  "Ca",  0.060, 0.190,
  "Mg",  0.020, 0.060,
  "SO4", 0.010, 0.030,
  "Cl",  0.010, 0.030,
  "DOC", 0.050, 0.160,
  "DSi", 0.200, 0.600
)

dl_tbl <- ccal_dl %>%
  dplyr::mutate(DL_mgL = if (USE_ML) ML_mgL else MDL_mgL) %>%
  dplyr::select(variable, DL_mgL)

# =============================================================================
# Helper: join DLs + flag zero replacements (½ MDL)
# =============================================================================
apply_dl_half_for_zeros <- function(df, dl_lookup = dl_tbl) {
  df %>%
    dplyr::left_join(dl_lookup, by = "variable") %>%
    dplyr::mutate(
      ReplacedZero_MDLhalf = dplyr::if_else(!is.na(value) & value == 0 & !is.na(DL_mgL), TRUE, FALSE),
      value = dplyr::if_else(ReplacedZero_MDLhalf, 0.5 * DL_mgL, value)
    ) %>%
    dplyr::select(-DL_mgL)
}

# =======================================================================
# 1) CF00201 = MASTER (discrete) + Q_cms from MEAN_LPS
# =======================================================================
cf001_file <- first_existing(
  file.path(raw_dir, c("CF00201_v7.csv", "CF00201_v6.csv")),
  "CF00201 chemistry file"
)

cf001_raw <- readr::read_csv(cf001_file, show_col_types = FALSE) %>%
  dplyr::mutate(
    Date        = as.Date(DATE_TIME),
    Stream_Name = as.character(SITECODE),
    Q_cms       = suppressWarnings(as.numeric(MEAN_LPS)) * 0.001
  )

chem_cf001 <- cf001_raw %>%
  dplyr::select(-dplyr::any_of(c("STCODE","ENTITY","SITECODE","WATERYEAR","DATE_TIME","LABNO","TYPE",
                                 "INTERVAL","Q_AREA_CM","QCODE","PVOL","PVOLCODE","ANCA","ANCACODE",
                                 "MEAN_LPS"))) %>%
  dplyr::select(-dplyr::ends_with("CODE")) %>%
  tidyr::pivot_longer(
    cols      = -c(Stream_Name, Date, Q_cms),
    names_to  = "variable",
    values_to = "value"
  ) %>%
  dplyr::filter(!is.na(value), !is.na(Date)) %>%
  dplyr::filter(variable %in% names(variable_rename_map)) %>%
  dplyr::mutate(
    variable = dplyr::recode(variable, !!!variable_rename_map),
    value    = suppressWarnings(as.numeric(value)),
    units    = "mg/L",
    source   = "CF00201"
  ) %>%
  dplyr::select(Stream_Name, Date, variable, value, Q_cms, units, source) %>%
  apply_dl_half_for_zeros()

cq_master <- chem_cf001 %>%
  dplyr::filter(variable %in% keep_vars) %>%
  dplyr::mutate(
    Year  = lubridate::year(Date),
    Month = lubridate::month(Date)
  ) %>%
  dplyr::arrange(Stream_Name, variable, Date)

readr::write_csv(cq_master, file.path(output_dir, "HJA_CQ_master.csv"))

# =======================================================================
# 2) CF00202 = MONTHLY
# =======================================================================
cf002_file <- file.path(raw_dir, "CF00202_v5.csv")
cf002_monthly <- tibble()
if (file.exists(cf002_file)) {
  cf002_monthly <- readr::read_csv(cf002_file, show_col_types = FALSE) %>%
    dplyr::mutate(
      .YEAR   = as.integer(if ("YEAR"  %in% names(.)) YEAR  else NA_integer_),
      .MONTH  = as.integer(if ("MONTH" %in% names(.)) MONTH else NA_integer_),
      MONTH_STR = dplyr::if_else(!is.na(.YEAR) & !is.na(.MONTH),
                                 sprintf("%04d-%02d", .YEAR, .MONTH), NA_character_),
      Date        = as.Date(paste0(MONTH_STR, "-01")),
      Stream_Name = as.character(dplyr::coalesce(SITECODE, NA_character_))
    ) %>%
    dplyr::select(-dplyr::any_of(c("STCODE","ENTITY","SITECODE","WATERYEAR","YEAR","MONTH","TYPE",
                                   "MEAN_LPS","Q_AREA_MO","QCODE_MO"))) %>%
    dplyr::select(-dplyr::ends_with("CODE_MO")) %>%
    {
      helper_cols <- c("Stream_Name","Date",".YEAR",".MONTH","MONTH_STR")
      pivot_cols  <- setdiff(names(.), helper_cols[helper_cols %in% names(.)])
      tidyr::pivot_longer(., cols = tidyselect::all_of(pivot_cols),
                          names_to = "variable", values_to = "value")
    } %>%
    dplyr::mutate(
      variable = stringr::str_remove(variable, "_MO$"),
      variable = dplyr::recode(variable, !!!variable_rename_map, .default = NA_character_),
      value    = suppressWarnings(as.numeric(value)),
      units    = "mg/L",
      source   = "CF00202"
    ) %>%
    dplyr::filter(!is.na(Stream_Name), !is.na(Date), !is.na(value), variable %in% keep_vars) %>%
    dplyr::mutate(
      Year  = lubridate::year(Date),
      Month = lubridate::month(Date),
      MONTH = format(Date, "%Y-%m")
    ) %>%
    apply_dl_half_for_zeros() %>%
    dplyr::relocate(Stream_Name, Date, Year, Month, MONTH, variable, value, ReplacedZero_MDLhalf, units, source)
}

# Fill months missing from the official monthly product using the newer
# discrete chemistry file. CF00202 remains primary wherever it has a value.
cf002_monthly_official <- cf002_monthly %>%
  dplyr::mutate(n_discrete_samples = NA_integer_)

cf001_monthly_supplement <- cq_master %>%
  dplyr::filter(!is.na(value), variable %in% keep_vars) %>%
  dplyr::group_by(Stream_Name, variable, Year, Month) %>%
  dplyr::summarise(
    value = mean(value, na.rm = TRUE),
    ReplacedZero_MDLhalf = any(ReplacedZero_MDLhalf, na.rm = TRUE),
    n_discrete_samples = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Date = as.Date(sprintf("%04d-%02d-01", Year, Month)),
    MONTH = format(Date, "%Y-%m"),
    units = "mg/L",
    source = "CF00201_monthly_mean"
  ) %>%
  dplyr::select(Stream_Name, Date, Year, Month, MONTH, variable, value,
                ReplacedZero_MDLhalf, units, source, n_discrete_samples)

if (nrow(cf002_monthly_official) > 0) {
  cf001_monthly_supplement <- cf001_monthly_supplement %>%
    dplyr::anti_join(
      cf002_monthly_official %>%
        dplyr::distinct(Stream_Name, variable, Year, Month),
      by = c("Stream_Name", "variable", "Year", "Month")
    )
}

cf002_monthly <- dplyr::bind_rows(cf002_monthly_official, cf001_monthly_supplement) %>%
  dplyr::arrange(Stream_Name, variable, Date)

message("Monthly chemistry inputs:")
message("  - Official CF00202 rows: ", nrow(cf002_monthly_official))
message("  - CF00201-derived monthly rows added where CF00202 is absent: ", nrow(cf001_monthly_supplement))
message("  - Monthly table year range: ", min(cf002_monthly$Year, na.rm = TRUE), "-", max(cf002_monthly$Year, na.rm = TRUE))

# =======================================================================
# 4) Monthly Averages
# =======================================================================

# First, make sure cf002_monthly has Year & Month (you already did this above),
# so we just add water_year here and then compute three products:
#   (a) climatological monthly means (old behavior)
#   (b) water-year-specific monthly means (for WY-based analyses)
#   (c) calendar-year-specific monthly means (for DTW clustering)

cf002_monthly <- cf002_monthly %>%
  dplyr::mutate(
    water_year = dplyr::if_else(Month >= 10, Year + 1L, Year),
    calendar_year = Year  # calendar year is just Year
  )

# ---- (a) CLIMATOLOGICAL monthly means across all years (what you already had) ----
cf002_monthlyave_wide <- cf002_monthly %>%
  dplyr::filter(variable %in% keep_vars) %>%
  dplyr::group_by(Stream_Name, variable, Month) %>%
  dplyr::summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(Month_Name = month.abb[Month]) %>%
  dplyr::select(Stream_Name, variable, Month_Name, mean_value) %>%
  tidyr::pivot_wider(names_from = Month_Name, values_from = mean_value) %>%
  dplyr::select(Stream_Name, variable, Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec)

readr::write_csv(
  cf002_monthlyave_wide,
  file.path(output_dir, "HJA_CF002_monthly_means.csv")
)

# ---- (b) WATER-YEAR-SPECIFIC monthly means (for WY-based analyses) ----
cf002_monthlyave_wide_wy <- cf002_monthly %>%
  dplyr::filter(variable %in% keep_vars) %>%
  dplyr::group_by(Stream_Name, variable, water_year, Month) %>%
  dplyr::summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from  = Month,
    values_from = mean_value
  ) %>%
  # ensure columns are in a nice order: Stream, var, water_year, months 1–12
  dplyr::arrange(Stream_Name, variable, water_year) %>%
  dplyr::select(Stream_Name, variable, water_year, `1`:`12`)

readr::write_csv(
  cf002_monthlyave_wide_wy,
  file.path(output_dir, "HJA_CF002_monthly_means_byWY.csv")
)

# ---- (c) CALENDAR-YEAR-SPECIFIC monthly means (for DTW clustering) ----
# This is the CORRECT way to assign clusters to calendar years
cf002_monthlyave_wide_cy <- cf002_monthly %>%
  dplyr::filter(variable %in% keep_vars) %>%
  dplyr::group_by(Stream_Name, variable, calendar_year, Month) %>%
  dplyr::summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from  = Month,
    values_from = mean_value
  ) %>%
  # ensure columns are in a nice order: Stream, var, calendar_year, months 1–12
  dplyr::arrange(Stream_Name, variable, calendar_year) %>%
  dplyr::select(Stream_Name, variable, calendar_year, any_of(as.character(1:12)))

readr::write_csv(
  cf002_monthlyave_wide_cy,
  file.path(output_dir, "HJA_CF002_monthly_means_byCY.csv")
)

message("Created monthly means files:")
message("  - HJA_CF002_monthly_means.csv (climatological)")
message("  - HJA_CF002_monthly_means_byWY.csv (by water year)")
message("  - HJA_CF002_monthly_means_byCY.csv (by calendar year - for clustering)")

# =======================================================================
# 5) Chemistry sampling counts per Month × Year × Site (Bar Plot)
#     (Saved under data/exploratory_plots/conc_sample_bias)
# =======================================================================
if (write_prep_diagnostic_plots) {
  chem_counts <- cq_master %>%
    dplyr::mutate(
      Date     = as.Date(Date),
      Year     = lubridate::year(Date),
      MonthNum = lubridate::month(Date)
    ) %>%
    dplyr::distinct(Stream_Name, Year, MonthNum, Date) %>%
    dplyr::group_by(Stream_Name, Year, MonthNum) %>%
    dplyr::summarise(n_samples = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(Stream_Name, Year) %>%
    tidyr::complete(MonthNum = 1:12, fill = list(n_samples = 0)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(MonthName = factor(month.abb[MonthNum], levels = month.abb))

  sites <- sort(unique(chem_counts$Stream_Name))

  for (s in sites) {
    df_s <- chem_counts %>% dplyr::filter(Stream_Name == s)

    p_site <- ggplot(df_s, aes(x = MonthName, y = n_samples)) +
      geom_col(fill = "grey35") +
      facet_wrap(~ Year, ncol = 6) +
      labs(
        x = "Month",
        y = "Sampling dates"
      ) +
      theme_bw(base_size = 11) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "grey92", color = NA),
        panel.grid.minor = element_blank()
      )

    out_file <- file.path(
      plot_dir,
      paste0("Chem_Samples_MonthByYear_", gsub("[^A-Za-z0-9_]+", "_", s), ".png")
    )
    ggsave(out_file, p_site, width = 14, height = 8, dpi = 300)
  }

  # =======================================================================
  # 6) Supplemental FIGURE: DISCHARGE QUANTILE TIMING
  # =======================================================================

  dq_output_path <- file.path(base_box_dir, "exploratory_plots", "01_data_prep", "discharge_analysis")
  if (!dir.exists(dq_output_path)) dir.create(dq_output_path, recursive = TRUE)

  hf_file <- first_existing(
    file.path(raw_dir, c("HF00402_v14.csv", "HF00402_v15.csv")),
    "HF00402 discharge file"
  )

  Q_raw <- readr::read_csv(hf_file, show_col_types = FALSE) %>%
    dplyr::transmute(
      Stream_Name = as.character(SITECODE),
      date        = parse_hja_date(DATE),
      Qcms        = as.numeric(MEAN_Q) * 0.028316846592   # cfs -> m3/s
    ) %>%
    dplyr::filter(!is.na(Stream_Name), !is.na(date), is.finite(Qcms)) %>%
    dplyr::mutate(
      Stream_Name = dplyr::case_when(
        Stream_Name == "GSLOOK_FULL" ~ "GSLOOK",
        TRUE ~ Stream_Name
      ),
      month       = lubridate::month(date),
      water_year  = dplyr::if_else(month >= 10, lubridate::year(date) + 1L,
                                   lubridate::year(date)) %>% as.integer()
    )

  quantile_timing <- Q_raw %>%
    dplyr::group_by(Stream_Name, water_year) %>%
    dplyr::filter(dplyr::n() >= 300) %>%
    dplyr::mutate(
      Q_25  = stats::quantile(Qcms, 0.25, na.rm = TRUE),
      Q_50  = stats::quantile(Qcms, 0.50, na.rm = TRUE),
      Q_75  = stats::quantile(Qcms, 0.75, na.rm = TRUE),
      Q_max = max(Qcms, na.rm = TRUE),
      Q_category = dplyr::case_when(
        Qcms >= Q_75 ~ "75-100",
        Qcms >= Q_50 ~ "50-75",
        Qcms >= Q_25 ~ "25-50",
        TRUE ~ "0-25"
      )
    ) %>%
    dplyr::ungroup()

  quantile_by_month <- quantile_timing %>%
    dplyr::group_by(month, Q_category) %>%
    dplyr::summarise(n_days = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(month) %>%
    dplyr::mutate(
      pct_of_month = 100 * n_days / sum(n_days),
      month_name   = month.abb[month]
    ) %>%
    dplyr::ungroup()

  month_labels <- c("J","F","M","A","M","J","J","A","S","O","N","D")

  p_discharge_quantiles <- ggplot(quantile_by_month,
                                  aes(x = month_name, y = Q_category, fill = pct_of_month)) +
    geom_tile(color = "white", linewidth = 1) +
    scale_fill_gradient2(
      low = "#F0F4F7", mid = "#9CAEC4", high = "#4B6C91",
      midpoint = 50, limits = c(0, 100), breaks = seq(0, 100, 25),
      name = "% of Days"
    ) +
    scale_x_discrete(limits = month.abb, labels = month_labels) +
    labs(x = "Month", y = "Discharge Quantile") +
    theme_bw(base_size = 20) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.title  = element_text(size = 22),
      legend.position = "right",
      legend.title = element_text(size = 20),
      legend.text  = element_text(size = 18),
      legend.key.width  = grid::unit(20, "pt"),
      legend.key.height = grid::unit(30, "pt"),
      plot.margin = margin(10, 10, 10, 10)
    )

  ggsave(file.path(dq_output_path, "discharge_quantiles_by_month.png"),
         p_discharge_quantiles, width = 12, height = 8, dpi = 300, bg = "white")
}
