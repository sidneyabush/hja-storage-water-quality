# =============================================================================
# VISUALIZE HYDROLOGIC SEASONS
# =============================================================================

library(tidyverse)
library(lubridate)
library(patchwork)

rm(list = ls())
try(source(file.path("/Users/sidneybush/Documents/GitHub/hja-water-quality", "00_helpers", "plot_theme_set.R")), silent = TRUE)

base_dir    <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality"
data_path   <- file.path(base_dir, "data")
output_dir  <- file.path(base_dir, "outputs")
output_path_seasons <- file.path(base_dir, "exploratory_plots", "02_exploration", "2a_hydro_seasons")

if (!dir.exists(output_path_seasons)) dir.create(output_path_seasons, recursive = TRUE)

# -----------------------------------------------------------------------------
# Load season data
# -----------------------------------------------------------------------------
Q_with_seasons <- read_csv(
  file.path(output_dir, "HJA_daily_Q_with_seasons.csv"),
  show_col_types = FALSE
)

season_bounds <- read_csv(
  file.path(output_dir, "season_boundaries.csv"),
  show_col_types = FALSE
)

# Derive year classification from season_boundaries peak_Q_value
# Classify years as Wet/Normal/Dry based on tertiles of peak discharge
year_class <- season_bounds %>%
  filter(!is.na(peak_Q_value)) %>%
  mutate(
    peak_Q_tertile = ntile(peak_Q_value, 3),
    year_class = case_when(
      peak_Q_tertile == 1 ~ "Dry",
      peak_Q_tertile == 2 ~ "Normal",
      peak_Q_tertile == 3 ~ "Wet",
      TRUE ~ "Unknown"
    )
  ) %>%
  select(water_year, year_class, peak_Q_value)

# -----------------------------------------------------------------------------
# Load C–Q master (for sample dates)
# -----------------------------------------------------------------------------
cq_master <- read_csv(
  file.path(output_dir, "HJA_CQ_master.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    date       = as.Date(Date),
    water_year = if_else(month(date) >= 10, year(date) + 1L, year(date))
  )

Q_plot_sites <- Q_with_seasons %>%
  left_join(
    year_class %>% select(water_year, year_class, peak_Q_value),
    by = "water_year"
  )

# Calculate median for median plots
Q_plot_median <- Q_with_seasons %>%
  group_by(date, water_year, hydrologic_season) %>%
  summarise(
    Q_median = median(Qcms, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    year_class %>% select(water_year, year_class, peak_Q_value),
    by = "water_year"
  )

# -----------------------------------------------------------------------------
# Helper: compute contiguous Wet / Dry spans from daily labels
# -----------------------------------------------------------------------------
get_season_spans <- function(df) {
  df %>%
    arrange(date) %>%
    mutate(
      season = hydrologic_season,
      season = if_else(is.na(season), "None", season),
      grp = cumsum(season != dplyr::lag(season, default = first(season)))
    ) %>%
    group_by(season, grp) %>%
    summarise(
      xmin = min(date),
      xmax = max(date),
      .groups = "drop"
    ) %>%
    filter(season %in% c("Wet", "Dry")) %>%
    mutate(season = factor(season, levels = c("Dry", "Wet")))
}

# -----------------------------------------------------------------------------
# Function: Plot water year - MEDIAN discharge (with sample dates)
# -----------------------------------------------------------------------------
plot_water_year_median <- function(year, Q_data, bounds_data) {
  
  # Get data from Jul of previous year to Oct of next year
  start_date <- as.Date(paste0(year - 1, "-07-01"))
  end_date   <- as.Date(paste0(year,     "-10-31"))
  
  year_Q <- Q_data %>%
    filter(date >= start_date, date <= end_date) %>%
    mutate(is_focal_year = water_year == year)
  
  year_bounds <- bounds_data %>% filter(water_year == year)
  
  if (nrow(year_Q) == 0 || nrow(year_bounds) == 0) {
    return(NULL)
  }
  
  # Year classification
  year_type <- unique(year_Q$year_class)
  if (length(year_type) == 0) year_type <- "Unknown"
  
  peak_Q <- unique(year_Q$peak_Q_value)
  if (length(peak_Q) == 0) peak_Q <- NA
  
  # Season boundaries for vertical lines (from algorithm)
  boundaries <- tibble(
    date = c(
      year_bounds$wet_start_date,
      year_bounds$wet_end_date
    ),
    label = c("Wet Start / Dry End", "Wet End / Dry Start"),
    line_color = c("darkblue", "darkblue")
  ) %>%
    filter(!is.na(date))
  
  # ---------------------------------------------------------------------------
  # Build Wet / Dry shading from daily hydrologic_season labels
  # ---------------------------------------------------------------------------
  spans <- get_season_spans(year_Q)
  
  dry_periods <- spans %>%
    filter(season == "Dry") %>%
    select(xmin, xmax)
  
  wet_period <- spans %>%
    filter(season == "Wet") %>%
    select(xmin, xmax)
  
  # ---------------------------------------------------------------------------
  # Sample dates (all sites, distinct dates)
  # ---------------------------------------------------------------------------
  sample_dates <- cq_master %>%
    filter(
      water_year == year,
      date >= start_date, date <= end_date
    ) %>%
    distinct(date)
  
  # Plot median discharge
  p <- ggplot(year_Q, aes(x = date, y = Q_median)) +
    # Shade dry seasons in light brown
    geom_rect(
      data = dry_periods,
      aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf),
      fill = "#D2B48C", alpha = 0.2, inherit.aes = FALSE
    ) +
    # Shade wet seasons in light blue
    geom_rect(
      data = wet_period,
      aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf),
      fill = "lightblue", alpha = 0.3, inherit.aes = FALSE
    ) +
    # Discharge line
    geom_line(linewidth = 0.7, color = "black") +
    # Sample dates rug (bottom)
    geom_rug(
      data = sample_dates,
      aes(x = date),
      inherit.aes = FALSE,
      sides = "b",
      alpha = 0.5
    ) +
    # Season boundaries (solid dark blue lines)
    geom_vline(
      data = boundaries,
      aes(xintercept = date),
      color = "darkblue",
      linetype = "solid", linewidth = 0.8
    ) +
    scale_y_log10() +
    scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
    labs(
      x = "Date",
      y = "Discharge (mm/d, log scale)"
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12),
      axis.title  = element_text(size = 13),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

# -----------------------------------------------------------------------------
# Function: Plot water year - ALL SITES (with sample dates per site)
# -----------------------------------------------------------------------------
plot_water_year_sites <- function(year, Q_data, bounds_data) {
  
  # Get data from Jul of previous year to Oct of next year
  start_date <- as.Date(paste0(year - 1, "-07-01"))
  end_date   <- as.Date(paste0(year,     "-10-31"))
  
  # Water year boundaries for data completeness
  wy_start <- as.Date(paste0(year - 1, "-10-01"))
  wy_end   <- as.Date(paste0(year,     "-09-30"))
  
  # Filter to sites with full water year of data (>= 300 days)
  # Exclude GSWSMA and GSWSMF (only include GSWSMC from MA/MC/MF sites)
  sites_with_full_data <- Q_data %>%
    filter(
      water_year == year,
      date >= wy_start, date <= wy_end,
      !Stream_Name %in% c("GSWSMA", "GSWSMF")
    ) %>%
    group_by(Stream_Name) %>%
    summarise(n_days = n(), .groups = "drop") %>%
    filter(n_days >= 300) %>%
    pull(Stream_Name)
  
  year_Q <- Q_data %>%
    filter(
      date >= start_date, date <= end_date,
      Stream_Name %in% sites_with_full_data
    ) %>%
    mutate(is_focal_year = water_year == year)
  
  year_bounds <- bounds_data %>% filter(water_year == year)
  
  if (nrow(year_Q) == 0 || nrow(year_bounds) == 0 || length(sites_with_full_data) == 0) {
    return(NULL)
  }
  
  # Year classification
  year_type <- unique(year_Q$year_class)
  if (length(year_type) == 0) year_type <- "Unknown"
  
  peak_Q <- unique(year_Q$peak_Q_value)
  if (length(peak_Q) == 0) peak_Q <- NA
  
  # Season boundaries for vertical lines
  boundaries <- tibble(
    date = c(
      year_bounds$wet_start_date,
      year_bounds$peak_Q_date,
      year_bounds$wet_end_date
    ),
    label = c("Wet Start", "Peak Q", "Wet End"),
    color = c("darkgreen", "blue", "orange")
  ) %>%
    filter(!is.na(date))
  
  # ---------------------------------------------------------------------------
  # Build Wet / Dry shading from daily hydrologic_season labels
  # Use basin-median series for spans so rectangles match across facets
  # ---------------------------------------------------------------------------
  year_Q_median <- year_Q %>%
    group_by(date, water_year, hydrologic_season) %>%
    summarise(Qcms = median(Qcms, na.rm = TRUE), .groups = "drop")
  
  spans <- get_season_spans(year_Q_median)
  
  dry_periods <- spans %>%
    filter(season == "Dry") %>%
    select(xmin, xmax)
  
  wet_period <- spans %>%
    filter(season == "Wet") %>%
    select(xmin, xmax)
  
  # ---------------------------------------------------------------------------
  # Sample dates per site (for rug)
  # ---------------------------------------------------------------------------
  sample_dates_sites <- cq_master %>%
    filter(
      water_year == year,
      date >= start_date, date <= end_date,
      Stream_Name %in% sites_with_full_data
    ) %>%
    distinct(Stream_Name, date)
  
  # Plot with facets for each site
  p <- ggplot(year_Q, aes(x = date, y = Qcms)) +
    # Shade dry season periods in light brown
    geom_rect(
      data = dry_periods,
      aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf),
      fill = "#D2B48C", alpha = 0.2, inherit.aes = FALSE
    ) +
    # Shade wet season periods in light blue
    geom_rect(
      data = wet_period,
      aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf),
      fill = "lightblue", alpha = 0.3, inherit.aes = FALSE
    ) +
    # Discharge line for each site (dark gray to avoid clash with shading)
    geom_line(linewidth = 0.5, color = "darkgray") +
    # Sample dates rug per site (bottom)
    geom_rug(
      data = sample_dates_sites,
      aes(x = date),
      inherit.aes = FALSE,
      sides = "b",
      alpha = 0.5
    ) +
    # Season boundaries
    geom_vline(
      data = boundaries,
      aes(xintercept = date, linetype = label),
      color = "black", linewidth = 0.6
    ) +
    scale_linetype_manual(
      values = c("Wet Start" = "dashed", "Peak Q" = "solid", "Wet End" = "dashed"),
      name = "Season Boundaries"
    ) +
    scale_y_log10() +
    scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
    facet_wrap(~ Stream_Name, ncol = 3, scales = "free_y") +
    labs(
      x = "Date",
      y = "Discharge (mm/d, log scale)"
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10),
      axis.title  = element_text(size = 12),
      strip.text = element_text(size = 12, face = "bold"),
      strip.background = element_rect(fill = "white"),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

# -----------------------------------------------------------------------------
# Create output directories
# -----------------------------------------------------------------------------
output_path_median    <- file.path(output_path_seasons, "median_behavior")
output_path_all_sites <- file.path(output_path_seasons, "individual_site_behavior")

if (!dir.exists(output_path_median))    dir.create(output_path_median, recursive = TRUE)
if (!dir.exists(output_path_all_sites)) dir.create(output_path_all_sites, recursive = TRUE)

# -----------------------------------------------------------------------------
# Create plots for all water years
# -----------------------------------------------------------------------------
years_to_plot <- season_bounds %>%
  filter(!is.na(wet_start_date)) %>%
  pull(water_year) %>%
  sort()

# MEDIAN plots
year_plots_median <- map(years_to_plot, ~ plot_water_year_median(.x, Q_plot_median, season_bounds))
year_plots_median <- compact(year_plots_median)

for (i in seq_along(year_plots_median)) {
  year <- years_to_plot[i]
  ggsave(
    file.path(output_path_median, paste0("season_validation_WY", year, "_MEDIAN.png")),
    year_plots_median[[i]],
    width = 12, height = 5, dpi = 300
  )
}

# ALL SITES plots
year_plots_sites <- map(years_to_plot, ~ plot_water_year_sites(.x, Q_plot_sites, season_bounds))
year_plots_sites <- compact(year_plots_sites)

for (i in seq_along(year_plots_sites)) {
  year <- years_to_plot[i]
  
  wy_start <- as.Date(paste0(year - 1, "-10-01"))
  wy_end   <- as.Date(paste0(year,     "-09-30"))
  
  n_sites_year <- Q_plot_sites %>%
    filter(
      water_year == year,
      date >= wy_start, date <= wy_end,
      !Stream_Name %in% c("GSWSMA", "GSWSMF")
    ) %>%
    group_by(Stream_Name) %>%
    summarise(n_days = n(), .groups = "drop") %>%
    filter(n_days >= 300) %>%
    nrow()
  
  n_rows      <- ceiling(n_sites_year / 3)
  plot_height <- max(6, 3 * n_rows)
  plot_width  <- 16
  
  ggsave(
    file.path(output_path_all_sites, paste0("season_validation_WY", year, "_ALL_SITES.png")),
    year_plots_sites[[i]],
    width = plot_width, height = plot_height, dpi = 300
  )
}

# -----------------------------------------------------------------------------
# Multi-panel overview plots by year class
# -----------------------------------------------------------------------------
years_by_class <- Q_plot_median %>%
  distinct(water_year, year_class) %>%
  arrange(water_year)

# Very Dry
very_dry_years <- years_by_class %>%
  filter(grepl("Very Dry", year_class)) %>%
  pull(water_year)

if (length(very_dry_years) > 0) {
  
  for (year in very_dry_years) {
    p <- plot_water_year_median(year, Q_plot_median, season_bounds)
    if (!is.null(p)) {
      ggsave(
        file.path(output_path_median, paste0("season_validation_WY", year, "_VERY_DRY_MEDIAN.png")),
        p,
        width = 12, height = 5, dpi = 300
      )
    }
  }
  
  for (year in very_dry_years) {
    p <- plot_water_year_sites(year, Q_plot_sites, season_bounds)
    if (!is.null(p)) {
      wy_start <- as.Date(paste0(year - 1, "-10-01"))
      wy_end   <- as.Date(paste0(year,     "-09-30"))
      n_sites_year <- Q_plot_sites %>%
        filter(
          water_year == year,
          date >= wy_start, date <= wy_end,
          !Stream_Name %in% c("GSWSMA", "GSWSMF")
        ) %>%
        group_by(Stream_Name) %>%
        summarise(n_days = n(), .groups = "drop") %>%
        filter(n_days >= 300) %>%
        nrow()
      n_rows      <- ceiling(n_sites_year / 3)
      plot_height <- max(6, 3 * n_rows)
      
      ggsave(
        file.path(output_path_all_sites, paste0("season_validation_WY", year, "_VERY_DRY_ALL_SITES.png")),
        p,
        width = 16, height = plot_height, dpi = 300
      )
    }
  }
}

# Very Wet
very_wet_years <- years_by_class %>%
  filter(grepl("Very Wet", year_class)) %>%
  pull(water_year)

if (length(very_wet_years) > 0) {
  
  for (year in very_wet_years) {
    p <- plot_water_year_median(year, Q_plot_median, season_bounds)
    if (!is.null(p)) {
      ggsave(
        file.path(output_path_median, paste0("season_validation_WY", year, "_VERY_WET_MEDIAN.png")),
        p,
        width = 12, height = 5, dpi = 300
      )
    }
  }
  
  for (year in very_wet_years) {
    p <- plot_water_year_sites(year, Q_plot_sites, season_bounds)
    if (!is.null(p)) {
      wy_start <- as.Date(paste0(year - 1, "-10-01"))
      wy_end   <- as.Date(paste0(year,     "-09-30"))
      n_sites_year <- Q_plot_sites %>%
        filter(
          water_year == year,
          date >= wy_start, date <= wy_end,
          !Stream_Name %in% c("GSWSMA", "GSWSMF")
        ) %>%
        group_by(Stream_Name) %>%
        summarise(n_days = n(), .groups = "drop") %>%
        filter(n_days >= 300) %>%
        nrow()
      n_rows      <- ceiling(n_sites_year / 3)
      plot_height <- max(6, 3 * n_rows)
      
      ggsave(
        file.path(output_path_all_sites, paste0("season_validation_WY", year, "_VERY_WET_ALL_SITES.png")),
        p,
        width = 16, height = plot_height, dpi = 300
      )
    }
  }
}

# Normal (sample 6 years)
normal_years <- years_by_class %>%
  filter(grepl("Normal", year_class)) %>%
  pull(water_year)

if (length(normal_years) >= 6) {
  
  normal_sample <- normal_years[round(seq(1, length(normal_years), length.out = 6))]
  
  for (year in normal_sample) {
    p <- plot_water_year_median(year, Q_plot_median, season_bounds)
    if (!is.null(p)) {
      ggsave(
        file.path(output_path_median, paste0("season_validation_WY", year, "_NORMAL_MEDIAN.png")),
        p,
        width = 12, height = 5, dpi = 300
      )
    }
  }
  
  for (year in normal_sample) {
    p <- plot_water_year_sites(year, Q_plot_sites, season_bounds)
    if (!is.null(p)) {
      wy_start <- as.Date(paste0(year - 1, "-10-01"))
      wy_end   <- as.Date(paste0(year,     "-09-30"))
      n_sites_year <- Q_plot_sites %>%
        filter(
          water_year == year,
          date >= wy_start, date <= wy_end,
          !Stream_Name %in% c("GSWSMA", "GSWSMF")
        ) %>%
        group_by(Stream_Name) %>%
        summarise(n_days = n(), .groups = "drop") %>%
        filter(n_days >= 300) %>%
        nrow()
      n_rows      <- ceiling(n_sites_year / 3)
      plot_height <- max(6, 3 * n_rows)
      
      ggsave(
        file.path(output_path_all_sites, paste0("season_validation_WY", year, "_NORMAL_ALL_SITES.png")),
        p,
        width = 16, height = plot_height, dpi = 300
      )
    }
  }
}

# -----------------------------------------------------------------------------
# Summary statistics
# -----------------------------------------------------------------------------
season_summary <- Q_plot_median %>%
  group_by(water_year, year_class, hydrologic_season) %>%
  summarise(n_days = n(), .groups = "drop") %>%
  filter(hydrologic_season == "Wet") %>%
  group_by(year_class) %>%
  summarise(
    n_years       = n(),
    mean_wet_days = round(mean(n_days, na.rm = TRUE)),
    sd_wet_days   = round(sd(n_days, na.rm = TRUE)),
    min_wet_days  = min(n_days, na.rm = TRUE),
    max_wet_days  = max(n_days, na.rm = TRUE),
    .groups = "drop"
  )

print(season_summary)
