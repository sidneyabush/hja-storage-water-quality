## ...existing code...
# 5_exploratory_analysis.R
# =============================================================================
# EXPLORATORY ANALYSIS: CQ BEHAVIOR
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

rm(list = ls())

# =============================================================================
# SETUP: Source helpers and configure paths
# =============================================================================

repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
try(source(file.path(repo_dir, "00_helpers", "plot_prefs.R")), silent = TRUE)
try(source(file.path(repo_dir, "00_helpers", "workflow_config.R")), silent = TRUE)
try(source(file.path(repo_dir, "00_helpers", "plot_theme_set.R")), silent = TRUE)

# =============================================================================
# SETUP: Paths & Data Loading
# =============================================================================

base_dir    <- "/Users/sidneybush/Library/CloudStorage/Box-Box"
project_dir <- file.path(base_dir, "Sidney_Bush", "HJA_Water_Quality")
out_dir     <- file.path(project_dir, "outputs")

# Create output directory structure (organized by step)
fig_base    <- file.path(project_dir, "exploratory_plots", "02_exploration", "2c_CQ_behavior_synchrony")
fig_dirs <- list(
  cq_behavior = file.path(fig_base, "cq_behavior"),
  storage     = file.path(fig_base, "storage"),
  synchrony   = file.path(fig_base, "synchrony")
)

invisible(lapply(fig_dirs, function(d) dir.create(d, showWarnings = FALSE, recursive = TRUE)))

# Helper function: save plots to appropriate directory
save_plot <- function(p, filename, scale = "cq_behavior", width = 8, height = 6) {
  dir <- fig_dirs[[scale]]
  if (is.null(dir)) stop("Unknown scale: ", scale)
  ggplot2::ggsave(
    file.path(dir, filename),
    plot = p,
    width = width,
    height = height,
    dpi = 300
  )
}

# Load all data
# Helper function to standardize stream names (GSWSMC → GSMACK)
standardize_stream_ids <- function(df) {
  df %>%
    mutate(across(where(is.character),
                  ~ case_when(. %in% c("GSWSMC", "GSWSMC_FULL") ~ "GSMACK", TRUE ~ as.character(.))))
}

seasonal <- readr::read_csv(
  file.path(out_dir, "HJA_master_seasonal.csv"),
  show_col_types = FALSE
) %>%
  standardize_stream_ids() %>%
  filter(!is.na(solute), !is.na(Stream_Name)) %>%
  filter(solute %in% solute_order) %>%
  mutate(
    hydrologic_season = factor(hydrologic_season, levels = c("Wet", "Dry")),
    Stream_Name       = factor(Stream_Name, levels = intersect(site_order, unique(Stream_Name))),
    solute            = factor(solute, levels = intersect(solute_order, unique(solute)))
  )

annual <- readr::read_csv(
  file.path(out_dir, "HJA_master_annual.csv"),
  show_col_types = FALSE
) %>%
  standardize_stream_ids() %>%
  filter(!is.na(solute), !is.na(Stream_Name)) %>%
  filter(solute %in% solute_order) %>%
  mutate(
    Stream_Name = factor(Stream_Name, levels = intersect(site_order, unique(Stream_Name))),
    solute      = factor(solute, levels = intersect(solute_order, unique(solute)))
  )

sync_annual <- readr::read_csv(
  file.path(out_dir, "HJA_composite_synchrony_annual.csv"),
  show_col_types = FALSE
) %>%
  standardize_stream_ids()

site_means <- readr::read_csv(
  file.path(out_dir, "HJA_master_site_means.csv"),
  show_col_types = FALSE
) %>%
  standardize_stream_ids() %>%
  filter(!is.na(solute), !is.na(Stream_Name)) %>%
  filter(solute %in% solute_order) %>%
  mutate(
    Stream_Name = factor(Stream_Name, levels = intersect(site_order, unique(Stream_Name))),
    solute      = factor(solute, levels = intersect(solute_order, unique(solute)))
  )

# Static catchment characteristics (available columns)
static_cols <- c(
  "Area_km2", "Elevation_mean_m", "Slope_mean",
  "Lava1_per", "Lava2_per", "Ash_Per", "Pyro_per",
  "Age", "Harvest", "Landslide_Total",
  "DR_Overall", "FYw_final"
)
static_cols <- static_cols[static_cols %in% names(site_means)]

# Join data
annual_full <- annual %>%
  left_join(sync_annual, by = c("Stream_Name", "solute", "water_year"))

seasonal_joined <- seasonal %>%
  left_join(
    site_means %>% select(Stream_Name, solute, all_of(static_cols)) %>% distinct(),
    by = c("Stream_Name", "solute")
  )

annual_joined <- annual_full %>%
  left_join(
    site_means %>% select(Stream_Name, solute, all_of(static_cols)) %>% distinct(),
    by = c("Stream_Name", "solute")
  )

# =============================================================================
# LOAD ROLLING WINDOWS FOR CQ BEHAVIOR PROPORTIONS
# =============================================================================
# The prop_enrich/prop_dilute/prop_chemostat columns may not exist in seasonal data
# Compute them from rolling windows if needed

rolling_windows <- readr::read_csv(
  file.path(out_dir, "HJA_master_rolling_windows.csv"),
  show_col_types = FALSE
) %>%
  standardize_stream_ids() %>%
  filter(!is.na(solute), !is.na(Stream_Name), !is.na(cq_behavior)) %>%
  filter(solute %in% solute_order) %>%
  mutate(
    Stream_Name = factor(Stream_Name, levels = intersect(site_order, unique(Stream_Name))),
    solute = factor(solute, levels = intersect(solute_order, unique(solute))),
    hydrologic_season = factor(hydrologic_season, levels = c("Wet", "Dry"))
  )

# Compute CQ behavior proportions at seasonal level
cq_behavior_seasonal <- rolling_windows %>%
  filter(!is.na(hydrologic_season)) %>%
  group_by(Stream_Name, solute, hydrologic_season) %>%
  summarise(
    n_windows = n(),
    prop_mobilizing = sum(cq_behavior == "mobilizing", na.rm = TRUE) / n(),
    prop_diluting = sum(cq_behavior == "diluting", na.rm = TRUE) / n(),
    prop_chemostatic = sum(cq_behavior == "chemostatic", na.rm = TRUE) / n(),
    .groups = "drop"
  )

# Compute CQ behavior proportions at site level (all seasons combined)
cq_behavior_site <- rolling_windows %>%
  group_by(Stream_Name, solute) %>%
  summarise(
    n_windows = n(),
    prop_mobilizing = sum(cq_behavior == "mobilizing", na.rm = TRUE) / n(),
    prop_diluting = sum(cq_behavior == "diluting", na.rm = TRUE) / n(),
    prop_chemostatic = sum(cq_behavior == "chemostatic", na.rm = TRUE) / n(),
    .groups = "drop"
  )

message("Data loaded successfully.")

# =============================================================================
# PLOTTING THEMES & HELPER FUNCTIONS
# =============================================================================

theme_seasonal_plot <- theme_bw(base_size = 12) +
  theme_hja() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right", axis.line.x = element_line(color = "black", linewidth = 0.4))

theme_annual_plot <- theme_bw(base_size = 12) +
  theme_hja() +
  theme(legend.position = "right", axis.line.x = element_line(color = "black", linewidth = 0.4))

theme_site_plot <- theme_bw(base_size = 12) +
  theme_hja() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right", axis.line.x = element_line(color = "black", linewidth = 0.4))

# Generic scatter plot helper (works for all scales)
plot_xy <- function(df, x_var, y_var,
                    facet_var = "solute",
                    color_var = "Stream_Name",
                    shape_var = NULL,
                    theme    = theme_seasonal_plot,
                    scale    = "seasonal",
                    title_prefix = "") {
  
  needed <- c(x_var, y_var, facet_var, color_var)
  if (!is.null(shape_var)) needed <- c(needed, shape_var)
  
  if (!all(needed %in% names(df))) {
    warning("Skipping x=", x_var, " y=", y_var, " (missing columns)")
    return(invisible(NULL))
  }
  
  aes_map <- ggplot2::aes(
    x = .data[[x_var]],
    y = .data[[y_var]],
    color = .data[[color_var]]
  )
  
  if (!is.null(shape_var)) {
    aes_map <- ggplot2::aes(
      x = .data[[x_var]],
      y = .data[[y_var]],
      color = .data[[color_var]],
      shape = .data[[shape_var]]
    )
  }
  
  df_flagged <- flag_outlet_stream(df, color_var)
  p <- df_flagged %>%
    ggplot(aes_map)
  if ("is_outlet" %in% names(df_flagged) && any(df_flagged$is_outlet, na.rm = TRUE)) {
    non_outlet <- df_flagged[is.na(df_flagged$is_outlet) | !df_flagged$is_outlet, , drop = FALSE]
    outlet_df <- df_flagged[df_flagged$is_outlet %in% TRUE, , drop = FALSE]
    if (nrow(non_outlet) > 0) {
      p <- p + geom_point(data = non_outlet, alpha = 0.6, size = 2)
    }
    if (nrow(outlet_df) > 0) {
      p <- p + geom_point(
        data = outlet_df,
        alpha = 0.95,
        size = 3,
        shape = 17,
        stroke = 0.6
      )
    }
  } else {
    p <- p + geom_point(alpha = 0.6, size = 2)
  }
  p <- p +
    facet_wrap(as.formula(paste("~", facet_var)), scales = "free") +
    scale_color_site_cq() +
    labs(
      x = x_var,
      y = y_var,
      color = "Site",
      title = paste0(title_prefix, y_var, " vs ", x_var)
    ) +
    theme +
    theme(axis.line.x = element_line(color = "black", linewidth = 0.4))
  
  if (!is.null(shape_var)) {
    p <- p + labs(shape = "Season")
  }
  
  fname <- paste0(
    tolower(gsub("_", "", y_var)), "_vs_", tolower(gsub("_", "", x_var)), ".png"
  )
  
  save_plot(p, fname, scale = scale)
  invisible(p)
}

# =============================================================================
# SCALE 1: SEASONAL (90-day scale)
# =============================================================================
# Questions:
#   - How do CQ slopes vary seasonally across sites?
#   - What is the composition of CQ behaviors (mobilizing/diluting/chemostatic)?
#   - How variable is within-site solute synchrony by season?
#   - How do storage metrics explain seasonal variability?


# 1.1 Basic seasonal summaries

# CQ slopes by site & season
if ("cq_slope" %in% names(seasonal)) {
  p <- seasonal %>%
    filter(!is.na(cq_slope)) %>%
    ggplot(aes(x = Stream_Name, y = cq_slope, fill = hydrologic_season)) +
    geom_boxplot(outlier.alpha = 0.4, alpha = 0.7) +
    facet_wrap(~ solute, scales = "free_y") +
    scale_fill_season() +
    labs(
      x = "Site", y = "Window C–Q slope",
      title = "Seasonal CQ slopes by site and solute"
    ) +
    theme_seasonal_plot +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.line.x = element_line(color = "black", linewidth = 0.4))
  save_plot(p, "01_seasonal_cqslope_by_site_solute.png", scale = "cq_behavior", width = 12, height = 8)
}

# CQ behavior composition (mobilizing/diluting/chemostatic) - using computed proportions
if (nrow(cq_behavior_seasonal) > 0) {
    p <- cq_behavior_seasonal %>%
      apply_factor_orders() %>%
      pivot_longer(cols = starts_with("prop_"), names_to = "behavior", values_to = "prop") %>%
      mutate(behavior = str_remove(behavior, "prop_"),
             behavior = factor(behavior, levels = cq_behavior_order)) %>%
      ggplot(aes(x = Stream_Name, y = prop, fill = behavior)) +
      geom_col(position = "fill") +
      facet_grid(hydrologic_season ~ solute, switch = "x") +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      scale_fill_cq_behavior() +
      labs(
        x = "Site", y = "Fraction of 90-day windows",
        title = "Seasonal CQ behavior composition by site and solute"
      ) +
      theme_seasonal_plot +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "bottom",
        axis.title.x = element_text(),
        
        # REMOVE SILVER/GREY BETWEEN FACETS ⬇️
        panel.background = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        panel.spacing.x = unit(0, "points"),   # no vertical white gaps
        panel.spacing.y = unit(4, "points"),   # keep a little separation horizontally
        
        axis.line.x = element_line(color = "black", linewidth = 0.4)
      )
    
    save_plot(p, "02_seasonal_cq_behavior_composition.png", scale = "cq_behavior", width = 12, height = 9)
}

# =============================================================================
# ADDITIONAL CQ BEHAVIOR FIGURES (standalone, no synchrony)
# =============================================================================

# 1. CQ slope distributions by solute (violin + boxplot)
if ("cq_slope" %in% names(seasonal)) {
  p <- seasonal %>%
    apply_factor_orders() %>%
    ggplot(aes(x = solute, y = cq_slope, fill = solute)) +
    geom_violin(alpha = 0.5, scale = "width") +
    geom_boxplot(width = 0.2, outlier.alpha = 0.3) +
    geom_hline(yintercept = 0.1, linetype = "dashed", color = "#CCCCCC", alpha = 0.9) +
    geom_hline(yintercept = -0.1, linetype = "dashed", color = "#CCCCCC", alpha = 0.9) +
    scale_fill_solute() +
    labs(
      x = "Solute", y = "C–Q slope",
      title = "Distribution of C–Q slopes across all sites and seasons",
      subtitle = "Gray dashed lines = ±0.1"
    ) +
    theme_hja() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p, "03_cqslope_distribution_by_solute.png", scale = "cq_behavior", width = 10, height = 6)
}

# 2. CQ slope heatmap: site × solute (mean values)
if ("cq_slope" %in% names(seasonal)) {
  slope_matrix <- seasonal %>%
    group_by(Stream_Name, solute) %>%
    summarise(mean_slope = mean(cq_slope, na.rm = TRUE), .groups = "drop")
  
  p <- slope_matrix %>%
    ggplot(aes(x = solute, y = Stream_Name, fill = mean_slope)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         name = "Mean C–Q slope") +
    labs(
      x = "Solute", y = "Site",
      title = "Mean C–Q slope by site and solute",
      subtitle = "Blue = dilution, Red = mobilization, White = chemostatic"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p, "04_cqslope_heatmap_site_solute.png", scale = "cq_behavior", width = 10, height = 6)
}

# 3. CQ behavior classification by site (overall, not by season) - using computed proportions
if (nrow(cq_behavior_site) > 0) {
  behavior_overall <- cq_behavior_site %>%
    apply_factor_orders() %>%
    pivot_longer(cols = starts_with("prop_"), 
                 names_to = "behavior", values_to = "prop") %>%
    mutate(behavior = str_remove(behavior, "prop_"),
           behavior = factor(behavior, levels = cq_behavior_order))
  
  p <- behavior_overall %>%
    ggplot(aes(x = solute, y = prop, fill = behavior)) +
    geom_col(position = "stack") +
    facet_wrap(~ Stream_Name, switch = "x") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_fill_cq_behavior() +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      x = "Solute", y = "Proportion",
      title = "CQ behavior composition by site (all seasons combined)"
    ) +
    theme_hja() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8), axis.title.x = element_text(), axis.line.x = element_line(color = "black", linewidth = 0.4))
  save_plot(p, "05_cq_behavior_by_site_stacked.png", scale = "cq_behavior", width = 12, height = 8)
}

# 4. Dominant CQ behavior by site × solute - using computed proportions
if (nrow(cq_behavior_site) > 0) {
  dominant_behavior <- cq_behavior_site %>%
    apply_factor_orders() %>%
    mutate(dominant = case_when(
      prop_mobilizing >= prop_diluting & prop_mobilizing >= prop_chemostatic ~ "mobilizing",
      prop_diluting >= prop_mobilizing & prop_diluting >= prop_chemostatic ~ "diluting",
      TRUE ~ "chemostatic"
    ),
    dominant = factor(dominant, levels = cq_behavior_order))
  
  p <- dominant_behavior %>%
    ggplot(aes(x = solute, y = Stream_Name, fill = dominant)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_cq_behavior(name = "Dominant Behavior") +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    labs(
      x = "Solute", y = "Site",
      title = "Dominant CQ behavior by site and solute"
    ) +
    theme_hja() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_text(), axis.line.x = element_line(color = "black", linewidth = 0.4))
  save_plot(p, "06_dominant_cq_behavior_heatmap.png", scale = "cq_behavior", width = 10, height = 6)
}

# 5. CVc/CVq ratio (if available) - chemodynamic vs chemostatic indicator
if ("cvc_cvq_ratio" %in% names(seasonal)) {
  p <- seasonal %>%
    ggplot(aes(x = solute, y = cvc_cvq_ratio, fill = hydrologic_season)) +
    geom_boxplot(outlier.alpha = 0.3) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
    scale_y_log10() +
    labs(
      x = "Solute", y = "CVc/CVq ratio (log scale)", fill = "Season",
      title = "CVc/CVq ratio by solute and season",
      subtitle = "Above 1 = chemodynamic, Below 1 = chemostatic"
    ) +
    theme_bw(base_size = 12) +
    theme_hja() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line.x = element_line(color = "black", linewidth = 0.4))
  save_plot(p, "07_cvcq_ratio_by_solute_season.png", scale = "cq_behavior", width = 10, height = 6)
}

# 6. CQ slope temporal variability by site
if ("cq_slope" %in% names(seasonal) && "water_year" %in% names(seasonal)) {
  p <- seasonal %>%
    ggplot(aes(x = factor(water_year), y = cq_slope, color = hydrologic_season)) +
    geom_point(alpha = 0.5, position = position_jitter(width = 0.2)) +
    stat_summary(fun = mean, geom = "line", aes(group = hydrologic_season), linewidth = 1) +
    facet_grid(solute ~ Stream_Name, scales = "free_y") +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    labs(
      x = "Water Year", y = "C–Q slope", color = "Season",
      title = "CQ slope temporal patterns by site and solute"
    ) +
    theme_bw(base_size = 9) +
    theme_hja() +
    scale_color_manual(values = c(Wet = "#56B4E9", Dry = "#E69F00"), name = "Season") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6), strip.text = element_text(size = 7), axis.line.x = element_line(color = "black", linewidth = 0.4))
  save_plot(p, "08_cqslope_temporal_by_site_solute.png", scale = "cq_behavior", width = 14, height = 12)
}

# 7. Wet vs Dry season CQ slope comparison (paired)
if ("cq_slope" %in% names(seasonal) && "hydrologic_season" %in% names(seasonal)) {
  wet_dry_compare <- seasonal %>%
    select(Stream_Name, solute, water_year, hydrologic_season, cq_slope) %>%
    pivot_wider(names_from = hydrologic_season, values_from = cq_slope, values_fn = mean) %>%
    filter(!is.na(Wet), !is.na(Dry))
  
  if (nrow(wet_dry_compare) > 0) {
    p <- wet_dry_compare %>%
      ggplot(aes(x = Wet, y = Dry, color = solute)) +
      geom_point(alpha = 0.6) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
      facet_wrap(~ Stream_Name) +
      labs(
        x = "Wet season C–Q slope", y = "Dry season C–Q slope", color = "Solute",
        title = "Wet vs Dry season CQ slopes by site",
        subtitle = "Points above diagonal = steeper slope in dry season"
      ) +
      scale_color_solute() +
      theme_hja() +
      theme(axis.line.x = element_line(color = "black", linewidth = 0.4))
    save_plot(p, "09_cqslope_wet_vs_dry_scatter.png", scale = "cq_behavior", width = 12, height = 8)

    # === CUSTOM: Seasonal Mixing Regime Shifts (solutes on y-axis, boxplot, gray dashed lines at ±0.1, no grid, no border) ===
    wet_dry_compare_long <- wet_dry_compare %>%
      mutate(seasonal_shift = Wet - Dry)
    p2 <- wet_dry_compare_long %>%
      apply_factor_orders() %>%
      ggplot(aes(y = solute, x = seasonal_shift)) +
      geom_boxplot(fill = "#B3CDE3", color = "black", outlier.alpha = 0.6) +
      geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
      labs(
        title = "Seasonal Mixing Regime Shifts",
        subtitle = "Positive = more dilution in wet season",
        x = "Seasonal Shift (Wet - Dry CQ Slope)",
        y = "Solute"
      ) +
      theme_hja()
    save_plot(p2, "custom_seasonal_mixing_regime_shifts.png", scale = "cq_behavior", width = 9, height = 7)

    # Custom CQ shift plot (site-faceted)
    p3 <- wet_dry_compare_long %>%
      apply_factor_orders() %>%
      ggplot(aes(y = solute, x = seasonal_shift)) +
      geom_boxplot(fill = "#B3CDE3", color = "#4682B4", outlier.alpha = 0.6) +
      geom_vline(xintercept = 0.1, linetype = "dashed", color = "#CCCCCC", alpha = 0.9) +
      geom_vline(xintercept = -0.1, linetype = "dashed", color = "#CCCCCC", alpha = 0.9) +
      facet_wrap(~ Stream_Name) +
      labs(
        title = "Seasonal Mixing Regime Shifts by Site",
        subtitle = "Positive = more dilution in wet season",
        x = "Seasonal Shift (Wet - Dry CQ Slope)",
        y = "Solute"
      ) +
      theme_hja()
    save_plot(p3, "custom_seasonal_mixing_regime_shifts_by_site.png", scale = "cq_behavior", width = 14, height = 10)
  }
}

# 8. CQ slope variability summary (CV or SD by site)
if ("cq_slope" %in% names(seasonal)) {
  slope_variability <- seasonal %>%
    group_by(Stream_Name, solute) %>%
    summarise(
      mean_slope = mean(cq_slope, na.rm = TRUE),
      sd_slope = sd(cq_slope, na.rm = TRUE),
      cv_slope = sd_slope / abs(mean_slope),
      .groups = "drop"
    )
  
  p <- slope_variability %>%
    ggplot(aes(x = solute, y = Stream_Name, fill = sd_slope)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(option = "plasma", name = "SD of C–Q slope") +
    labs(
      x = "Solute", y = "Site",
      title = "Temporal variability in CQ slopes (SD)",
      subtitle = "Higher values = more variable CQ behavior over time"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p, "10_cqslope_variability_heatmap.png", scale = "cq_behavior", width = 10, height = 6)
}

# Within-site solute synchrony
if ("prop_pair_sync" %in% names(seasonal)) {
  p <- seasonal %>%
    ggplot(aes(x = Stream_Name, y = prop_pair_sync, color = hydrologic_season, group = hydrologic_season)) +
    geom_point(position = position_dodge(width = 0.4), size = 2) +
    facet_wrap(~ solute) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      x = "Site", y = "Within-site solute synchrony (% pairs)",
      color = "Season",
      title = "Seasonal within-site solute synchrony"
    ) +
    theme_seasonal_plot
  save_plot(p, "03_seasonal_within_site_pair_sync.png", scale = "synchrony", width = 12, height = 8)
}

# 1.2 Relationships with storage & hydro metrics
cq_vars_seas <- c("cq_slope", "prop_enrich", "prop_dilute", "prop_chemostat")
cq_vars_seas <- cq_vars_seas[cq_vars_seas %in% names(seasonal_joined)]

sync_vars_seas <- c("prop_pair_sync", "conc_sync_allpairs", "cqslope_sync_allpairs")
sync_vars_seas <- sync_vars_seas[sync_vars_seas %in% names(seasonal_joined)]

primary_storage_seas <- intersect(PRIMARY_STORAGE_METRIC, names(seasonal_joined))
supp_storage_seas <- setdiff(
  intersect(SUPPLEMENTAL_STORAGE_METRICS, names(seasonal_joined)),
  primary_storage_seas
)

# CQ vs storage
for (y_var in cq_vars_seas) {
  if (length(primary_storage_seas) == 1) {
    plot_xy(seasonal_joined, primary_storage_seas, y_var,
            facet_var = "solute", color_var = "Stream_Name",
            shape_var = "hydrologic_season",
            theme = theme_seasonal_plot, scale = "cq_behavior",
            title_prefix = "SEASONAL (Primary): ")
  }
  if (length(supp_storage_seas) > 0) {
    for (x_var in supp_storage_seas) {
      plot_xy(seasonal_joined, x_var, y_var,
              facet_var = "solute", color_var = "Stream_Name",
              shape_var = "hydrologic_season",
              theme = theme_seasonal_plot, scale = "storage",
              title_prefix = "SEASONAL (Supplemental): ")
    }
  }
}

# Synchrony vs storage
for (y_var in sync_vars_seas) {
  if (length(primary_storage_seas) == 1) {
    plot_xy(seasonal_joined, primary_storage_seas, y_var,
            facet_var = "solute", color_var = "Stream_Name",
            shape_var = "hydrologic_season",
            theme = theme_seasonal_plot, scale = "synchrony",
            title_prefix = "SEASONAL (Primary): ")
  }
  if (length(supp_storage_seas) > 0) {
    for (x_var in supp_storage_seas) {
      plot_xy(seasonal_joined, x_var, y_var,
              facet_var = "solute", color_var = "Stream_Name",
              shape_var = "hydrologic_season",
              theme = theme_seasonal_plot, scale = "storage",
              title_prefix = "SEASONAL (Supplemental): ")
    }
  }
}

# 1.3 Relationships with catchment characteristics
for (y_var in cq_vars_seas) {
  for (x_var in static_cols[1:4]) {  # Limit to avoid explosion of plots
    plot_xy(seasonal_joined, x_var, y_var,
            facet_var = "solute", color_var = "Stream_Name",
            shape_var = "hydrologic_season",
            theme = theme_seasonal_plot, scale = "cq_behavior",
            title_prefix = "SEASONAL: ")
  }
}

# =============================================================================
# SCALE 2: ANNUAL (year-to-year variation)
# =============================================================================
# Questions:
#   - How do annual CQ slopes and synchrony vary by year?
#   - Do all sites show consistent patterns in synchrony across years?
#   - How do annual storage ranges relate to synchrony?
#   - Which sites show synchrony, and which are asynchronous?


# 2.1 Time series of key synchrony metrics

if ("cqslope_sync_allpairs" %in% names(annual_joined)) {
  p <- annual_joined %>%
    filter(!is.na(cqslope_sync_allpairs)) %>%
    ggplot(aes(x = water_year, y = cqslope_sync_allpairs, color = solute)) +
    geom_line(alpha = 0.7, linewidth = 0.8) +
    geom_point(size = 1.5) +
    facet_wrap(~ Stream_Name) +
    labs(
      x = "Water Year", y = get_sync_label("cqslope_sync_allpairs"),
      color = "Solute",
      title = "Annual absolute CQ-slope synchrony by site and solute"
    ) +
    theme_annual_plot
  save_plot(p, "01_annual_ts_cqslope_sync.png", scale = "synchrony", width = 12, height = 8)
}

if ("conc_sync_allpairs" %in% names(annual_joined)) {
  p <- annual_joined %>%
    filter(!is.na(conc_sync_allpairs)) %>%
    ggplot(aes(x = water_year, y = conc_sync_allpairs, color = solute)) +
    geom_line(alpha = 0.7, linewidth = 0.8) +
    geom_point(size = 1.5) +
    facet_wrap(~ Stream_Name) +
    labs(
      x = "Water Year", y = get_sync_label("conc_sync_allpairs"),
      color = "Solute",
      title = "Annual absolute concentration synchrony by site and solute"
    ) +
    theme_annual_plot
  save_plot(p, "02_annual_ts_conc_sync.png", scale = "synchrony", width = 12, height = 8)
}

# 2.2 Annual relationships: CQ & synchrony vs storage

sync_vars_ann <- c("conc_sync_allpairs", "cqslope_sync_allpairs",
                   "conc_sync_outlet", "cqslope_sync_outlet", "prop_sync")
sync_vars_ann <- sync_vars_ann[sync_vars_ann %in% names(annual_joined)]

cq_vars_ann <- c("cq_slope", "mean_C", "median_C")
cq_vars_ann <- cq_vars_ann[cq_vars_ann %in% names(annual_joined)]

primary_storage_ann <- intersect(PRIMARY_STORAGE_METRIC, names(annual_joined))
supp_storage_ann <- setdiff(
  intersect(SUPPLEMENTAL_STORAGE_METRICS, names(annual_joined)),
  primary_storage_ann
)

for (y_var in sync_vars_ann) {
  if (length(primary_storage_ann) == 1) {
    plot_xy(annual_joined, primary_storage_ann, y_var,
            facet_var = "solute", color_var = "Stream_Name",
            theme = theme_annual_plot, scale = "synchrony",
            title_prefix = "ANNUAL (Primary): ")
  }
  if (length(supp_storage_ann) > 0) {
    for (x_var in supp_storage_ann) {
      plot_xy(annual_joined, x_var, y_var,
              facet_var = "solute", color_var = "Stream_Name",
              theme = theme_annual_plot, scale = "storage",
              title_prefix = "ANNUAL (Supplemental): ")
    }
  }
}

for (y_var in cq_vars_ann) {
  if (length(primary_storage_ann) == 1) {
    plot_xy(annual_joined, primary_storage_ann, y_var,
            facet_var = "solute", color_var = "Stream_Name",
            theme = theme_annual_plot, scale = "cq_behavior",
            title_prefix = "ANNUAL (Primary): ")
  }
  if (length(supp_storage_ann) > 0) {
    for (x_var in supp_storage_ann) {
      plot_xy(annual_joined, x_var, y_var,
              facet_var = "solute", color_var = "Stream_Name",
              theme = theme_annual_plot, scale = "storage",
              title_prefix = "ANNUAL (Supplemental): ")
    }
  }
}

# 2.3 Abbott vs Wymore synchrony comparison
if (all(c("cqslope_sync_allpairs", "prop_sync") %in% names(annual_joined))) {
  p <- annual_joined %>%
    filter(!is.na(cqslope_sync_allpairs), !is.na(prop_sync)) %>%
    ggplot(aes(x = prop_sync, y = cqslope_sync_allpairs, color = solute)) +
    geom_point(alpha = 0.7, size = 2) +
    facet_wrap(~ Stream_Name) +
    labs(
      x = "Wymore sync fraction", y = get_sync_label("cqslope_sync_allpairs"),
      color = "Solute",
      title = "Annual synchrony: Abbott vs Wymore methods"
    ) +
    theme_annual_plot
  save_plot(p, "10_annual_abbott_vs_wymore.png", scale = "synchrony", width = 12, height = 8)
}

# =============================================================================
# SCALE 3: SITE-MEANS (long-term cross-site patterns)
# =============================================================================
# Questions:
#   - Do all sites have similar long-term synchrony?
#   - How do storage ranges (Q_dS_range, WB_dS_range) explain synchrony?
#   - Which sites are most/least synchronous?
#   - How do catchment characteristics (area, forest type) relate to synchrony?

# 3.1 Basic site-level comparison

# Cross-stream synchrony by site
if ("conc_sync_allpairs" %in% names(site_means)) {
  p <- site_means %>%
    ggplot(aes(x = Stream_Name, y = conc_sync_allpairs, fill = solute)) +
    geom_col(position = "dodge") +
    labs(
      x = "Site", y = get_sync_label("conc_sync_allpairs"),
      fill = "Solute",
      title = "Cross-site variation in absolute concentration synchrony"
    ) +
    theme_site_plot
  save_plot(p, "01_site_conc_sync_comparison.png", scale = "synchrony", width = 10, height = 6)
}

if ("cqslope_sync_allpairs" %in% names(site_means)) {
  p <- site_means %>%
    ggplot(aes(x = Stream_Name, y = cqslope_sync_allpairs, fill = solute)) +
    geom_col(position = "dodge") +
    labs(
      x = "Site", y = get_sync_label("cqslope_sync_allpairs"),
      fill = "Solute",
      title = "Cross-site variation in absolute CQ-slope synchrony"
    ) +
    theme_site_plot
  save_plot(p, "02_site_cqslope_sync_comparison.png", scale = "synchrony", width = 10, height = 6)
}

# 3.2 Synchrony vs storage metrics (KEY QUESTION)
sync_vars_site <- c("conc_sync_allpairs", "cqslope_sync_allpairs",
                    "conc_sync_outlet", "cqslope_sync_outlet", "prop_pair_sync")
sync_vars_site <- sync_vars_site[sync_vars_site %in% names(site_means)]

primary_storage_site <- intersect(PRIMARY_STORAGE_METRIC, names(site_means))
supp_storage_site <- setdiff(
  intersect(SUPPLEMENTAL_STORAGE_METRICS, names(site_means)),
  primary_storage_site
)

for (y_var in sync_vars_site) {
  if (length(primary_storage_site) == 1) {
    plot_xy(site_means, primary_storage_site, y_var,
            facet_var = "solute", color_var = "Stream_Name",
            theme = theme_site_plot, scale = "synchrony",
            title_prefix = "SITE-MEAN (Primary): ")
  }
  if (length(supp_storage_site) > 0) {
    for (x_var in supp_storage_site) {
      plot_xy(site_means, x_var, y_var,
              facet_var = "solute", color_var = "Stream_Name",
              theme = theme_site_plot, scale = "storage",
              title_prefix = "SITE-MEAN (Supplemental): ")
    }
  }
}

# 3.3 Synchrony vs catchment characteristics

for (y_var in sync_vars_site) {
  for (x_var in static_cols[1:6]) {  # Sample to avoid explosion
    plot_xy(site_means, x_var, y_var,
            facet_var = "solute", color_var = "Stream_Name",
            theme = theme_site_plot, scale = "synchrony",
            title_prefix = "SITE-MEAN: ")
  }
}

# 3.4 CQ behavior vs storage/catchment

cq_vars_site <- c("cq_slope", "prop_enrich", "prop_dilute", "prop_chemostatic")
cq_vars_site <- cq_vars_site[cq_vars_site %in% names(site_means)]

for (y_var in cq_vars_site) {
  if (length(primary_storage_site) == 1) {
    plot_xy(site_means, primary_storage_site, y_var,
            facet_var = "solute", color_var = "Stream_Name",
            theme = theme_site_plot, scale = "cq_behavior",
            title_prefix = "SITE-MEAN (Primary): ")
  }
  if (length(supp_storage_site) > 0) {
    for (x_var in supp_storage_site) {
      plot_xy(site_means, x_var, y_var,
              facet_var = "solute", color_var = "Stream_Name",
              theme = theme_site_plot, scale = "storage",
              title_prefix = "SITE-MEAN (Supplemental): ")
    }
  }
}

# Calculate wet_dry_compare and generate custom CQ shift plots
if (all(c("Stream_Name", "solute", "water_year", "hydrologic_season", "cq_slope") %in% names(seasonal))) {
  wet_dry_compare <- seasonal %>%
    select(Stream_Name, solute, water_year, hydrologic_season, cq_slope) %>%
    pivot_wider(names_from = hydrologic_season, values_from = cq_slope, values_fn = mean) %>%
    filter(!is.na(Wet), !is.na(Dry))
  print("wet_dry_compare head:")
  print(head(wet_dry_compare))

  # Custom CQ shift plot (solutes on y-axis)
  wet_dry_compare_long <- wet_dry_compare %>%
    mutate(seasonal_shift = Wet - Dry)
  p2 <- wet_dry_compare_long %>%
    apply_factor_orders() %>%
    ggplot(aes(y = solute, x = seasonal_shift)) +
    geom_boxplot(fill = "#B3CDE3", color = "black", outlier.alpha = 0.6) +
    geom_vline(xintercept = 0.1, linetype = "dashed", color = "#CCCCCC", alpha = 0.9) +
    geom_vline(xintercept = -0.1, linetype = "dashed", color = "#CCCCCC", alpha = 0.9) +
    labs(
      title = "Seasonal Mixing Regime Shifts",
      subtitle = "Positive = more dilution in wet season",
      x = "Seasonal Shift (Wet - Dry CQ Slope)",
      y = "Solute"
    ) +
    theme_hja()
  save_plot(p2, "custom_seasonal_mixing_regime_shifts.png", scale = "cq_behavior", width = 9, height = 7)

  # Custom CQ shift plot (site-faceted)
  p3 <- wet_dry_compare_long %>%
    apply_factor_orders() %>%
    ggplot(aes(y = solute, x = seasonal_shift)) +
    geom_boxplot(fill = "#B3CDE3", color = "black", outlier.alpha = 0.6) +
    geom_vline(xintercept = 0.1, linetype = "dashed", color = "#CCCCCC", alpha = 0.9) +
    geom_vline(xintercept = -0.1, linetype = "dashed", color = "#CCCCCC", alpha = 0.9) +
    facet_wrap(~ Stream_Name) +
    labs(
      title = "Seasonal Mixing Regime Shifts by Site",
      subtitle = "Positive = more dilution in wet season",
      x = "Seasonal Shift (Wet - Dry CQ Slope)",
      y = "Solute"
    ) +
    theme_hja()
  save_plot(p3, "custom_seasonal_mixing_regime_shifts_by_site.png", scale = "cq_behavior", width = 14, height = 10)
}
