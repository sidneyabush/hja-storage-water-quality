# =============================================================================
# STEP 03c: COMPREHENSIVE CORRELATION ANALYSIS
# =============================================================================
# Goal: Identify key relationships between:
#   - Catchment characteristics (slope, elevation, harvest, geology)
#   - Hydrologic metrics (RBI, MTT, FYw, DR, storage)
#   - CQ behavior (slope, synchrony)
#   - Cluster membership and stability
#
# Analyze at multiple scales:
#   - Site-level (long-term means)
#   - Annual (water year aggregates)
#   - Seasonal (wet vs dry)
#   - Window-level (90-day rolling)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(corrplot)
  library(patchwork)
})

rm(list = ls())

# Source helpers
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
try(source(file.path(repo_dir, "00_helpers", "plot_prefs.R")), silent = TRUE)
try(source(file.path(repo_dir, "00_helpers", "workflow_config.R")), silent = TRUE)

# Paths
base_dir   <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality"
output_dir <- file.path(base_dir, "outputs")
plot_dir   <- file.path(base_dir, "exploratory_plots", "03_stats", "3c_correlations")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

fig_dir <- file.path(base_dir, "exploratory_plots", "02_exploration", "3c_correlations")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Theme
theme_clean <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(panel.grid = element_blank(), strip.background = element_blank())
}

message("\n=== LOADING DATA ===\n")

# Load all datasets
catchment <- read_csv(file.path(output_dir, "Catchment_site_characteristics.csv"), show_col_types = FALSE)
site_data <- read_csv(file.path(output_dir, "HJA_exploratory_site.csv"), show_col_types = FALSE)
mega <- read_csv(file.path(output_dir, "HJA_mega_90d_windows_CQ_hydro_static.csv"), show_col_types = FALSE)
clusters_wy <- read_csv(file.path(output_dir, "ClusterStreams_allSolutes_byWaterYear.csv"), show_col_types = FALSE)
clusters_modal <- read_csv(file.path(output_dir, "ClusterStreams_allSolutes_modalClusters.csv"), show_col_types = FALSE)
stability <- read_csv(file.path(output_dir, "ClusterStreams_stability_metrics.csv"), show_col_types = FALSE)
season_bounds <- tryCatch(read_csv(file.path(output_dir, "season_boundaries.csv"), show_col_types = FALSE), error = function(e) NULL)

primary_storage_metric <- if (exists("PRIMARY_STORAGE_METRIC")) PRIMARY_STORAGE_METRIC else "WB_dS_range_mm"
storage_candidates <- unique(c(primary_storage_metric, "WB_dS_range_mm", "Q_dS_range_mm"))

choose_storage_col <- function(df) {
  matches <- storage_candidates[storage_candidates %in% names(df)]
  if (length(matches) == 0) return(NA_character_)
  matches[[1]]
}

storage_col_site <- choose_storage_col(site_data)
storage_col_mega <- choose_storage_col(mega)
storage_col_clusters <- choose_storage_col(clusters_wy)

storage_label_site <- get_storage_label(ifelse(is.na(storage_col_site), primary_storage_metric, storage_col_site))
storage_label_mega <- get_storage_label(ifelse(is.na(storage_col_mega), primary_storage_metric, storage_col_mega))

# #############################################################################
# SCALE 1: SITE-LEVEL CORRELATIONS
# #############################################################################
message("\n=== SCALE 1: SITE-LEVEL CORRELATIONS ===\n")

# Build site-level summary
site_summary <- site_data %>%
  group_by(Stream_Name) %>%
  summarise(
    mean_cq_slope = mean(cqslope_sync_allpairs, na.rm = TRUE),
    mean_RBI = mean(RBI, na.rm = TRUE),
    mean_storage = if (!is.na(storage_col_site)) mean(.data[[storage_col_site]], na.rm = TRUE) else NA_real_,
    mean_DR = mean(DR_Overall, na.rm = TRUE),
    mean_MTT = mean(MTT_final, na.rm = TRUE),
    mean_FYw = mean(FYw_final, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(catchment, by = "Stream_Name") %>%
  left_join(
    stability %>% group_by(Stream_Name) %>% summarise(mean_stability = mean(stability, na.rm = TRUE), .groups = "drop"),
    by = "Stream_Name"
  )

# Define variable groups
catchment_vars <- c("Area_km2", "Elevation_mean_m", "Slope_mean", "Harvest", "Age", "Lava1_per", "Ash_Per")
hydro_vars <- c("mean_RBI", "mean_storage", "mean_DR", "mean_MTT", "mean_FYw", "DR_Overall", "MTT_overall", "Fyw_overall")
cq_vars <- c("mean_cq_slope", "mean_stability")

# Filter to available variables
all_vars <- c(catchment_vars, hydro_vars, cq_vars)
all_vars <- all_vars[all_vars %in% names(site_summary)]
all_vars <- all_vars[!duplicated(all_vars)]

# Create correlation matrix
site_corr_data <- site_summary %>%
  select(all_of(all_vars)) %>%
  drop_na()

if (nrow(site_corr_data) >= 5) {
  site_corr <- cor(site_corr_data, use = "pairwise.complete.obs")
  
  # Save correlation plot
  png(file.path(plot_dir, "01_site_level_correlations.png"), width = 12, height = 10, units = "in", res = 300)
  corrplot(site_corr, method = "color", type = "full",
           addCoef.col = "black", number.cex = 0.7,
           tl.col = "black", tl.srt = 45, tl.cex = 0.8,
           title = "Site-Level Correlations (n = 8-9 sites)", mar = c(0,0,2,0))
  dev.off()
  
  # Export correlation table
  site_corr_df <- as.data.frame(site_corr) %>%
    rownames_to_column("Variable1") %>%
    pivot_longer(-Variable1, names_to = "Variable2", values_to = "r") %>%
    filter(Variable1 != Variable2)
  
  write_csv(site_corr_df, file.path(output_dir, "03_stats/site_level_correlations.csv"))
  
  message("  Site-level correlations saved\n")
}

# #############################################################################
# SCALE 2: WINDOW-LEVEL CORRELATIONS (Rolling 90-day)
# #############################################################################
message("\n=== SCALE 2: WINDOW-LEVEL CORRELATIONS ===\n")

# Key variables from mega dataset
window_vars <- c("cq_slope.x", storage_col_mega, "RBI", "DR_Overall", "MTT_final", "FYw_final")
window_vars <- window_vars[window_vars %in% names(mega)]

if (length(window_vars) >= 3) {
  # Overall window correlations
  window_corr_data <- mega %>%
    select(all_of(window_vars)) %>%
    drop_na()
  
  if (nrow(window_corr_data) > 100) {
    window_corr <- cor(window_corr_data, use = "pairwise.complete.obs")
    
    png(file.path(plot_dir, "02_window_level_correlations_overall.png"), width = 10, height = 10, units = "in", res = 300)
    corrplot(window_corr, method = "color", type = "upper",
             addCoef.col = "black", number.cex = 0.9,
             tl.col = "black", tl.srt = 45,
             title = paste0("Window-Level Correlations (n = ", nrow(window_corr_data), ")"), mar = c(0,0,2,0))
    dev.off()
  }
  
  # By season
  if ("hydrologic_season" %in% names(mega)) {
    for (season in c("wet", "dry")) {
      s_data <- mega %>%
        filter(hydrologic_season == season) %>%
        select(all_of(window_vars)) %>%
        drop_na()
      
      if (nrow(s_data) > 100) {
        s_corr <- cor(s_data, use = "pairwise.complete.obs")
        
        png(file.path(plot_dir, paste0("02_window_level_correlations_", season, ".png")), 
            width = 800, height = 800, res = 100)
        corrplot(s_corr, method = "color", type = "upper",
                 addCoef.col = "black", number.cex = 0.9,
                 tl.col = "black", tl.srt = 45,
                 title = paste0(tools::toTitleCase(season), " Season (n = ", nrow(s_data), ")"), mar = c(0,0,2,0))
        dev.off()
      }
    }
  }
  
  message("  Window-level correlations saved\n")
}

# #############################################################################
# SCALE 3: ANNUAL CORRELATIONS
# #############################################################################
message("\n=== SCALE 3: ANNUAL CORRELATIONS ===\n")

# Aggregate mega to annual
annual_summary <- mega %>%
  group_by(Stream_Name, solute, water_year) %>%
  summarise(
    mean_cq = mean(cq_slope.x, na.rm = TRUE),
    mean_storage = if (!is.na(storage_col_mega)) mean(.data[[storage_col_mega]], na.rm = TRUE) else NA_real_,
    mean_RBI = mean(RBI, na.rm = TRUE),
    pct_sync = mean(cq_sync == "sync", na.rm = TRUE) * 100,
    .groups = "drop"
  )

# Join with cluster info
annual_cluster <- annual_summary %>%
  left_join(
    clusters_wy %>% select(Stream_Name, chemical, water_year, Cluster_climRef),
    by = c("Stream_Name", "solute" = "chemical", "water_year")
  ) %>%
  rename(Cluster = Cluster_climRef)

# Join with season timing if available
if (!is.null(season_bounds) && "wet_start_doy" %in% names(season_bounds)) {
  annual_cluster <- annual_cluster %>%
    left_join(season_bounds %>% select(water_year, wet_start_doy, wet_duration), by = "water_year")
}

# Compute annual correlations
annual_vars <- c("mean_cq", "mean_storage", "mean_RBI", "pct_sync")
if ("wet_start_doy" %in% names(annual_cluster)) {
  annual_vars <- c(annual_vars, "wet_start_doy", "wet_duration")
}

annual_corr_data <- annual_cluster %>%
  select(any_of(annual_vars)) %>%
  drop_na()

if (nrow(annual_corr_data) > 50) {
  annual_corr <- cor(annual_corr_data, use = "pairwise.complete.obs")
  
  png(file.path(plot_dir, "03_annual_correlations.png"), width = 10, height = 10, units = "in", res = 300)
  corrplot(annual_corr, method = "color", type = "upper",
           addCoef.col = "black", number.cex = 0.9,
           tl.col = "black", tl.srt = 45,
           title = paste0("Annual Correlations (n = ", nrow(annual_corr_data), ")"), mar = c(0,0,2,0))
  dev.off()
  
  message("  Annual correlations saved\n")
}

# #############################################################################
# KEY RELATIONSHIP SCATTERPLOTS
# #############################################################################
message("\n=== KEY RELATIONSHIP SCATTERPLOTS ===\n")

# Storage vs CQ slope (faceted by site)
if (!is.na(storage_col_mega) && all(c(storage_col_mega, "cq_slope.x", "Stream_Name") %in% names(mega))) {
  storage_cq_data <- mega %>%
    filter(is.finite(.data[[storage_col_mega]]), is.finite(cq_slope.x)) %>%
    apply_factor_orders()
  if (nrow(storage_cq_data) > 0) {
    p1 <- storage_cq_data %>%
      ggplot(aes(x = .data[[storage_col_mega]], y = cq_slope.x)) +
      geom_point(alpha = 0.25, size = 1.2, color = "#1F4E79") +
      geom_smooth(method = "lm", se = TRUE, color = "#B3365B", linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dotted", color = "grey60", linewidth = 0.6) +
      geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed", color = "grey55", linewidth = 0.6) +
      facet_wrap(~ Stream_Name, ncol = 3) +
      labs(
        x = storage_label_mega,
        y = get_label("cq_slope"),
        title = "CQ slope vs storage divergence by site",
        subtitle = "Each point = 90-day window; shaded ribbon is 95% CI for site-level linear fit",
        caption = "Horizontal dashed lines at ±0.1 highlight mobilizing vs diluting thresholds."
      ) +
      theme_hja(14)
    if (all(storage_cq_data[[storage_col_mega]] > 0, na.rm = TRUE)) {
      p1 <- p1 + scale_x_log10()
    }
    ggsave(file.path(plot_dir, "04_storage_vs_cq_slope.png"), p1, width = 18, height = 10, dpi = 300)
  }
}

# RBI vs CQ slope (per site)
if (all(c("RBI", "cq_slope.x", "Stream_Name") %in% names(mega))) {
  rbi_cq_data <- mega %>%
    filter(is.finite(RBI), is.finite(cq_slope.x)) %>%
    apply_factor_orders()
  if (nrow(rbi_cq_data) > 0) {
    p2 <- rbi_cq_data %>%
      ggplot(aes(x = RBI, y = cq_slope.x)) +
      geom_point(alpha = 0.25, size = 1.2, color = "#155E75") +
      geom_smooth(method = "lm", se = TRUE, color = "#B3365B", linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dotted", color = "grey60", linewidth = 0.6) +
      geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed", color = "grey55", linewidth = 0.6) +
      facet_wrap(~ Stream_Name, ncol = 3) +
      labs(
        x = get_label("RBI"),
        y = get_label("cq_slope"),
        title = "CQ slope vs flashiness (RBI) by site",
        subtitle = "Scatter shows 90-day windows; ribbon = 95% CI for site-level fit",
        caption = "Dashed bands mark ±0.1 thresholds for CQ slope."
      ) +
      theme_hja(14)
    ggsave(file.path(plot_dir, "04_rbi_vs_cq_slope.png"), p2, width = 13, height = 10, dpi = 300)
  }
}

# Storage vs synchrony (annual aggregates)
if (all(c("mean_storage", "pct_sync", "Stream_Name") %in% names(annual_cluster))) {
  storage_sync_data <- annual_cluster %>%
    mutate(
      sync_fraction = pct_sync / 100,
      Cluster = factor(as.character(Cluster), levels = cluster_levels)
    ) %>%
    filter(is.finite(mean_storage), is.finite(sync_fraction)) %>%
    apply_factor_orders()
  if (nrow(storage_sync_data) > 0) {
    p3 <- storage_sync_data %>%
      ggplot(aes(x = mean_storage, y = sync_fraction, color = Cluster)) +
      geom_point(size = 2.2, alpha = 0.8) +
      geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "#4B5563", linewidth = 0.7) +
      facet_wrap(~ Stream_Name, ncol = 3) +
      scale_color_cluster() +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
      labs(
        x = storage_label_mega,
        y = "% Windows synchronous (Q1 + Q3)",
        color = "Cluster",
        title = "Annual storage vs synchrony",
        subtitle = "Points represent site-solute water years; color reflects cluster assignment",
        caption = "Synchrony calculated as share of windows in Q1 or Q3 for each site-solute-year."
      ) +
      theme_hja(14)
    if (all(storage_sync_data$mean_storage > 0, na.rm = TRUE)) {
      p3 <- p3 + scale_x_log10()
    }
    ggsave(file.path(plot_dir, "04_storage_vs_synchrony.png"), p3, width = 13, height = 10, dpi = 300)
  }
}

# Catchment characteristics vs CQ (site level)
if (nrow(site_summary) >= 5) {
  # Elevation vs CQ
  if ("Elevation_mean_m" %in% names(site_summary) && "mean_cq_slope" %in% names(site_summary)) {
    p4 <- ggplot(site_summary, aes(x = Elevation_mean_m, y = mean_cq_slope)) +
      geom_point(size = 4, alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE, color = "#2166AC") +
      geom_text(aes(label = Stream_Name), vjust = -1, size = 3) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(x = "Mean elevation (m)", y = "Mean CQ slope") +
      theme_clean()
    
    ggsave(file.path(plot_dir, "05_elevation_vs_cq.png"), p4, width = 10, height = 7, dpi = 300)
  }
  
  # Slope vs CQ
  if ("Slope_mean" %in% names(site_summary) && "mean_cq_slope" %in% names(site_summary)) {
    p5 <- ggplot(site_summary, aes(x = Slope_mean, y = mean_cq_slope)) +
      geom_point(size = 4, alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE, color = "#2166AC") +
      geom_text(aes(label = Stream_Name), vjust = -1, size = 3) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(x = "Mean basin slope (degrees)", y = "Mean CQ slope") +
      theme_clean()
    
    ggsave(file.path(plot_dir, "05_slope_vs_cq.png"), p5, width = 10, height = 7, dpi = 300)
  }
  
  # Harvest vs stability
  if ("Harvest" %in% names(site_summary) && "mean_stability" %in% names(site_summary)) {
    p6 <- ggplot(site_summary, aes(x = Harvest, y = mean_stability)) +
      geom_point(size = 4, alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE, color = "#2166AC") +
      geom_text(aes(label = Stream_Name), vjust = -1, size = 3) +
      labs(x = "Harvest (%)", y = "Mean cluster stability") +
      theme_clean()
    
    ggsave(file.path(plot_dir, "05_harvest_vs_stability.png"), p6, width = 10, height = 7, dpi = 300)
  }
}

message("  Scatterplots saved\n")

# #############################################################################
# SUMMARY STATISTICS
# #############################################################################
message("\n=== SUMMARY STATISTICS ===\n")

# Export key correlations summary
if (exists("site_corr")) {
  # Find strongest correlations
  strong_corrs <- site_corr_df %>%
    mutate(abs_r = abs(r)) %>%
    filter(abs_r > 0.5) %>%
    arrange(desc(abs_r)) %>%
    head(20)
  
  message("\nTop site-level correlations (|r| > 0.5):\n")
  print(strong_corrs, n = 20)
  
  write_csv(strong_corrs, file.path(output_dir, "03_stats/top_site_correlations.csv"))
}

message("\n╔════════════════════════════════════════════════════════════════╗\n")
message("║  CORRELATION ANALYSIS COMPLETE                                ║\n")
message("╚════════════════════════════════════════════════════════════════╝\n\n")
message("Outputs saved to:\n")
message("  Plots: ", plot_dir, "\n")
message("  Tables: ", file.path(output_dir, "03_stats"), "\n\n")
