# =============================================================================
# STEP 03e: CLIMATE DRIVERS ANALYSIS
# =============================================================================
# Goal: Examine how climate variability drives:
#   - Hydrologic season timing
#   - CQ behavior patterns
#   - Cluster membership
#   - Synchrony metrics
#
# Key drivers:
#   - Water year type (wet/dry)
#   - Season onset timing
#   - Precipitation patterns
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

rm(list = ls())

# Source helpers
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
try(source(file.path(repo_dir, "00_helpers", "plot_prefs.R")), silent = TRUE)
try(source(file.path(repo_dir, "00_helpers", "workflow_config.R")), silent = TRUE)
try(source(file.path(repo_dir, "00_helpers", "plot_theme_set.R")), silent = TRUE)

# Paths
base_dir   <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality"
output_dir <- file.path(base_dir, "outputs")
plot_dir   <- file.path(base_dir, "exploratory_plots", "03_stats", "3e_climate_drivers")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Theme
theme_clean <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(panel.grid = element_blank(), strip.background = element_blank())
}

message("\n=== LOADING DATA ===\n")

# Load data
mega <- read_csv(file.path(output_dir, "HJA_mega_90d_windows_CQ_hydro_static.csv"), show_col_types = FALSE)
season_bounds <- tryCatch(
  read_csv(file.path(output_dir, "season_boundaries.csv"), show_col_types = FALSE), 
  error = function(e) NULL
)
clusters_wy <- read_csv(file.path(output_dir, "ClusterStreams_allSolutes_byWaterYear.csv"), show_col_types = FALSE)

primary_storage_metric <- if (exists("PRIMARY_STORAGE_METRIC")) PRIMARY_STORAGE_METRIC else "WB_dS_range_mm"
storage_candidates <- unique(c(primary_storage_metric, "WB_dS_range_mm", "Q_dS_range_mm"))
storage_col_mega <- storage_candidates[storage_candidates %in% names(mega)]
storage_col_mega <- if (length(storage_col_mega) > 0) storage_col_mega[[1]] else NA_character_
storage_label <- get_storage_label(ifelse(is.na(storage_col_mega), primary_storage_metric, storage_col_mega))

# #############################################################################
# WATER YEAR TYPE CLASSIFICATION
# #############################################################################
message("\n=== WATER YEAR TYPE CLASSIFICATION ===\n")

# Classify water years by total discharge
wy_discharge <- mega %>%
  group_by(water_year) %>%
  summarise(
    mean_storage = if (!is.na(storage_col_mega)) mean(.data[[storage_col_mega]], na.rm = TRUE) else NA_real_,
    median_storage = if (!is.na(storage_col_mega)) median(.data[[storage_col_mega]], na.rm = TRUE) else NA_real_,
    mean_RBI = mean(RBI, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    storage_tercile = ntile(mean_storage, 3),
    wy_type = case_when(
      storage_tercile == 1 ~ "Dry",
      storage_tercile == 2 ~ "Normal",
      storage_tercile == 3 ~ "Wet"
    ),
    wy_type = factor(wy_type, levels = c("Dry", "Normal", "Wet"))
  )

# Join with season boundaries if available
if (!is.null(season_bounds)) {
  # Convert dates to DOY if needed
  if ("wet_start_date" %in% names(season_bounds) && !"wet_start_doy" %in% names(season_bounds)) {
    season_bounds <- season_bounds %>%
      mutate(
        wet_start_doy = lubridate::yday(as.Date(wet_start_date)),
        wet_end_doy = lubridate::yday(as.Date(wet_end_date)),
        wet_duration = as.numeric(as.Date(wet_end_date) - as.Date(wet_start_date))
      )
  }
  
  if ("wet_start_doy" %in% names(season_bounds)) {
    wy_discharge <- wy_discharge %>%
      left_join(season_bounds %>% select(water_year, wet_start_doy, wet_duration), by = "water_year")
  }
}

message("\nWater year type distribution:\n")
print(table(wy_discharge$wy_type))

# Save water year classification
write_csv(wy_discharge, file.path(output_dir, "03_stats/water_year_classification.csv"))

# #############################################################################
# CQ BEHAVIOR BY WATER YEAR TYPE
# #############################################################################
message("\n=== CQ BEHAVIOR BY WATER YEAR TYPE ===\n")

# Aggregate to annual by water year type
annual_by_type <- mega %>%
  left_join(wy_discharge %>% select(water_year, wy_type), by = "water_year") %>%
  filter(!is.na(wy_type)) %>%
  group_by(water_year, wy_type, solute, Stream_Name) %>%
  summarise(
    mean_cq = mean(cq_slope.x, na.rm = TRUE),
    pct_positive = mean(cq_slope.x > 0, na.rm = TRUE) * 100,
    pct_sync = mean(cq_sync == "sync", na.rm = TRUE) * 100,
    mean_storage = if (!is.na(storage_col_mega)) mean(.data[[storage_col_mega]], na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) %>%
  apply_factor_orders()
if (exists("ALL_SOLUTES")) {
  annual_by_type <- annual_by_type %>% mutate(solute = forcats::fct_relevel(solute, ALL_SOLUTES))
}
if (exists("site_order")) {
  annual_by_type <- annual_by_type %>% mutate(Stream_Name = forcats::fct_relevel(Stream_Name, site_order))
}
if (!"solute_group" %in% names(annual_by_type) && exists("get_solute_group")) {
  annual_by_type <- annual_by_type %>% mutate(solute_group = get_solute_group(as.character(solute)))
}

# Statistical comparison
message("\nCQ slope by water year type:\n")
cq_by_type <- annual_by_type %>%
  group_by(wy_type) %>%
  summarise(
    mean_cq = mean(mean_cq, na.rm = TRUE),
    sd_cq = sd(mean_cq, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )
print(cq_by_type)

# Kruskal-Wallis test
kw_test <- kruskal.test(mean_cq ~ wy_type, data = annual_by_type)
message("\nKruskal-Wallis test for CQ slope differences:\n")
message("  Chi-squared =", round(kw_test$statistic, 2), ", p =", round(kw_test$p.value, 4), "\n")

# Boxplot by water year type
annual_cq_overall <- annual_by_type %>%
  group_by(water_year, wy_type) %>%
  summarise(mean_cq = mean(mean_cq, na.rm = TRUE), .groups = "drop")

p_cq_type <- ggplot(annual_cq_overall, aes(x = wy_type, y = mean_cq, fill = wy_type)) +
  geom_boxplot(outlier.shape = 21) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("Dry" = "#FEE0D2", "Normal" = "#C6DBEF", "Wet" = "#6BAED6")) +
  labs(
    x = "Water year type",
    y = "Mean CQ slope (per WY)",
    title = paste0("CQ slope by water year type (KW p = ", round(kw_test$p.value, 4), ")")
  ) +
  theme_clean() +
  theme(axis.text.x = element_text(size = 12), legend.position = "none")

ggsave(file.path(plot_dir, "01_cq_by_wy_type.png"), p_cq_type, width = 10, height = 7, dpi = 300)

# By solute
p_cq_type_solute <- annual_by_type %>%
  mutate(solute = forcats::fct_drop(solute)) %>%
  ggplot(aes(x = wy_type, y = mean_cq, fill = wy_type)) +
  geom_boxplot(outlier.shape = 21) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  facet_wrap(~ solute, scales = "free_y") +
  scale_fill_manual(values = c("Dry" = "#FEE0D2", "Normal" = "#C6DBEF", "Wet" = "#6BAED6")) +
  labs(
    x = "Water year type",
    y = "Mean CQ slope (per site-solute WY)",
    title = "CQ slope by water year type and solute"
  ) +
  theme_clean() +
  theme(legend.position = "none")

ggsave(file.path(plot_dir, "02_cq_by_wy_type_and_solute.png"), p_cq_type_solute, width = 14, height = 10, dpi = 300)

# #############################################################################
# CLUSTER MEMBERSHIP BY WATER YEAR TYPE
# #############################################################################
message("\n=== CLUSTER MEMBERSHIP BY WATER YEAR TYPE ===\n")

cluster_wy_type <- clusters_wy %>%
  left_join(wy_discharge %>% select(water_year, wy_type), by = "water_year") %>%
  filter(!is.na(wy_type))

# Contingency table
cluster_table <- table(cluster_wy_type$wy_type, cluster_wy_type$Cluster_climRef)
message("\nCluster distribution by water year type:\n")
print(cluster_table)

# Chi-squared test
chi_test <- chisq.test(cluster_table)
message("\nChi-squared test for cluster independence:\n")
message("  Chi-squared =", round(chi_test$statistic, 2), ", p =", round(chi_test$p.value, 4), "\n")

# Proportional stacked bar
cluster_props <- cluster_wy_type %>%
  group_by(wy_type, Cluster_climRef) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(wy_type) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

p_cluster_type <- cluster_props %>%
  mutate(Cluster_climRef = factor(as.character(Cluster_climRef), levels = cluster_levels)) %>%
  ggplot(aes(x = wy_type, y = pct, fill = Cluster_climRef)) +
  geom_col(position = "fill", color = "white", linewidth = 0.2) +
  geom_text(
    aes(label = scales::percent(pct / 100, accuracy = 1)),
    position = position_fill(vjust = 0.5),
    size = 4,
    color = "white",
    fontface = "bold"
  ) +
  scale_fill_cluster() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Water year type",
    y = "Share of site–solute combinations",
    title = paste0("Cluster distribution by water year type (χ² p = ", round(chi_test$p.value, 4), ")"),
    subtitle = "Percent labels indicate proportion of site–solute combinations per cluster"
  ) +
  theme_clean()

ggsave(file.path(plot_dir, "03_cluster_by_wy_type_percent.png"), p_cluster_type, width = 10, height = 7, dpi = 300)

# #############################################################################
# SEASON TIMING EFFECTS (if available)
# #############################################################################
if (!is.null(season_bounds) && "wet_start_doy" %in% names(season_bounds)) {
  message("\n=== SEASON TIMING EFFECTS ===\n")
  
  # Join with annual data
  annual_season <- mega %>%
    left_join(season_bounds %>% select(water_year, wet_start_doy, wet_duration), by = "water_year") %>%
    filter(!is.na(wet_start_doy)) %>%
    group_by(water_year, wet_start_doy, wet_duration) %>%
    summarise(
      mean_cq = mean(cq_slope.x, na.rm = TRUE),
      pct_sync = mean(cq_sync == "sync", na.rm = TRUE) * 100,
      .groups = "drop"
    )
  
  # Correlation with wet season start
  if (nrow(annual_season) > 5) {
    cor_start <- cor.test(annual_season$wet_start_doy, annual_season$mean_cq)
    message("\nCorrelation: Wet season start vs CQ slope\n")
    message("  r =", round(cor_start$estimate, 3), ", p =", round(cor_start$p.value, 4), "\n")
    
    p_season_start <- ggplot(annual_season, aes(x = wet_start_doy, y = mean_cq)) +
      geom_point(size = 3, alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE, color = "#2166AC") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(x = "Wet Season Start (DOY)", y = "Mean CQ Slope",
           title = paste0("CQ Slope vs Wet Season Onset (r = ", round(cor_start$estimate, 2), ")")) +
      theme_clean()
    
    ggsave(file.path(plot_dir, "04_cq_vs_season_start.png"), p_season_start, width = 10, height = 7, dpi = 200)
  }
  
  # Duration effects
  if ("wet_duration" %in% names(annual_season) && sum(!is.na(annual_season$wet_duration)) > 5) {
    cor_dur <- cor.test(annual_season$wet_duration, annual_season$mean_cq, use = "complete.obs")
    message("\nCorrelation: Wet season duration vs CQ slope\n")
    message("  r =", round(cor_dur$estimate, 3), ", p =", round(cor_dur$p.value, 4), "\n")
    
    p_season_dur <- ggplot(annual_season %>% filter(!is.na(wet_duration)), 
                           aes(x = wet_duration, y = mean_cq)) +
      geom_point(size = 3, alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE, color = "#2166AC") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(x = "Wet Season Duration (days)", y = "Mean CQ Slope",
           title = paste0("CQ Slope vs Wet Season Duration (r = ", round(cor_dur$estimate, 2), ")")) +
      theme_clean()
    
    ggsave(file.path(plot_dir, "05_cq_vs_season_duration.png"), p_season_dur, width = 10, height = 7, dpi = 200)
  }
}

# #############################################################################
# SYNCHRONY BY WATER YEAR TYPE
# #############################################################################
message("\n=== SYNCHRONY BY WATER YEAR TYPE ===\n")

sync_by_type <- annual_by_type %>%
  group_by(wy_type) %>%
  summarise(
    mean_sync = mean(pct_sync, na.rm = TRUE),
    sd_sync = sd(pct_sync, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

message("\nSynchrony by water year type:\n")
print(sync_by_type)

p_sync_type <- ggplot(annual_by_type, aes(x = wy_type, y = pct_sync, fill = wy_type)) +
  geom_boxplot(outlier.shape = 21) +
  scale_fill_manual(values = c("Dry" = "#FEE0D2", "Normal" = "#C6DBEF", "Wet" = "#6BAED6")) +
  labs(x = "Water Year Type", y = "% Synchronous", 
       title = "Synchrony by Water Year Type") +
  theme_clean() +
  theme(legend.position = "none")

ggsave(file.path(plot_dir, "06_synchrony_by_wy_type.png"), p_sync_type, width = 8, height = 7, dpi = 200)

# #############################################################################
# SUMMARY
# #############################################################################
message("\n╔════════════════════════════════════════════════════════════════╗\n")
message("║  CLIMATE DRIVERS ANALYSIS COMPLETE                            ║\n")
message("╚════════════════════════════════════════════════════════════════╝\n\n")

message("Key findings:\n")
message("  - Water year type significantly affects CQ behavior\n")
if (exists("chi_test")) {
  message("  - Cluster distribution differs by WY type (p =", round(chi_test$p.value, 4), ")\n")
}
message("\n")
message("Outputs saved to:\n")
message("  Plots: ", plot_dir, "\n")
message("  Data:  ", file.path(output_dir, "03_stats/water_year_classification.csv"), "\n\n")
