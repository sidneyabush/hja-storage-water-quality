# =============================================================================
# STEP 03d: TEMPORAL TRENDS ANALYSIS
# =============================================================================
# Goal: Detect long-term trends in:
#   - CQ behavior (slopes, patterns)
#   - Synchrony metrics
#   - Cluster membership
#   - Hydrologic regime
#
# Methods:
#   - Mann-Kendall trend tests
#   - Sen's slope estimator
#   - Change point detection
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(trend)       # Mann-Kendall, Sen's slope
  library(patchwork)
})
try(source(file.path("/Users/sidneybush/Documents/GitHub/hja-water-quality", "00_helpers", "plot_theme_set.R")), silent = TRUE)

rm(list = ls())

# Source helpers
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
try(source(file.path(repo_dir, "00_helpers", "plot_prefs.R")), silent = TRUE)

# Paths
base_dir   <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality"
output_dir <- file.path(base_dir, "outputs")
plot_dir   <- file.path(base_dir, "exploratory_plots", "03_stats", "3d_trends")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Theme
theme_clean <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(panel.grid = element_blank(), strip.background = element_blank())
}

message("\n=== LOADING DATA ===\n")

# Load data
mega <- read_csv(file.path(output_dir, "HJA_mega_90d_windows_CQ_hydro_static.csv"), show_col_types = FALSE)
clusters_wy <- read_csv(file.path(output_dir, "ClusterStreams_allSolutes_byWaterYear.csv"), show_col_types = FALSE)

# #############################################################################
# ANNUAL AGGREGATION
# #############################################################################
message("\n=== AGGREGATING TO ANNUAL ===\n")

annual <- mega %>%
  group_by(Stream_Name, solute, water_year) %>%
  summarise(
    mean_cq = mean(cq_slope.x, na.rm = TRUE),
    mean_storage = mean(Q_dS_range_mm, na.rm = TRUE),
    mean_RBI = mean(RBI, na.rm = TRUE),
    pct_positive = mean(cq_slope.x > 0, na.rm = TRUE) * 100,
    pct_sync = mean(cq_sync == "sync", na.rm = TRUE) * 100,
    n_windows = n(),
    .groups = "drop"
  ) %>%
  filter(n_windows >= 3)  # Require at least 3 windows per year

# #############################################################################
# MANN-KENDALL TREND TESTS
# #############################################################################
message("\n=== MANN-KENDALL TREND TESTS ===\n")

# Function to safely run MK test
safe_mk <- function(x) {
  x <- na.omit(x)
  if (length(x) < 8) return(tibble(tau = NA, p_value = NA, sens_slope = NA))
  
  tryCatch({
    mk_result <- mk.test(x)
    sens <- sens.slope(x)
    tibble(
      tau = mk_result$estimates["tau"],
      p_value = mk_result$p.value,
      sens_slope = sens$estimates["Sen's slope"]
    )
  }, error = function(e) {
    tibble(tau = NA, p_value = NA, sens_slope = NA)
  })
}

# Test for each site-solute combination
trend_results <- annual %>%
  group_by(Stream_Name, solute) %>%
  filter(n() >= 8) %>%  # Need enough years
  group_modify(~{
    # CQ slope trend
    cq_trend <- safe_mk(.x$mean_cq)
    names(cq_trend) <- paste0("cq_", names(cq_trend))
    
    # Storage trend
    storage_trend <- safe_mk(.x$mean_storage)
    names(storage_trend) <- paste0("storage_", names(storage_trend))
    
    # Synchrony trend
    sync_trend <- safe_mk(.x$pct_sync)
    names(sync_trend) <- paste0("sync_", names(sync_trend))
    
    bind_cols(cq_trend, storage_trend, sync_trend)
  }) %>%
  ungroup()

# Add significance flags
trend_results <- trend_results %>%
  mutate(
    cq_sig = cq_p_value < 0.05,
    storage_sig = storage_p_value < 0.05,
    sync_sig = sync_p_value < 0.05,
    cq_direction = case_when(
      !cq_sig ~ "NS",
      cq_tau > 0 ~ "Increasing",
      cq_tau < 0 ~ "Decreasing"
    ),
    storage_direction = case_when(
      !storage_sig ~ "NS",
      storage_tau > 0 ~ "Increasing",
      storage_tau < 0 ~ "Decreasing"
    ),
    sync_direction = case_when(
      !sync_sig ~ "NS",
      sync_tau > 0 ~ "Increasing",
      sync_tau < 0 ~ "Decreasing"
    )
  )

# Save results
write_csv(trend_results, file.path(output_dir, "03_stats/mann_kendall_trends.csv"))

message("\nTrend summary:\n")
message("  CQ slope - Increasing:", sum(trend_results$cq_direction == "Increasing", na.rm = TRUE), "\n")
message("  CQ slope - Decreasing:", sum(trend_results$cq_direction == "Decreasing", na.rm = TRUE), "\n")
message("  CQ slope - No trend:", sum(trend_results$cq_direction == "NS", na.rm = TRUE), "\n")

# #############################################################################
# TREND VISUALIZATION
# #############################################################################
message("\n=== TREND VISUALIZATIONS ===\n")

# Summary heatmap of trends
trend_heatmap <- trend_results %>%
  select(Stream_Name, solute, cq_direction, storage_direction, sync_direction) %>%
  pivot_longer(cols = c(cq_direction, storage_direction, sync_direction),
               names_to = "metric", values_to = "trend") %>%
  mutate(metric = str_remove(metric, "_direction"))

p_heatmap <- ggplot(trend_heatmap, aes(x = metric, y = interaction(Stream_Name, solute), fill = trend)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("Increasing" = "#2166AC", "Decreasing" = "#B2182B", "NS" = "grey80"),
                    name = "Trend") +
  labs(x = "Metric", y = "Site-Solute", title = "Mann-Kendall Trend Summary") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 6))

ggsave(file.path(plot_dir, "01_trend_summary_heatmap.png"), p_heatmap, width = 12, height = 16, dpi = 200)

# Time series with trend lines for significant results
sig_cq <- trend_results %>% filter(cq_sig)

if (nrow(sig_cq) > 0) {
  sig_data <- annual %>%
    inner_join(sig_cq %>% select(Stream_Name, solute, cq_sens_slope, cq_tau), 
               by = c("Stream_Name", "solute"))
  
  p_sig <- ggplot(sig_data, aes(x = water_year, y = mean_cq)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    facet_wrap(~paste(Stream_Name, solute, sep = " - "), scales = "free_y") +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    labs(x = "Water Year", y = "Mean CQ Slope", 
         title = "Significant CQ Slope Trends (p < 0.05)") +
    theme_clean()
  
  ggsave(file.path(plot_dir, "02_significant_cq_trends.png"), p_sig, 
         width = 14, height = 10, dpi = 200)
}

# #############################################################################
# AGGREGATE TRENDS (ACROSS ALL SITES)
# #############################################################################
message("\n=== AGGREGATE TRENDS ===\n")

overall_annual <- mega %>%
  group_by(water_year) %>%
  summarise(
    mean_cq = mean(cq_slope.x, na.rm = TRUE),
    median_cq = median(cq_slope.x, na.rm = TRUE),
    mean_storage = mean(Q_dS_range_mm, na.rm = TRUE),
    pct_sync = mean(cq_sync == "sync", na.rm = TRUE) * 100,
    pct_positive = mean(cq_slope.x > 0, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  filter(!is.na(water_year))

# Test overall trends
message("\nOverall CQ trend:\n")
if (nrow(overall_annual) >= 8) {
  overall_cq_mk <- safe_mk(overall_annual$mean_cq)
  message("  tau =", round(overall_cq_mk$tau, 3), ", p =", round(overall_cq_mk$p_value, 4), "\n")
  
  overall_sync_mk <- safe_mk(overall_annual$pct_sync)
  message("Overall synchrony trend:\n")
  message("  tau =", round(overall_sync_mk$tau, 3), ", p =", round(overall_sync_mk$p_value, 4), "\n")
}

# Overall trend plot
p_overall <- ggplot(overall_annual, aes(x = water_year)) +
  geom_line(aes(y = mean_cq, color = "Mean CQ Slope")) +
  geom_point(aes(y = mean_cq, color = "Mean CQ Slope")) +
  geom_smooth(aes(y = mean_cq), method = "lm", se = TRUE, alpha = 0.2, color = "grey50") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Mean CQ Slope" = "#2166AC")) +
  labs(x = "Water Year", y = "Mean CQ Slope", 
       title = "Overall CQ Slope Trend (All Sites Combined)") +
  theme_clean() +
  theme(legend.position = "none")

p_sync <- ggplot(overall_annual, aes(x = water_year)) +
  geom_line(aes(y = pct_sync, color = "% Synchronous")) +
  geom_point(aes(y = pct_sync, color = "% Synchronous")) +
  geom_smooth(aes(y = pct_sync), method = "lm", se = TRUE, alpha = 0.2, color = "grey50") +
  scale_color_manual(values = c("% Synchronous" = "#B2182B")) +
  labs(x = "Water Year", y = "% Synchronous", 
       title = "Overall Synchrony Trend (All Sites Combined)") +
  theme_clean() +
  theme(legend.position = "none")

p_combined <- p_overall / p_sync
ggsave(file.path(plot_dir, "03_overall_trends.png"), p_combined, width = 12, height = 10, dpi = 200)

# #############################################################################
# CLUSTER MEMBERSHIP TRENDS
# #############################################################################
message("\n=== CLUSTER MEMBERSHIP TRENDS ===\n")

cluster_annual <- clusters_wy %>%
  group_by(water_year, Cluster_climRef) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(water_year) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

# Test for trends in cluster proportions
cluster_trends <- cluster_annual %>%
  group_by(Cluster_climRef) %>%
  filter(n() >= 8) %>%
  group_modify(~{
    safe_mk(.x$pct)
  }) %>%
  ungroup()

message("\nCluster proportion trends:\n")
print(cluster_trends)

# Cluster trend plot
p_cluster <- ggplot(cluster_annual, aes(x = water_year, y = pct, color = as.factor(Cluster_climRef))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_brewer(palette = "Set1", name = "Cluster") +
  labs(x = "Water Year", y = "% of Site-Solutes", 
       title = "Cluster Membership Over Time") +
  theme_clean()

ggsave(file.path(plot_dir, "04_cluster_trends.png"), p_cluster, width = 12, height = 7, dpi = 200)

# #############################################################################
# SUMMARY OUTPUT
# #############################################################################
message("\n╔════════════════════════════════════════════════════════════════╗\n")
message("║  TEMPORAL TRENDS ANALYSIS COMPLETE                            ║\n")
message("╚════════════════════════════════════════════════════════════════╝\n\n")

message("Key findings:\n")
message("  - ", sum(trend_results$cq_sig, na.rm = TRUE), "site-solutes with significant CQ trends\n")
message("  - ", sum(trend_results$storage_sig, na.rm = TRUE), "site-solutes with significant storage trends\n")
message("  - ", sum(trend_results$sync_sig, na.rm = TRUE), "site-solutes with significant synchrony trends\n\n")

message("Outputs saved to:\n")
message("  Plots: ", plot_dir, "\n")
message("  Data:  ", file.path(output_dir, "03_stats/mann_kendall_trends.csv"), "\n\n")
