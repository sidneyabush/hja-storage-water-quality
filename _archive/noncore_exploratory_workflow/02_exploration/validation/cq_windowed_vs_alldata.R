#!/usr/bin/env Rscript
# =============================================================================
# Compare Windowed vs All-Data CQ Slopes
# =============================================================================
# Purpose: Compare CQ slopes calculated from:
#   1. Rolling windows (75/150 days) averaged to site level
#   2. All data (no windowing) at site level
#
# This helps determine if the windowing approach introduces bias or if
# both methods give similar site-level CQ relationships.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

rm(list = ls())

# Paths
repo_dir      <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
base_dir      <- "/Users/sidneybush/Library/CloudStorage/Box-Box"
project_dir   <- file.path(base_dir, "Sidney_Bush", "HJA_Water_Quality")
out_dir       <- file.path(project_dir, "outputs")
fig_dir       <- file.path(project_dir, "exploratory_plots", "02_exploration", "validation")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Source plot preferences for consistent colors
source(file.path(repo_dir, "00_helpers", "plot_prefs.R"))

message("=== COMPARING WINDOWED VS ALL-DATA CQ SLOPES ===\n")

# =============================================================================
# LOAD DATA
# =============================================================================

# All-data CQ slopes (just calculated)
cq_all <- readr::read_csv(file.path(out_dir, "HJA_CQ_all_data_site_means.csv"),
                          show_col_types = FALSE)

# Windowed CQ slopes averaged to site level
# Need to load the clean site means or calculate from master
site_means <- readr::read_csv(file.path(out_dir, "HJA_clean_site_means.csv"),
                              show_col_types = FALSE)

# If cq_slope not in site_means, try master
if (!"cq_slope" %in% names(site_means)) {
  message("  CQ slopes not in clean site means, calculating from master...")

  master <- readr::read_csv(file.path(out_dir, "HJA_master_site_means.csv"),
                            show_col_types = FALSE)

  if ("cq_slope" %in% names(master)) {
    site_means <- master %>%
      select(Stream_Name, solute, cq_slope, cq_CVc_CVq)
  } else {
    stop("cq_slope not found in site means files. Run data prep first.")
  }
}

# Clean windowed data
cq_windowed <- site_means %>%
  select(Stream_Name, solute,
         cq_slope_windowed = cq_slope,
         cq_CVc_CVq_windowed = cq_CVc_CVq) %>%
  filter(!is.na(cq_slope_windowed))

message("  Loaded windowed CQ slopes: ", nrow(cq_windowed), " site-solute combinations")
message("  Loaded all-data CQ slopes: ", nrow(cq_all), " site-solute combinations")

# =============================================================================
# MERGE AND COMPARE
# =============================================================================

comparison <- cq_all %>%
  select(Stream_Name, solute,
         cq_slope_all, cq_CVc_CVq_all, cq_behavior_all,
         n_obs, n_years, cq_r2_all) %>%
  inner_join(cq_windowed, by = c("Stream_Name", "solute"))

message("\n  Matched ", nrow(comparison), " site-solute combinations for comparison")

# Calculate differences and convert solute to factor with proper ordering
comparison <- comparison %>%
  mutate(
    slope_diff = cq_slope_windowed - cq_slope_all,
    slope_pct_diff = 100 * (cq_slope_windowed - cq_slope_all) / abs(cq_slope_all),
    CVratio_diff = cq_CVc_CVq_windowed - cq_CVc_CVq_all,
    # Use plot_prefs solute ordering
    solute = factor(solute, levels = names(solute_colors))
  )

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

message("\n=== SUMMARY STATISTICS ===")
message("\nCQ Slope Comparison:")
message("  Correlation (r): ", round(cor(comparison$cq_slope_windowed,
                                          comparison$cq_slope_all,
                                          use = "complete.obs"), 3))
message("  Mean difference (windowed - all): ",
        round(mean(comparison$slope_diff, na.rm = TRUE), 4))
message("  Median difference: ",
        round(median(comparison$slope_diff, na.rm = TRUE), 4))
message("  Mean absolute difference: ",
        round(mean(abs(comparison$slope_diff), na.rm = TRUE), 4))
message("  Mean % difference: ",
        round(mean(abs(comparison$slope_pct_diff), na.rm = TRUE), 1), "%")

message("\nCV Ratio Comparison:")
message("  Correlation (r): ", round(cor(comparison$cq_CVc_CVq_windowed,
                                          comparison$cq_CVc_CVq_all,
                                          use = "complete.obs"), 3))
message("  Mean difference: ",
        round(mean(comparison$CVratio_diff, na.rm = TRUE), 4))

# Agreement in behavior classification
behavior_comparison <- comparison %>%
  mutate(
    behavior_windowed = case_when(
      abs(cq_slope_windowed) < 0.05 ~ "chemostatic",
      cq_slope_windowed > 0.05 ~ "mobilization",
      cq_slope_windowed < -0.05 ~ "dilution"
    ),
    behavior_match = behavior_windowed == cq_behavior_all
  )

message("\nBehavior Classification Agreement:")
message("  Matching classifications: ",
        sum(behavior_comparison$behavior_match, na.rm = TRUE), "/",
        sum(!is.na(behavior_comparison$behavior_match)),
        " (", round(100 * mean(behavior_comparison$behavior_match, na.rm = TRUE), 1), "%)")

# =============================================================================
# VISUALIZATIONS
# =============================================================================

message("\n=== CREATING COMPARISON PLOTS ===")

# Plot 1: CQ slope windowed vs all-data
p1 <- ggplot(comparison, aes(x = cq_slope_all, y = cq_slope_windowed)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(color = solute), alpha = 0.7, size = 2) +
  scale_color_solute() +
  labs(
                                        comparison$cq_slope_all,
                                        use = "complete.obs"), 3)),
    x = "CQ Slope (All Data)",
    y = "CQ Slope (Windowed Average)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

# Plot 2: Difference vs all-data slope
p2 <- ggplot(comparison, aes(x = cq_slope_all, y = slope_diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(color = solute), alpha = 0.7, size = 2) +
  scale_color_solute() +
  labs(
                     round(mean(comparison$slope_diff, na.rm = TRUE), 4)),
    x = "CQ Slope (All Data)",
    y = "Difference (Windowed - All Data)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

# Plot 3: CV ratio comparison
p3 <- ggplot(comparison, aes(x = cq_CVc_CVq_all, y = cq_CVc_CVq_windowed)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(color = solute), alpha = 0.7, size = 2) +
  scale_color_solute() +
  labs(
                                        comparison$cq_CVc_CVq_all,
                                        use = "complete.obs"), 3)),
    x = "CV(C)/CV(Q) (All Data)",
    y = "CV(C)/CV(Q) (Windowed Average)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

# Plot 4: Distribution of differences
p4 <- ggplot(comparison, aes(x = slope_diff)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = median(comparison$slope_diff, na.rm = TRUE),
             linetype = "solid", color = "darkblue") +
  labs(
    x = "Difference (Windowed - All Data)",
    y = "Count"
  ) +
  theme_minimal(base_size = 12)

# Combine plots
combined <- (p1 | p2) / (p3 | p4)

ggsave(file.path(fig_dir, "CQ_windowed_vs_alldata_comparison.png"),
       combined, width = 14, height = 10, dpi = 300, bg = "white")

message("  Saved: CQ_windowed_vs_alldata_comparison.png")

# =============================================================================
# SAVE COMPARISON TABLE
# =============================================================================

readr::write_csv(comparison, file.path(out_dir, "CQ_windowed_vs_alldata_comparison.csv"))

message("\n  Saved: CQ_windowed_vs_alldata_comparison.csv")

# =============================================================================
# RECOMMENDATIONS
# =============================================================================

message("\n=== RECOMMENDATIONS ===")

corr <- cor(comparison$cq_slope_windowed, comparison$cq_slope_all, use = "complete.obs")
mean_abs_diff <- mean(abs(comparison$slope_diff), na.rm = TRUE)
agreement <- mean(behavior_comparison$behavior_match, na.rm = TRUE)

if (corr > 0.95 && mean_abs_diff < 0.1 && agreement > 0.9) {
  message("\n✓ HIGH AGREEMENT: Windowed and all-data approaches give very similar results.")
  message("  → Recommendation: Stick with windowed approach for consistency.")
  message("    Windowing is essential for temporal synchrony analysis.")
} else if (corr > 0.85 && mean_abs_diff < 0.2) {
  message("\n⚠ MODERATE AGREEMENT: Some differences exist but overall correlation is strong.")
  message("  → Recommendation: Use windowed data for temporal analyses (synchrony).")
  message("    Consider using all-data CQ slopes for site characterization.")
} else {
  message("\n✗ LOW AGREEMENT: Substantial differences between approaches.")
  message("  → Recommendation: Investigate why differences are large.")
  message("    May need to refine window size or minimum observation thresholds.")
}

message("\n=== COMPARISON COMPLETE ===")
