#!/usr/bin/env Rscript
# =============================================================================
# Compare Univariate vs Multivariate Synchrony
# =============================================================================
# Purpose: Compare synchrony metrics calculated using:
#   1. UNIVARIATE: Each solute analyzed separately, then averaged to site level
#   2. MULTIVARIATE: All solutes analyzed together as a chemical profile
#
# RESEARCH QUESTION:
#   Does multivariate synchrony add information beyond averaging univariate metrics?
#   Or do they give similar site-level characterizations?
#
# METHODS:
#   - Abbott: Concentration synchrony (univariate) vs chemical profile synchrony (multivariate)
#   - Wymore: Individual solute CQ synchrony (univariate) vs multi-solute CQ synchrony (multivariate)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
})

rm(list = ls())

# Paths
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
base_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box"
project_dir <- file.path(base_dir, "Sidney_Bush", "HJA_Water_Quality")
out_dir <- file.path(project_dir, "outputs")
fig_dir <- file.path(project_dir, "exploratory_plots", "02_exploration", "synchrony_comparisons")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Source plot preferences
source(file.path(repo_dir, "00_helpers", "plot_prefs.R"))

message("=== UNIVARIATE vs MULTIVARIATE SYNCHRONY ===\n")

# =============================================================================
# LOAD DATA
# =============================================================================

# Univariate synchrony (site-solute level, one solute at a time)
sync_uni <- readr::read_csv(file.path(out_dir, "HJA_composite_synchrony.csv"),
                            show_col_types = FALSE)

# Multivariate synchrony (site level, all solutes together)
sync_mv <- readr::read_csv(file.path(out_dir, "HJA_composite_synchrony_multivariate.csv"),
                           show_col_types = FALSE)

message("Loaded univariate: ", nrow(sync_uni), " site-solute combinations")
message("Loaded multivariate: ", nrow(sync_mv), " sites\n")

# =============================================================================
# PREPARE UNIVARIATE DATA (AVERAGE ACROSS SOLUTES)
# =============================================================================

# Average univariate metrics across all solutes to get site-level
uni_site_all <- sync_uni %>%
  group_by(Stream_Name) %>%
  summarise(
    # Abbott concentration synchrony (averaged across solutes)
    conc_sync_uni_all = mean(conc_sync_allpairs, na.rm = TRUE),
    # Wymore CQ synchrony (averaged across solutes)
    wymore_sync_uni_all = mean(wymore_crosssite_allpairs, na.rm = TRUE),
    # Wymore Q1 proportion (averaged across solutes)
    wymore_Q1_proportion = mean(wymore_Q1_proportion, na.rm = TRUE),
    # Wymore Q3 proportion (averaged across solutes)
    wymore_Q3_proportion = mean(wymore_Q3_proportion, na.rm = TRUE),
    n_solutes_all = n(),
    .groups = "drop"
  )

# Average by solute group (using 3-way grouping)
uni_site_biogenic <- sync_uni %>%
  filter(solute %in% BIOGENIC_SOLUTES) %>%
  group_by(Stream_Name) %>%
  summarise(
    conc_sync_uni_bio = mean(conc_sync_allpairs, na.rm = TRUE),
    wymore_sync_uni_bio = mean(wymore_crosssite_allpairs, na.rm = TRUE),
    n_solutes_bio = n(),
    .groups = "drop"
  )

uni_site_geogenic <- sync_uni %>%
  filter(solute %in% GEOGENIC_SOLUTES) %>%
  group_by(Stream_Name) %>%
  summarise(
    conc_sync_uni_geo = mean(conc_sync_allpairs, na.rm = TRUE),
    wymore_sync_uni_geo = mean(wymore_crosssite_allpairs, na.rm = TRUE),
    n_solutes_geo = n(),
    .groups = "drop"
  )

# Combine
uni_site <- uni_site_all %>%
  left_join(uni_site_biogenic, by = "Stream_Name") %>%
  left_join(uni_site_geogenic, by = "Stream_Name")

message("Averaged univariate metrics to ", nrow(uni_site), " sites")

# =============================================================================
# MERGE UNIVARIATE AND MULTIVARIATE
# =============================================================================

comparison <- sync_mv %>%
  left_join(uni_site, by = "Stream_Name") %>%
  mutate(Stream_Name = factor(Stream_Name, levels = site_order))

message("Merged data: ", nrow(comparison), " sites for comparison\n")

# =============================================================================
# ABBOTT COMPARISON: UNIVARIATE VS MULTIVARIATE
# =============================================================================
message("=== ABBOTT CONCENTRATION SYNCHRONY ===")

# All solutes
abbott_all <- comparison %>%
  select(Stream_Name, concentration_mv_all_sync_allpairs, conc_sync_uni_all) %>%
  filter(!is.na(concentration_mv_all_sync_allpairs), !is.na(conc_sync_uni_all))

cor_abbott_all <- cor(abbott_all$concentration_mv_all_sync_allpairs,
                      abbott_all$conc_sync_uni_all, use = "complete.obs")
message("All solutes: r = ", round(cor_abbott_all, 3))

# Biogenic
abbott_bio <- comparison %>%
  select(Stream_Name, concentration_mv_biogenic_sync_allpairs, conc_sync_uni_bio) %>%
  filter(!is.na(concentration_mv_biogenic_sync_allpairs), !is.na(conc_sync_uni_bio))

cor_abbott_bio <- cor(abbott_bio$concentration_mv_biogenic_sync_allpairs,
                      abbott_bio$conc_sync_uni_bio, use = "complete.obs")
message("Biogenic: r = ", round(cor_abbott_bio, 3))

# Geogenic
abbott_geo <- comparison %>%
  select(Stream_Name, concentration_mv_geogenic_sync_allpairs, conc_sync_uni_geo) %>%
  filter(!is.na(concentration_mv_geogenic_sync_allpairs), !is.na(conc_sync_uni_geo))

cor_abbott_geo <- cor(abbott_geo$concentration_mv_geogenic_sync_allpairs,
                      abbott_geo$conc_sync_uni_geo, use = "complete.obs")
message("Geogenic: r = ", round(cor_abbott_geo, 3))

# =============================================================================
# WYMORE COMPARISON: UNIVARIATE VS MULTIVARIATE
# =============================================================================
message("\n=== WYMORE CQ SYNCHRONY ===")

# All solutes
wymore_all <- comparison %>%
  select(Stream_Name, wymore_mv_all_prop_sync_allpairs, wymore_sync_uni_all) %>%
  filter(!is.na(wymore_mv_all_prop_sync_allpairs), !is.na(wymore_sync_uni_all))

cor_wymore_all <- cor(wymore_all$wymore_mv_all_prop_sync_allpairs,
                     wymore_all$wymore_sync_uni_all, use = "complete.obs")
message("All solutes: r = ", round(cor_wymore_all, 3))

# Biogenic
wymore_bio <- comparison %>%
  select(Stream_Name, wymore_mv_biogenic_prop_sync_allpairs, wymore_sync_uni_bio) %>%
  filter(!is.na(wymore_mv_biogenic_prop_sync_allpairs), !is.na(wymore_sync_uni_bio))

cor_wymore_bio <- cor(wymore_bio$wymore_mv_biogenic_prop_sync_allpairs,
                     wymore_bio$wymore_sync_uni_bio, use = "complete.obs")
message("Biogenic: r = ", round(cor_wymore_bio, 3))

# Geogenic
wymore_geo <- comparison %>%
  select(Stream_Name, wymore_mv_geogenic_prop_sync_allpairs, wymore_sync_uni_geo) %>%
  filter(!is.na(wymore_mv_geogenic_prop_sync_allpairs), !is.na(wymore_sync_uni_geo))

cor_wymore_geo <- cor(wymore_geo$wymore_mv_geogenic_prop_sync_allpairs,
                     wymore_geo$wymore_sync_uni_geo, use = "complete.obs")
message("Geogenic: r = ", round(cor_wymore_geo, 3))

# =============================================================================
# PLOTS
# =============================================================================
message("\n=== CREATING PLOTS ===\n")

# Abbott All Solutes
p_abbott_all <- ggplot(abbott_all,
                       aes(x = conc_sync_uni_all, y = concentration_mv_all_sync_allpairs)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(color = Stream_Name), size = 4, alpha = 0.8) +
  geom_text_repel(aes(label = Stream_Name), size = 2.5, max.overlaps = 20, seed = 42) +
  scale_color_site() +
  labs(
    title = "ALL SOLUTES (n=11)",
    subtitle = "Abbott Concentration Synchrony",
    x = "Univariate (averaged across solutes)",
    y = "Multivariate (all solutes together)"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 9))

# Wymore All Solutes
p_wymore_all <- ggplot(wymore_all,
                       aes(x = wymore_sync_uni_all, y = wymore_mv_all_prop_sync_allpairs)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(color = Stream_Name), size = 4, alpha = 0.8) +
  geom_text_repel(aes(label = Stream_Name), size = 2.5, max.overlaps = 20, seed = 42) +
  scale_color_site() +
  labs(
    title = "ALL SOLUTES (n=11)",
    subtitle = "Wymore CQ Behavior Synchrony",
    x = "Univariate (averaged across solutes)",
    y = "Multivariate (proportion of all windows)"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 9))

# Abbott Biogenic
p_abbott_bio <- ggplot(abbott_bio,
                       aes(x = conc_sync_uni_bio, y = concentration_mv_biogenic_sync_allpairs)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(color = Stream_Name), size = 4, alpha = 0.8) +
  geom_text_repel(aes(label = Stream_Name), size = 2.5, max.overlaps = 20, seed = 42) +
  scale_color_site() +
  labs(
    title = "BIOGENIC (n=4)",
    subtitle = "NO3, PO4, NH3, DOC",
    x = "Univariate (averaged)",
    y = "Multivariate"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 8, face = "italic"))

# Wymore Biogenic
p_wymore_bio <- ggplot(wymore_bio,
                       aes(x = wymore_sync_uni_bio, y = wymore_mv_biogenic_prop_sync_allpairs)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(color = Stream_Name), size = 4, alpha = 0.8) +
  geom_text_repel(aes(label = Stream_Name), size = 2.5, max.overlaps = 20, seed = 42) +
  scale_color_site() +
  labs(
    title = "BIOGENIC (n=4)",
    subtitle = "NO3, PO4, NH3, DOC",
    x = "Univariate (averaged)",
    y = "Multivariate (proportion of all windows)"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 8, face = "italic"))

# Abbott Geogenic
p_abbott_geo <- ggplot(abbott_geo,
                       aes(x = conc_sync_uni_geo, y = concentration_mv_geogenic_sync_allpairs)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(color = Stream_Name), size = 4, alpha = 0.8) +
  geom_text_repel(aes(label = Stream_Name), size = 2.5, max.overlaps = 20, seed = 42) +
  scale_color_site() +
  labs(
    title = "GEOGENIC (n=6)",
    subtitle = "Ca, Mg, Na, K, Cl, SO4",
    x = "Univariate (averaged)",
    y = "Multivariate"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 8, face = "italic"))

# Wymore Geogenic
p_wymore_geo <- ggplot(wymore_geo,
                       aes(x = wymore_sync_uni_geo, y = wymore_mv_geogenic_prop_sync_allpairs)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(color = Stream_Name), size = 4, alpha = 0.8) +
  geom_text_repel(aes(label = Stream_Name), size = 2.5, max.overlaps = 20, seed = 42) +
  scale_color_site() +
  labs(
    title = "GEOGENIC (n=6)",
    subtitle = "Ca, Mg, Na, K, Cl, SO4",
    x = "Univariate (averaged)",
    y = "Multivariate (proportion of all windows)"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 8, face = "italic")) +
  guides(color = guide_legend(nrow = 2, title = "Site"))

# Combine plots
combined <- (p_abbott_all | p_wymore_all) /
            (p_abbott_bio | p_wymore_bio) /
            (p_abbott_geo | p_wymore_geo) +
  plot_layout()

ggsave(file.path(fig_dir, "univariate_vs_multivariate_comparison.png"),
       combined, width = 12, height = 14, dpi = 300, bg = "white")

message("✓ Saved: univariate_vs_multivariate_comparison.png")

# =============================================================================
# INTERPRETATION
# =============================================================================
message("\n=== INTERPRETATION ===")

# Threshold for "high agreement"
high_agreement_threshold <- 0.90
moderate_agreement_threshold <- 0.70

interpret_correlation <- function(r, name) {
  if (r > high_agreement_threshold) {
    message(sprintf("  %s: r=%.3f - HIGH AGREEMENT", name, r))
    message("    → Multivariate does NOT add much information beyond averaging univariate")
  } else if (r > moderate_agreement_threshold) {
    message(sprintf("  %s: r=%.3f - MODERATE AGREEMENT", name, r))
    message("    → Multivariate adds SOME unique information")
  } else {
    message(sprintf("  %s: r=%.3f - LOW AGREEMENT", name, r))
    message("    → Multivariate captures DIFFERENT patterns than univariate average")
  }
}

message("\nAbbott (Concentration Synchrony):")
interpret_correlation(cor_abbott_all, "All solutes")
interpret_correlation(cor_abbott_bio, "Biogenic")
interpret_correlation(cor_abbott_geo, "Geogenic")

message("\nWymore (CQ Synchrony):")
interpret_correlation(cor_wymore_all, "All solutes")
interpret_correlation(cor_wymore_bio, "Biogenic")
interpret_correlation(cor_wymore_geo, "Geogenic")

message("\n=== KEY INSIGHTS ===")
message("If correlations are HIGH (r > 0.90):")
message("  → Averaging univariate metrics is sufficient for site characterization")
message("  → Multivariate approach confirms univariate patterns but adds little new info")
message("\nIf correlations are LOW (r < 0.70):")
message("  → Multivariate captures emergent patterns not visible in individual solutes")
message("  → Chemical profile synchrony differs from average solute-by-solute synchrony")
message("  → Multivariate approach reveals site-level chemical coherence")

message("\n=== COMPLETE ===\n")

# =============================================================================
# WYMORE CQ QUADRANT PROPORTIONS: Q1 (Mobilizing) vs Q3 (Diluting)
# =============================================================================
message("\n=== CREATING WYMORE QUADRANT PROPORTION FIGURES ===\n")

# Prepare data for Q1 and Q3 proportions (all solutes)
q1_all <- comparison %>%
  select(Stream_Name, wymore_mv_all_prop_Q1_allpairs, wymore_Q1_proportion) %>%
  filter(!is.na(wymore_mv_all_prop_Q1_allpairs), !is.na(wymore_Q1_proportion))

q3_all <- comparison %>%
  select(Stream_Name, wymore_mv_all_prop_Q3_allpairs, wymore_Q3_proportion) %>%
  filter(!is.na(wymore_mv_all_prop_Q3_allpairs), !is.na(wymore_Q3_proportion))

cor_q1 <- cor(q1_all$wymore_mv_all_prop_Q1_allpairs, q1_all$wymore_Q1_proportion, use = "complete.obs")
cor_q3 <- cor(q3_all$wymore_mv_all_prop_Q3_allpairs, q3_all$wymore_Q3_proportion, use = "complete.obs")

message("Q1 (Mobilizing) correlation: r = ", round(cor_q1, 3))
message("Q3 (Diluting) correlation: r = ", round(cor_q3, 3))

# Plot 1: Overall synchrony (reuse wymore_all from above)
p_sync <- ggplot(wymore_all,
                 aes(x = wymore_sync_uni_all, y = wymore_mv_all_prop_sync_allpairs)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(color = Stream_Name), size = 4, alpha = 0.8) +
  geom_text_repel(aes(label = Stream_Name), size = 2.5, max.overlaps = 20, seed = 42) +
  scale_color_site() +
  labs(
    x = "Univariate (averaged across solutes)",
    y = "Multivariate (proportion synchronized)"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

# Plot 2: Q1 Proportion (Mobilizing)
p_q1 <- ggplot(q1_all,
               aes(x = wymore_Q1_proportion, y = wymore_mv_all_prop_Q1_allpairs)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(color = Stream_Name), size = 4, alpha = 0.8) +
  geom_text_repel(aes(label = Stream_Name), size = 2.5, max.overlaps = 20, seed = 42) +
  scale_color_site() +
  labs(
    x = "Univariate (averaged across solutes)",
    y = "Multivariate (proportion in Q1)"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

# Plot 3: Q3 Proportion (Diluting)
p_q3 <- ggplot(q3_all,
               aes(x = wymore_Q3_proportion, y = wymore_mv_all_prop_Q3_allpairs)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(color = Stream_Name), size = 4, alpha = 0.8) +
  geom_text_repel(aes(label = Stream_Name), size = 2.5, max.overlaps = 20, seed = 42) +
  scale_color_site() +
  labs(
    x = "Univariate (averaged across solutes)",
    y = "Multivariate (proportion in Q3)"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2, title = "Site"))

# Combine three rows
combined_quadrants <- p_sync / p_q1 / p_q3 +
  plot_layout()

ggsave(file.path(fig_dir, "wymore_quadrant_proportions_comparison.png"),
       combined_quadrants, width = 10, height = 13, dpi = 300, bg = "white")

message("✓ Saved: wymore_quadrant_proportions_comparison.png\n")
