#!/usr/bin/env Rscript
# =============================================================================
# Plot Abbott vs Wymore MULTIVARIATE Synchrony
# =============================================================================
# Purpose: Compare multivariate Abbott concentration synchrony with multivariate
#          Wymore CQ slope synchrony across sites
#
# MULTIVARIATE APPROACH:
#   - Abbott MV: All solutes compared together as a chemical profile
#   - Wymore MV: Proportion of solutes synchronized (Q1 or Q3) across site pairs
#
# INTERPRETATION:
#   - Abbott MV sync: Do sites have similar CHEMICAL PROFILES over time?
#   - Wymore MV sync: Do sites mobilize MULTIPLE solutes the same way?
#
# This complements univariate analysis by asking if sites synchronize in
# their overall chemical behavior, not just individual solute patterns.
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

# Source plot preferences for consistent colors
source(file.path(repo_dir, "00_helpers", "plot_prefs.R"))

message("=== MULTIVARIATE: ABBOTT vs WYMORE SYNCHRONY ===\n")

# =============================================================================
# LOAD DATA
# =============================================================================

# Load multivariate synchrony (site-level, averaged across years)
sync_mv <- readr::read_csv(file.path(out_dir, "HJA_composite_synchrony_multivariate.csv"),
                           show_col_types = FALSE)

message("Loaded ", nrow(sync_mv), " site-level multivariate synchrony metrics\n")

# Check what columns exist
message("Available columns:")
message(paste(names(sync_mv), collapse = ", "))

# =============================================================================
# PREPARE DATA FOR PLOTTING
# =============================================================================

# Abbott multivariate synchrony columns: concentration_mv_*_sync_allpairs
# Wymore multivariate synchrony columns: wymore_mv_*_prop_sync_allpairs

# Extract Q1/Q3 proportions to calculate dominant quadrant
wymore_q1 <- sync_mv %>%
  select(Stream_Name, matches("wymore_mv.*prop_Q1_allpairs")) %>%
  pivot_longer(
    cols = -Stream_Name,
    names_to = "metric",
    values_to = "prop_Q1"
  ) %>%
  mutate(
    solute_group = case_when(
      str_detect(metric, "_all_") ~ "All Solutes (n=11)",
      str_detect(metric, "_nutrient_") ~ "Nutrient (n=1)",
      str_detect(metric, "_biogenic_") ~ "Biogenic (n=4)",
      str_detect(metric, "_geogenic_") ~ "Geogenic (n=6)",
      TRUE ~ "Unknown"
    )
  ) %>%
  select(Stream_Name, solute_group, prop_Q1) %>%
  filter(solute_group != "Unknown")

wymore_q3 <- sync_mv %>%
  select(Stream_Name, matches("wymore_mv.*prop_Q3_allpairs")) %>%
  pivot_longer(
    cols = -Stream_Name,
    names_to = "metric",
    values_to = "prop_Q3"
  ) %>%
  mutate(
    solute_group = case_when(
      str_detect(metric, "_all_") ~ "All Solutes (n=11)",
      str_detect(metric, "_nutrient_") ~ "Nutrient (n=1)",
      str_detect(metric, "_biogenic_") ~ "Biogenic (n=4)",
      str_detect(metric, "_geogenic_") ~ "Geogenic (n=6)",
      TRUE ~ "Unknown"
    )
  ) %>%
  select(Stream_Name, solute_group, prop_Q3) %>%
  filter(solute_group != "Unknown")

wymore_quadrants <- wymore_q1 %>%
  left_join(wymore_q3, by = c("Stream_Name", "solute_group")) %>%
  mutate(
    dominant_quadrant = case_when(
      prop_Q1 > prop_Q3 ~ "Dual Mobilizing (Q1)",
      prop_Q3 > prop_Q1 ~ "Dual Diluting (Q3)",
      TRUE ~ "Mixed"
    )
  ) %>%
  select(Stream_Name, solute_group, dominant_quadrant)

sync_plot <- sync_mv %>%
  mutate(Stream_Name = factor(Stream_Name, levels = site_order)) %>%
  select(Stream_Name,
         # Abbott multivariate
         matches("concentration_mv.*sync_allpairs"),
         # Wymore multivariate
         matches("wymore_mv.*prop_sync_allpairs")) %>%
  # Reshape to long format
  pivot_longer(
    cols = -Stream_Name,
    names_to = "metric",
    values_to = "sync_value"
  ) %>%
  mutate(
    # Extract method (Abbott or Wymore) and solute group (all, bio, geo)
    method = case_when(
      str_detect(metric, "concentration_mv") ~ "Abbott",
      str_detect(metric, "wymore_mv") ~ "Wymore",
      TRUE ~ "Unknown"
    ),
    solute_group = case_when(
      str_detect(metric, "_all_") ~ "All Solutes (n=11)",
      str_detect(metric, "_nutrient_") ~ "Nutrient (n=1)",
      str_detect(metric, "_biogenic_") ~ "Biogenic (n=4)",
      str_detect(metric, "_geogenic_") ~ "Geogenic (n=6)",
      TRUE ~ "Unknown"
    ),
    solute_group = factor(solute_group,
                         levels = c("All Solutes (n=11)", "Geogenic (n=6)", "Biogenic (n=4)", "Nutrient (n=1)"))
  ) %>%
  filter(!is.na(sync_value), method != "Unknown", solute_group != "Unknown")

# Reshape to wide format for correlation plot
sync_wide <- sync_plot %>%
  select(Stream_Name, method, solute_group, sync_value) %>%
  pivot_wider(names_from = method, values_from = sync_value) %>%
  filter(!is.na(Abbott), !is.na(Wymore)) %>%
  # Add dominant quadrant information
  left_join(wymore_quadrants, by = c("Stream_Name", "solute_group"))

message("\nPrepared ", nrow(sync_wide), " site-group combinations for plotting")

# Calculate correlations by solute group
correlations <- sync_wide %>%
  group_by(solute_group) %>%
  summarise(
    cor_r = cor(Abbott, Wymore, use = "complete.obs"),
    n = n(),
    .groups = "drop"
  )

message("\nCorrelations by solute group:")
for (i in seq_len(nrow(correlations))) {
  message("  ", correlations$solute_group[i], ": r = ",
          round(correlations$cor_r[i], 3), " (n=", correlations$n[i], ")")
}

# =============================================================================
# PLOT 1: SCATTER BY SOLUTE GROUP
# =============================================================================
message("\nCreating scatter plot by solute group...")

p1 <- ggplot(sync_wide, aes(x = Abbott, y = Wymore)) +
  # Add reference lines
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray60", alpha = 0.4) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray60", alpha = 0.4) +

  # Points colored by site, shaped by dominant CQ behavior
  geom_point(aes(color = Stream_Name, shape = dominant_quadrant), size = 4, alpha = 0.8, stroke = 1) +

  # Add text labels for stream names
  geom_text_repel(aes(label = Stream_Name), size = 2.5, max.overlaps = 20, seed = 42) +

  # Use plot_prefs site colors
  scale_color_site() +

  # Shape scale for Q1 vs Q3
  scale_shape_manual(
    values = c("Dual Mobilizing (Q1)" = 16,  # filled circle
               "Dual Diluting (Q3)" = 17,     # filled triangle
               "Mixed" = 15),                  # filled square
    name = "CQ Behavior",
    na.value = 4  # open circle for NA
  ) +

  # Facet by solute group
  facet_wrap(~ solute_group, ncol = 3) +

  labs(
    x = "Abbott Multivariate Concentration Synchrony\n(Do sites have similar chemical profiles?)",
    y = "Wymore Multivariate CQ Synchrony (proportion of all windows)\nHigher = greater proportion of windows with dual mobilizing (Q1) or dual diluting (Q3)",
    legend.position = "bottom",
    legend.box = "vertical",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "gray95", color = NA)
  ) +
  guides(color = guide_legend(nrow = 2, order = 1),
         shape = guide_legend(order = 2))

ggsave(file.path(fig_dir, "abbott_vs_wymore_multivariate_by_group.png"),
       p1, width = 14, height = 6, dpi = 300, bg = "white")

message("✓ Saved: abbott_vs_wymore_multivariate_by_group.png")

# =============================================================================
# PLOT 2: CORRELATION BY SOLUTE GROUP
# =============================================================================
message("\nCreating correlation bar plot...")

# Create a data frame with correlations
correlation_df <- correlations %>%
  mutate(
    cor_strength = case_when(
      abs(cor_r) > 0.7 ~ "Strong",
      abs(cor_r) > 0.4 ~ "Moderate",
      TRUE ~ "Weak"
    ),
    cor_direction = ifelse(cor_r > 0, "Positive", "Negative")
  )

p2 <- ggplot(correlation_df, aes(x = solute_group, y = cor_r)) +
  # Add reference line at 0
  geom_hline(yintercept = 0, linetype = "solid", color = "gray50", linewidth = 0.5) +
  geom_hline(yintercept = c(-0.7, 0.7), linetype = "dashed", color = "gray70", linewidth = 0.3) +

  # Bars colored by correlation strength
  geom_col(aes(fill = cor_r), width = 0.6, alpha = 0.8) +

  # Add correlation values as text
  geom_text(aes(label = sprintf("r = %.3f", cor_r),
                y = cor_r + ifelse(cor_r > 0, 0.05, -0.05)),
            size = 4, fontface = "bold") +

  scale_fill_gradient2(
    low = "#E07A5F", mid = "gray90", high = "#2E86AB",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Correlation"
  ) +

  labs(
    x = "Solute Group",
    y = "Pearson Correlation (r)",
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(fig_dir, "abbott_vs_wymore_multivariate_correlations.png"),
       p2, width = 10, height = 7, dpi = 300, bg = "white")

message("✓ Saved: abbott_vs_wymore_multivariate_correlations.png")

# =============================================================================
# OUTLET-SPECIFIC MULTIVARIATE SYNCHRONY
# =============================================================================
message("\n=== OUTLET-SPECIFIC MULTIVARIATE SYNCHRONY (with GSLOOK) ===")

# Extract outlet-specific Q1/Q3 proportions
wymore_q1_outlet <- sync_mv %>%
  select(Stream_Name, matches("wymore_mv.*prop_Q1_outlet")) %>%
  pivot_longer(
    cols = -Stream_Name,
    names_to = "metric",
    values_to = "prop_Q1"
  ) %>%
  mutate(
    solute_group = case_when(
      str_detect(metric, "_all_") ~ "All Solutes (n=11)",
      str_detect(metric, "_nutrient_") ~ "Nutrient (n=1)",
      str_detect(metric, "_biogenic_") ~ "Biogenic (n=4)",
      str_detect(metric, "_geogenic_") ~ "Geogenic (n=6)",
      TRUE ~ "Unknown"
    )
  ) %>%
  select(Stream_Name, solute_group, prop_Q1) %>%
  filter(solute_group != "Unknown")

wymore_q3_outlet <- sync_mv %>%
  select(Stream_Name, matches("wymore_mv.*prop_Q3_outlet")) %>%
  pivot_longer(
    cols = -Stream_Name,
    names_to = "metric",
    values_to = "prop_Q3"
  ) %>%
  mutate(
    solute_group = case_when(
      str_detect(metric, "_all_") ~ "All Solutes (n=11)",
      str_detect(metric, "_nutrient_") ~ "Nutrient (n=1)",
      str_detect(metric, "_biogenic_") ~ "Biogenic (n=4)",
      str_detect(metric, "_geogenic_") ~ "Geogenic (n=6)",
      TRUE ~ "Unknown"
    )
  ) %>%
  select(Stream_Name, solute_group, prop_Q3) %>%
  filter(solute_group != "Unknown")

wymore_quadrants_outlet <- wymore_q1_outlet %>%
  left_join(wymore_q3_outlet, by = c("Stream_Name", "solute_group")) %>%
  mutate(
    dominant_quadrant = case_when(
      prop_Q1 > prop_Q3 ~ "Dual Mobilizing (Q1)",
      prop_Q3 > prop_Q1 ~ "Dual Diluting (Q3)",
      TRUE ~ "Mixed"
    )
  ) %>%
  select(Stream_Name, solute_group, dominant_quadrant)

# Prepare outlet synchrony data for plotting
sync_outlet_plot <- sync_mv %>%
  mutate(Stream_Name = factor(Stream_Name, levels = site_order)) %>%
  select(Stream_Name,
         # Abbott multivariate outlet
         matches("concentration_mv.*sync_outlet"),
         # Wymore multivariate outlet
         matches("wymore_mv.*prop_sync_outlet")) %>%
  # Reshape to long format
  pivot_longer(
    cols = -Stream_Name,
    names_to = "metric",
    values_to = "sync_value"
  ) %>%
  mutate(
    # Extract method (Abbott or Wymore) and solute group (all, bio, geo)
    method = case_when(
      str_detect(metric, "concentration_mv") ~ "Abbott",
      str_detect(metric, "wymore_mv") ~ "Wymore",
      TRUE ~ "Unknown"
    ),
    solute_group = case_when(
      str_detect(metric, "_all_") ~ "All Solutes (n=11)",
      str_detect(metric, "_nutrient_") ~ "Nutrient (n=1)",
      str_detect(metric, "_biogenic_") ~ "Biogenic (n=4)",
      str_detect(metric, "_geogenic_") ~ "Geogenic (n=6)",
      TRUE ~ "Unknown"
    ),
    solute_group = factor(solute_group,
                         levels = c("All Solutes (n=11)", "Geogenic (n=6)", "Biogenic (n=4)", "Nutrient (n=1)"))
  ) %>%
  filter(!is.na(sync_value), method != "Unknown", solute_group != "Unknown")

# Reshape to wide format for correlation plot
sync_outlet_wide <- sync_outlet_plot %>%
  select(Stream_Name, method, solute_group, sync_value) %>%
  pivot_wider(names_from = method, values_from = sync_value) %>%
  filter(!is.na(Abbott), !is.na(Wymore)) %>%
  # Add dominant quadrant information
  left_join(wymore_quadrants_outlet, by = c("Stream_Name", "solute_group"))

message("Prepared ", nrow(sync_outlet_wide), " site-group combinations for outlet plotting")

# Calculate correlations by solute group for outlet
correlations_outlet <- sync_outlet_wide %>%
  group_by(solute_group) %>%
  summarise(
    cor_r = cor(Abbott, Wymore, use = "complete.obs"),
    n = n(),
    .groups = "drop"
  )

message("\nOutlet correlations by solute group:")
for (i in seq_len(nrow(correlations_outlet))) {
  message("  ", correlations_outlet$solute_group[i], ": r = ",
          round(correlations_outlet$cor_r[i], 3), " (n=", correlations_outlet$n[i], ")")
}

# =============================================================================
# PLOT 3: OUTLET SCATTER BY SOLUTE GROUP
# =============================================================================
message("\nCreating outlet scatter plot by solute group...")

p3 <- ggplot(sync_outlet_wide, aes(x = Abbott, y = Wymore)) +
  # Add reference lines
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray60", alpha = 0.4) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray60", alpha = 0.4) +

  # Points colored by site, shaped by dominant CQ behavior
  geom_point(aes(color = Stream_Name, shape = dominant_quadrant), size = 4, alpha = 0.8, stroke = 1) +

  # Add text labels for stream names
  geom_text_repel(aes(label = Stream_Name), size = 2.5, max.overlaps = 20, seed = 42) +

  # Use plot_prefs site colors
  scale_color_site() +

  # Shape scale for Q1 vs Q3
  scale_shape_manual(
    values = c("Dual Mobilizing (Q1)" = 16,  # filled circle
               "Dual Diluting (Q3)" = 17,     # filled triangle
               "Mixed" = 15),                  # filled square
    name = "CQ Behavior",
    na.value = 4  # open circle for NA
  ) +

  # Facet by solute group
  facet_wrap(~ solute_group, ncol = 3) +

  labs(
    x = "Abbott Multivariate Concentration Synchrony (with outlet)\n(Does site have similar chemical profile to GSLOOK?)",
    y = "Wymore Multivariate CQ Synchrony (proportion of all windows)\nHigher = greater proportion of windows with dual mobilizing (Q1) or dual diluting (Q3)",
    legend.position = "bottom",
    legend.box = "vertical",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "gray95", color = NA)
  ) +
  guides(color = guide_legend(nrow = 2, order = 1),
         shape = guide_legend(order = 2))

ggsave(file.path(fig_dir, "abbott_vs_wymore_multivariate_outlet_by_group.png"),
       p3, width = 14, height = 6, dpi = 300, bg = "white")

message("✓ Saved: abbott_vs_wymore_multivariate_outlet_by_group.png")

# =============================================================================
# PLOT 4: OUTLET CORRELATION BY SOLUTE GROUP
# =============================================================================
message("\nCreating outlet correlation bar plot...")

# Create a data frame with outlet correlations
correlation_outlet_df <- correlations_outlet %>%
  mutate(
    cor_strength = case_when(
      abs(cor_r) > 0.7 ~ "Strong",
      abs(cor_r) > 0.4 ~ "Moderate",
      TRUE ~ "Weak"
    ),
    cor_direction = ifelse(cor_r > 0, "Positive", "Negative")
  )

p4 <- ggplot(correlation_outlet_df, aes(x = solute_group, y = cor_r)) +
  # Add reference line at 0
  geom_hline(yintercept = 0, linetype = "solid", color = "gray50", linewidth = 0.5) +
  geom_hline(yintercept = c(-0.7, 0.7), linetype = "dashed", color = "gray70", linewidth = 0.3) +

  # Bars colored by correlation strength
  geom_col(aes(fill = cor_r), width = 0.6, alpha = 0.8) +

  # Add correlation values as text
  geom_text(aes(label = sprintf("r = %.3f", cor_r),
                y = cor_r + ifelse(cor_r > 0, 0.05, -0.05)),
            size = 4, fontface = "bold") +

  scale_fill_gradient2(
    low = "#E07A5F", mid = "gray90", high = "#2E86AB",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Correlation"
  ) +

  labs(
    x = "Solute Group",
    y = "Pearson Correlation (r)",
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(fig_dir, "abbott_vs_wymore_multivariate_outlet_correlations.png"),
       p4, width = 10, height = 7, dpi = 300, bg = "white")

message("✓ Saved: abbott_vs_wymore_multivariate_outlet_correlations.png")

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================
message("\n=== SUMMARY STATISTICS ===")

summary_stats <- sync_plot %>%
  group_by(method, solute_group) %>%
  summarise(
    mean_sync = mean(sync_value, na.rm = TRUE),
    sd_sync = sd(sync_value, na.rm = TRUE),
    min_sync = min(sync_value, na.rm = TRUE),
    max_sync = max(sync_value, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

message("\nSynchrony by Method and Solute Group:")
for (i in seq_len(nrow(summary_stats))) {
  message(sprintf("  %s - %s: %.3f ± %.3f (range: %.3f - %.3f, n=%d)",
                 summary_stats$method[i],
                 summary_stats$solute_group[i],
                 summary_stats$mean_sync[i],
                 summary_stats$sd_sync[i],
                 summary_stats$min_sync[i],
                 summary_stats$max_sync[i],
                 summary_stats$n[i]))
}

message("\n=== INTERPRETATION ===")
message("Abbott Multivariate: Measures whether sites have synchronized CHEMICAL PROFILES")
message("  - Higher = sites show similar multi-solute patterns over time")
message("  - Captures overall chemical similarity, not just individual solutes")
message("\nWymore Multivariate: Measures proportion of solutes mobilized similarly")
message("  - Higher = more solutes show synchronized CQ behavior (Q1 or Q3)")
message("  - Captures consistency of mobilization mechanisms across chemistry")

# Compare biogenic vs geogenic synchrony
bio_abbott <- sync_plot %>% filter(method == "Abbott", solute_group == "Biogenic (n=4)") %>% pull(sync_value)
geo_abbott <- sync_plot %>% filter(method == "Abbott", solute_group == "Geogenic (n=6)") %>% pull(sync_value)
bio_wymore <- sync_plot %>% filter(method == "Wymore", solute_group == "Biogenic (n=4)") %>% pull(sync_value)
geo_wymore <- sync_plot %>% filter(method == "Wymore", solute_group == "Geogenic (n=6)") %>% pull(sync_value)

message("\nBio vs Geo comparison:")
message(sprintf("  Abbott: Bio mean = %.3f, Geo mean = %.3f (Δ = %.3f)",
               mean(bio_abbott, na.rm = TRUE),
               mean(geo_abbott, na.rm = TRUE),
               mean(bio_abbott, na.rm = TRUE) - mean(geo_abbott, na.rm = TRUE)))
message(sprintf("  Wymore: Bio mean = %.3f, Geo mean = %.3f (Δ = %.3f)",
               mean(bio_wymore, na.rm = TRUE),
               mean(geo_wymore, na.rm = TRUE),
               mean(bio_wymore, na.rm = TRUE) - mean(geo_wymore, na.rm = TRUE)))

message("\n=== COMPLETE ===\n")
