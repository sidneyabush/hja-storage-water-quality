#!/usr/bin/env Rscript
# =============================================================================
# Plot Abbott Concentration Sync vs Wymore CQ Slope Sync
# =============================================================================
# Purpose: Visualize the relationship between concentration synchrony and
#          CQ behavior synchrony across sites
#
# INTERPRETATION:
#   - Abbott conc sync: Do sites have HIGH/LOW concentrations together?
#   - Wymore CQ sync: Do sites MOBILIZE solutes the same way (proportion)
#
# QUADRANTS:
#   High Abbott, High Wymore: Sites synchronize in BOTH conc AND mobilization
#   High Abbott, Low Wymore: Sites synchronize in conc but DIFFERENT mobilization
#   Low Abbott, High Wymore: Sites have DIFFERENT concs but SIMILAR mobilization
#   Low Abbott, Low Wymore: Sites are INDEPENDENT in both
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
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

message("=== ABBOTT vs WYMORE SYNCHRONY ===\n")

# Load composite synchrony
sync_data <- readr::read_csv(file.path(out_dir, "HJA_composite_synchrony.csv"),
                             show_col_types = FALSE) %>%
  filter(!is.na(conc_sync_allpairs), !is.na(wymore_crosssite_allpairs)) %>%
  mutate(
    # Classify dominant CQ behavior (Q1 = dual mobilizing, Q3 = dual diluting)
    dominant_behavior = case_when(
      wymore_Q1_proportion > wymore_Q3_proportion ~ "Dual Mobilizing (Q1)",
      wymore_Q3_proportion > wymore_Q1_proportion ~ "Dual Diluting (Q3)",
      TRUE ~ "Mixed"
    ),
    # Convert solute and site to factors with proper ordering for plot_prefs
    solute = factor(solute, levels = names(solute_colors)),
    Stream_Name = factor(Stream_Name, levels = site_order)
  )

message("Loaded ", nrow(sync_data), " site-solute combinations\n")

# Summary statistics
message("Abbott Concentration Sync (allpairs):")
message("  Range: ", round(min(sync_data$conc_sync_allpairs, na.rm = TRUE), 3), " - ",
        round(max(sync_data$conc_sync_allpairs, na.rm = TRUE), 3))
message("  Mean: ", round(mean(sync_data$conc_sync_allpairs, na.rm = TRUE), 3))

message("\nWymore CQ Slope Sync (proportion synchronized):")
message("  Range: ", round(min(sync_data$wymore_crosssite_allpairs, na.rm = TRUE), 3), " - ",
        round(max(sync_data$wymore_crosssite_allpairs, na.rm = TRUE), 3))
message("  Mean: ", round(mean(sync_data$wymore_crosssite_allpairs, na.rm = TRUE), 3))

# Correlation
corr <- cor(sync_data$conc_sync_allpairs, sync_data$wymore_crosssite_allpairs,
            use = "complete.obs")
message("\nCorrelation: r = ", round(corr, 3))

# =============================================================================
# PLOT 1: FACET BY SOLUTE, COLOR BY SITE
# =============================================================================
message("\nCreating Plot 1: Facet by solute, color by site...")

p1 <- ggplot(sync_data, aes(x = conc_sync_allpairs, y = wymore_crosssite_allpairs)) +
  # Add quadrant dividers at overall medians
  geom_vline(xintercept = median(sync_data$conc_sync_allpairs, na.rm = TRUE),
             linetype = "dashed", color = "gray60", alpha = 0.4, linewidth = 0.5) +
  geom_hline(yintercept = median(sync_data$wymore_crosssite_allpairs, na.rm = TRUE),
             linetype = "dashed", color = "gray60", alpha = 0.4, linewidth = 0.5) +

  # Points colored by site, shaped by dominant CQ behavior
  geom_point(aes(color = Stream_Name, shape = dominant_behavior),
             size = 3.5, alpha = 0.8, stroke = 1) +

  # Add text labels for stream names
  geom_text_repel(aes(label = Stream_Name), size = 2.5, max.overlaps = 20, seed = 42) +

  # Use plot_prefs site color scale
  scale_color_site() +

  # Shape scale for Q1 vs Q3
  scale_shape_manual(
    values = c("Dual Mobilizing (Q1)" = 16,  # filled circle
               "Dual Diluting (Q3)" = 17,     # filled triangle
               "Mixed" = 15),                  # filled square
    name = "CQ Behavior"
  ) +

  # Facet by solute
  facet_wrap(~ solute, ncol = 4) +

  labs(
    x = "Abbott Concentration Synchrony (all pairs)\nHigher = concentrations vary together across sites",
    y = "Wymore CQ Slope Synchrony (proportion of all windows)\nHigher = greater proportion of windows with dual mobilizing (Q1) or dual diluting (Q3)",
    legend.position = "bottom",
    legend.box = "vertical",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "gray95", color = NA)
  )

ggsave(file.path(fig_dir, "abbott_vs_wymore_synchrony_by_solute.png"),
       p1, width = 14, height = 10, dpi = 300, bg = "white")

message("✓ Plot 1 saved: abbott_vs_wymore_synchrony_by_solute.png")

# =============================================================================
# PLOT 2: FACET BY SITE, COLOR BY SOLUTE
# =============================================================================
message("\nCreating Plot 2: Facet by site, color by solute...")

p2 <- ggplot(sync_data, aes(x = conc_sync_allpairs, y = wymore_crosssite_allpairs)) +
  # Add quadrant dividers at overall medians
  geom_vline(xintercept = median(sync_data$conc_sync_allpairs, na.rm = TRUE),
             linetype = "dashed", color = "gray60", alpha = 0.4, linewidth = 0.5) +
  geom_hline(yintercept = median(sync_data$wymore_crosssite_allpairs, na.rm = TRUE),
             linetype = "dashed", color = "gray60", alpha = 0.4, linewidth = 0.5) +

  # Points colored by solute, shaped by dominant CQ behavior
  geom_point(aes(color = solute, shape = dominant_behavior),
             size = 3.5, alpha = 0.8, stroke = 1) +

  # Add text labels for solute names
  geom_text_repel(aes(label = solute), size = 2.5, max.overlaps = 20, seed = 42) +

  # Use plot_prefs solute color scale
  scale_color_solute() +

  # Shape scale for Q1 vs Q3
  scale_shape_manual(
    values = c("Dual Mobilizing (Q1)" = 16,  # filled circle
               "Dual Diluting (Q3)" = 17,     # filled triangle
               "Mixed" = 15),                  # filled square
    name = "CQ Behavior"
  ) +

  # Facet by site
  facet_wrap(~ Stream_Name, ncol = 4) +

  labs(
    x = "Abbott Concentration Synchrony (all pairs)\nHigher = concentrations vary together across sites",
    y = "Wymore CQ Slope Synchrony (proportion of all windows)\nHigher = greater proportion of windows with dual mobilizing (Q1) or dual diluting (Q3)",
    legend.position = "bottom",
    legend.box = "vertical",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "gray95", color = NA)
  )

ggsave(file.path(fig_dir, "abbott_vs_wymore_synchrony_by_site.png"),
       p2, width = 14, height = 10, dpi = 300, bg = "white")

message("✓ Plot 2 saved: abbott_vs_wymore_synchrony_by_site.png")

# =============================================================================
# OUTLET-SPECIFIC SYNCHRONY PLOTS (GSLOOK only)
# =============================================================================
message("\n=== OUTLET-SPECIFIC SYNCHRONY (with GSLOOK) ===")

# Prepare outlet-specific data
sync_outlet <- sync_data %>%
  filter(!is.na(conc_sync_outlet), !is.na(wymore_crosssite_outlet)) %>%
  mutate(
    # Classify dominant CQ behavior for outlet pairs
    dominant_behavior_outlet = case_when(
      wymore_Q1_proportion_outlet > wymore_Q3_proportion_outlet ~ "Dual Mobilizing (Q1)",
      wymore_Q3_proportion_outlet > wymore_Q1_proportion_outlet ~ "Dual Diluting (Q3)",
      TRUE ~ "Mixed"
    )
  )

message("Loaded ", nrow(sync_outlet), " site-solute combinations for outlet analysis\n")

# Correlation for outlet pairs
corr_outlet <- cor(sync_outlet$conc_sync_outlet, sync_outlet$wymore_crosssite_outlet,
                   use = "complete.obs")
message("Outlet correlation: r = ", round(corr_outlet, 3))

# =============================================================================
# PLOT 3: OUTLET SYNCHRONY - FACET BY SOLUTE
# =============================================================================
message("\nCreating Plot 3: Outlet synchrony faceted by solute...")

p3 <- ggplot(sync_outlet, aes(x = conc_sync_outlet, y = wymore_crosssite_outlet)) +
  # Add quadrant dividers at overall medians
  geom_vline(xintercept = median(sync_outlet$conc_sync_outlet, na.rm = TRUE),
             linetype = "dashed", color = "gray60", alpha = 0.4, linewidth = 0.5) +
  geom_hline(yintercept = median(sync_outlet$wymore_crosssite_outlet, na.rm = TRUE),
             linetype = "dashed", color = "gray60", alpha = 0.4, linewidth = 0.5) +

  # Points colored by site, shaped by dominant CQ behavior
  geom_point(aes(color = Stream_Name, shape = dominant_behavior_outlet),
             size = 3.5, alpha = 0.8, stroke = 1) +

  # Add text labels for stream names
  geom_text_repel(aes(label = Stream_Name), size = 2.5, max.overlaps = 20, seed = 42) +

  # Use plot_prefs site color scale
  scale_color_site() +

  # Shape scale for Q1 vs Q3
  scale_shape_manual(
    values = c("Dual Mobilizing (Q1)" = 16,  # filled circle
               "Dual Diluting (Q3)" = 17,     # filled triangle
               "Mixed" = 15),                  # filled square
    name = "CQ Behavior"
  ) +

  # Facet by solute
  facet_wrap(~ solute, ncol = 4) +

  labs(
    x = "Abbott Concentration Synchrony (with outlet)\nHigher = concentrations vary together with GSLOOK",
    y = "Wymore CQ Slope Synchrony (proportion of all windows)\nHigher = greater proportion of windows with dual mobilizing (Q1) or dual diluting (Q3)",
    legend.position = "bottom",
    legend.box = "vertical",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "gray95", color = NA)
  )

ggsave(file.path(fig_dir, "abbott_vs_wymore_synchrony_outlet_by_solute.png"),
       p3, width = 14, height = 10, dpi = 300, bg = "white")

message("✓ Plot 3 saved: abbott_vs_wymore_synchrony_outlet_by_solute.png")

# =============================================================================
# PLOT 4: OUTLET SYNCHRONY - FACET BY SITE
# =============================================================================
message("\nCreating Plot 4: Outlet synchrony faceted by site...")

p4 <- ggplot(sync_outlet, aes(x = conc_sync_outlet, y = wymore_crosssite_outlet)) +
  # Add quadrant dividers at overall medians
  geom_vline(xintercept = median(sync_outlet$conc_sync_outlet, na.rm = TRUE),
             linetype = "dashed", color = "gray60", alpha = 0.4, linewidth = 0.5) +
  geom_hline(yintercept = median(sync_outlet$wymore_crosssite_outlet, na.rm = TRUE),
             linetype = "dashed", color = "gray60", alpha = 0.4, linewidth = 0.5) +

  # Points colored by solute, shaped by dominant CQ behavior
  geom_point(aes(color = solute, shape = dominant_behavior_outlet),
             size = 3.5, alpha = 0.8, stroke = 1) +

  # Add text labels for solute names
  geom_text_repel(aes(label = solute), size = 2.5, max.overlaps = 20, seed = 42) +

  # Use plot_prefs solute color scale
  scale_color_solute() +

  # Shape scale for Q1 vs Q3
  scale_shape_manual(
    values = c("Dual Mobilizing (Q1)" = 16,  # filled circle
               "Dual Diluting (Q3)" = 17,     # filled triangle
               "Mixed" = 15),                  # filled square
    name = "CQ Behavior"
  ) +

  # Facet by site
  facet_wrap(~ Stream_Name, ncol = 4) +

  labs(
    x = "Abbott Concentration Synchrony (with outlet)\nHigher = concentrations vary together with GSLOOK",
    y = "Wymore CQ Slope Synchrony (proportion of all windows)\nHigher = greater proportion of windows with dual mobilizing (Q1) or dual diluting (Q3)",
    legend.position = "bottom",
    legend.box = "vertical",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "gray95", color = NA)
  )

ggsave(file.path(fig_dir, "abbott_vs_wymore_synchrony_outlet_by_site.png"),
       p4, width = 14, height = 10, dpi = 300, bg = "white")

message("✓ Plot 4 saved: abbott_vs_wymore_synchrony_outlet_by_site.png")

# Classify into quadrants
sync_data_classified <- sync_data %>%
  mutate(
    abbott_cat = ifelse(conc_sync_allpairs > median(conc_sync_allpairs, na.rm = TRUE),
                        "High Abbott", "Low Abbott"),
    wymore_cat = ifelse(wymore_crosssite_allpairs > median(wymore_crosssite_allpairs, na.rm = TRUE),
                        "High Wymore", "Low Wymore"),
    quadrant = paste(abbott_cat, wymore_cat, sep = " + ")
  )

# Summary by quadrant
message("\n=== QUADRANT DISTRIBUTION ===")
quadrant_summary <- sync_data_classified %>%
  count(quadrant) %>%
  arrange(desc(n))

for (i in seq_len(nrow(quadrant_summary))) {
  message("  ", quadrant_summary$quadrant[i], ": ", quadrant_summary$n[i],
          " (", round(100 * quadrant_summary$n[i] / sum(quadrant_summary$n), 1), "%)")
}

message("\n=== INTERPRETATION ===")
message("High Abbott + High Wymore: Sites track together in BOTH concentration")
message("                           AND mobilization behavior (most synchronized)")
message("\nHigh Abbott + Low Wymore: Concentrations sync but DIFFERENT mechanisms")
message("                          (e.g., both high, but one dilutes, one enriches)")
message("\nLow Abbott + High Wymore: Different concentrations but SIMILAR mobilization")
message("                          (same CQ slopes despite different absolute concs)")
message("\nLow Abbott + Low Wymore: Independent behavior in both")
message("                         (each site responds differently)")
message("\n=== Q1 vs Q3 BEHAVIOR ===")
message("Dual Mobilizing (Q1): Sites both show POSITIVE CQ slopes (enrich with Q)")
message("Dual Diluting (Q3): Sites both show NEGATIVE CQ slopes (dilute with Q)")
message("Mixed: Sites alternate between Q1 and Q3 across different time windows")

message("\n=== COMPLETE ===\n")
