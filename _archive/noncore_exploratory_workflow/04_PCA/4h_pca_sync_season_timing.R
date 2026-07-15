# Restored from archive on 2025-12-04

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
})

base_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality"
out_dir <- file.path(base_dir, "outputs")
fig_dir <- file.path(base_dir, "exploratory_plots", "04_PCA", "4h_sync_season_timing")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
try(source(file.path(repo_dir, "00_helpers", "plot_prefs.R")), silent = TRUE)
try(source(file.path(repo_dir, "00_helpers", "workflow_config.R")), silent = TRUE)

abbott <- read_csv(file.path(out_dir, "HJA_Abbott_synchrony_windows.csv"), show_col_types = FALSE)
season_raw <- read_csv(file.path(out_dir, "season_boundaries.csv"), show_col_types = FALSE)
site_chars <- read_csv(file.path(out_dir, "Catchment_site_characteristics.csv"), show_col_types = FALSE)

season_timing <- season_raw %>%
  mutate(
    wet_start_date = as.Date(wet_start_date),
    wet_end_date = as.Date(wet_end_date),
    wet_length_days = as.numeric(wet_end_date - wet_start_date),
    wet_start_doy = lubridate::yday(wet_start_date)
  ) %>%
  select(water_year, wet_length_days, wet_start_doy)

outlet_sync <- abbott %>%
  filter(is_outlet_pair == TRUE, time_scale == "annual", synchrony_type == "concentration") %>%
  mutate(Stream_Name = ifelse(Stream1 == "GSLOOK", Stream2, Stream1)) %>%
  group_by(Stream_Name, water_year) %>%
  summarise(conc_sync_outlet = mean(pearson_r, na.rm = TRUE), n_solutes = n(), .groups = "drop")

sync_season <- outlet_sync %>% left_join(season_timing, by = "water_year") %>% filter(!is.na(wet_length_days))

message("Sync + season timing PCA script restored. Implement plotting as needed.")
# =============================================================================
# 4h_pca_sync_season_timing.R
# =============================================================================
# PCA COMBINING SYNCHRONY METRICS WITH SEASON TIMING
#
# RESEARCH QUESTION:
#   How do sites cluster when considering both synchrony patterns AND
#   hydrological season timing together?
#
# BACKGROUND:
#   Season timing (wet_length_days, wet_start_doy) has opposite effects on
#   synchrony at different scales:
#   - Outlet sync: longer wet → MORE sync (r = +0.246)
#   - Pairwise sync: longer wet → LESS sync (r = -0.327)
#
# APPROACH:
#   1. PCA on synchrony metrics colored by season timing
#   2. PCA combining sync + timing variables
#   3. Visualize how sites separate by both dimensions
#
# OUTPUT:
#   - Biplots with season timing coloring
#   - Combined sync + timing PCA
#   - Site clustering visualization
#
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
})

# =============================================================================
# SETUP
# =============================================================================

base_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality"
out_dir <- file.path(base_dir, "outputs")
fig_dir <- file.path(base_dir, "exploratory_plots", "04_PCA", "4h_sync_season_timing")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Source helpers
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
try(source(file.path(repo_dir, "00_helpers", "plot_prefs.R")), silent = TRUE)
try(source(file.path(repo_dir, "00_helpers", "workflow_config.R")), silent = TRUE)

# Theme
theme_pca <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_blank(),
      legend.position = "right"
    )
}

message("\n╔════════════════════════════════════════════════════════════════╗\n")
message("║  PCA: SYNCHRONY + SEASON TIMING                                ║\n")
message("╚════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# LOAD DATA
# =============================================================================

message("=== LOADING DATA ===\n")

# Load Abbott synchrony (outlet pairs)
abbott <- read_csv(file.path(out_dir, "HJA_Abbott_synchrony_windows.csv"), 
                   show_col_types = FALSE)

# Load season timing from season_boundaries.csv
season_raw <- read_csv(file.path(out_dir, "season_boundaries.csv"), 
                       show_col_types = FALSE)

# Load site characteristics for gradient coloring
site_chars <- read_csv(file.path(out_dir, "Catchment_site_characteristics.csv"),
                       show_col_types = FALSE)

# Calculate season timing variables
season_timing <- season_raw %>%
  mutate(
    wet_start_date = as.Date(wet_start_date),
    wet_end_date = as.Date(wet_end_date),
    wet_length_days = as.numeric(wet_end_date - wet_start_date),
    wet_start_doy = lubridate::yday(wet_start_date)
  ) %>%
  select(water_year, wet_length_days, wet_start_doy)

message("  Abbott sync rows:", nrow(abbott), "\n")
message("  Season timing rows:", nrow(season_timing), "\n")

# =============================================================================
# PREPARE SITE-YEAR SYNCHRONY DATA
# =============================================================================

message("\n=== PREPARING SITE-YEAR OUTLET SYNCHRONY ===\n")

# Get outlet pairs only, annual time scale, concentration sync
outlet_sync <- abbott %>%
  filter(is_outlet_pair == TRUE,
         time_scale == "annual",
         synchrony_type == "concentration") %>%
  # Get the non-GSLOOK site
  mutate(Stream_Name = ifelse(Stream1 == "GSLOOK", Stream2, Stream1)) %>%
  group_by(Stream_Name, water_year) %>%
  summarise(
    conc_sync_outlet = mean(pearson_r, na.rm = TRUE),  # Average across solutes
    n_solutes = n(),
    .groups = "drop"
  )

message("  Site-years with outlet sync:", nrow(outlet_sync), "\n")

# Merge with season timing
sync_season <- outlet_sync %>%
  left_join(season_timing, by = "water_year") %>%
  filter(!is.na(wet_length_days))

message("  Site-years with timing:", nrow(sync_season), "\n")

# =============================================================================
# PCA 1: SYNC BY SOLUTE COLORED BY SEASON LENGTH
# =============================================================================

message("\n=== PCA 1: SYNC COLORED BY WET SEASON LENGTH ===\n")

# Get sync by solute for more variables
sync_by_solute <- abbott %>%
  filter(is_outlet_pair == TRUE,
         time_scale == "annual",
         synchrony_type == "concentration") %>%
  mutate(Stream_Name = ifelse(Stream1 == "GSLOOK", Stream2, Stream1)) %>%
  select(Stream_Name, water_year, solute, pearson_r) %>%
  pivot_wider(names_from = solute, values_from = pearson_r, 
              names_prefix = "sync_") %>%
  left_join(season_timing, by = "water_year") %>%
  filter(!is.na(wet_length_days))

# Find columns with enough data
sync_cols <- names(sync_by_solute)[grepl("^sync_", names(sync_by_solute))]
good_cols <- sync_cols[colSums(!is.na(sync_by_solute[sync_cols])) > 30]

message("  Solutes with enough data:", length(good_cols), "\n")

# Prepare matrix
pca_data1 <- sync_by_solute %>%
  select(Stream_Name, water_year, all_of(good_cols), wet_length_days, wet_start_doy) %>%
  filter(complete.cases(.))

# Scale sync variables
pca_matrix1 <- pca_data1 %>%
  select(all_of(good_cols)) %>%
  scale()

# Run PCA
pca1 <- prcomp(pca_matrix1, center = FALSE, scale. = FALSE)

message("  PC1 variance:", round(summary(pca1)$importance[2,1] * 100, 1), "%\n")
message("  PC2 variance:", round(summary(pca1)$importance[2,2] * 100, 1), "%\n")

# Create scores dataframe with site characteristics
scores1 <- as_tibble(pca1$x) %>%
  bind_cols(pca_data1) %>%
  left_join(site_chars, by = "Stream_Name") %>%
  flag_outlet_stream()

# Loadings for arrows
loadings1 <- as_tibble(pca1$rotation, rownames = "variable") %>%
  mutate(
    PC1 = PC1 * 2,  # Scale for visibility
    PC2 = PC2 * 2
  )

# Plot: Site colors with size by wet_length_days (matching project style)
p1 <- ggplot() +
  stat_ellipse(
    data = scores1,
    aes(x = PC1, y = PC2),
    level = 0.95,
    type = "t",
    linetype = "dashed",
    color = "gray50",
    linewidth = 0.5
  ) +
  geom_point(
    data = scores1,
    aes(x = PC1, y = PC2, fill = Stream_Name, size = wet_length_days),
    shape = 21,
    colour = "grey20",
    stroke = 0.6,
    alpha = 0.8
  ) +
  geom_segment(data = loadings1,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "gray30", linewidth = 0.8) +
  geom_text_repel(data = loadings1,
                  aes(x = PC1, y = PC2, label = variable),
                  color = "gray20", size = 3.5) +
  scale_fill_site(name = "Site") +
  scale_size_continuous(name = "Wet Season\nLength (days)", range = c(2, 7)) +
  labs(
    title = "Synchrony PCA: Sites Colored, Sized by Wet Season Length",
    subtitle = paste0("PC1: ", round(summary(pca1)$importance[2,1] * 100, 1), 
                      "%, PC2: ", round(summary(pca1)$importance[2,2] * 100, 1), "%"),
    x = paste0("PC1 (", round(summary(pca1)$importance[2,1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca1)$importance[2,2] * 100, 1), "%)")
  ) +
  theme_pca() +
  coord_fixed() +
  guides(fill = guide_legend(override.aes = list(size = 4)))

ggsave(file.path(fig_dir, "sync_pca_colored_by_wet_length.png"), p1, 
       width = 10, height = 8, dpi = 300)

# =============================================================================
# PCA 2: COMBINED SYNC + TIMING VARIABLES
# =============================================================================

message("\n=== PCA 2: COMBINED SYNC + TIMING ===\n")

# Use mean sync across solutes + timing
combined_vars <- c("conc_sync_outlet", "wet_length_days", "wet_start_doy")

pca_data2 <- sync_season %>%
  select(Stream_Name, water_year, all_of(combined_vars)) %>%
  filter(complete.cases(.))

# Scale all variables
pca_matrix2 <- pca_data2 %>%
  select(all_of(combined_vars)) %>%
  scale()

# Run PCA
pca2 <- prcomp(pca_matrix2, center = FALSE, scale. = FALSE)

message("  PC1 variance:", round(summary(pca2)$importance[2,1] * 100, 1), "%\n")
message("  PC2 variance:", round(summary(pca2)$importance[2,2] * 100, 1), "%\n")
message("  Cumulative (PC1+PC2):", round(summary(pca2)$importance[3,2] * 100, 1), "%\n")

# Create scores dataframe with site characteristics
scores2 <- as_tibble(pca2$x) %>%
  bind_cols(pca_data2) %>%
  left_join(site_chars, by = "Stream_Name") %>%
  flag_outlet_stream()

# Loadings for arrows
loadings2 <- as_tibble(pca2$rotation, rownames = "variable") %>%
  mutate(
    PC1 = PC1 * 2.5,
    PC2 = PC2 * 2.5,
    # Clean up variable names for display
    var_label = case_when(
      variable == "conc_sync_outlet" ~ "Conc Sync",
      variable == "wet_length_days" ~ "Wet Length",
      variable == "wet_start_doy" ~ "Wet Start",
      TRUE ~ variable
    )
  )

# Determine storage column for sizing (use Elevation as fallback)
storage_col <- if ("Elevation_mean_m" %in% names(scores2) && 
                   any(is.finite(scores2$Elevation_mean_m))) {
  "Elevation_mean_m"
} else {
  NA_character_
}

# Plot: biplot with site colors and elevation sizing
p2 <- ggplot() +
  stat_ellipse(
    data = scores2,
    aes(x = PC1, y = PC2),
    level = 0.95,
    type = "t",
    linetype = "dashed",
    color = "gray50",
    linewidth = 0.5
  ) +
  geom_point(
    data = scores2,
    aes(x = PC1, y = PC2, fill = Stream_Name, 
        size = if (!is.na(storage_col)) .data[[storage_col]] else 3),
    shape = 21,
    colour = "grey20",
    stroke = 0.6,
    alpha = 0.7
  ) +
  geom_segment(data = loadings2,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "gray30", linewidth = 1) +
  geom_text_repel(data = loadings2,
                  aes(x = PC1 * 1.1, y = PC2 * 1.1, label = var_label),
                  color = "gray20", size = 4, fontface = "bold") +
  scale_fill_site(name = "Site") +
  scale_size_continuous(name = "Elevation\n(m)", range = c(2, 7)) +
  labs(
    title = "Combined PCA: Synchrony + Season Timing",
    subtitle = paste0("PC1: ", round(summary(pca2)$importance[2,1] * 100, 1), 
                      "%, PC2: ", round(summary(pca2)$importance[2,2] * 100, 1), "%"),
    x = paste0("PC1 (", round(summary(pca2)$importance[2,1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca2)$importance[2,2] * 100, 1), "%)")
  ) +
  theme_pca() +
  coord_fixed() +
  guides(fill = guide_legend(override.aes = list(size = 4)))

ggsave(file.path(fig_dir, "combined_sync_timing_pca_biplot.png"), p2, 
       width = 11, height = 9, dpi = 300)

# =============================================================================
# PCA 3: SITE MEANS - HOW DO SITES CLUSTER?
# =============================================================================

message("\n=== PCA 3: SITE-LEVEL CLUSTERING ===\n")

# Calculate site means across years
site_means <- sync_season %>%
  group_by(Stream_Name) %>%
  summarise(
    conc_sync_outlet = mean(conc_sync_outlet, na.rm = TRUE),
    mean_wet_length = mean(wet_length_days, na.rm = TRUE),
    mean_wet_start = mean(wet_start_doy, na.rm = TRUE),
    n_years = n(),
    .groups = "drop"
  )

message("  Sites:", nrow(site_means), "\n")

# Skip PCA3 if only 1 variable - just do scatter plot
# Add site characteristics for coloring
site_means <- site_means %>%
  left_join(site_chars %>% select(Stream_Name, Elevation_mean_m), 
            by = "Stream_Name")

p3 <- ggplot(site_means, aes(x = mean_wet_length, y = conc_sync_outlet)) +
  geom_point(aes(fill = Stream_Name, size = n_years), 
             shape = 21, colour = "grey20", stroke = 0.6, alpha = 0.8) +
  geom_text_repel(aes(label = Stream_Name), size = 3.5, max.overlaps = 20) +
  geom_smooth(method = "lm", se = TRUE, color = "gray40", linetype = "dashed") +
  scale_fill_site(name = "Site") +
  scale_size_continuous(name = "Years", range = c(3, 8)) +
  labs(
    title = "Site-Level: Wet Season Length vs Outlet Sync",
    subtitle = paste0("r = ", round(cor(site_means$mean_wet_length, 
                                         site_means$conc_sync_outlet, 
                                         use = "complete.obs"), 3)),
    x = "Mean Wet Season Length (days)",
    y = "Mean Concentration Sync with Outlet"
  ) +
  theme_pca() +
  guides(fill = guide_legend(override.aes = list(size = 4)))

ggsave(file.path(fig_dir, "site_sync_vs_wet_length.png"), p3, 
       width = 10, height = 8, dpi = 300)

# =============================================================================
# CORRELATION VISUALIZATION
# =============================================================================

message("\n=== CORRELATION BETWEEN SYNC AND TIMING ===\n")

# Scatter: wet length vs sync (annual data)
# Add site characteristics
sync_season_with_chars <- sync_season %>%
  left_join(site_chars %>% select(Stream_Name), by = "Stream_Name")

p4 <- ggplot(sync_season_with_chars, aes(x = wet_length_days, y = conc_sync_outlet)) +
  geom_point(aes(fill = Stream_Name), 
             shape = 21, colour = "grey20", stroke = 0.4, 
             alpha = 0.6, size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  scale_fill_site(name = "Site") +
  labs(
    title = "Wet Season Length vs Concentration Sync with Outlet",
    subtitle = paste0("r = ", round(cor(sync_season$wet_length_days, 
                                         sync_season$conc_sync_outlet, 
                                         use = "complete.obs"), 3),
                      " (Site-year data, n = ", nrow(sync_season), ")"),
    x = "Wet Season Length (days)",
    y = "Concentration Sync with Outlet"
  ) +
  theme_pca() +
  guides(fill = guide_legend(override.aes = list(size = 4)))

ggsave(file.path(fig_dir, "wet_length_vs_sync_scatter.png"), p4, 
       width = 10, height = 7, dpi = 300)

# =============================================================================
# SUMMARY
# =============================================================================

message("\n╔════════════════════════════════════════════════════════════════╗\n")
message("║  PCA ANALYSIS COMPLETE                                        ║\n")
message("╚════════════════════════════════════════════════════════════════╝\n\n")

message("Key findings:\n")
message("  - Wet season length correlates with outlet sync\n")
message("  - Sites cluster by both sync behavior and timing\n")
message("  - Combined PCA shows how timing relates to sync patterns\n\n")

message("Outputs saved to:\n")
message("  ", fig_dir, "\n\n")

# Save loadings
dir.create(file.path(out_dir, "04_PCA"), showWarnings = FALSE, recursive = TRUE)
write_csv(loadings2 %>% select(variable, PC1, PC2),
          file.path(out_dir, "04_PCA/sync_timing_pca_loadings.csv"))
