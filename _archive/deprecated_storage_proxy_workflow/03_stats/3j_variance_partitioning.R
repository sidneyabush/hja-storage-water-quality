# =============================================================================
# STEP 03j: VARIANCE PARTITIONING ANALYSIS
# =============================================================================
# Purpose: Decompose synchrony variance into contributions from:
#   1. Storage dynamics (Q_dS, WB_dS, RBI, flashiness)
#   2. Geology (Lava%, Ash%, Pyroclastic%)
#   3. Topography (Elevation, Slope, Area)
#   4. Spatial proximity (distance between sites)
#
# Method: Variation partitioning (Borcard et al. 1992) via vegan::varpart()
#
# This analysis quantifies:
#   - How much variance in synchrony is explained by each factor?
#   - What are the unique vs. shared contributions?
#   - What is the relative importance of storage vs. other drivers?
#
# Research Questions Addressed:
#   RQ1: Does dynamic storage explain shared chemical behavior?
#   RQ2: Does storage variability weaken synchrony?
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)      # For varpart()
  library(geosphere)  # For distance calculations
  library(patchwork)
})

rm(list = ls())

# Source helpers
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
try(source(file.path(repo_dir, "00_helpers", "plot_prefs.R")), silent = TRUE)
try(source(file.path(repo_dir, "00_helpers", "workflow_config.R")), silent = TRUE)

# Paths
paths <- get_project_paths()
out_dir <- paths$out_dir
fig_dir <- file.path(paths$fig_root, "03_stats", "3j_variance_partition")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "03_stats"), showWarnings = FALSE, recursive = TRUE)

message("\n╔════════════════════════════════════════════════════════════════╗\n")
message("║  VARIANCE PARTITIONING: Storage vs. Geology vs. Topography   ║\n")
message("╚════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# LOAD DATA
# =============================================================================
message("=== LOADING DATA ===\n")

# Site-level data with all metrics
site_means <- read_csv(file.path(out_dir, "HJA_master_site_means.csv"),
                       show_col_types = FALSE) %>%
  mutate(
    Stream_Name = factor(Stream_Name, levels = site_order),
    solute = factor(solute, levels = solute_order)
  )

# Catchment characteristics (for coordinates)
catchment <- read_csv(file.path(out_dir, "Catchment_site_characteristics.csv"),
                     show_col_types = FALSE)

# Composite synchrony
composite_sync <- read_csv(file.path(out_dir, "HJA_composite_synchrony.csv"),
                           show_col_types = FALSE) %>%
  mutate(
    Stream_Name = factor(Stream_Name, levels = site_order),
    solute = factor(solute, levels = solute_order)
  )

# Merge
data_full <- site_means %>%
  left_join(composite_sync, by = c("Stream_Name", "solute"))

message("  Data:", nrow(data_full), "site-solute combinations\n")

# Identify available columns
storage_metric <- intersect(PRIMARY_STORAGE_METRIC, names(data_full))
if (length(storage_metric) == 0) {
  storage_metric <- intersect(c("Q_dS_range_mm", "WB_dS_range_mm"), names(data_full))
}
storage_metric <- storage_metric[[1]]
message("  Primary storage metric:", storage_metric, "\n")

# =============================================================================
# PREPARE VARIABLE GROUPS
# =============================================================================
message("\n=== PREPARING VARIABLE GROUPS ===\n")

# Storage dynamics variables (including DS_drawdown)
storage_vars <- c(storage_metric, "RBI", "FDC_slope_5_95", "Q05_mm_d",
                  "DS_drawdown_mean_mm", "DS_drawdown_sd_mm", "DS_drawdown_range_mm")
storage_vars <- intersect(storage_vars, names(data_full))
message("  Storage variables (", length(storage_vars), "):", paste(storage_vars, collapse = ", "), "\n")

# Geology variables
geology_vars <- c("Lava1_per", "Lava2_per", "Ash_Per", "Pyro_per")
geology_vars <- intersect(geology_vars, names(data_full))
message("  Geology variables (", length(geology_vars), "):", paste(geology_vars, collapse = ", "), "\n")

# Topography variables
topo_vars <- c("Area_km2", "Elevation_mean_m", "Slope_mean")
topo_vars <- intersect(topo_vars, names(data_full))
message("  Topography variables (", length(topo_vars), "):", paste(topo_vars, collapse = ", "), "\n")

# =============================================================================
# COMPUTE SPATIAL DISTANCE MATRIX
# =============================================================================
message("\n=== COMPUTING SPATIAL DISTANCES ===\n")

# Get site coordinates (if available)
if (all(c("Latitude", "Longitude") %in% names(catchment))) {

  site_coords <- catchment %>%
    select(Stream_Name, Longitude, Latitude) %>%
    filter(Stream_Name %in% levels(data_full$Stream_Name))

  # Compute pairwise distances (km)
  coords_matrix <- as.matrix(site_coords[, c("Longitude", "Latitude")])
  rownames(coords_matrix) <- site_coords$Stream_Name

  dist_matrix <- distm(coords_matrix, fun = distHaversine) / 1000  # Convert to km
  rownames(dist_matrix) <- site_coords$Stream_Name
  colnames(dist_matrix) <- site_coords$Stream_Name

  # Create distance variable for each observation (mean distance to all other sites)
  data_full <- data_full %>%
    rowwise() %>%
    mutate(
      mean_distance = if (as.character(Stream_Name) %in% rownames(dist_matrix)) {
        mean(dist_matrix[as.character(Stream_Name), ], na.rm = TRUE)
      } else {
        NA_real_
      }
    ) %>%
    ungroup()

  spatial_vars <- "mean_distance"
  message("  Spatial distance computed for", n_distinct(data_full$Stream_Name), "sites\n")

} else {
  message("  Warning: No coordinates available; skipping spatial distance\n")
  spatial_vars <- character(0)
}

# =============================================================================
# VARIANCE PARTITIONING: 3-WAY
# =============================================================================
message("\n=== 3-WAY VARIANCE PARTITIONING ===\n")

# Response variables (synchrony metrics)
response_vars <- c("conc_sync_allpairs", "cqslope_sync_allpairs",
                   "wymore_crosssite_allpairs")

# Check for .y suffix (in case of duplicate columns from join)
if (any(grepl("\\.y$", names(data_full)))) {
  response_vars <- paste0(response_vars, ".y")
}

response_vars <- intersect(response_vars, names(data_full))

message("  Response variables:", paste(response_vars, collapse = ", "), "\n")

# Run variance partitioning for each response
varpart_results <- list()

for (resp in response_vars) {

  message("\n--- Response:", resp, "---\n")

  # Filter complete cases
  data_vp <- data_full %>%
    select(all_of(c(resp, storage_vars, geology_vars, topo_vars))) %>%
    drop_na()

  if (nrow(data_vp) < 10) {
    message("  Insufficient data (n =", nrow(data_vp), "); skipping\n")
    next
  }

  message("  N =", nrow(data_vp), "complete observations\n")

  # Prepare matrices
  Y <- data_vp[[resp]]
  X_storage <- as.matrix(data_vp[, storage_vars, drop = FALSE])
  X_geology <- as.matrix(data_vp[, geology_vars, drop = FALSE])
  X_topo <- as.matrix(data_vp[, topo_vars, drop = FALSE])

  # Check for sufficient variation
  if (ncol(X_storage) < 1 || ncol(X_geology) < 1 || ncol(X_topo) < 1) {
    message("  Insufficient variable groups; skipping\n")
    next
  }

  # Variance partitioning (3-way)
  vp <- varpart(Y, X_storage, X_geology, X_topo)

  message("\n  Variance partitioning results:\n")
  print(vp)

  # Extract fractions (individual components)
  fractions <- as.data.frame(vp$part$indfract)
  fractions$fraction <- rownames(fractions)
  fractions$response <- resp

  # Also extract total R² from the full model
  total_r2 <- vp$part$fract["[a+b+c+d+e+f+g] = All", "Adj.R.square"]
  fractions_with_total <- rbind(
    fractions,
    data.frame(
      Df = NA,
      R.square = NA,
      Adj.R.square = total_r2,
      Testable = NA,
      fraction = "[a+b+c]",
      response = resp
    )
  )

  varpart_results[[resp]] <- fractions_with_total

  # Plot variance partitioning
  png(file.path(fig_dir, paste0("varpart_", resp, ".png")),
      width = 10, height = 8, units = "in", res = 300)
  plot(vp, digits = 2,
       Xnames = c("Storage", "Geology", "Topography"),
       bg = c(cluster_colors[["1"]], cluster_colors[["2"]], cluster_colors[["3"]]))
  title(main = paste("Variance Partitioning:", get_sync_label(resp)),
        cex.main = 1.2, font.main = 2)
  dev.off()

  message("  Plot saved:", paste0("varpart_", resp, ".png"), "\n")
}

# Combine results
if (length(varpart_results) > 0) {
  varpart_combined <- bind_rows(varpart_results)
  write_csv(varpart_combined, file.path(out_dir, "03_stats/variance_partitioning_results.csv"))
}

# =============================================================================
# SUMMARY TABLE
# =============================================================================
message("\n=== CREATING SUMMARY TABLE ===\n")

if (length(varpart_results) > 0) {

  # Extract unique variance explained by each factor
  summary_table <- varpart_combined %>%
    filter(fraction %in% c("[a]", "[b]", "[c]", "[a+b+c]")) %>%
    mutate(
      component = case_when(
        fraction == "[a]" ~ "Storage (unique)",
        fraction == "[b]" ~ "Geology (unique)",
        fraction == "[c]" ~ "Topography (unique)",
        fraction == "[a+b+c]" ~ "Total explained"
      )
    ) %>%
    select(response, component, Adj.R.square) %>%
    pivot_wider(names_from = component, values_from = Adj.R.square) %>%
    mutate(
      Unexplained = 1 - `Total explained`,
      Storage_pct = `Storage (unique)` / `Total explained` * 100,
      Geology_pct = `Geology (unique)` / `Total explained` * 100,
      Topography_pct = `Topography (unique)` / `Total explained` * 100
    )

  message("\nVariance partitioning summary:\n")
  print(summary_table, width = Inf)

  write_csv(summary_table, file.path(out_dir, "03_stats/variance_partitioning_summary.csv"))

  # Visualization: Bar plot of unique contributions
  plot_data <- summary_table %>%
    select(response, Storage_pct, Geology_pct, Topography_pct) %>%
    pivot_longer(cols = c(Storage_pct, Geology_pct, Topography_pct),
                 names_to = "component", values_to = "percent") %>%
    mutate(
      component = str_remove(component, "_pct"),
      component = factor(component, levels = c("Storage", "Geology", "Topography"))
    )

  p_varpart <- ggplot(plot_data, aes(x = response, y = percent, fill = component)) +
    geom_col(position = "stack", color = "white", linewidth = 0.5) +
    scale_fill_manual(
      values = c("Storage" = cluster_colors[["1"]],
                 "Geology" = cluster_colors[["2"]],
                 "Topography" = cluster_colors[["3"]]),
      name = "Component"
    ) +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    labs(
      x = "Synchrony Metric",
      y = "% of Explained Variance",
      title = "Variance Partitioning: Relative Importance of Drivers",
      subtitle = "Unique contributions of Storage, Geology, and Topography to synchrony",
      caption = "Percentages are relative to total explained variance (excludes shared components)"
    ) +
    theme_hja() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    legend_bottom()

  ggsave(file.path(fig_dir, "variance_partition_summary.png"), p_varpart,
         width = 11, height = 8, dpi = 300)
}

# =============================================================================
# 4-WAY PARTITIONING (if spatial data available)
# =============================================================================

if (length(spatial_vars) > 0 && "mean_distance" %in% names(data_full)) {

  message("\n=== 4-WAY VARIANCE PARTITIONING (including Spatial) ===\n")

  for (resp in response_vars) {

    message("\n--- Response:", resp, "---\n")

    # Filter complete cases
    data_vp4 <- data_full %>%
      select(all_of(c(resp, storage_vars, geology_vars, topo_vars, spatial_vars))) %>%
      drop_na()

    if (nrow(data_vp4) < 10) {
      message("  Insufficient data; skipping\n")
      next
    }

    message("  N =", nrow(data_vp4), "complete observations\n")

    # Prepare matrices
    Y <- data_vp4[[resp]]
    X_storage <- as.matrix(data_vp4[, storage_vars, drop = FALSE])
    X_geology <- as.matrix(data_vp4[, geology_vars, drop = FALSE])
    X_topo <- as.matrix(data_vp4[, topo_vars, drop = FALSE])
    X_spatial <- as.matrix(data_vp4[, spatial_vars, drop = FALSE])

    # 4-way variance partitioning
    vp4 <- varpart(Y, X_storage, X_geology, X_topo, X_spatial)

    message("\n  4-way variance partitioning:\n")
    print(vp4)

    # Plot
    png(file.path(fig_dir, paste0("varpart_4way_", resp, ".png")),
        width = 12, height = 10, units = "in", res = 300)
    plot(vp4, digits = 2,
         Xnames = c("Storage", "Geology", "Topography", "Distance"),
         bg = c(cluster_colors[["1"]], cluster_colors[["2"]],
                cluster_colors[["3"]], cluster_colors[["4"]]))
    title(main = paste("4-Way Variance Partitioning:", get_sync_label(resp)),
          cex.main = 1.2, font.main = 2)
    dev.off()

    message("  4-way plot saved\n")
  }
}

# =============================================================================
# INTERPRETATION GUIDE
# =============================================================================
message("\n╔════════════════════════════════════════════════════════════════╗\n")
message("║  INTERPRETATION GUIDE                                         ║\n")
message("╚════════════════════════════════════════════════════════════════╝\n\n")

message("Variance partitioning components:\n\n")

message("[a] = Unique variance explained by Storage\n")
message("  → How much synchrony variance is ONLY explained by storage?\n\n")

message("[b] = Unique variance explained by Geology\n")
message("  → How much is ONLY explained by geology?\n\n")

message("[c] = Unique variance explained by Topography\n")
message("  → How much is ONLY explained by topography?\n\n")

message("[a+b] = Shared variance (Storage + Geology)\n")
message("  → Variance that could be attributed to either (confounded)\n\n")

message("[a+b+c] = Total explained variance\n")
message("  → R² for the full model with all predictors\n\n")

message("Residual = Unexplained variance\n")
message("  → Variation not explained by any of the measured factors\n\n")

message("Relative importance:\n")
message("  Storage % = [a] / [a+b+c] × 100\n")
message("  → What % of the explained variance is uniquely due to storage?\n\n")

# =============================================================================
# SUMMARY
# =============================================================================
message("\n╔════════════════════════════════════════════════════════════════╗\n")
message("║  ANALYSIS COMPLETE                                            ║\n")
message("╚════════════════════════════════════════════════════════════════╝\n\n")

if (exists("summary_table")) {
  message("KEY FINDINGS:\n\n")

  for (i in 1:nrow(summary_table)) {
    resp <- summary_table$response[i]
    storage_pct <- round(summary_table$Storage_pct[i], 1)
    geology_pct <- round(summary_table$Geology_pct[i], 1)
    topo_pct <- round(summary_table$Topography_pct[i], 1)
    total_r2 <- round(summary_table$`Total explained`[i], 3)

    message(resp, ":\n")
    message("  Total R² =", total_r2, "\n")
    message("  Storage:", storage_pct, "% of explained variance\n")
    message("  Geology:", geology_pct, "%\n")
    message("  Topography:", topo_pct, "%\n\n")
  }
}

message("Outputs saved to:\n")
message("  Tables:", file.path(out_dir, "03_stats/"), "\n")
message("  Figures:", fig_dir, "\n\n")

message("Key files:\n")
message("  - variance_partitioning_results.csv (full results)\n")
message("  - variance_partitioning_summary.csv (summary table)\n")
message("  - varpart_*.png (Venn diagrams)\n\n")
