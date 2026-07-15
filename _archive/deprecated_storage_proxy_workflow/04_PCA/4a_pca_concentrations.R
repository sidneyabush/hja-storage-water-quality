# Restored from archive on 2025-12-04

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
})
try(source(file.path("/Users/sidneybush/Documents/GitHub/hja-water-quality", "00_helpers", "plot_theme_set.R")), silent = TRUE)

base_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality"
out_dir <- file.path(base_dir, "outputs")
fig_dir <- file.path(base_dir, "exploratory_plots", "04_PCA", "4a_concentrations")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
try(source(file.path(repo_dir, "00_helpers", "plot_prefs.R")), silent = TRUE)
try(source(file.path(repo_dir, "00_helpers", "workflow_config.R")), silent = TRUE)

if (!exists("flag_outlet_stream")) {
  flag_outlet_stream <- function(df) {
    df %>%
      mutate(
        is_outlet = grepl("GSLOOK|Lookout", Stream_Name, ignore.case = TRUE),
        outlet_marker = if_else(is_outlet, "Outlet", "Headwater")
      )
  }
}

cq_master <- read_csv(file.path(out_dir, "HJA_CQ_master.csv"), show_col_types = FALSE)
if ("variable" %in% names(cq_master) && !"solute" %in% names(cq_master)) {
  cq_master <- cq_master %>% rename(solute = variable)
}
if ("value" %in% names(cq_master) && !"Concentration" %in% names(cq_master)) {
  cq_master <- cq_master %>% rename(Concentration = value)
}
cq_master <- cq_master %>%
  filter(!is.na(Stream_Name) & Stream_Name != "",
         !is.na(solute) & solute != "")

conc_means <- cq_master %>%
  group_by(Stream_Name, solute) %>%
  summarize(
    mean_conc = mean(Concentration, na.rm = TRUE),
    sd_conc = sd(Concentration, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(mean_conc) & n_obs >= 10)

site_chars <- read_csv(file.path(out_dir, "HJA_clean_site_means.csv"), show_col_types = FALSE) %>%
  filter(!is.na(Stream_Name) & Stream_Name != "") %>%
  select(Stream_Name, any_of(c(
    "RBI", "FDC_slope_5_95", "WB_dS_range_mm", "Q_dS_range_mm",
    "Area_km2", "Elevation_mean_m", "Slope_mean", "MTT_final"
  ))) %>%
  distinct(Stream_Name, .keep_all = TRUE)

conc_wide <- conc_means %>%
  select(Stream_Name, solute, mean_conc) %>%
  pivot_wider(names_from = solute, values_from = mean_conc, names_prefix = "conc_")

site_names <- conc_wide$Stream_Name
conc_mat <- conc_wide %>% select(-Stream_Name) %>% as.matrix()
rownames(conc_mat) <- site_names

row_na_frac <- rowMeans(is.na(conc_mat))
keep_rows <- row_na_frac <= 0.5
conc_mat <- conc_mat[keep_rows, , drop = FALSE]
site_names <- site_names[keep_rows]

col_var <- apply(conc_mat, 2, function(x) var(x, na.rm = TRUE))
keep_cols <- !is.na(col_var) & col_var > 0
conc_mat <- conc_mat[, keep_cols, drop = FALSE]

if (nrow(conc_mat) >= 3 && ncol(conc_mat) >= 2) {
  cc <- complete.cases(conc_mat)
  if (sum(cc) < 3) {
    for (j in 1:ncol(conc_mat)) {
      conc_mat[is.na(conc_mat[, j]), j] <- median(conc_mat[, j], na.rm = TRUE)
    }
    cc <- rep(TRUE, nrow(conc_mat))
  }
  conc_mat_pca <- conc_mat[cc, , drop = FALSE]
  sites_pca <- site_names[cc]
  conc_scaled <- scale(conc_mat_pca, center = TRUE, scale = TRUE)
  pca_result <- prcomp(conc_scaled, center = FALSE, scale. = FALSE)
  var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100

  scores <- as_tibble(pca_result$x[, 1:min(3, ncol(pca_result$x))]) %>%
    mutate(Stream_Name = sites_pca) %>%
    left_join(site_chars, by = "Stream_Name") %>%
    flag_outlet_stream()
  scores <- apply_factor_orders(scores)

  loadings <- as.data.frame(pca_result$rotation[, 1:min(3, ncol(pca_result$rotation))]) %>%
    rownames_to_column("solute") %>%
    as_tibble() %>%
    mutate(solute = gsub("conc_", "", solute))

  loading_scale <- max(abs(scores$PC1), abs(scores$PC2)) / max(abs(loadings$PC1), abs(loadings$PC2)) * 0.8
  loadings_scaled <- loadings %>% mutate(PC1_plot = PC1 * loading_scale, PC2_plot = PC2 * loading_scale)

  storage_candidates <- c("Q_dS_range_mm", "WB_dS_range_mm")
  storage_available <- purrr::keep(storage_candidates, ~ .x %in% names(scores) && any(is.finite(scores[[.x]])))
  if (length(storage_available) == 0) storage_available <- NA_character_

  build_biplot <- function(storage_col) {
    has_storage <- !is.na(storage_col)
    storage_label <- if (has_storage) get_storage_label(storage_col) else "Storage metric unavailable"
    point_mapping <- if (has_storage) {
      aes(x = PC1, y = PC2, fill = Stream_Name, size = .data[[storage_col]])
    } else {
      aes(x = PC1, y = PC2, fill = Stream_Name)
    }
    p <- ggplot() +
      stat_ellipse(data = scores, aes(x = PC1, y = PC2), level = 0.95, type = "t", linetype = "dashed", color = "gray50", linewidth = 0.5) +
      geom_point(data = scores, mapping = point_mapping, shape = 21, colour = "grey20", stroke = 0.6, alpha = 0.9) +
      geom_text_repel(data = scores, aes(x = PC1, y = PC2, label = Stream_Name), size = 3, max.overlaps = 20, segment.alpha = 0.3) +
      geom_segment(data = loadings_scaled, aes(x = 0, y = 0, xend = PC1_plot, yend = PC2_plot), arrow = arrow(length = unit(0.02, "npc")), color = "gray40", linewidth = 0.8) +
      geom_text_repel(data = loadings_scaled, aes(x = PC1_plot, y = PC2_plot, label = solute), size = 3.5, color = "gray30", fontface = "bold", max.overlaps = 20) +
      scale_fill_site(name = "Site") +
      labs(title = storage_label, x = paste0("PC1 (", round(var_explained[1], 1), "%)"), y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
      theme_bw(base_size = 12) + theme(legend.position = "right") + coord_fixed() +
      guides(fill = guide_legend(override.aes = list(shape = 21, colour = "grey20")))
    if (has_storage) {
      p <- p + scale_size_continuous(name = storage_label, range = c(2.5, 8), breaks = scales::pretty_breaks(n = 4))
    } else {
      p <- p + guides(size = "none")
    }
    p
  }

  plot_list <- purrr::map(storage_available, build_biplot)
  combined_plot <- if (length(plot_list) == 1) plot_list[[1]] else patchwork::wrap_plots(plotlist = plot_list, nrow = 1)
  combined_plot <- combined_plot + patchwork::plot_annotation(title = "PCA: Geochemical Fingerprint (Mean Concentrations)", subtitle = "Panels distinguish dynamic storage metrics (Q-dS vs WB-dS)", caption = "Dashed ellipse: 95% confidence interval; points outlined for visibility")
  panel_width <- if (length(plot_list) > 1) 14 else 11
  ggsave(file.path(fig_dir, "01_pca_concentrations_biplot.png"), combined_plot, width = panel_width, height = 9, dpi = 300)

  scree_df <- tibble(PC = 1:length(var_explained), Variance = var_explained, Cumulative = cumsum(var_explained))
  p_scree <- ggplot(scree_df, aes(x = PC)) +
    geom_col(aes(y = Variance), fill = "steelblue", alpha = 0.7) +
    geom_line(aes(y = Cumulative), color = "red", linewidth = 1) +
    geom_point(aes(y = Cumulative), color = "red", size = 2) +
    geom_hline(yintercept = 80, linetype = "dashed", color = "gray50") +
    labs(title = "Scree Plot: Concentration PCA", x = "Principal Component", y = "Variance Explained (%)") +
    theme_bw(base_size = 12)
  ggsave(file.path(fig_dir, "03_scree_plot.png"), p_scree, width = 8, height = 6, dpi = 300)

  write_csv(scores, file.path(fig_dir, "pca_concentration_scores.csv"))
  write_csv(loadings, file.path(fig_dir, "pca_concentration_loadings.csv"))
  write_csv(scree_df, file.path(fig_dir, "pca_concentration_variance.csv"))
}
# =============================================================================
# 4a_pca_concentrations.R
# =============================================================================
# PCA ON MEAN SOLUTE CONCENTRATIONS - "GEOCHEMICAL FINGERPRINT"
#
# RESEARCH QUESTION:
#   How do sites cluster based on their overall geochemical fingerprint?
#   Do sites with similar geology/weathering have similar chemistry?
#
# APPROACH:
#   - Create wide matrix: sites × mean concentrations for each solute
#   - Standardize (z-score) to give equal weight to all solutes
#   - Run PCA and visualize biplot with 95% confidence ellipses
#   - Color by catchment characteristics (storage, geology, etc.)
#
# OUTPUT:
#   - Biplot with site scores and solute loadings
#   - 95% confidence ellipses showing spread
#   - Storage-colored version
#   - Scree plot
#
# =============================================================================

# rm(list = ls())  # Commented out when run from master script

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
fig_dir <- file.path(base_dir, "exploratory_plots", "04_PCA", "4a_concentrations")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Source helpers
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
try(source(file.path(repo_dir, "00_helpers", "plot_prefs.R")), silent = TRUE)
try(source(file.path(repo_dir, "00_helpers", "workflow_config.R")), silent = TRUE)

# Define outlet stream for shape coding
if (!exists("flag_outlet_stream")) {
  flag_outlet_stream <- function(df) {
    df %>%
      mutate(
        is_outlet = grepl("GSLOOK|Lookout", Stream_Name, ignore.case = TRUE),
        outlet_marker = if_else(is_outlet, "Outlet", "Headwater")
      )
  }
}

# =============================================================================
# LOAD DATA
# =============================================================================

message("Loading CQ master data for concentrations...")

cq_master <- read_csv(file.path(out_dir, "HJA_CQ_master.csv"), show_col_types = FALSE)

# Identify the correct column names
# CQ_master uses 'variable' for solute and 'value' for concentration
if ("variable" %in% names(cq_master) && !"solute" %in% names(cq_master)) {
  cq_master <- cq_master %>% rename(solute = variable)
}
if ("value" %in% names(cq_master) && !"Concentration" %in% names(cq_master)) {
  cq_master <- cq_master %>% rename(Concentration = value)
}

# Filter out NA Stream_Name and solute
cq_master <- cq_master %>%
  filter(!is.na(Stream_Name) & Stream_Name != "",
         !is.na(solute) & solute != "")

# Calculate mean concentration per site × solute
conc_means <- cq_master %>%
  group_by(Stream_Name, solute) %>%
  summarize(
    mean_conc = mean(Concentration, na.rm = TRUE),
    sd_conc = sd(Concentration, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(mean_conc) & n_obs >= 10)  # Require at least 10 observations

# Get catchment characteristics for coloring
site_chars <- read_csv(file.path(out_dir, "HJA_clean_site_means.csv"), show_col_types = FALSE) %>%
  filter(!is.na(Stream_Name) & Stream_Name != "") %>%
  select(Stream_Name, any_of(c(
    "RBI", "FDC_slope_5_95", "WB_dS_range_mm", "Q_dS_range_mm",
    "Area_km2", "Elevation_mean_m", "Slope_mean", "MTT_final"
  ))) %>%
  distinct(Stream_Name, .keep_all = TRUE)

# =============================================================================
# CREATE WIDE MATRIX FOR PCA
# =============================================================================

message("Creating concentration matrix...")

conc_wide <- conc_means %>%
  select(Stream_Name, solute, mean_conc) %>%
  pivot_wider(
    names_from = solute,
    values_from = mean_conc,
    names_prefix = "conc_"
  )

# Store site names
site_names <- conc_wide$Stream_Name

# Extract numeric matrix
conc_mat <- conc_wide %>% 
  select(-Stream_Name) %>%
  as.matrix()
rownames(conc_mat) <- site_names

# QC: Remove sites with >50% missing solutes
row_na_frac <- rowMeans(is.na(conc_mat))
keep_rows <- row_na_frac <= 0.5
conc_mat <- conc_mat[keep_rows, , drop = FALSE]
site_names <- site_names[keep_rows]

# Remove solutes with no variance or all NA
col_var <- apply(conc_mat, 2, function(x) var(x, na.rm = TRUE))
keep_cols <- !is.na(col_var) & col_var > 0
conc_mat <- conc_mat[, keep_cols, drop = FALSE]

message("Matrix dimensions: ", nrow(conc_mat), " sites × ", ncol(conc_mat), " solutes")

# =============================================================================
# RUN PCA
# =============================================================================

if (nrow(conc_mat) >= 3 && ncol(conc_mat) >= 2) {
  
  # Handle missing values: use complete cases or impute
  cc <- complete.cases(conc_mat)
  if (sum(cc) < 3) {
    message("Not enough complete cases. Imputing missing values with column medians...")
    for (j in 1:ncol(conc_mat)) {
      conc_mat[is.na(conc_mat[, j]), j] <- median(conc_mat[, j], na.rm = TRUE)
    }
    cc <- rep(TRUE, nrow(conc_mat))
  }
  
  conc_mat_pca <- conc_mat[cc, , drop = FALSE]
  sites_pca <- site_names[cc]
  
  # Scale (z-score standardization)
  conc_scaled <- scale(conc_mat_pca, center = TRUE, scale = TRUE)
  
  # Run PCA
  pca_result <- prcomp(conc_scaled, center = FALSE, scale. = FALSE)
  
  # Variance explained
  var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100
  
  message("\nVariance explained:")
  message("  PC1: ", round(var_explained[1], 1), "%")
  message("  PC2: ", round(var_explained[2], 1), "%")
  if (length(var_explained) >= 3) message("  PC3: ", round(var_explained[3], 1), "%")
  
  # ==========================================================================
  # EXTRACT SCORES AND LOADINGS
  # ==========================================================================
  
  scores <- as_tibble(pca_result$x[, 1:min(3, ncol(pca_result$x))]) %>%
    mutate(Stream_Name = sites_pca) %>%
    left_join(site_chars, by = "Stream_Name") %>%
    flag_outlet_stream()
  scores <- apply_factor_orders(scores)
  
  loadings <- as.data.frame(pca_result$rotation[, 1:min(3, ncol(pca_result$rotation))]) %>%
    rownames_to_column("solute") %>%
    as_tibble() %>%
    mutate(solute = gsub("conc_", "", solute))
  
  # Scale loadings for biplot visualization
  loading_scale <- max(abs(scores$PC1), abs(scores$PC2)) / max(abs(loadings$PC1), abs(loadings$PC2)) * 0.8
  loadings_scaled <- loadings %>%
    mutate(
      PC1_plot = PC1 * loading_scale,
      PC2_plot = PC2 * loading_scale
    )
  
  # ==========================================================================
  # BIPLOT: SITES + SOLUTE LOADINGS WITH CONFIDENCE ELLIPSE
  # Color by site, size by dynamic storage
  # ==========================================================================
  
  storage_candidates <- c("Q_dS_range_mm", "WB_dS_range_mm")
  storage_available <- purrr::keep(
    storage_candidates,
    ~ .x %in% names(scores) && any(is.finite(scores[[.x]]))
  )
  if (length(storage_available) == 0) {
    storage_available <- NA_character_
  }

  build_biplot <- function(storage_col) {
    has_storage <- !is.na(storage_col)
    storage_label <- if (has_storage) get_storage_label(storage_col) else "Storage metric unavailable"
    point_mapping <- if (has_storage) {
      aes(x = PC1, y = PC2, fill = Stream_Name, size = .data[[storage_col]])
    } else {
      aes(x = PC1, y = PC2, fill = Stream_Name)
    }

    p <- ggplot() +
      stat_ellipse(
        data = scores,
        aes(x = PC1, y = PC2),
        level = 0.95,
        type = "t",
        linetype = "dashed",
        color = "gray50",
        linewidth = 0.5
      ) +
      geom_point(
        data = scores,
        mapping = point_mapping,
        shape = 21,
        colour = "grey20",
        stroke = 0.6,
        alpha = 0.9
      ) +
      geom_text_repel(
        data = scores,
        aes(x = PC1, y = PC2, label = Stream_Name),
        size = 3,
        max.overlaps = 20,
        segment.alpha = 0.3
      ) +
      geom_segment(
        data = loadings_scaled,
        aes(x = 0, y = 0, xend = PC1_plot, yend = PC2_plot),
        arrow = arrow(length = unit(0.02, "npc")),
        color = "gray40",
        linewidth = 0.8
      ) +
      geom_text_repel(
        data = loadings_scaled,
        aes(x = PC1_plot, y = PC2_plot, label = solute),
        size = 3.5,
        color = "gray30",
        fontface = "bold",
        max.overlaps = 20
      ) +
      scale_fill_site(name = "Site") +
      labs(
        title = storage_label,
        x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
        y = paste0("PC2 (", round(var_explained[2], 1), "%)")
      ) +
      theme_bw(base_size = 12) +
      theme(legend.position = "right") +
      coord_fixed() +
      guides(
        fill = guide_legend(override.aes = list(shape = 21, colour = "grey20"))
      )

    if (has_storage) {
      p <- p + scale_size_continuous(
        name = storage_label,
        range = c(2.5, 8),
        breaks = scales::pretty_breaks(n = 4)
      )
    } else {
      p <- p + guides(size = "none")
    }

    p
  }

  plot_list <- purrr::map(storage_available, build_biplot)
  combined_plot <- if (length(plot_list) == 1) {
    plot_list[[1]]
  } else {
    patchwork::wrap_plots(plotlist = plot_list, nrow = 1)
  }

  combined_plot <- combined_plot + patchwork::plot_annotation(
    title = "PCA: Geochemical Fingerprint (Mean Concentrations)",
    subtitle = "Panels distinguish dynamic storage metrics (Q-dS vs WB-dS)",
    caption = "Dashed ellipse: 95% confidence interval; points outlined for visibility"
  )

  panel_width <- if (length(plot_list) > 1) 14 else 11
  ggsave(
    file.path(fig_dir, "01_pca_concentrations_biplot.png"),
    combined_plot,
    width = panel_width,
    height = 9,
    dpi = 300
  )
  
  # ==========================================================================
  # BIPLOT WITH SITE ELLIPSES (individual ellipses not possible with n=1 per site)
  # Alternative: show convex hulls or skip if single point per site
  # ==========================================================================
  
  # Note: With only 1 observation per site, we can't draw per-site ellipses
  # The overall ellipse in the biplot above serves this purpose
  
  # ==========================================================================
  # SCREE PLOT
  # ==========================================================================
  
  scree_df <- tibble(
    PC = 1:length(var_explained),
    Variance = var_explained,
    Cumulative = cumsum(var_explained)
  )
  
  p_scree <- ggplot(scree_df, aes(x = PC)) +
    geom_col(aes(y = Variance), fill = "steelblue", alpha = 0.7) +
    geom_line(aes(y = Cumulative), color = "red", linewidth = 1) +
    geom_point(aes(y = Cumulative), color = "red", size = 2) +
    geom_hline(yintercept = 80, linetype = "dashed", color = "gray50") +
    labs(
      title = "Scree Plot: Concentration PCA",
      x = "Principal Component",
      y = "Variance Explained (%)"
    ) +
    theme_bw(base_size = 12)
  
  ggsave(file.path(fig_dir, "03_scree_plot.png"), p_scree,
         width = 8, height = 6, dpi = 300)
  
  # ==========================================================================
  # SAVE RESULTS
  # ==========================================================================
  
  write_csv(scores, file.path(fig_dir, "pca_concentration_scores.csv"))
  write_csv(loadings, file.path(fig_dir, "pca_concentration_loadings.csv"))
  write_csv(scree_df, file.path(fig_dir, "pca_concentration_variance.csv"))
  
  message("\n✓ Concentration PCA complete!")
  message("  Figures saved to: ", fig_dir)
  
} else {
  message("Insufficient data for PCA (need ≥3 sites and ≥2 solutes)")
}
