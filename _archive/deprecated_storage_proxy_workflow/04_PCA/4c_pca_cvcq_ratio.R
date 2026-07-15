# Restored from archive on 2025-12-04

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
})
try(source(file.path("/Users/sidneybush/Documents/GitHub/hja-water-quality", "00_helpers", "plot_theme_set.R")), silent = TRUE)

base_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality"
out_dir <- file.path(base_dir, "outputs")
fig_dir <- file.path(base_dir, "exploratory_plots", "04_PCA", "4c_cvcq_ratio")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
try(source(file.path(repo_dir, "00_helpers", "plot_prefs.R")), silent = TRUE)
try(source(file.path(repo_dir, "00_helpers", "workflow_config.R")), silent = TRUE)

if (!exists("flag_outlet_stream")) {
  flag_outlet_stream <- function(df) {
    df %>% mutate(
      is_outlet = grepl("GSLOOK|Lookout", Stream_Name, ignore.case = TRUE),
      outlet_marker = if_else(is_outlet, "Outlet", "Headwater")
    )
  }
}

site_means <- read_csv(file.path(out_dir, "HJA_clean_site_means.csv"), show_col_types = FALSE)
if (!"cq_CVc_CVq" %in% names(site_means)) stop("cq_CVc_CVq column not found in HJA_clean_site_means.csv")

cvcq_wide <- site_means %>%
  select(Stream_Name, solute, cq_CVc_CVq) %>%
  filter(!is.na(cq_CVc_CVq) & is.finite(cq_CVc_CVq)) %>%
  group_by(Stream_Name, solute) %>%
  summarize(cq_CVc_CVq = mean(cq_CVc_CVq, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = solute, values_from = cq_CVc_CVq, names_prefix = "cvcq_")

site_names <- cvcq_wide$Stream_Name
cvcq_mat <- cvcq_wide %>% select(-Stream_Name) %>% as.matrix()
rownames(cvcq_mat) <- site_names

row_na_frac <- rowMeans(is.na(cvcq_mat))
keep_rows <- row_na_frac <= 0.5
cvcq_mat <- cvcq_mat[keep_rows, , drop = FALSE]
site_names <- site_names[keep_rows]

col_var <- apply(cvcq_mat, 2, function(x) var(x, na.rm = TRUE))
keep_cols <- !is.na(col_var) & col_var > 0
cvcq_mat <- cvcq_mat[, keep_cols, drop = FALSE]

site_chars <- site_means %>%
  select(Stream_Name, any_of(c("RBI", "FDC_slope_5_95", "WB_dS_range_mm", "Q_dS_range_mm", "Area_km2", "Elevation_mean_m", "Slope_mean", "MTT_final"))) %>%
  distinct(Stream_Name, .keep_all = TRUE)

if (nrow(cvcq_mat) >= 3 && ncol(cvcq_mat) >= 2) {
  cvcq_mat_log <- log10(cvcq_mat + 0.01)
  cc <- complete.cases(cvcq_mat_log)
  if (sum(cc) < 3) {
    for (j in 1:ncol(cvcq_mat_log)) cvcq_mat_log[is.na(cvcq_mat_log[, j]), j] <- median(cvcq_mat_log[, j], na.rm = TRUE)
    cc <- rep(TRUE, nrow(cvcq_mat_log))
  }
  cvcq_mat_pca <- cvcq_mat_log[cc, , drop = FALSE]
  sites_pca <- site_names[cc]
  cvcq_scaled <- scale(cvcq_mat_pca, center = TRUE, scale = TRUE)
  pca_result <- prcomp(cvcq_scaled, center = FALSE, scale. = FALSE)
  var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100

  scores <- as_tibble(pca_result$x[, 1:min(3, ncol(pca_result$x))]) %>%
    mutate(Stream_Name = sites_pca) %>%
    left_join(site_chars, by = "Stream_Name") %>%
    flag_outlet_stream()
  scores <- apply_factor_orders(scores)

  loadings <- as.data.frame(pca_result$rotation[, 1:min(3, ncol(pca_result$rotation))]) %>%
    rownames_to_column("solute") %>% as_tibble() %>% mutate(solute = gsub("cvcq_", "", solute))

  loading_scale <- max(abs(scores$PC1), abs(scores$PC2)) / max(abs(loadings$PC1), abs(loadings$PC2)) * 0.8
  loadings_scaled <- loadings %>% mutate(PC1_plot = PC1 * loading_scale, PC2_plot = PC2 * loading_scale)

  storage_candidates <- c("Q_dS_range_mm", "WB_dS_range_mm")
  storage_available <- purrr::keep(storage_candidates, ~ .x %in% names(scores) && any(is.finite(scores[[.x]])))
  if (length(storage_available) == 0) storage_available <- NA_character_

  build_biplot <- function(storage_col) {
    has_storage <- !is.na(storage_col)
    storage_label <- if (has_storage) get_storage_label(storage_col) else "Storage metric unavailable"
    point_mapping <- if (has_storage) aes(x = PC1, y = PC2, fill = Stream_Name, size = .data[[storage_col]]) else aes(x = PC1, y = PC2, fill = Stream_Name)
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
    if (has_storage) p <- p + scale_size_continuous(name = storage_label, range = c(2.5, 8), breaks = scales::pretty_breaks(n = 4)) else p <- p + guides(size = "none")
    p
  }

  plot_list <- purrr::map(storage_available, build_biplot)
  combined_plot <- if (length(plot_list) == 1) plot_list[[1]] else patchwork::wrap_plots(plotlist = plot_list, nrow = 1)
  combined_plot <- combined_plot + patchwork::plot_annotation(title = "PCA: Buffering Fingerprint (CVc/CVq Ratio)", subtitle = "Panels distinguish dynamic storage metrics (Q-dS vs WB-dS)", caption = "Dashed ellipse: 95% CI; positive loadings = more chemodynamic; points outlined for visibility")
  panel_width <- if (length(plot_list) > 1) 14 else 11
  ggsave(file.path(fig_dir, "01_pca_cvcq_biplot.png"), combined_plot, width = panel_width, height = 9, dpi = 300)

  mean_cvcq <- cvcq_wide %>%
    pivot_longer(cols = starts_with("cvcq_"), names_to = "solute", values_to = "cvcq") %>%
    group_by(Stream_Name) %>%
    summarize(mean_cvcq = mean(cvcq, na.rm = TRUE), n_chemostatic = sum(cvcq < 1, na.rm = TRUE), n_chemodynamic = sum(cvcq > 1, na.rm = TRUE), .groups = "drop") %>%
    mutate(behavior_type = case_when(mean_cvcq < 0.5 ~ "Strongly Chemostatic", mean_cvcq < 1.0 ~ "Chemostatic", mean_cvcq < 2.0 ~ "Chemodynamic", TRUE ~ "Strongly Chemodynamic"))

  storage_col_class <- storage_available[!is.na(storage_available)][1]
  if (length(storage_col_class) == 0) storage_col_class <- NULL
  has_storage_class <- !is.null(storage_col_class)

  scores_class <- scores %>% left_join(mean_cvcq, by = "Stream_Name")
  if (has_storage_class) scores_class <- scores_class %>% mutate(storage_value = .data[[storage_col_class]])

  p_class <- ggplot(scores_class, aes(x = PC1, y = PC2)) +
    stat_ellipse(level = 0.95, type = "t", linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Strongly Chemostatic" = "#1a9850", "Chemostatic" = "#91cf60", "Chemodynamic" = "#fc8d59", "Strongly Chemodynamic" = "#d73027"), name = "Overall Behavior") +
    labs(title = "PCA: Sites Classified by Chemostatic Behavior", x = paste0("PC1 (", round(var_explained[1], 1), "%)"), y = paste0("PC2 (", round(var_explained[2], 1), "%)"), caption = if (has_storage_class) "Point size = dynamic storage" else NULL) +
    theme_bw(base_size = 12) + theme(legend.position = "right") + coord_fixed()
  if (has_storage_class) {
    p_class <- p_class + geom_point(aes(color = behavior_type, size = storage_value), alpha = 0.8) + scale_size_continuous(name = "Dynamic Storage\n(mm)", range = c(2, 8))
  } else {
    p_class <- p_class + geom_point(aes(color = behavior_type), alpha = 0.8, size = 3) + guides(size = "none")
  }
  p_class <- p_class + geom_text_repel(aes(label = Stream_Name), size = 3, max.overlaps = 20)
  ggsave(file.path(fig_dir, "03_pca_cvcq_behavior_class.png"), p_class, width = 10, height = 8, dpi = 300)

  scree_df <- tibble(PC = 1:length(var_explained), Variance = var_explained, Cumulative = cumsum(var_explained))
  p_scree <- ggplot(scree_df, aes(x = PC)) +
    geom_col(aes(y = Variance), fill = "purple4", alpha = 0.7) +
    geom_line(aes(y = Cumulative), color = "red", linewidth = 1) +
    geom_point(aes(y = Cumulative), color = "red", size = 2) +
    geom_hline(yintercept = 80, linetype = "dashed", color = "gray50") +
    labs(title = "Scree Plot: CVc/CVq PCA", x = "Principal Component", y = "Variance Explained (%)") +
    theme_bw(base_size = 12)
  ggsave(file.path(fig_dir, "04_scree_plot.png"), p_scree, width = 8, height = 6, dpi = 300)

  write_csv(scores_class, file.path(fig_dir, "pca_cvcq_scores.csv"))
  write_csv(loadings, file.path(fig_dir, "pca_cvcq_loadings.csv"))
  write_csv(scree_df, file.path(fig_dir, "pca_cvcq_variance.csv"))
}
# =============================================================================
# 4c_pca_cvcq_ratio.R
# =============================================================================
# PCA ON CVc/CVq RATIO - "BUFFERING FINGERPRINT"
#
# RESEARCH QUESTION:
#   How do sites cluster based on their variability fingerprint?
#   Which sites are chemostatic vs chemodynamic across solutes?
#
# INTERPRETATION:
#   - CVc/CVq > 1: Chemodynamic (concentration varies more than flow)
#   - CVc/CVq < 1: Chemostatic (chemistry is buffered)
#   - CVc/CVq ≈ 1: Concentration and flow vary similarly
#
#   Sites clustering together have similar patterns of buffering vs
#   responsiveness across their solute suite.
#
# OUTPUT:
#   - Biplot with 95% confidence ellipse
#   - Storage-colored version
#   - Behavior classification plot
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
fig_dir <- file.path(base_dir, "exploratory_plots", "04_PCA", "4c_cvcq_ratio")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Source helpers
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

# =============================================================================
# LOAD DATA
# =============================================================================

message("Loading site-level CVc/CVq data...")

site_means <- read_csv(file.path(out_dir, "HJA_clean_site_means.csv"), show_col_types = FALSE)

# Check for required column
if (!"cq_CVc_CVq" %in% names(site_means)) {
  stop("cq_CVc_CVq column not found in HJA_clean_site_means.csv")
}

# =============================================================================
# CREATE WIDE MATRIX: Sites × Solute CVc/CVq
# =============================================================================

message("Creating CVc/CVq matrix...")

cvcq_wide <- site_means %>%
  select(Stream_Name, solute, cq_CVc_CVq) %>%
  filter(!is.na(cq_CVc_CVq) & is.finite(cq_CVc_CVq)) %>%
  group_by(Stream_Name, solute) %>%
  summarize(cq_CVc_CVq = mean(cq_CVc_CVq, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = solute,
    values_from = cq_CVc_CVq,
    names_prefix = "cvcq_"
  )

# Store site names
site_names <- cvcq_wide$Stream_Name

# Extract numeric matrix
cvcq_mat <- cvcq_wide %>%
  select(-Stream_Name) %>%
  as.matrix()
rownames(cvcq_mat) <- site_names

# QC: Remove sites with >50% missing
row_na_frac <- rowMeans(is.na(cvcq_mat))
keep_rows <- row_na_frac <= 0.5
cvcq_mat <- cvcq_mat[keep_rows, , drop = FALSE]
site_names <- site_names[keep_rows]

# Remove solutes with no variance
col_var <- apply(cvcq_mat, 2, function(x) var(x, na.rm = TRUE))
keep_cols <- !is.na(col_var) & col_var > 0
cvcq_mat <- cvcq_mat[, keep_cols, drop = FALSE]

message("Matrix dimensions: ", nrow(cvcq_mat), " sites × ", ncol(cvcq_mat), " solutes")

# Get catchment characteristics
site_chars <- site_means %>%
  select(Stream_Name, any_of(c(
    "RBI", "FDC_slope_5_95", "WB_dS_range_mm", "Q_dS_range_mm",
    "Area_km2", "Elevation_mean_m", "Slope_mean", "MTT_final"
  ))) %>%
  distinct(Stream_Name, .keep_all = TRUE)

# =============================================================================
# RUN PCA (on log-transformed data since ratios can be skewed)
# =============================================================================

if (nrow(cvcq_mat) >= 3 && ncol(cvcq_mat) >= 2) {
  
  # Log-transform (add small constant to avoid log(0))
  cvcq_mat_log <- log10(cvcq_mat + 0.01)
  
  # Handle missing values
  cc <- complete.cases(cvcq_mat_log)
  if (sum(cc) < 3) {
    message("Not enough complete cases. Imputing missing values...")
    for (j in 1:ncol(cvcq_mat_log)) {
      cvcq_mat_log[is.na(cvcq_mat_log[, j]), j] <- median(cvcq_mat_log[, j], na.rm = TRUE)
    }
    cc <- rep(TRUE, nrow(cvcq_mat_log))
  }
  
  cvcq_mat_pca <- cvcq_mat_log[cc, , drop = FALSE]
  sites_pca <- site_names[cc]
  
  # Scale
  cvcq_scaled <- scale(cvcq_mat_pca, center = TRUE, scale = TRUE)
  
  # Run PCA
  pca_result <- prcomp(cvcq_scaled, center = FALSE, scale. = FALSE)
  
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
    mutate(solute = gsub("cvcq_", "", solute))
  
  # Scale loadings for biplot
  loading_scale <- max(abs(scores$PC1), abs(scores$PC2)) / 
                   max(abs(loadings$PC1), abs(loadings$PC2)) * 0.8
  loadings_scaled <- loadings %>%
    mutate(
      PC1_plot = PC1 * loading_scale,
      PC2_plot = PC2 * loading_scale
    )
  
  # ==========================================================================
  # BIPLOT WITH CONFIDENCE ELLIPSE
  # Color by site, size by dynamic storage
  # ==========================================================================
  
  # Determine storage column to use for sizing
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
      guides(fill = guide_legend(override.aes = list(shape = 21, colour = "grey20")))

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
    title = "PCA: Buffering Fingerprint (CVc/CVq Ratio)",
    subtitle = "Panels distinguish dynamic storage metrics (Q-dS vs WB-dS)",
    caption = "Dashed ellipse: 95% CI; positive loadings = more chemodynamic; points outlined for visibility"
  )

  panel_width <- if (length(plot_list) > 1) 14 else 11
  ggsave(
    file.path(fig_dir, "01_pca_cvcq_biplot.png"),
    combined_plot,
    width = panel_width,
    height = 9,
    dpi = 300
  )
  
  
  # ==========================================================================
  # CHEMOSTATIC VS CHEMODYNAMIC CLASSIFICATION PLOT
  # ==========================================================================
  
  # Calculate mean CVc/CVq across solutes for each site
  mean_cvcq <- cvcq_wide %>%
    pivot_longer(cols = starts_with("cvcq_"), names_to = "solute", values_to = "cvcq") %>%
    group_by(Stream_Name) %>%
    summarize(
      mean_cvcq = mean(cvcq, na.rm = TRUE),
      n_chemostatic = sum(cvcq < 1, na.rm = TRUE),
      n_chemodynamic = sum(cvcq > 1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      behavior_type = case_when(
        mean_cvcq < 0.5 ~ "Strongly Chemostatic",
        mean_cvcq < 1.0 ~ "Chemostatic",
        mean_cvcq < 2.0 ~ "Chemodynamic",
        TRUE ~ "Strongly Chemodynamic"
      )
    )
  
  storage_col_class <- storage_available[!is.na(storage_available)][1]
  if (length(storage_col_class) == 0) storage_col_class <- NULL
  has_storage_class <- !is.null(storage_col_class)

  scores_class <- scores %>%
    left_join(mean_cvcq, by = "Stream_Name")
  if (has_storage_class) {
    scores_class <- scores_class %>%
      mutate(storage_value = .data[[storage_col_class]])
  }

  p_class <- ggplot(scores_class, aes(x = PC1, y = PC2)) +
    stat_ellipse(level = 0.95, type = "t", linetype = "dashed", color = "gray50") +
    scale_color_manual(
      values = c(
        "Strongly Chemostatic" = "#1a9850",
        "Chemostatic" = "#91cf60",
        "Chemodynamic" = "#fc8d59",
        "Strongly Chemodynamic" = "#d73027"
      ),
      name = "Overall Behavior"
    ) +
    labs(
      title = "PCA: Sites Classified by Chemostatic Behavior",
      x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
      y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
      caption = if (has_storage_class) "Point size = dynamic storage" else NULL
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right") +
    coord_fixed()
  
  if (has_storage_class) {
    p_class <- p_class +
      geom_point(aes(color = behavior_type, size = storage_value), alpha = 0.8) +
      scale_size_continuous(name = "Dynamic Storage\n(mm)", range = c(2, 8))
  } else {
    p_class <- p_class +
      geom_point(aes(color = behavior_type), alpha = 0.8, size = 3) +
      guides(size = "none")
  }
  
  p_class <- p_class +
    geom_text_repel(aes(label = Stream_Name), size = 3, max.overlaps = 20)
  
  ggsave(file.path(fig_dir, "03_pca_cvcq_behavior_class.png"), p_class,
         width = 10, height = 8, dpi = 300)
  
  # ==========================================================================
  # SCREE PLOT
  # ==========================================================================
  
  scree_df <- tibble(
    PC = 1:length(var_explained),
    Variance = var_explained,
    Cumulative = cumsum(var_explained)
  )
  
  p_scree <- ggplot(scree_df, aes(x = PC)) +
    geom_col(aes(y = Variance), fill = "purple4", alpha = 0.7) +
    geom_line(aes(y = Cumulative), color = "red", linewidth = 1) +
    geom_point(aes(y = Cumulative), color = "red", size = 2) +
    geom_hline(yintercept = 80, linetype = "dashed", color = "gray50") +
    labs(
      title = "Scree Plot: CVc/CVq PCA",
      x = "Principal Component",
      y = "Variance Explained (%)"
    ) +
    theme_bw(base_size = 12)
  
  ggsave(file.path(fig_dir, "04_scree_plot.png"), p_scree,
         width = 8, height = 6, dpi = 300)
  
  # ==========================================================================
  # SAVE RESULTS
  # ==========================================================================
  
  write_csv(scores_class, file.path(fig_dir, "pca_cvcq_scores.csv"))
  write_csv(loadings, file.path(fig_dir, "pca_cvcq_loadings.csv"))
  write_csv(scree_df, file.path(fig_dir, "pca_cvcq_variance.csv"))
  
  message("\n✓ CVc/CVq PCA complete!")
  message("  Figures saved to: ", fig_dir)
  
} else {
  message("Insufficient data for PCA")
}
