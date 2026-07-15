# Restored from archive on 2025-12-04

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
})
try(source(file.path("/Users/sidneybush/Documents/GitHub/hja-water-quality", "00_helpers", "plot_theme_set.R")), silent = TRUE)

base_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality"
out_dir <- file.path(base_dir, "outputs")
fig_dir <- file.path(base_dir, "exploratory_plots", "04_PCA", "4d_synchrony_abbott")
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

sync_data <- read_csv(file.path(out_dir, "HJA_composite_synchrony.csv"), show_col_types = FALSE)
sync_cols <- c("conc_sync_allpairs", "cqslope_sync_allpairs", "conc_sync_outlet", "cqslope_sync_outlet")
available_sync <- intersect(sync_cols, names(sync_data))
if (length(available_sync) == 0) stop("No synchrony columns found in data")

site_chars <- read_csv(file.path(out_dir, "HJA_clean_site_means.csv"), show_col_types = FALSE) %>%
  select(Stream_Name, any_of(c("RBI", "FDC_slope_5_95", "WB_dS_range_mm", "Q_dS_range_mm", "Area_km2", "Elevation_mean_m", "Slope_mean", "MTT_final"))) %>%
  distinct(Stream_Name, .keep_all = TRUE)

sync_site <- sync_data %>%
  group_by(Stream_Name) %>%
  summarize(across(all_of(available_sync), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  filter(!is.na(Stream_Name))

sync_mat <- sync_site %>% select(all_of(available_sync)) %>% as.matrix()
rownames(sync_mat) <- sync_site$Stream_Name
keep_rows <- rowSums(!is.na(sync_mat)) > 0
sync_mat <- sync_mat[keep_rows, , drop = FALSE]
site_names <- sync_site$Stream_Name[keep_rows]
col_var <- apply(sync_mat, 2, function(x) var(x, na.rm = TRUE))
keep_cols <- !is.na(col_var) & col_var > 0
sync_mat <- sync_mat[, keep_cols, drop = FALSE]

if (nrow(sync_mat) >= 3 && ncol(sync_mat) >= 2) {
  cc <- complete.cases(sync_mat)
  if (sum(cc) < 3) {
    for (j in 1:ncol(sync_mat)) sync_mat[is.na(sync_mat[, j]), j] <- median(sync_mat[, j], na.rm = TRUE)
    cc <- rep(TRUE, nrow(sync_mat))
  }
  sync_mat_pca <- sync_mat[cc, , drop = FALSE]
  sites_pca <- site_names[cc]
  sync_scaled <- scale(sync_mat_pca, center = TRUE, scale = TRUE)
  pca_result <- prcomp(sync_scaled, center = FALSE, scale. = FALSE)
  var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100

  scores <- as_tibble(pca_result$x[, 1:min(2, ncol(pca_result$x))]) %>%
    mutate(Stream_Name = sites_pca) %>%
    left_join(site_chars, by = "Stream_Name") %>%
    flag_outlet_stream() %>%
    apply_factor_orders()

  loadings <- as.data.frame(pca_result$rotation[, 1:min(2, ncol(pca_result$rotation))]) %>%
    rownames_to_column("metric") %>% as_tibble()

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
      geom_text_repel(data = loadings_scaled, aes(x = PC1_plot, y = PC2_plot, label = metric), size = 3.5, color = "gray30", fontface = "bold", max.overlaps = 20) +
      scale_fill_site(name = "Site") +
      labs(title = storage_label, x = paste0("PC1 (", round(var_explained[1], 1), "%)"), y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
      theme_bw(base_size = 12) + theme(legend.position = "right") + coord_fixed() +
      guides(fill = guide_legend(override.aes = list(shape = 21, colour = "grey20")))
    if (has_storage) p <- p + scale_size_continuous(name = storage_label, range = c(2.5, 8), breaks = scales::pretty_breaks(n = 4)) else p <- p + guides(size = "none")
    p
  }

  plot_list <- purrr::map(storage_available, build_biplot)
  combined_plot <- if (length(plot_list) == 1) plot_list[[1]] else patchwork::wrap_plots(plotlist = plot_list, nrow = 1)
  combined_plot <- combined_plot + patchwork::plot_annotation(title = "PCA: Abbott Synchrony Fingerprint", subtitle = "Panels distinguish dynamic storage metrics (Q-dS vs WB-dS)", caption = "Dashed ellipse: 95% CI; loadings show synchrony metrics; points outlined for visibility")
  panel_width <- if (length(plot_list) > 1) 14 else 11
  ggsave(file.path(fig_dir, "01_pca_sync_site_biplot.png"), combined_plot, width = panel_width, height = 9, dpi = 300)

  write_csv(scores, file.path(fig_dir, "pca_sync_site_scores.csv"))
  write_csv(loadings, file.path(fig_dir, "pca_sync_site_loadings.csv"))
}
# =============================================================================
# 4d_pca_synchrony_abbott.R
# =============================================================================
# PCA ON ABBOTT SYNCHRONY METRICS - "SPATIAL COHERENCE FINGERPRINT"
#
# RESEARCH QUESTION:
#   How do sites cluster based on spatial coherence (synchrony)?
#   Which sites behave "in sync" with others vs have local control?
#
# UPDATED DECEMBER 2025:
#   - Focus on OUTLET synchrony (sync with GSLOOK), not all-pairs
#   - Abbott CQ-slope sync deprecated; use Wymore CQ-slope sync instead
#   - Key finding: Catchment characteristics don't predict outlet synchrony well
#     (but DO predict CQ behavior itself — see 3m script)
#
# METRICS:
#   - conc_sync_outlet: Absolute concentration synchrony with the outlet (GSLOOK)
#   - cqslope_sync_outlet: Wymore CQ-slope synchrony with the outlet
#   - conc_sync_allpairs: (secondary) All-pairs concentration sync
#
# INTERPRETATION:
#   - High outlet sync = site tracks watershed-scale dynamics
#   - Low outlet sync = site has local/decoupled behavior
#
# OUTPUT:
#   - Biplot with 95% confidence ellipse
#   - Storage-colored version
#   - Solute-specific synchrony PCA
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
fig_dir <- file.path(base_dir, "exploratory_plots", "04_PCA", "4d_synchrony_abbott")
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

message("Loading synchrony data...")

sync_data <- read_csv(file.path(out_dir, "HJA_composite_synchrony.csv"), show_col_types = FALSE)

# Check available synchrony columns
sync_cols <- c("conc_sync_allpairs", "cqslope_sync_allpairs", 
               "conc_sync_outlet", "cqslope_sync_outlet")
available_sync <- intersect(sync_cols, names(sync_data))

if (length(available_sync) == 0) {
  stop("No synchrony columns found in data")
}

message("Available synchrony metrics: ", paste(available_sync, collapse = ", "))

# Get site characteristics for coloring
site_chars <- read_csv(file.path(out_dir, "HJA_clean_site_means.csv"), show_col_types = FALSE) %>%
  select(Stream_Name, any_of(c(
    "RBI", "FDC_slope_5_95", "WB_dS_range_mm", "Q_dS_range_mm",
    "Area_km2", "Elevation_mean_m", "Slope_mean", "MTT_final"
  ))) %>%
  distinct(Stream_Name, .keep_all = TRUE)

# =============================================================================
# APPROACH 1: PCA ON SITE-LEVEL SYNCHRONY (averaging across solutes)
# =============================================================================

message("\n=== Site-level Synchrony PCA ===")

sync_site <- sync_data %>%
  group_by(Stream_Name) %>%
  summarize(
    across(all_of(available_sync), ~mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  filter(!is.na(Stream_Name))

# Create matrix
sync_mat <- sync_site %>%
  select(all_of(available_sync)) %>%
  as.matrix()
rownames(sync_mat) <- sync_site$Stream_Name

# Remove rows with all NA
keep_rows <- rowSums(!is.na(sync_mat)) > 0
sync_mat <- sync_mat[keep_rows, , drop = FALSE]
site_names <- sync_site$Stream_Name[keep_rows]

# Remove columns with no variance
col_var <- apply(sync_mat, 2, function(x) var(x, na.rm = TRUE))
keep_cols <- !is.na(col_var) & col_var > 0
sync_mat <- sync_mat[, keep_cols, drop = FALSE]

if (nrow(sync_mat) >= 3 && ncol(sync_mat) >= 2) {
  
  # Handle missing values
  cc <- complete.cases(sync_mat)
  if (sum(cc) < 3) {
    for (j in 1:ncol(sync_mat)) {
      sync_mat[is.na(sync_mat[, j]), j] <- median(sync_mat[, j], na.rm = TRUE)
    }
    cc <- rep(TRUE, nrow(sync_mat))
  }
  
  sync_mat_pca <- sync_mat[cc, , drop = FALSE]
  sites_pca <- site_names[cc]
  
  # Scale
  sync_scaled <- scale(sync_mat_pca, center = TRUE, scale = TRUE)
  
  # Run PCA
  pca_result <- prcomp(sync_scaled, center = FALSE, scale. = FALSE)
  
  var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100
  
  message("Variance explained: PC1 = ", round(var_explained[1], 1), 
          "%, PC2 = ", round(var_explained[2], 1), "%")
  
  # Scores and loadings
  scores <- as_tibble(pca_result$x[, 1:min(2, ncol(pca_result$x))]) %>%
    mutate(Stream_Name = sites_pca) %>%
    left_join(site_chars, by = "Stream_Name") %>%
    flag_outlet_stream() %>%
    apply_factor_orders()
  
  loadings <- as.data.frame(pca_result$rotation[, 1:min(2, ncol(pca_result$rotation))]) %>%
    rownames_to_column("metric") %>%
    as_tibble()
  
  # Scale loadings for biplot
  loading_scale <- max(abs(scores$PC1), abs(scores$PC2)) / 
                   max(abs(loadings$PC1), abs(loadings$PC2)) * 0.8
  loadings_scaled <- loadings %>%
    mutate(
      PC1_plot = PC1 * loading_scale,
      PC2_plot = PC2 * loading_scale
    )
  
  # Biplot with confidence ellipse - color by site, size by storage
  
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
        aes(x = PC1_plot, y = PC2_plot, label = metric),
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
    title = "PCA: Abbott Synchrony Fingerprint",
    subtitle = "Panels distinguish dynamic storage metrics (Q-dS vs WB-dS)",
    caption = "Dashed ellipse: 95% CI; loadings show synchrony metrics; points outlined for visibility"
  )

  panel_width <- if (length(plot_list) > 1) 14 else 11
  ggsave(
    file.path(fig_dir, "01_pca_sync_site_biplot.png"),
    combined_plot,
    width = panel_width,
    height = 9,
    dpi = 300
  )
  
  # Save results
  write_csv(scores, file.path(fig_dir, "pca_sync_site_scores.csv"))
  write_csv(loadings, file.path(fig_dir, "pca_sync_site_loadings.csv"))
}

# =============================================================================
# APPROACH 2: PCA ON SYNCHRONY BY SOLUTE
# =============================================================================

message("\n=== Synchrony by Solute PCA ===")

# Create wide matrix: sites × solute synchrony
if ("conc_sync_allpairs" %in% names(sync_data)) {
  
  sync_wide <- sync_data %>%
    select(Stream_Name, solute, conc_sync_allpairs) %>%
    filter(!is.na(conc_sync_allpairs)) %>%
    pivot_wider(
      names_from = solute,
      values_from = conc_sync_allpairs,
      names_prefix = "sync_"
    )
  
  site_names <- sync_wide$Stream_Name
  sync_sol_mat <- sync_wide %>%
    select(-Stream_Name) %>%
    as.matrix()
  rownames(sync_sol_mat) <- site_names
  
  # QC
  row_na_frac <- rowMeans(is.na(sync_sol_mat))
  keep_rows <- row_na_frac <= 0.5
  sync_sol_mat <- sync_sol_mat[keep_rows, , drop = FALSE]
  site_names <- site_names[keep_rows]
  
  col_var <- apply(sync_sol_mat, 2, function(x) var(x, na.rm = TRUE))
  keep_cols <- !is.na(col_var) & col_var > 0
  sync_sol_mat <- sync_sol_mat[, keep_cols, drop = FALSE]
  
  if (nrow(sync_sol_mat) >= 3 && ncol(sync_sol_mat) >= 2) {
    
    # Handle NAs
    cc <- complete.cases(sync_sol_mat)
    if (sum(cc) < 3) {
      for (j in 1:ncol(sync_sol_mat)) {
        sync_sol_mat[is.na(sync_sol_mat[, j]), j] <- median(sync_sol_mat[, j], na.rm = TRUE)
      }
      cc <- rep(TRUE, nrow(sync_sol_mat))
    }
    
    sync_sol_pca <- sync_sol_mat[cc, , drop = FALSE]
    sites_pca <- site_names[cc]
    
    sync_sol_scaled <- scale(sync_sol_pca, center = TRUE, scale = TRUE)
    pca_sol <- prcomp(sync_sol_scaled, center = FALSE, scale. = FALSE)
    
    var_exp_sol <- (pca_sol$sdev^2 / sum(pca_sol$sdev^2)) * 100
    
    scores_sol <- as_tibble(pca_sol$x[, 1:min(2, ncol(pca_sol$x))]) %>%
      mutate(Stream_Name = sites_pca) %>%
      left_join(site_chars, by = "Stream_Name") %>%
      flag_outlet_stream()
    
    loadings_sol <- as.data.frame(pca_sol$rotation[, 1:min(2, ncol(pca_sol$rotation))]) %>%
      rownames_to_column("solute") %>%
      as_tibble() %>%
      mutate(solute = gsub("sync_", "", solute))
    
    loading_scale <- max(abs(scores_sol$PC1), abs(scores_sol$PC2)) / 
                     max(abs(loadings_sol$PC1), abs(loadings_sol$PC2)) * 0.8
    loadings_sol_scaled <- loadings_sol %>%
      mutate(
        PC1_plot = PC1 * loading_scale,
        PC2_plot = PC2 * loading_scale
      )
    
    storage_col_sol <- if ("Q_dS_range_mm" %in% names(scores_sol) && any(is.finite(scores_sol$Q_dS_range_mm))) {
      "Q_dS_range_mm"
    } else if ("WB_dS_range_mm" %in% names(scores_sol) && any(is.finite(scores_sol$WB_dS_range_mm))) {
      "WB_dS_range_mm"
    } else {
      NULL
    }

    if (!is.null(storage_col_sol)) {
      scores_sol <- scores_sol %>% mutate(storage_value = .data[[storage_col_sol]])
    }

    p_sol_biplot <- ggplot() +
      stat_ellipse(
        data = scores_sol,
        aes(x = PC1, y = PC2),
        level = 0.95, type = "t",
        linetype = "dashed", color = "gray50", linewidth = 0.5
      ) +
      geom_segment(
        data = loadings_sol_scaled,
        aes(x = 0, y = 0, xend = PC1_plot, yend = PC2_plot),
        arrow = arrow(length = unit(0.02, "npc")),
        color = "gray40", linewidth = 0.8
      ) +
      geom_text_repel(
        data = loadings_sol_scaled,
        aes(x = PC1_plot, y = PC2_plot, label = solute),
        size = 3.5, color = "gray30", fontface = "bold",
        max.overlaps = 20
      ) +
      scale_color_site() +
      labs(
        title = "PCA: Synchrony Fingerprint by Solute",
        subtitle = "Which solutes drive spatial coherence patterns?",
        x = paste0("PC1 (", round(var_exp_sol[1], 1), "%)"),
        y = paste0("PC2 (", round(var_exp_sol[2], 1), "%)"),
        caption = if (!is.null(storage_col_sol)) "Dashed ellipse: 95% CI\nPoint size = dynamic storage" else "Dashed ellipse: 95% CI"
      ) +
      theme_bw(base_size = 12) +
      theme(legend.position = "right") +
      coord_fixed() +
      guides(color = guide_legend(override.aes = list(size = 4)))

    if (!is.null(storage_col_sol)) {
      p_sol_biplot <- p_sol_biplot +
        geom_point(
          data = scores_sol,
          aes(x = PC1, y = PC2, color = Stream_Name, size = storage_value),
          alpha = 0.8
        ) +
        scale_size_continuous(name = "Dynamic Storage\n(mm)", range = c(2, 8))
    } else {
      p_sol_biplot <- p_sol_biplot +
        geom_point(
          data = scores_sol,
          aes(x = PC1, y = PC2, color = Stream_Name),
          size = 3,
          alpha = 0.8
        ) +
        guides(size = "none")
    }

    p_sol_biplot <- p_sol_biplot +
      geom_text_repel(
        data = scores_sol,
        aes(x = PC1, y = PC2, label = Stream_Name),
        size = 3, max.overlaps = 20, segment.alpha = 0.3
      )
    
    ggsave(file.path(fig_dir, "03_pca_sync_by_solute_biplot.png"), p_sol_biplot,
           width = 10, height = 9, dpi = 300)
    
    write_csv(scores_sol, file.path(fig_dir, "pca_sync_solute_scores.csv"))
    write_csv(loadings_sol, file.path(fig_dir, "pca_sync_solute_loadings.csv"))
  }
}

# =============================================================================
# SYNCHRONY INTERPRETATION PLOTS
# =============================================================================

message("\n=== Synchrony Summary Plots ===")

# Correlation between concentration and C-Q slope synchrony
if (all(c("conc_sync_allpairs", "cqslope_sync_allpairs") %in% names(sync_site))) {
  
  sync_corr <- sync_site %>%
    left_join(site_chars, by = "Stream_Name") %>%
    flag_outlet_stream()
  
  p_corr <- ggplot(sync_corr, aes(x = conc_sync_allpairs, y = cqslope_sync_allpairs)) +
    geom_point(aes(color = Stream_Name), size = 3, alpha = 0.8) +
    geom_text_repel(aes(label = Stream_Name), size = 3, max.overlaps = 15) +
    geom_smooth(method = "lm", se = TRUE, color = "gray50", alpha = 0.2) +
    scale_color_site() +
    labs(
      title = "Absolute synchrony tradeoffs",
      subtitle = "Do sites synchronous in concentrations also synchronize in C-Q behavior?",
      x = get_sync_label("conc_sync_allpairs"),
      y = get_sync_label("cqslope_sync_allpairs")
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")
  
  ggsave(file.path(fig_dir, "04_conc_vs_cqslope_sync.png"), p_corr,
         width = 9, height = 8, dpi = 300)
}

message("\n✓ Abbott Synchrony PCA complete!")
message("  Figures saved to: ", fig_dir)
