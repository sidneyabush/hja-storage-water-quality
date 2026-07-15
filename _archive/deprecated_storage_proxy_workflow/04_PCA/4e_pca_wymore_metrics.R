# Restored from archive on 2025-12-04

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
})
try(source(file.path("/Users/sidneybush/Documents/GitHub/hja-water-quality", "00_helpers", "plot_theme_set.R")), silent = TRUE)

base_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality"
out_dir <- file.path(base_dir, "outputs")
fig_dir <- file.path(base_dir, "exploratory_plots", "04_PCA", "4e_wymore_metrics")
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
wymore_cols <- c("wymore_cvcq_consistency", "wymore_crosssite_allpairs", "wymore_crosssite_outlet")
available_wymore <- intersect(wymore_cols, names(sync_data))
if (length(available_wymore) == 0) {
  site_means <- read_csv(file.path(out_dir, "HJA_clean_site_means.csv"), show_col_types = FALSE)
  available_wymore <- intersect(wymore_cols, names(site_means))
  if (length(available_wymore) > 0) sync_data <- site_means
}
if (length(available_wymore) == 0) stop("No Wymore metrics found in data files")

site_chars <- read_csv(file.path(out_dir, "HJA_clean_site_means.csv"), show_col_types = FALSE) %>%
  select(Stream_Name, any_of(c("RBI", "FDC_slope_5_95", "WB_dS_range_mm", "Q_dS_range_mm", "Area_km2", "Elevation_mean_m", "MTT_final", "cq_CVc_CVq"))) %>%
  distinct(Stream_Name, .keep_all = TRUE)

wymore_site <- sync_data %>%
  group_by(Stream_Name) %>%
  summarize(across(any_of(available_wymore), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  filter(!is.na(Stream_Name))

if (length(available_wymore) >= 2) {
  wymore_mat <- wymore_site %>% select(all_of(available_wymore)) %>% as.matrix()
  rownames(wymore_mat) <- wymore_site$Stream_Name
  keep_rows <- rowSums(!is.na(wymore_mat)) > 0
  wymore_mat <- wymore_mat[keep_rows, , drop = FALSE]
  site_names <- wymore_site$Stream_Name[keep_rows]
  col_var <- apply(wymore_mat, 2, function(x) var(x, na.rm = TRUE))
  keep_cols <- !is.na(col_var) & col_var > 0
  wymore_mat <- wymore_mat[, keep_cols, drop = FALSE]
  if (nrow(wymore_mat) >= 3 && ncol(wymore_mat) >= 2) {
    cc <- complete.cases(wymore_mat)
    if (sum(cc) < 3) {
      for (j in 1:ncol(wymore_mat)) wymore_mat[is.na(wymore_mat[, j]), j] <- median(wymore_mat[, j], na.rm = TRUE)
      cc <- rep(TRUE, nrow(wymore_mat))
    }
    wymore_pca <- wymore_mat[cc, , drop = FALSE]
    sites_pca <- site_names[cc]
    wymore_scaled <- scale(wymore_pca, center = TRUE, scale = TRUE)
    pca_result <- prcomp(wymore_scaled, center = FALSE, scale. = FALSE)
    var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100
    scores <- as_tibble(pca_result$x[, 1:min(2, ncol(pca_result$x))]) %>%
      mutate(Stream_Name = sites_pca) %>%
      left_join(site_chars, by = "Stream_Name") %>%
      flag_outlet_stream() %>%
      apply_factor_orders()
    loadings <- as.data.frame(pca_result$rotation[, 1:min(2, ncol(pca_result$rotation))]) %>% rownames_to_column("metric") %>% as_tibble()
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
    combined_plot <- combined_plot + patchwork::plot_annotation(title = "PCA: Wymore Consistency Fingerprint", subtitle = "Panels distinguish dynamic storage metrics (Q-dS vs WB-dS)", caption = "Dashed ellipse: 95% CI; loadings show Wymore consistency metrics; points outlined for visibility")
    panel_width <- if (length(plot_list) > 1) 14 else 11
    ggsave(file.path(fig_dir, "01_pca_wymore_site_biplot.png"), combined_plot, width = panel_width, height = 9, dpi = 300)
    write_csv(scores, file.path(fig_dir, "pca_wymore_scores.csv"))
    write_csv(loadings, file.path(fig_dir, "pca_wymore_loadings.csv"))
  }
}
# =============================================================================
# 4e_pca_wymore_metrics.R
# =============================================================================
# PCA ON WYMORE CROSS-SITE CONSISTENCY METRICS - "STABILITY FINGERPRINT"
#
# RESEARCH QUESTION:
#   How consistent are CVc/CVq patterns across sites?
#   Do some sites show stable solute behavior while others are variable?
#
# METRICS:
#   - wymore_cvcq_consistency: Within-site CVc/CVq consistency across time
#   - wymore_crosssite_allpairs: CVc/CVq similarity to other sites
#   - wymore_crosssite_outlet: CVc/CVq similarity to outlet
#
# INTERPRETATION:
#   High consistency = stable chemostatic/chemodynamic behavior
#   Low consistency = variable behavior across time/conditions
#
# OUTPUT:
#   - Biplot with 95% confidence ellipse
#   - Solute-specific consistency PCA
#   - CVc/CVq vs consistency summary plot
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
fig_dir <- file.path(base_dir, "exploratory_plots", "04_PCA", "4e_wymore_metrics")
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

message("Loading Wymore metrics...")

sync_data <- read_csv(file.path(out_dir, "HJA_composite_synchrony.csv"), show_col_types = FALSE)

# Check available Wymore columns
wymore_cols <- c("wymore_cvcq_consistency", "wymore_crosssite_allpairs", "wymore_crosssite_outlet")
available_wymore <- intersect(wymore_cols, names(sync_data))

if (length(available_wymore) == 0) {
  message("No Wymore columns found. Checking site_means...")
  site_means <- read_csv(file.path(out_dir, "HJA_clean_site_means.csv"), show_col_types = FALSE)
  available_wymore <- intersect(wymore_cols, names(site_means))
  if (length(available_wymore) > 0) {
    sync_data <- site_means
  }
}

if (length(available_wymore) == 0) {
  stop("No Wymore metrics found in data files")
}

message("Available Wymore metrics: ", paste(available_wymore, collapse = ", "))

# Get site characteristics
site_chars <- read_csv(file.path(out_dir, "HJA_clean_site_means.csv"), show_col_types = FALSE) %>%
  select(Stream_Name, any_of(c(
    "RBI", "FDC_slope_5_95", "WB_dS_range_mm", "Q_dS_range_mm",
    "Area_km2", "Elevation_mean_m", "MTT_final", "cq_CVc_CVq"
  ))) %>%
  distinct(Stream_Name, .keep_all = TRUE)

# =============================================================================
# SITE-LEVEL WYMORE METRICS
# =============================================================================

message("\n=== Site-level Wymore PCA ===")

wymore_site <- sync_data %>%
  group_by(Stream_Name) %>%
  summarize(
    across(any_of(available_wymore), ~mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  filter(!is.na(Stream_Name))

# If only one metric, can't do PCA - do summary plots instead
if (length(available_wymore) < 2) {
  
  message("Only one Wymore metric available. Creating summary plots...")
  
  wymore_plot <- wymore_site %>%
    left_join(site_chars, by = "Stream_Name") %>%
    flag_outlet_stream()
  
  metric_name <- available_wymore[1]
  
  p_dist <- ggplot(wymore_plot, aes_string(x = "reorder(Stream_Name, -", metric_name, ")", y = metric_name)) +
    geom_col(aes(fill = outlet_marker)) +
    coord_flip() +
    scale_fill_manual(values = c("Headwater" = "steelblue", "Outlet" = "coral")) +
    labs(
      title = paste("Distribution of", metric_name),
      x = "Stream",
      y = metric_name,
      fill = "Stream Type"
    ) +
    theme_bw(base_size = 12)
  
  ggsave(file.path(fig_dir, "01_wymore_distribution.png"), p_dist,
         width = 10, height = 8, dpi = 300)
  
} else {
  
  # Create matrix
  wymore_mat <- wymore_site %>%
    select(all_of(available_wymore)) %>%
    as.matrix()
  rownames(wymore_mat) <- wymore_site$Stream_Name
  
  # QC
  keep_rows <- rowSums(!is.na(wymore_mat)) > 0
  wymore_mat <- wymore_mat[keep_rows, , drop = FALSE]
  site_names <- wymore_site$Stream_Name[keep_rows]
  
  col_var <- apply(wymore_mat, 2, function(x) var(x, na.rm = TRUE))
  keep_cols <- !is.na(col_var) & col_var > 0
  wymore_mat <- wymore_mat[, keep_cols, drop = FALSE]
  
  if (nrow(wymore_mat) >= 3 && ncol(wymore_mat) >= 2) {
    
    cc <- complete.cases(wymore_mat)
    if (sum(cc) < 3) {
      for (j in 1:ncol(wymore_mat)) {
        wymore_mat[is.na(wymore_mat[, j]), j] <- median(wymore_mat[, j], na.rm = TRUE)
      }
      cc <- rep(TRUE, nrow(wymore_mat))
    }
    
    wymore_pca <- wymore_mat[cc, , drop = FALSE]
    sites_pca <- site_names[cc]
    
    wymore_scaled <- scale(wymore_pca, center = TRUE, scale = TRUE)
    pca_result <- prcomp(wymore_scaled, center = FALSE, scale. = FALSE)
    
    var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100
    
    scores <- as_tibble(pca_result$x[, 1:min(2, ncol(pca_result$x))]) %>%
      mutate(Stream_Name = sites_pca) %>%
      left_join(site_chars, by = "Stream_Name") %>%
      flag_outlet_stream() %>%
      apply_factor_orders()
    
    loadings <- as.data.frame(pca_result$rotation[, 1:min(2, ncol(pca_result$rotation))]) %>%
      rownames_to_column("metric") %>%
      as_tibble()
    
    loading_scale <- max(abs(scores$PC1), abs(scores$PC2)) / 
                     max(abs(loadings$PC1), abs(loadings$PC2)) * 0.8
    loadings_scaled <- loadings %>%
      mutate(
        PC1_plot = PC1 * loading_scale,
        PC2_plot = PC2 * loading_scale
      )
    
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
      title = "PCA: Wymore Consistency Fingerprint",
      subtitle = "Panels distinguish dynamic storage metrics (Q-dS vs WB-dS)",
      caption = "Dashed ellipse: 95% CI; loadings show Wymore consistency metrics; points outlined for visibility"
    )

    panel_width <- if (length(plot_list) > 1) 14 else 11
    ggsave(
      file.path(fig_dir, "01_pca_wymore_site_biplot.png"),
      combined_plot,
      width = panel_width,
      height = 9,
      dpi = 300
    )
    
    write_csv(scores, file.path(fig_dir, "pca_wymore_scores.csv"))
    write_csv(loadings, file.path(fig_dir, "pca_wymore_loadings.csv"))
  }
}

# =============================================================================
# WYMORE BY SOLUTE
# =============================================================================

message("\n=== Wymore by Solute ===")

if ("wymore_cvcq_consistency" %in% names(sync_data) && "solute" %in% names(sync_data)) {
  
  wymore_sol <- sync_data %>%
    select(Stream_Name, solute, wymore_cvcq_consistency) %>%
    filter(!is.na(wymore_cvcq_consistency)) %>%
    pivot_wider(
      names_from = solute,
      values_from = wymore_cvcq_consistency,
      names_prefix = "wymore_"
    )
  
  site_names <- wymore_sol$Stream_Name
  wymore_sol_mat <- wymore_sol %>%
    select(-Stream_Name) %>%
    as.matrix()
  rownames(wymore_sol_mat) <- site_names
  
  # QC
  row_na_frac <- rowMeans(is.na(wymore_sol_mat))
  keep_rows <- row_na_frac <= 0.5
  wymore_sol_mat <- wymore_sol_mat[keep_rows, , drop = FALSE]
  site_names <- site_names[keep_rows]
  
  col_var <- apply(wymore_sol_mat, 2, function(x) var(x, na.rm = TRUE))
  keep_cols <- !is.na(col_var) & col_var > 0
  wymore_sol_mat <- wymore_sol_mat[, keep_cols, drop = FALSE]
  
  if (nrow(wymore_sol_mat) >= 3 && ncol(wymore_sol_mat) >= 2) {
    
    cc <- complete.cases(wymore_sol_mat)
    if (sum(cc) < 3) {
      for (j in 1:ncol(wymore_sol_mat)) {
        wymore_sol_mat[is.na(wymore_sol_mat[, j]), j] <- median(wymore_sol_mat[, j], na.rm = TRUE)
      }
      cc <- rep(TRUE, nrow(wymore_sol_mat))
    }
    
    wymore_sol_pca <- wymore_sol_mat[cc, , drop = FALSE]
    sites_pca <- site_names[cc]
    
    wymore_sol_scaled <- scale(wymore_sol_pca, center = TRUE, scale = TRUE)
    pca_sol <- prcomp(wymore_sol_scaled, center = FALSE, scale. = FALSE)
    
    var_exp <- (pca_sol$sdev^2 / sum(pca_sol$sdev^2)) * 100
    
    scores_sol <- as_tibble(pca_sol$x[, 1:min(2, ncol(pca_sol$x))]) %>%
      mutate(Stream_Name = sites_pca) %>%
      left_join(site_chars, by = "Stream_Name") %>%
      flag_outlet_stream()
    
    loadings_sol <- as.data.frame(pca_sol$rotation[, 1:min(2, ncol(pca_sol$rotation))]) %>%
      rownames_to_column("solute") %>%
      as_tibble() %>%
      mutate(solute = gsub("wymore_", "", solute))
    
    loading_scale <- max(abs(scores_sol$PC1), abs(scores_sol$PC2)) / 
                     max(abs(loadings_sol$PC1), abs(loadings_sol$PC2)) * 0.8
    loadings_sol_scaled <- loadings_sol %>%
      mutate(
        PC1_plot = PC1 * loading_scale,
        PC2_plot = PC2 * loading_scale
      )
    
    # Determine storage column for this scores dataframe
    storage_col_sol <- if ("Q_dS_range_mm" %in% names(scores_sol) && any(!is.na(scores_sol$Q_dS_range_mm))) {
      "Q_dS_range_mm"
    } else if ("WB_dS_range_mm" %in% names(scores_sol) && any(!is.na(scores_sol$WB_dS_range_mm))) {
      "WB_dS_range_mm"
    } else {
      NULL
    }
    
    p_sol <- ggplot() +
      # 95% confidence ellipse
      stat_ellipse(
        data = scores_sol,
        aes(x = PC1, y = PC2),
        level = 0.95, type = "t",
        linetype = "dashed", color = "gray50", linewidth = 0.5
      ) +
      geom_point(
        data = scores_sol,
        aes(x = PC1, y = PC2, color = Stream_Name,
            size = if (!is.null(storage_col_sol)) .data[[storage_col_sol]] else NULL),
        alpha = 0.8
      ) +
      geom_text_repel(
        data = scores_sol,
        aes(x = PC1, y = PC2, label = Stream_Name),
        size = 3, max.overlaps = 20
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
        size = 3.5, color = "gray30", fontface = "bold"
      ) +
      scale_color_site() +
      scale_size_continuous(name = "Dynamic Storage\n(mm)", range = c(2, 8)) +
      labs(
        title = "PCA: Stability Fingerprint by Solute",
        subtitle = "Which solutes have consistent behavior across sites?",
        x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
        y = paste0("PC2 (", round(var_exp[2], 1), "%)"),
        caption = "Dashed ellipse: 95% CI\nPoint size = dynamic storage"
      ) +
      theme_bw(base_size = 12) +
      theme(legend.position = "right") +
      coord_fixed() +
      guides(color = guide_legend(override.aes = list(size = 4)))
    
    ggsave(file.path(fig_dir, "02_pca_wymore_by_solute.png"), p_sol,
           width = 10, height = 9, dpi = 300)
    
    write_csv(scores_sol, file.path(fig_dir, "pca_wymore_solute_scores.csv"))
    write_csv(loadings_sol, file.path(fig_dir, "pca_wymore_solute_loadings.csv"))
  }
}

# =============================================================================
# SUMMARY PLOT: CONSISTENCY VS CVc/CVq
# =============================================================================

message("\n=== Summary Plots ===")

if ("wymore_cvcq_consistency" %in% names(sync_data)) {
  
  summary_df <- sync_data %>%
    group_by(Stream_Name) %>%
    summarize(
      mean_consistency = mean(wymore_cvcq_consistency, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(
      site_chars %>% select(Stream_Name, cq_CVc_CVq, WB_dS_range_mm),
      by = "Stream_Name"
    ) %>%
    flag_outlet_stream()
  
  if ("cq_CVc_CVq" %in% names(summary_df) && any(!is.na(summary_df$cq_CVc_CVq))) {
    
    # Calculate mean CVc/CVq per site from site_means
    cvcq_means <- read_csv(file.path(out_dir, "HJA_clean_site_means.csv"), show_col_types = FALSE) %>%
      group_by(Stream_Name) %>%
      summarize(mean_cvcq = mean(cq_CVc_CVq, na.rm = TRUE), .groups = "drop")
    
    summary_df <- summary_df %>%
      left_join(cvcq_means, by = "Stream_Name")
    
    p_summary <- ggplot(summary_df, aes(x = mean_cvcq, y = mean_consistency)) +
      geom_point(aes(color = Stream_Name), size = 3, alpha = 0.8) +
      geom_text_repel(aes(label = Stream_Name), size = 3, max.overlaps = 15) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
      geom_smooth(method = "lm", se = TRUE, color = "gray50", alpha = 0.2) +
      scale_color_site() +
      labs(
        title = "CVc/CVq vs Behavioral Consistency",
        subtitle = "Do chemostatic sites have more consistent behavior?",
        x = "Mean CVc/CVq (log scale)",
        y = "Wymore Consistency",
        caption = "Vertical line: CVc/CVq = 1 (boundary between chemostatic/chemodynamic)"
      ) +
      scale_x_log10() +
      theme_bw(base_size = 12) +
      theme(legend.position = "none")
    
    ggsave(file.path(fig_dir, "03_cvcq_vs_consistency.png"), p_summary,
           width = 9, height = 8, dpi = 300)
  }
}

message("\n✓ Wymore PCA complete!")
message("  Figures saved to: ", fig_dir)
