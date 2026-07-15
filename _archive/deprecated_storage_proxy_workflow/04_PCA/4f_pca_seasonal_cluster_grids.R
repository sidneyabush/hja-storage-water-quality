# Restored from archive on 2025-12-04

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
})
try(source(file.path("/Users/sidneybush/Documents/GitHub/hja-water-quality", "00_helpers", "plot_theme_set.R")), silent = TRUE)

base_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality"
out_dir <- file.path(base_dir, "outputs")
fig_dir <- file.path(base_dir, "exploratory_plots", "04_PCA", "4f_seasonal_cluster_grids")
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

run_subset_pca <- function(data_wide, site_col = "Stream_Name", prefix = "", title = "", arrow_color = "gray40", site_chars = NULL) {
  site_names <- data_wide[[site_col]]
  data_mat <- data_wide %>% select(-all_of(site_col)) %>% as.matrix()
  rownames(data_mat) <- site_names
  row_na_frac <- rowMeans(is.na(data_mat))
  keep_rows <- row_na_frac <= 0.5
  if (sum(keep_rows) < 3) return(list(plot = NULL, scores = NULL, message = "Too few valid sites"))
  data_mat <- data_mat[keep_rows, , drop = FALSE]
  site_names <- site_names[keep_rows]
  col_var <- apply(data_mat, 2, function(x) var(x, na.rm = TRUE))
  keep_cols <- !is.na(col_var) & col_var > 0
  if (sum(keep_cols) < 2) return(list(plot = NULL, scores = NULL, message = "Too few valid variables"))
  data_mat <- data_mat[, keep_cols, drop = FALSE]
  for (j in 1:ncol(data_mat)) data_mat[is.na(data_mat[, j]), j] <- median(data_mat[, j], na.rm = TRUE)
  data_scaled <- scale(data_mat, center = TRUE, scale = TRUE)
  pca_result <- prcomp(data_scaled, center = FALSE, scale. = FALSE)
  var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100
  scores <- as_tibble(pca_result$x[, 1:min(2, ncol(pca_result$x))]) %>% mutate(Stream_Name = site_names)
  if (!is.null(site_chars)) scores <- scores %>% left_join(site_chars, by = "Stream_Name")
  scores <- scores %>% flag_outlet_stream()
  scores <- apply_factor_orders(scores)
  loadings <- as.data.frame(pca_result$rotation[, 1:min(2, ncol(pca_result$rotation))]) %>% rownames_to_column("variable") %>% as_tibble() %>% mutate(variable = gsub(paste0("^", prefix, "_?"), "", variable))
  loading_scale <- max(abs(scores$PC1), abs(scores$PC2)) / max(abs(loadings$PC1), abs(loadings$PC2)) * 0.7
  loadings_scaled <- loadings %>% mutate(PC1_plot = PC1 * loading_scale, PC2_plot = PC2 * loading_scale)
  storage_candidates <- c("Q_dS_range_mm", "WB_dS_range_mm")
  storage_available <- purrr::keep(storage_candidates, ~ .x %in% names(scores) && any(is.finite(scores[[.x]])))
  if (length(storage_available) == 0) storage_available <- NA_character_
  build_panel <- function(storage_col) {
    has_storage <- !is.na(storage_col)
    storage_label <- if (has_storage) get_storage_label(storage_col) else "Storage metric unavailable"
    point_mapping <- if (has_storage) aes(x = PC1, y = PC2, fill = Stream_Name, size = .data[[storage_col]]) else aes(x = PC1, y = PC2, fill = Stream_Name)
    panel <- ggplot() +
      stat_ellipse(data = scores, aes(x = PC1, y = PC2), level = 0.95, type = "t", linetype = "dashed", color = "gray50", linewidth = 0.4) +
      geom_point(data = scores, mapping = point_mapping, shape = 21, colour = "grey20", stroke = 0.5, alpha = 0.9) +
      geom_text_repel(data = scores, aes(x = PC1, y = PC2, label = Stream_Name), size = 2.5, max.overlaps = 15, segment.alpha = 0.3) +
      geom_segment(data = loadings_scaled, aes(x = 0, y = 0, xend = PC1_plot, yend = PC2_plot), arrow = arrow(length = unit(0.015, "npc")), color = "gray40", linewidth = 0.6) +
      geom_text_repel(data = loadings_scaled, aes(x = PC1_plot, y = PC2_plot, label = variable), size = 2.5, color = "gray30", fontface = "bold", max.overlaps = 15) +
      scale_fill_site(name = "Site") +
      labs(title = storage_label, x = paste0("PC1 (", round(var_explained[1], 1), "%)"), y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
      theme_bw(base_size = 10) + theme(legend.position = "none", plot.title = element_text(size = 10, face = "bold")) + coord_fixed()
    if (has_storage) panel <- panel + scale_size_continuous(name = storage_label, range = c(1.5, 5), breaks = scales::pretty_breaks(n = 3)) else panel <- panel + guides(size = "none")
    panel
  }
  plot_list <- purrr::map(storage_available, build_panel)
  combined_plot <- if (length(plot_list) == 1) plot_list[[1]] else patchwork::wrap_plots(plotlist = plot_list, nrow = 1)
  combined_plot <- combined_plot + patchwork::plot_annotation(title = title)
  return(list(plot = combined_plot, scores = scores, var_explained = var_explained))
}

site_means <- read_csv(file.path(out_dir, "HJA_clean_site_means.csv"), show_col_types = FALSE) %>%
  filter(!is.na(Stream_Name) & Stream_Name != "", !is.na(solute) & solute != "")
cq_master <- read_csv(file.path(out_dir, "HJA_CQ_master.csv"), show_col_types = FALSE) %>%
  filter(!is.na(Stream_Name) & Stream_Name != "")
if ("variable" %in% names(cq_master) && !"solute" %in% names(cq_master)) cq_master <- cq_master %>% rename(solute = variable)
if ("value" %in% names(cq_master) && !"Concentration" %in% names(cq_master)) cq_master <- cq_master %>% rename(Concentration = value)
site_chars <- site_means %>% select(Stream_Name, any_of(c("RBI", "WB_dS_range_mm", "Q_dS_range_mm"))) %>% distinct(Stream_Name, .keep_all = TRUE)
has_season <- "hydrologic_season" %in% names(site_means) || "hydrologic_season" %in% names(cq_master)
has_cluster <- "Cluster_mode_wy" %in% names(site_means)

if (has_season) {
  if ("hydrologic_season" %in% names(site_means)) {
    season_data <- site_means
  } else {
    season_info <- cq_master %>% select(Stream_Name, solute, hydrologic_season) %>% distinct()
    season_data <- site_means %>% left_join(season_info, by = c("Stream_Name", "solute"))
  }
  seasons <- unique(season_data$hydrologic_season)
  seasons <- seasons[!is.na(seasons)]
  if (length(seasons) >= 2) {
    seasonal_plots <- list()
    for (szn in c("Wet", "Dry")) {
      if (szn %in% seasons) {
        slope_wide <- season_data %>% filter(hydrologic_season == szn) %>% select(Stream_Name, solute, cq_slope) %>% group_by(Stream_Name, solute) %>% summarize(cq_slope = mean(cq_slope, na.rm = TRUE), .groups = "drop") %>% pivot_wider(names_from = solute, values_from = cq_slope, names_prefix = "slope_")
        result <- run_subset_pca(slope_wide, site_col = "Stream_Name", prefix = "slope", title = paste(szn, "Season"), arrow_color = ifelse(szn == "Wet", "steelblue", "darkorange"), site_chars = site_chars)
        if (!is.null(result$plot)) seasonal_plots[[szn]] <- result$plot
      }
    }
    if (length(seasonal_plots) == 2) {
      p_grid <- seasonal_plots[["Wet"]] + seasonal_plots[["Dry"]] + plot_annotation(title = "Chemodynamic Fingerprint: Seasonal Comparison (C-Q Slopes)", subtitle = "How does chemistry-flow response differ between wet and dry seasons?", caption = "Dashed ellipse: 95% CI | Each panel run as separate PCA", theme = theme(plot.title = element_text(face = "bold", size = 12)))
      ggsave(file.path(fig_dir, "01_cqslope_seasonal_grid.png"), p_grid, width = 14, height = 7, dpi = 300)
    }
    seasonal_cvcq <- list()
    for (szn in c("Wet", "Dry")) {
      if (szn %in% seasons) {
        cvcq_wide <- season_data %>% filter(hydrologic_season == szn) %>% select(Stream_Name, solute, cq_CVc_CVq) %>% filter(!is.na(cq_CVc_CVq) & is.finite(cq_CVc_CVq)) %>% group_by(Stream_Name, solute) %>% summarize(cq_CVc_CVq = mean(cq_CVc_CVq, na.rm = TRUE), .groups = "drop") %>% mutate(cq_CVc_CVq = log10(cq_CVc_CVq + 0.01)) %>% pivot_wider(names_from = solute, values_from = cq_CVc_CVq, names_prefix = "cvcq_")
        result <- run_subset_pca(cvcq_wide, site_col = "Stream_Name", prefix = "cvcq", title = paste(szn, "Season"), arrow_color = ifelse(szn == "Wet", "purple4", "coral"), site_chars = site_chars)
        if (!is.null(result$plot)) seasonal_cvcq[[szn]] <- result$plot
      }
    }
    if (length(seasonal_cvcq) == 2) {
      p_grid <- seasonal_cvcq[["Wet"]] + seasonal_cvcq[["Dry"]] + plot_annotation(title = "Buffering Fingerprint: Seasonal Comparison (CVc/CVq)", subtitle = "Does chemostatic/chemodynamic behavior differ between seasons?", caption = "Dashed ellipse: 95% CI | Log10-transformed CVc/CVq", theme = theme(plot.title = element_text(face = "bold", size = 12)))
      ggsave(file.path(fig_dir, "02_cvcq_seasonal_grid.png"), p_grid, width = 14, height = 7, dpi = 300)
    }
  }
}

if (has_cluster) {
  clusters <- sort(unique(site_means$Cluster_mode_wy))
  clusters <- clusters[!is.na(clusters)]
  if (length(clusters) >= 2) {
    cluster_plots <- list()
    for (cl in clusters) {
      slope_wide <- site_means %>% filter(Cluster_mode_wy == cl) %>% select(Stream_Name, solute, cq_slope) %>% group_by(Stream_Name, solute) %>% summarize(cq_slope = mean(cq_slope, na.rm = TRUE), .groups = "drop") %>% pivot_wider(names_from = solute, values_from = cq_slope, names_prefix = "slope_")
      cluster_colors <- c("1" = "#CFA980", "2" = "#98B89F", "3" = "#5E8AA1", "4" = "#526B8E")
      result <- run_subset_pca(slope_wide, site_col = "Stream_Name", prefix = "slope", title = paste("Cluster", cl), arrow_color = cluster_colors[as.character(cl)], site_chars = site_chars)
      if (!is.null(result$plot)) cluster_plots[[as.character(cl)]] <- result$plot
    }
    if (length(cluster_plots) >= 2) {
      if (length(cluster_plots) == 2) {
        p_grid <- cluster_plots[[1]] + cluster_plots[[2]]
      } else if (length(cluster_plots) == 3) {
        p_grid <- (cluster_plots[[1]] + cluster_plots[[2]]) / (cluster_plots[[3]] + plot_spacer())
      } else {
        p_grid <- (cluster_plots[[1]] + cluster_plots[[2]]) / (cluster_plots[[3]] + cluster_plots[[4]])
      }
      p_grid <- p_grid + plot_annotation(title = "Chemodynamic Fingerprint: Cluster Comparison (C-Q Slopes)", subtitle = "How does chemistry-flow response vary across hydrologic clusters?", caption = "Dashed ellipse: 95% CI | Each panel run as separate PCA", theme = theme(plot.title = element_text(face = "bold", size = 12)))
      ggsave(file.path(fig_dir, "03_cqslope_cluster_grid.png"), p_grid, width = 14, height = 12, dpi = 300)
    }
    cluster_cvcq <- list()
    for (cl in clusters) {
      cvcq_wide <- site_means %>% filter(Cluster_mode_wy == cl) %>% select(Stream_Name, solute, cq_CVc_CVq) %>% filter(!is.na(cq_CVc_CVq) & is.finite(cq_CVc_CVq)) %>% group_by(Stream_Name, solute) %>% summarize(cq_CVc_CVq = mean(cq_CVc_CVq, na.rm = TRUE), .groups = "drop") %>% mutate(cq_CVc_CVq = log10(cq_CVc_CVq + 0.01)) %>% pivot_wider(names_from = solute, values_from = cq_CVc_CVq, names_prefix = "cvcq_")
      cluster_colors <- c("1" = "#CFA980", "2" = "#98B89F", "3" = "#5E8AA1", "4" = "#526B8E")
      result <- run_subset_pca(cvcq_wide, site_col = "Stream_Name", prefix = "cvcq", title = paste("Cluster", cl), arrow_color = cluster_colors[as.character(cl)], site_chars = site_chars)
      if (!is.null(result$plot)) cluster_cvcq[[as.character(cl)]] <- result$plot
    }
    if (length(cluster_cvcq) >= 2) {
      if (length(cluster_cvcq) == 2) {
        p_grid <- cluster_cvcq[[1]] + cluster_cvcq[[2]]
      } else if (length(cluster_cvcq) == 3) {
        p_grid <- (cluster_cvcq[[1]] + cluster_cvcq[[2]]) / (cluster_cvcq[[3]] + plot_spacer())
      } else {
        p_grid <- (cluster_cvcq[[1]] + cluster_cvcq[[2]]) / (cluster_cvcq[[3]] + cluster_cvcq[[4]])
      }
      p_grid <- p_grid + plot_annotation(title = "Buffering Fingerprint: Cluster Comparison (CVc/CVq)", subtitle = "Does chemostatic/chemodynamic behavior vary across hydrologic clusters?", caption = "Dashed ellipse: 95% CI | Log10-transformed CVc/CVq", theme = theme(plot.title = element_text(face = "bold", size = 12)))
      ggsave(file.path(fig_dir, "04_cvcq_cluster_grid.png"), p_grid, width = 14, height = 12, dpi = 300)
    }
  }
}
# =============================================================================
# 4f_pca_seasonal_cluster_grids.R
# =============================================================================
# PCAs SPLIT BY SEASON AND CLUSTER
#
# PURPOSE:
#   Create grid-layout PCA comparisons showing:
#   1. Wet vs Dry season fingerprints side-by-side
#   2. Cluster 1-4 fingerprints in 2×2 grid
#
# APPROACH:
#   - Run separate PCA for each subset (season or cluster)
#   - Use consistent axis scaling across panels for comparison
#   - Combine with patchwork for publication-ready grids
#
# OUTPUT:
#   - 2-panel seasonal comparison grids
#   - 4-panel cluster comparison grids
#   - For concentrations, C-Q slopes, and CVc/CVq
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
fig_dir <- file.path(base_dir, "exploratory_plots", "04_PCA", "4f_seasonal_cluster_grids")
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

# =============================================================================
# HELPER FUNCTION: Run PCA on subset and return plot
# =============================================================================

run_subset_pca <- function(data_wide, site_col = "Stream_Name", 
                           prefix = "", title = "", 
                           arrow_color = "gray40",
                           site_chars = NULL) {
  
  # Extract site names
  site_names <- data_wide[[site_col]]
  
  # Extract numeric matrix
  data_mat <- data_wide %>%
    select(-all_of(site_col)) %>%
    as.matrix()
  rownames(data_mat) <- site_names
  
  # Remove rows with all NA
  row_na_frac <- rowMeans(is.na(data_mat))
  keep_rows <- row_na_frac <= 0.5
  if (sum(keep_rows) < 3) {
    return(list(plot = NULL, scores = NULL, message = "Too few valid sites"))
  }
  
  data_mat <- data_mat[keep_rows, , drop = FALSE]
  site_names <- site_names[keep_rows]
  
  # Remove columns with no variance
  col_var <- apply(data_mat, 2, function(x) var(x, na.rm = TRUE))
  keep_cols <- !is.na(col_var) & col_var > 0
  if (sum(keep_cols) < 2) {
    return(list(plot = NULL, scores = NULL, message = "Too few valid variables"))
  }
  
  data_mat <- data_mat[, keep_cols, drop = FALSE]
  
  # Impute missing values
  for (j in 1:ncol(data_mat)) {
    data_mat[is.na(data_mat[, j]), j] <- median(data_mat[, j], na.rm = TRUE)
  }
  
  # Scale and run PCA
  data_scaled <- scale(data_mat, center = TRUE, scale = TRUE)
  pca_result <- prcomp(data_scaled, center = FALSE, scale. = FALSE)
  
  # Variance explained
  var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100
  
  # Scores
  scores <- as_tibble(pca_result$x[, 1:min(2, ncol(pca_result$x))]) %>%
    mutate(Stream_Name = site_names)
  
  if (!is.null(site_chars)) {
    scores <- scores %>% left_join(site_chars, by = "Stream_Name")
  }
  scores <- scores %>% flag_outlet_stream()
  scores <- apply_factor_orders(scores)
  
  # Loadings
  loadings <- as.data.frame(pca_result$rotation[, 1:min(2, ncol(pca_result$rotation))]) %>%
    rownames_to_column("variable") %>%
    as_tibble() %>%
    mutate(variable = gsub(paste0("^", prefix, "_?"), "", variable))
  
  # Scale loadings
  loading_scale <- max(abs(scores$PC1), abs(scores$PC2)) / 
                   max(abs(loadings$PC1), abs(loadings$PC2)) * 0.7
  loadings_scaled <- loadings %>%
    mutate(
      PC1_plot = PC1 * loading_scale,
      PC2_plot = PC2 * loading_scale
    )
  
  storage_candidates <- c("Q_dS_range_mm", "WB_dS_range_mm")
  storage_available <- purrr::keep(
    storage_candidates,
    ~ .x %in% names(scores) && any(is.finite(scores[[.x]]))
  )
  if (length(storage_available) == 0) {
    storage_available <- NA_character_
  }

  build_panel <- function(storage_col) {
    has_storage <- !is.na(storage_col)
    storage_label <- if (has_storage) get_storage_label(storage_col) else "Storage metric unavailable"
    point_mapping <- if (has_storage) {
      aes(x = PC1, y = PC2, fill = Stream_Name, size = .data[[storage_col]])
    } else {
      aes(x = PC1, y = PC2, fill = Stream_Name)
    }

    panel <- ggplot() +
      stat_ellipse(
        data = scores,
        aes(x = PC1, y = PC2),
        level = 0.95,
        type = "t",
        linetype = "dashed",
        color = "gray50",
        linewidth = 0.4
      ) +
      geom_point(
        data = scores,
        mapping = point_mapping,
        shape = 21,
        colour = "grey20",
        stroke = 0.5,
        alpha = 0.9
      ) +
      geom_text_repel(
        data = scores,
        aes(x = PC1, y = PC2, label = Stream_Name),
        size = 2.5,
        max.overlaps = 15,
        segment.alpha = 0.3
      ) +
      geom_segment(
        data = loadings_scaled,
        aes(x = 0, y = 0, xend = PC1_plot, yend = PC2_plot),
        arrow = arrow(length = unit(0.015, "npc")),
        color = "gray40",
        linewidth = 0.6
      ) +
      geom_text_repel(
        data = loadings_scaled,
        aes(x = PC1_plot, y = PC2_plot, label = variable),
        size = 2.5,
        color = "gray30",
        fontface = "bold",
        max.overlaps = 15
      ) +
      scale_fill_site(name = "Site") +
      labs(
        title = storage_label,
        x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
        y = paste0("PC2 (", round(var_explained[2], 1), "%)")
      ) +
      theme_bw(base_size = 10) +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 10, face = "bold")
      ) +
      coord_fixed()

    if (has_storage) {
      panel <- panel + scale_size_continuous(
        name = storage_label,
        range = c(1.5, 5),
        breaks = scales::pretty_breaks(n = 3)
      )
    } else {
      panel <- panel + guides(size = "none")
    }

    panel
  }

  plot_list <- purrr::map(storage_available, build_panel)
  combined_plot <- if (length(plot_list) == 1) {
    plot_list[[1]]
  } else {
    patchwork::wrap_plots(plotlist = plot_list, nrow = 1)
  }

  combined_plot <- combined_plot + patchwork::plot_annotation(title = title)
  
  return(list(plot = combined_plot, scores = scores, var_explained = var_explained))
}

# =============================================================================
# LOAD DATA
# =============================================================================

message("Loading data...")

site_means <- read_csv(file.path(out_dir, "HJA_clean_site_means.csv"), show_col_types = FALSE) %>%
  filter(!is.na(Stream_Name) & Stream_Name != "",
         !is.na(solute) & solute != "")

cq_master <- read_csv(file.path(out_dir, "HJA_CQ_master.csv"), show_col_types = FALSE) %>%
  filter(!is.na(Stream_Name) & Stream_Name != "")

# Rename columns if needed
if ("variable" %in% names(cq_master) && !"solute" %in% names(cq_master)) {
  cq_master <- cq_master %>% rename(solute = variable)
}
if ("value" %in% names(cq_master) && !"Concentration" %in% names(cq_master)) {
  cq_master <- cq_master %>% rename(Concentration = value)
}

# Get site characteristics
site_chars <- site_means %>%
  select(Stream_Name, any_of(c("RBI", "WB_dS_range_mm", "Q_dS_range_mm"))) %>%
  distinct(Stream_Name, .keep_all = TRUE)

# Check for hydrologic_season column
has_season <- "hydrologic_season" %in% names(site_means) || "hydrologic_season" %in% names(cq_master)
has_cluster <- "Cluster_mode_wy" %in% names(site_means)

message("Has season data: ", has_season)
message("Has cluster data: ", has_cluster)

# =============================================================================
# 1. SEASONAL PCAs (Wet vs Dry)
# =============================================================================

if (has_season) {
  message("\n=== Creating Seasonal PCA Grids ===")
  
  # Identify season column location
  if ("hydrologic_season" %in% names(site_means)) {
    season_data <- site_means
  } else {
    # Need to get season from cq_master and join
    season_info <- cq_master %>%
      select(Stream_Name, solute, hydrologic_season) %>%
      distinct()
    season_data <- site_means %>%
      left_join(season_info, by = c("Stream_Name", "solute"))
  }
  
  # Check what seasons we have
  seasons <- unique(season_data$hydrologic_season)
  seasons <- seasons[!is.na(seasons)]
  message("Seasons found: ", paste(seasons, collapse = ", "))
  
  if (length(seasons) >= 2) {
    
    # ----- C-Q SLOPES BY SEASON -----
    message("Creating C-Q slope seasonal comparison...")
    
    seasonal_plots <- list()
    
    for (szn in c("Wet", "Dry")) {
      if (szn %in% seasons) {
        
        slope_wide <- season_data %>%
          filter(hydrologic_season == szn) %>%
          select(Stream_Name, solute, cq_slope) %>%
          group_by(Stream_Name, solute) %>%
          summarize(cq_slope = mean(cq_slope, na.rm = TRUE), .groups = "drop") %>%
          pivot_wider(
            names_from = solute,
            values_from = cq_slope,
            names_prefix = "slope_"
          )
        
        result <- run_subset_pca(
          slope_wide, 
          site_col = "Stream_Name",
          prefix = "slope",
          title = paste(szn, "Season"),
          arrow_color = ifelse(szn == "Wet", "steelblue", "darkorange"),
          site_chars = site_chars
        )
        
        if (!is.null(result$plot)) {
          seasonal_plots[[szn]] <- result$plot
        }
      }
    }
    
    if (length(seasonal_plots) == 2) {
      p_grid <- seasonal_plots[["Wet"]] + seasonal_plots[["Dry"]] +
        plot_annotation(
          title = "Chemodynamic Fingerprint: Seasonal Comparison (C-Q Slopes)",
          subtitle = "How does chemistry-flow response differ between wet and dry seasons?",
          caption = "Dashed ellipse: 95% CI | Each panel run as separate PCA",
          theme = theme(plot.title = element_text(face = "bold", size = 12))
        )
      
      ggsave(file.path(fig_dir, "01_cqslope_seasonal_grid.png"), p_grid,
             width = 14, height = 7, dpi = 300)
      message("  Saved: 01_cqslope_seasonal_grid.png")
    }
    
    # ----- CVc/CVq BY SEASON -----
    message("Creating CVc/CVq seasonal comparison...")
    
    seasonal_cvcq <- list()
    
    for (szn in c("Wet", "Dry")) {
      if (szn %in% seasons) {
        
        cvcq_wide <- season_data %>%
          filter(hydrologic_season == szn) %>%
          select(Stream_Name, solute, cq_CVc_CVq) %>%
          filter(!is.na(cq_CVc_CVq) & is.finite(cq_CVc_CVq)) %>%
          group_by(Stream_Name, solute) %>%
          summarize(cq_CVc_CVq = mean(cq_CVc_CVq, na.rm = TRUE), .groups = "drop") %>%
          mutate(cq_CVc_CVq = log10(cq_CVc_CVq + 0.01)) %>%
          pivot_wider(
            names_from = solute,
            values_from = cq_CVc_CVq,
            names_prefix = "cvcq_"
          )
        
        result <- run_subset_pca(
          cvcq_wide,
          site_col = "Stream_Name",
          prefix = "cvcq",
          title = paste(szn, "Season"),
          arrow_color = ifelse(szn == "Wet", "purple4", "coral"),
          site_chars = site_chars
        )
        
        if (!is.null(result$plot)) {
          seasonal_cvcq[[szn]] <- result$plot
        }
      }
    }
    
    if (length(seasonal_cvcq) == 2) {
      p_grid <- seasonal_cvcq[["Wet"]] + seasonal_cvcq[["Dry"]] +
        plot_annotation(
          title = "Buffering Fingerprint: Seasonal Comparison (CVc/CVq)",
          subtitle = "Does chemostatic/chemodynamic behavior differ between seasons?",
          caption = "Dashed ellipse: 95% CI | Log10-transformed CVc/CVq",
          theme = theme(plot.title = element_text(face = "bold", size = 12))
        )
      
      ggsave(file.path(fig_dir, "02_cvcq_seasonal_grid.png"), p_grid,
             width = 14, height = 7, dpi = 300)
      message("  Saved: 02_cvcq_seasonal_grid.png")
    }
  }
} else {
  message("No hydrologic_season column found - skipping seasonal PCAs")
}

# =============================================================================
# 2. CLUSTER PCAs (Clusters 1-4)
# =============================================================================

if (has_cluster) {
  message("\n=== Creating Cluster PCA Grids ===")
  
  # Check what clusters we have
  clusters <- sort(unique(site_means$Cluster_mode_wy))
  clusters <- clusters[!is.na(clusters)]
  message("Clusters found: ", paste(clusters, collapse = ", "))
  
  if (length(clusters) >= 2) {
    
    # ----- C-Q SLOPES BY CLUSTER -----
    message("Creating C-Q slope cluster comparison...")
    
    cluster_plots <- list()
    
    for (cl in clusters) {
      
      slope_wide <- site_means %>%
        filter(Cluster_mode_wy == cl) %>%
        select(Stream_Name, solute, cq_slope) %>%
        group_by(Stream_Name, solute) %>%
        summarize(cq_slope = mean(cq_slope, na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(
          names_from = solute,
          values_from = cq_slope,
          names_prefix = "slope_"
        )
      
      # Cluster palette aligned with climatology reference colors
      cluster_colors <- c(
        "1" = "#CFA980",
        "2" = "#98B89F",
        "3" = "#5E8AA1",
        "4" = "#526B8E"
      )
      
      result <- run_subset_pca(
        slope_wide,
        site_col = "Stream_Name",
        prefix = "slope",
        title = paste("Cluster", cl),
        arrow_color = cluster_colors[as.character(cl)],
        site_chars = site_chars
      )
      
      if (!is.null(result$plot)) {
        cluster_plots[[as.character(cl)]] <- result$plot
      }
    }
    
    if (length(cluster_plots) >= 2) {
      # Arrange in grid
      if (length(cluster_plots) == 2) {
        p_grid <- cluster_plots[[1]] + cluster_plots[[2]]
      } else if (length(cluster_plots) == 3) {
        p_grid <- (cluster_plots[[1]] + cluster_plots[[2]]) / 
                  (cluster_plots[[3]] + plot_spacer())
      } else {
        p_grid <- (cluster_plots[[1]] + cluster_plots[[2]]) / 
                  (cluster_plots[[3]] + cluster_plots[[4]])
      }
      
      p_grid <- p_grid +
        plot_annotation(
          title = "Chemodynamic Fingerprint: Cluster Comparison (C-Q Slopes)",
          subtitle = "How does chemistry-flow response vary across hydrologic clusters?",
          caption = "Dashed ellipse: 95% CI | Each panel run as separate PCA",
          theme = theme(plot.title = element_text(face = "bold", size = 12))
        )
      
      ggsave(file.path(fig_dir, "03_cqslope_cluster_grid.png"), p_grid,
             width = 14, height = 12, dpi = 300)
      message("  Saved: 03_cqslope_cluster_grid.png")
    }
    
    # ----- CVc/CVq BY CLUSTER -----
    message("Creating CVc/CVq cluster comparison...")
    
    cluster_cvcq <- list()
    
    for (cl in clusters) {
      
      cvcq_wide <- site_means %>%
        filter(Cluster_mode_wy == cl) %>%
        select(Stream_Name, solute, cq_CVc_CVq) %>%
        filter(!is.na(cq_CVc_CVq) & is.finite(cq_CVc_CVq)) %>%
        group_by(Stream_Name, solute) %>%
        summarize(cq_CVc_CVq = mean(cq_CVc_CVq, na.rm = TRUE), .groups = "drop") %>%
        mutate(cq_CVc_CVq = log10(cq_CVc_CVq + 0.01)) %>%
        pivot_wider(
          names_from = solute,
          values_from = cq_CVc_CVq,
          names_prefix = "cvcq_"
        )
      
      # Cluster palette aligned with climatology reference colors
      cluster_colors <- c(
        "1" = "#CFA980",
        "2" = "#98B89F",
        "3" = "#5E8AA1",
        "4" = "#526B8E"
      )
      
      result <- run_subset_pca(
        cvcq_wide,
        site_col = "Stream_Name",
        prefix = "cvcq",
        title = paste("Cluster", cl),
        arrow_color = cluster_colors[as.character(cl)],
        site_chars = site_chars
      )
      
      if (!is.null(result$plot)) {
        cluster_cvcq[[as.character(cl)]] <- result$plot
      }
    }
    
    if (length(cluster_cvcq) >= 2) {
      if (length(cluster_cvcq) == 2) {
        p_grid <- cluster_cvcq[[1]] + cluster_cvcq[[2]]
      } else if (length(cluster_cvcq) == 3) {
        p_grid <- (cluster_cvcq[[1]] + cluster_cvcq[[2]]) / 
                  (cluster_cvcq[[3]] + plot_spacer())
      } else {
        p_grid <- (cluster_cvcq[[1]] + cluster_cvcq[[2]]) / 
                  (cluster_cvcq[[3]] + cluster_cvcq[[4]])
      }
      
      p_grid <- p_grid +
        plot_annotation(
          title = "Buffering Fingerprint: Cluster Comparison (CVc/CVq)",
          subtitle = "Does chemostatic/chemodynamic behavior vary across hydrologic clusters?",
          caption = "Dashed ellipse: 95% CI | Log10-transformed CVc/CVq",
          theme = theme(plot.title = element_text(face = "bold", size = 12))
        )
      
      ggsave(file.path(fig_dir, "04_cvcq_cluster_grid.png"), p_grid,
             width = 14, height = 12, dpi = 300)
      message("  Saved: 04_cvcq_cluster_grid.png")
    }
  }
} else {
  message("No Cluster_mode_wy column found - skipping cluster PCAs")
}

# =============================================================================
# SUMMARY
# =============================================================================

message("\n✓ Seasonal and Cluster PCA grids complete!")
message("  Figures saved to: ", fig_dir)
