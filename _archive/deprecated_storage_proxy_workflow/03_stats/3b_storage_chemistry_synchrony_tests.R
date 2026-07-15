# =============================================================================
# 3b_storage_chemistry_synchrony_tests.R
# =============================================================================
# Statistical analysis of storage-chemistry-synchrony relationships
# Key insight: Overall GAM R² ~ 0.01, but within-cluster R² ~ 0.70
# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(grid)
  library(patchwork)
})

if (!exists("%||%", inherits = TRUE)) {
  `%||%` <- function(x, y) if (!is.null(x)) x else y
}

# Source helpers and shared configuration
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
try(source(file.path(repo_dir, "00_helpers", "plot_prefs.R")), silent = TRUE)
try(source(file.path(repo_dir, "00_helpers", "workflow_config.R")), silent = TRUE)

paths   <- get_project_paths()
out_dir <- paths$out_dir
fig_dir <- file.path(paths$fig_root, "03_stats", "3b_storage_sync_tests")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# LOAD DATA
# =============================================================================
message("Loading data...")

seasonal <- read_csv(file.path(out_dir, "HJA_clean_seasonal.csv"), show_col_types = FALSE)
site_means <- read_csv(file.path(out_dir, "HJA_rolling_hydro_storage_90d_site_means.csv"), show_col_types = FALSE)
sync_long <- read_csv(file.path(out_dir, "HJA_composite_synchrony.csv"), show_col_types = FALSE)
clusters <- read_csv(file.path(out_dir, "ClusterStreams_allSolutes_byWaterYear.csv"), show_col_types = FALSE)

seasonal <- seasonal %>%
  mutate(
    Stream_Name = factor(Stream_Name, levels = site_order),
    solute = factor(solute, levels = solute_order),
    hydrologic_season = stringr::str_to_title(hydrologic_season)
  )

site_means <- site_means %>%
  mutate(Stream_Name = factor(Stream_Name, levels = site_order))

sync_long <- sync_long %>%
  mutate(
    Stream_Name = factor(Stream_Name, levels = site_order),
    solute = factor(solute, levels = solute_order)
  )

clusters <- clusters %>%
  mutate(Stream_Name = factor(Stream_Name, levels = site_order))

metric_priority <- c(PRIMARY_STORAGE_METRIC, "WB_dS_range_mm", "Q_dS_range_mm")
storage_metrics_to_plot <- unique(metric_priority[metric_priority %in% names(site_means)])

if (length(storage_metrics_to_plot) == 0) {
  stop("No recognized storage metrics found in site_means. Expected WB_dS_range_mm and/or Q_dS_range_mm.")
}

cluster_by_solute <- clusters %>%
  group_by(Stream_Name, chemical) %>%
  summarise(cluster = as.character(round(median(Cluster_mode, na.rm = TRUE))), .groups = "drop") %>%
  mutate(
    Stream_Name = factor(Stream_Name, levels = site_order),
    solute = factor(chemical, levels = solute_order),
    cluster = factor(cluster, levels = cluster_levels)
  ) %>%
  select(Stream_Name, solute, cluster)

message("Cluster distribution (Stream × Solute):")
print(table(cluster_by_solute$cluster, useNA = "ifany"))

storage_by_stream <- function(metric) {
  site_means %>%
    group_by(Stream_Name) %>%
    summarize(
      storage = mean(.data[[metric]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(Stream_Name = factor(Stream_Name, levels = site_order))
}

storage_by_stream_solute <- function(metric) {
  site_means %>%
    group_by(Stream_Name, solute) %>%
    summarize(
      storage = mean(.data[[metric]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Stream_Name = factor(Stream_Name, levels = site_order),
      solute = factor(solute, levels = solute_order)
    )
}

# Summarize synchrony by site  
# =============================================================================
# PLOT 1: Seasonal PCA of CQ slopes (Wet vs Dry side by side)
# =============================================================================
message("Creating seasonal PCA...")

make_seasonal_pca <- function(data, season_name, storage_df, legend_label, caption_label) {
  pca_data <- data %>%
    filter(hydrologic_season == season_name) %>%
    group_by(Stream_Name, solute) %>%
    summarize(cq_slope = mean(cq_slope, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = solute, values_from = cq_slope) %>%
    column_to_rownames("Stream_Name")
  
  # Remove cols with NA
  pca_data <- pca_data[, colSums(!is.na(pca_data)) == nrow(pca_data)]
  if (ncol(pca_data) < 2 || nrow(pca_data) < 3) return(NULL)
  
  pc <- prcomp(scale(pca_data), center = FALSE)
  var_exp <- round((pc$sdev^2 / sum(pc$sdev^2)) * 100, 1)
  
  scores <- tibble(Stream_Name = rownames(pca_data), PC1 = pc$x[,1], PC2 = pc$x[,2]) %>%
    left_join(storage_df, by = "Stream_Name") %>%
    flag_outlet_stream()
  non_outlet <- scores %>% filter(!is_outlet)
  outlet_pts <- scores %>% filter(is_outlet)
  
  loadings <- tibble(solute = colnames(pca_data), PC1 = pc$rotation[,1], PC2 = pc$rotation[,2])
  
  ggplot() +
    geom_point(data = non_outlet, aes(PC1, PC2, color = storage), size = 5) +
    {
      if (nrow(outlet_pts) > 0) {
        geom_point(
          data = outlet_pts,
          aes(PC1, PC2),
          color = "black",
          shape = 21,
          stroke = 1.2,
          size = 5.5,
          fill = "gold"
        )
      }
    } +
    geom_text_repel(data = scores, aes(PC1, PC2, label = Stream_Name), size = 4, fontface = "bold") +
    geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1*3, yend = PC2*3),
                 arrow = arrow(length = unit(0.02, "npc")), color = "grey40", linewidth = 0.7) +
    geom_text(data = loadings, aes(PC1*3.3, PC2*3.3, label = solute), size = 3.5, color = "grey30") +
    scale_color_viridis_c(name = legend_label, option = "plasma", na.value = "grey70") +
    labs(
      title = paste(season_name, "Season: Stream CQ Personalities"),
      subtitle = paste0("PC1 = ", var_exp[1], "%, PC2 = ", var_exp[2], "% variance explained"),
      x = paste0("PC1 (", var_exp[1], "%)"),
      y = paste0("PC2 (", var_exp[2], "%)"),
      caption = paste0("Points colored by ", caption_label, "; GSLOOK highlighted in gold. Arrows = solute CQ loadings.")
    ) +
    guides(color = guide_colorbar(
      title.position = "top",
      barwidth = 20,
      barheight = 0.6
    )) +
    theme_hja() +
    theme(
      legend.title = element_text(face = "bold"),
      legend.key.width = grid::unit(2, "cm"),
      legend.text = element_text(size = BASE_SIZE - 2)
    )
}

for (storage_metric in storage_metrics_to_plot) {
  storage_label_full <- get_storage_label(storage_metric)
  storage_label_short <- get_storage_label(storage_metric, short = TRUE)
  storage_summary_stream <- storage_by_stream(storage_metric)
  metric_suffix <- stringr::str_replace_all(storage_metric, "[^A-Za-z0-9]+", "_")
  display_name <- storage_label_short

  p_wet <- make_seasonal_pca(seasonal, "Wet", storage_summary_stream, storage_label_short, storage_label_full)
  p_dry <- make_seasonal_pca(seasonal, "Dry", storage_summary_stream, storage_label_short, storage_label_full)

  if (!is.null(p_wet) && !is.null(p_dry)) {
    p_combined <- (p_wet + p_dry + plot_layout(guides = "collect")) +
      plot_annotation(title = paste0("Seasonal CQ slope PCA (", display_name, ")"))
    p_combined <- p_combined & theme(legend.position = "bottom")
    file_name <- paste0("01_", metric_suffix, "_seasonal_pca_wet_vs_dry.png")
    save_plot(p_combined, file_name, fig_dir, width = 16, height = 8)
    message("Saved:", file_name)
  } else {
    warning("Skipping PCA composite for ", storage_metric, " due to insufficient data.")
  }
}

# =============================================================================
# PLOT 2: Cluster-stratified storage vs synchrony
# =============================================================================
message("Creating cluster-stratified analysis...")

cluster_stats_all <- list()
metric_r2_summary <- list()

for (storage_metric in storage_metrics_to_plot) {
  storage_label_full <- get_storage_label(storage_metric)
  storage_label_short <- get_storage_label(storage_metric, short = TRUE)
  metric_suffix <- stringr::str_replace_all(storage_metric, "[^A-Za-z0-9]+", "_")
  display_name <- storage_label_short

  storage_solute <- storage_by_stream_solute(storage_metric)

  site_data_metric <- storage_solute %>%
    left_join(sync_long, by = c("Stream_Name", "solute")) %>%
    left_join(cluster_by_solute, by = c("Stream_Name", "solute")) %>%
    drop_na(storage, cluster, conc_sync_allpairs) %>%
    mutate(cluster = droplevels(cluster)) %>%
    flag_outlet_stream()

  message("  ", display_name, ":", nrow(site_data_metric), " site–solute combinations")

  cluster_ids <- levels(droplevels(site_data_metric$cluster))
  if (length(cluster_ids) == 0) {
    warning("No cluster assignments available for ", storage_metric, ". Skipping panel.")
    next
  }

  plot_list <- list()
  cluster_stats_metric <- tibble()

  for (cl in cluster_ids) {
    cl_data <- site_data_metric %>% filter(cluster == cl)
    if (nrow(cl_data) < 3) next

    fit <- tryCatch(lm(conc_sync_allpairs ~ storage, data = cl_data), error = function(e) NULL)
    if (is.null(fit)) next

    summary_fit <- summary(fit)
    r2 <- round(summary_fit$r.squared %||% NA_real_, 3)
    p_val <- tryCatch(summary_fit$coefficients[2, 4], error = function(e) NA_real_)

    storage_range <- range(cl_data$storage, na.rm = TRUE)
    if (!all(is.finite(storage_range)) || storage_range[1] == storage_range[2]) {
      grid <- tibble(storage = cl_data$storage, fit = fitted(fit), lwr = fitted(fit), upr = fitted(fit))
    } else {
      grid <- tibble(storage = seq(storage_range[1], storage_range[2], length.out = 100))
      preds <- predict(fit, newdata = grid, interval = "confidence")
      grid <- bind_cols(grid, as_tibble(preds)) %>%
        rename(fit = fit, lwr = lwr, upr = upr)
    }

    cluster_stats_metric <- bind_rows(
      cluster_stats_metric,
      tibble(cluster = cl, r2 = r2, p = p_val, n = nrow(cl_data))
    )

    non_outlet <- cl_data %>% filter(!is_outlet)
    outlet_pts <- cl_data %>% filter(is_outlet)
    label_data <- cl_data %>%
      group_by(Stream_Name) %>%
      slice_max(abs(conc_sync_allpairs), n = 1, with_ties = FALSE) %>%
      ungroup()

    cluster_color <- cluster_colors[as.character(cl)] %||% "grey60"

    plot_list[[cl]] <- ggplot(cl_data, aes(storage, conc_sync_allpairs)) +
      {if (nrow(grid) > 0) geom_ribbon(data = grid, aes(x = storage, ymin = lwr, ymax = upr), inherit.aes = FALSE,
                                       alpha = 0.18, fill = cluster_color)} +
      {if (nrow(grid) > 0) geom_line(data = grid, aes(x = storage, y = fit), inherit.aes = FALSE,
                                     color = cluster_color, linewidth = 1)} +
      geom_point(data = non_outlet, aes(color = solute), size = 3.2, alpha = 0.85) +
      {
        if (nrow(outlet_pts) > 0) {
          geom_point(
            data = outlet_pts,
            aes(color = solute),
            size = 4,
            shape = OUTLET_SHAPE_TRIANGLE,
            stroke = 1
          )
        }
      } +
      geom_text_repel(
        data = label_data,
        aes(label = Stream_Name),
        size = 3,
        min.segment.length = 0,
        box.padding = 0.25
      ) +
      scale_color_solute() +
      labs(
        title = paste("Cluster", cl),
        subtitle = paste0("R² = ", r2, ", p = ", round(p_val, 3), ", n = ", nrow(cl_data)),
        x = storage_label_full,
        y = get_sync_label("conc_sync_allpairs")
      ) +
      theme_hja()
  }

  if (length(plot_list) == 0) {
    warning("No cluster panels generated for ", storage_metric)
    next
  }

  overall_fit <- tryCatch(lm(conc_sync_allpairs ~ storage, data = site_data_metric), error = function(e) NULL)
  overall_r2 <- if (!is.null(overall_fit)) round(summary(overall_fit)$r.squared %||% NA_real_, 3) else NA_real_
  max_r2 <- suppressWarnings(max(cluster_stats_metric$r2, na.rm = TRUE))
  if (!is.finite(max_r2)) max_r2 <- NA_real_

  p_final <- wrap_plots(plot_list, ncol = 2) +
    plot_annotation(
      title = paste0("Storage vs synchrony by cluster (", display_name, ")"),
      subtitle = paste0(
        "Overall R² = ", overall_r2,
        " • Max within-cluster R² = ", max_r2
      )
    ) +
    plot_layout(guides = "collect")
  p_final <- p_final & theme(legend.position = "bottom")

  file_name <- paste0("02_", metric_suffix, "_cluster_stratified_storage_sync.png")
  save_plot(p_final, file_name, fig_dir, width = 14, height = 10)
  message("Saved:", file_name)

  cluster_stats_all[[storage_metric]] <- cluster_stats_metric %>%
    mutate(storage_metric = storage_metric)
  metric_r2_summary[[storage_metric]] <- list(overall = overall_r2, max = max_r2)
}

if (length(cluster_stats_all) > 0) {
  cluster_stats_path <- file.path(out_dir, "03_stats", "cluster_stratified_stats.csv")
  dir.create(dirname(cluster_stats_path), showWarnings = FALSE, recursive = TRUE)
  cluster_stats_all %>% bind_rows() %>% write_csv(cluster_stats_path)
}

# =============================================================================
# PLOT 3: All 4 clusters side by side with CQ patterns
# =============================================================================
message("Creating cluster CQ pattern comparison...")

cq_by_cluster <- seasonal %>%
  left_join(cluster_by_solute, by = c("Stream_Name", "solute")) %>%
  drop_na(cluster, cq_slope)

p_cq_clusters <- ggplot(cq_by_cluster, aes(x = solute, y = cq_slope, fill = cluster)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~hydrologic_season, ncol = 1) +
  scale_fill_cluster() +
  labs(
    title = "CQ Slopes by Cluster and Season",
    subtitle = "Positive = flushing/mobilization, Negative = dilution, Zero = chemostatic",
    x = "Solute",
    y = "CQ Slope (log-log)",
    caption = "Each box shows distribution of CQ slopes within that cluster."
  ) +
  theme_hja() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p_cq_clusters, "03_cq_slopes_by_cluster_season.png", fig_dir, width = 12, height = 10)
message("Saved: 03_cq_slopes_by_cluster_season.png")

# =============================================================================
# SUMMARY
# =============================================================================
message("\n========== ANALYSIS COMPLETE ==========")
message(paste("Plots saved to:", fig_dir))
message("\nKey finding:")
if (length(metric_r2_summary) > 0) {
  for (metric in names(metric_r2_summary)) {
    label_short <- get_storage_label(metric, short = TRUE)
    vals <- metric_r2_summary[[metric]]
    if (!is.na(vals$overall)) {
      message(paste0("  ", label_short, " overall R² = ", vals$overall))
    }
    if (!is.na(vals$max)) {
      message(paste0("  ", label_short, " max within-cluster R² = ", vals$max))
    }
  }
} else {
  message("  Storage–synchrony model fits not available (insufficient data).")
}
