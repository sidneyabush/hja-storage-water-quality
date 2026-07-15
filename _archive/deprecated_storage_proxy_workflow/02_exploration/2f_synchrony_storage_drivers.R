# =============================================================================
# 6c_synchrony_storage_drivers.R
# =============================================================================
# WHAT DRIVES SYNCHRONY? Storage divergence and CQ regime analysis
#
# RESEARCH QUESTIONS:
#   1. When sites have divergent storage dynamics, does synchrony decrease?
#   2. Are sites in different CQ regimes also asynchronous?
#   3. Do these patterns differ seasonally vs annually?
#
# APPROACH:
#   Section 1: Storage divergence hypothesis
#     - Calculate inter-site storage SD (by year, by season)
#     - Relate to mean synchrony across that year/season
#     - GAM to detect nonlinearity
#
#   Section 2: CQ regime hypothesis
#     - Classify sites by CQ regime (mobilizing/diluting/chemostatic)
#     - Compare synchrony across regimes
#     - Question: do mixed-regime site-pairs show lower synchrony?
# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  quietly_loaded_mgcv <- requireNamespace("mgcv", quietly = TRUE)
})

# Provide safe fallback operator for script discovery
if (!exists("%||%", inherits = TRUE)) {
  `%||%` <- function(x, y) if (!is.null(x)) x else y
}

# Source helpers from repo root
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
try(source(file.path(repo_dir, "00_helpers", "plot_prefs.R")), silent = TRUE)
try(source(file.path(repo_dir, "00_helpers", "workflow_config.R")), silent = TRUE)

# =============================================================================
# SETUP
# =============================================================================

base_dir    <- "/Users/sidneybush/Library/CloudStorage/Box-Box"
project_dir <- file.path(base_dir, "Sidney_Bush", "HJA_Water_Quality")
out_dir     <- file.path(project_dir, "outputs")

fig_dir     <- file.path(project_dir, "exploratory_plots", "02_exploration", "2f_synchrony_drivers")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, filename, width = 8, height = 6) {
  ggplot2::ggsave(
    file.path(fig_dir, filename),
    plot = p, width = width, height = height, dpi = 300
  )
}

theme_sync <- theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    legend.position = "right"
  )

compute_lm_stats <- function(formula, data) {
  fit <- tryCatch(stats::lm(formula, data = data), error = function(e) NULL)
  if (is.null(fit)) {
    return(NULL)
  }
  smry <- tryCatch(summary(fit), error = function(e) NULL)
  if (is.null(smry)) {
    return(NULL)
  }
  coef_table <- smry$coefficients
  p_value <- if (!is.null(coef_table) && nrow(coef_table) >= 2 && ncol(coef_table) >= 4) {
    coef_table[2, 4]
  } else {
    NA_real_
  }
  list(r2 = smry$r.squared, p = p_value)
}

format_lm_label <- function(stats) {
  if (is.null(stats) || is.null(stats$r2) || is.null(stats$p)) {
    return(NULL)
  }
  sprintf("R² = %.2f (p = %.3f)", stats$r2, stats$p)
}

annotate_lm <- function(p, data, x_var, y_var, stats) {
  if (is.null(stats)) {
    return(p)
  }
  label <- format_lm_label(stats)
  if (is.null(label) || !all(c(x_var, y_var) %in% names(data))) {
    return(p)
  }
  x_vals <- data[[x_var]]
  y_vals <- data[[y_var]]
  if (!any(is.finite(x_vals)) || !any(is.finite(y_vals))) {
    return(p)
  }
  x_pos <- max(x_vals, na.rm = TRUE)
  y_pos <- max(y_vals, na.rm = TRUE)
  p + annotate(
    "text",
    x = x_pos,
    y = y_pos,
    label = label,
    hjust = 1.05,
    vjust = 1.1,
    fontface = "italic"
  )
}

# =============================================================================
# LOAD DATA
# =============================================================================

message("Loading data...")

seasonal <- readr::read_csv(
  file.path(out_dir, "HJA_clean_seasonal.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    hydrologic_season = factor(hydrologic_season, levels = c("Wet", "Dry")),
    Stream_Name = as.factor(Stream_Name),
    solute = as.factor(solute)
  ) %>%
  apply_factor_orders()

annual <- readr::read_csv(
  file.path(out_dir, "HJA_clean_annual.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    Stream_Name = as.factor(Stream_Name),
    solute = as.factor(solute)
  ) %>%
  apply_factor_orders()

sync_annual <- readr::read_csv(
  file.path(out_dir, "HJA_composite_synchrony_annual.csv"),
  show_col_types = FALSE
)

# Load multivariate synchrony (site-year level, not solute level)
sync_multivariate <- readr::read_csv(
  file.path(out_dir, "HJA_composite_synchrony_multivariate.csv"),
  show_col_types = FALSE
) %>%
  select(Stream_Name, water_year,
         concentration_mv_geogenic_sync_allpairs,
         concentration_mv_biogenic_sync_allpairs,
         concentration_mv_nutrient_sync_allpairs)

modal_clusters <- readr::read_csv(
  file.path(out_dir, "ClusterStreams_allSolutes_modalClusters.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    Cluster_mode = factor(Cluster_mode, levels = cluster_levels)
  ) %>%
  apply_factor_orders()

annual_joined <- annual %>%
  left_join(sync_annual, by = c("Stream_Name", "solute", "water_year")) %>%
  left_join(
    select(modal_clusters, Stream_Name, chemical, Cluster_mode),
    by = c("Stream_Name", "solute" = "chemical")
  ) %>%
  rename(Cluster = Cluster_mode) %>%
  mutate(
    Stream_Name = as.factor(Stream_Name),
    solute = as.factor(solute),
    solute_type = factor(categorize_solute(solute), levels = c("Geogenic", "Biogenic", "Nutrient")),
    Cluster = factor(Cluster, levels = cluster_levels)
  ) %>%
  apply_factor_orders()

primary_storage_metric <- intersect(PRIMARY_STORAGE_METRIC, names(annual_joined))
label_metric <- if (length(primary_storage_metric) == 1) primary_storage_metric else PRIMARY_STORAGE_METRIC
primary_label <- get_label(label_metric)
annual_divergence <- tibble()

# =============================================================================
# SECTION 1: STORAGE DIVERGENCE VS SYNCHRONY
# =============================================================================

message("\n=== STORAGE DIVERGENCE HYPOTHESIS ===\n")

# 1a. Annual: inter-site storage divergence vs synchrony

message("Computing annual storage divergence...")

if (length(primary_storage_metric) == 1 &&
    all(c("Stream_Name", "water_year", primary_storage_metric) %in% names(annual_joined)) &&
    any(c("conc_sync_allpairs", "cqslope_sync_allpairs") %in% names(annual_joined))) {
  metric <- primary_storage_metric[[1]]

  has_conc <- "conc_sync_allpairs" %in% names(annual_joined)
  has_cqslope <- "cqslope_sync_allpairs" %in% names(annual_joined)

  annual_divergence <- annual_joined %>%
    group_by(water_year) %>%
    summarise(
      storage_sd = sd(.data[[metric]], na.rm = TRUE),
      conc_sync_mean = if (has_conc) mean(conc_sync_allpairs, na.rm = TRUE) else NA_real_,
      cqslope_sync_mean = if (has_cqslope) mean(cqslope_sync_allpairs, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) %>%
    filter(!is.na(storage_sd))
  
  # Plot 1a: Conc synchrony vs storage divergence
  if (any(!is.na(annual_divergence$conc_sync_mean))) {
    df_plot <- annual_divergence %>% filter(!is.na(conc_sync_mean))
    lm_stats <- compute_lm_stats(conc_sync_mean ~ storage_sd, df_plot)
    lm_label <- format_lm_label(lm_stats)
    subtitle_text <- if (!is.null(lm_label)) {
      paste("Each point = one water year", lm_label, sep = "\n")
    } else {
      "Each point = one water year"
    }
    
    p <- df_plot %>%
      ggplot(aes(x = storage_sd, y = conc_sync_mean)) +
      geom_point(alpha = 0.8, size = 3) +
      geom_smooth(method = "lm", se = TRUE, color = "darkblue", fill = "lightblue") +
      labs(
        x = paste0("SD of ", primary_label, " across sites (annual)"),
        y = get_sync_label("conc_sync_allpairs"),
      ) +
      theme_sync

    p <- annotate_lm(p, df_plot, "storage_sd", "conc_sync_mean", lm_stats)
    
    save_plot(p, "01_annual_storage_divergence_conc_sync.png", width = 8, height = 6)
    
    corr <- cor(df_plot$storage_sd, df_plot$conc_sync_mean, use = "complete.obs")
    message("  Conc sync vs storage SD: r = ", round(corr, 3))
  }
  
  # Plot 1b: CQ-slope synchrony vs storage divergence
  if (any(!is.na(annual_divergence$cqslope_sync_mean))) {
    df_plot <- annual_divergence %>% filter(!is.na(cqslope_sync_mean))
    lm_stats <- compute_lm_stats(cqslope_sync_mean ~ storage_sd, df_plot)
    lm_label <- format_lm_label(lm_stats)
    subtitle_text <- if (!is.null(lm_label)) {
      paste("Each point = one water year", lm_label, sep = "\n")
    } else {
      "Each point = one water year"
    }
    
    p <- df_plot %>%
      ggplot(aes(x = storage_sd, y = cqslope_sync_mean)) +
      geom_point(alpha = 0.8, size = 3) +
      geom_smooth(method = "lm", se = TRUE, color = "darkblue", fill = "lightblue") +
      labs(
        x = paste0("SD of ", primary_label, " across sites (annual)"),
        y = get_sync_label("cqslope_sync_allpairs"),
      ) +
      theme_sync

    p <- annotate_lm(p, df_plot, "storage_sd", "cqslope_sync_mean", lm_stats)
    
    save_plot(p, "02_annual_storage_divergence_cqslope_sync.png", width = 8, height = 6)
    
    corr <- cor(df_plot$storage_sd, df_plot$cqslope_sync_mean, use = "complete.obs")
    message("  CQ-slope sync vs storage SD: r = ", round(corr, 3))
  }

  if (length(primary_storage_metric) == 1 &&
      "solute_type" %in% names(annual_joined) &&
      "conc_sync_allpairs" %in% names(annual_joined)) {
    metric <- primary_storage_metric[[1]]

    annual_solute <- annual_joined %>%
      filter(!is.na(solute_type)) %>%
      group_by(solute_type, water_year) %>%
      summarise(
        storage_sd = sd(.data[[metric]], na.rm = TRUE),
        conc_sync_mean = mean(conc_sync_allpairs, na.rm = TRUE),
        n_obs = sum(is.finite(.data[[metric]]) & is.finite(conc_sync_allpairs)),
        .groups = "drop"
      ) %>%
      filter(n_obs >= 2, is.finite(storage_sd), is.finite(conc_sync_mean))

    valid_solute <- annual_solute %>%
      group_by(solute_type) %>%
      summarise(n_points = n(), .groups = "drop") %>%
      filter(n_points >= 2)

    if (nrow(valid_solute) > 0) {
      annual_solute_filtered <- annual_solute %>%
        semi_join(valid_solute, by = "solute_type") %>%
        mutate(solute_type = droplevels(solute_type))

      type_stats <- annual_solute_filtered %>%
        group_by(solute_type) %>%
        summarise(
          lm_stats = list(compute_lm_stats(conc_sync_mean ~ storage_sd, cur_data())),
          .groups = "drop"
        ) %>%
        mutate(
          label_raw = vapply(
            lm_stats,
            function(stat) {
              label <- format_lm_label(stat)
              if (is.null(label)) "Fit unavailable" else label
            },
            character(1)
          ),
          subtitle_entry = paste(solute_type, label_raw, sep = " – ")
        )

      solute_levels <- levels(droplevels(annual_solute_filtered$solute_type))
      palette_vals <- solute_type_colors[solute_levels]
      subtitle_entries <- type_stats %>%
        filter(label_raw != "Fit unavailable") %>%
        pull(subtitle_entry)

      subtitle_text <- if (length(subtitle_entries) > 0) {
        paste(subtitle_entries, collapse = " | ")
      } else {
        "Linear fits unavailable (>= 2 obs required per group)"
      }

      annotation_df <- annual_solute_filtered %>%
        group_by(solute_type) %>%
        summarise(
          x = max(storage_sd, na.rm = TRUE),
          y = max(conc_sync_mean, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        left_join(select(type_stats, solute_type, annotation = label_raw), by = "solute_type") %>%
        filter(!is.na(annotation), annotation != "Fit unavailable")

      p <- annual_solute_filtered %>%
        ggplot(aes(x = storage_sd, y = conc_sync_mean, color = solute_type)) +
        geom_point(alpha = 0.8, size = 2.8) +
        geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
        facet_wrap(~ solute_type) +
        scale_color_manual(values = palette_vals) +
        labs(
          x = paste0("SD of ", primary_label, " across sites (annual)"),
          y = get_sync_label("conc_sync_allpairs"),
          caption = paste0(
            "n = ", nrow(annual_solute_filtered),
            " | Synchrony metric = conc_sync_allpairs (Abbott)"
          )
        ) +
        theme_sync +
        theme(legend.position = "none") +
        geom_text(
          data = annotation_df,
          aes(x = x, y = y, label = annotation, color = solute_type),
          inherit.aes = FALSE,
          hjust = 1.03,
          vjust = 1.1,
          fontface = "italic",
          show.legend = FALSE
        )

      save_plot(p, "01b_annual_storage_divergence_conc_sync_by_solute_type.png", width = 9, height = 7)

      message("  Annual solute-type faceted plot saved (n groups = ", nrow(type_stats), ")")
    }
  }

  if (length(primary_storage_metric) == 1 &&
      "Cluster" %in% names(annual_joined) &&
      "conc_sync_allpairs" %in% names(annual_joined)) {
    metric <- primary_storage_metric[[1]]

    annual_cluster <- annual_joined %>%
      filter(!is.na(Cluster)) %>%
      group_by(Cluster, water_year) %>%
      summarise(
        storage_sd = sd(.data[[metric]], na.rm = TRUE),
        conc_sync_mean = mean(conc_sync_allpairs, na.rm = TRUE),
        n_obs = sum(is.finite(.data[[metric]]) & is.finite(conc_sync_allpairs)),
        .groups = "drop"
      ) %>%
      filter(n_obs >= 2, is.finite(storage_sd), is.finite(conc_sync_mean))

    valid_cluster <- annual_cluster %>%
      group_by(Cluster) %>%
      summarise(n_points = n(), .groups = "drop") %>%
      filter(n_points >= 2)

    if (nrow(valid_cluster) > 0) {
      annual_cluster_filtered <- annual_cluster %>%
        semi_join(valid_cluster, by = "Cluster") %>%
        mutate(Cluster = droplevels(Cluster))

      cluster_stats <- annual_cluster_filtered %>%
        group_by(Cluster) %>%
        summarise(
          lm_stats = list(compute_lm_stats(conc_sync_mean ~ storage_sd, cur_data())),
          .groups = "drop"
        ) %>%
        mutate(
          label_raw = vapply(
            lm_stats,
            function(stat) {
              label <- format_lm_label(stat)
              if (is.null(label)) "Fit unavailable" else label
            },
            character(1)
          ),
          subtitle_entry = paste("Cluster", Cluster, label_raw, sep = " – ")
        )

      subtitle_entries <- cluster_stats %>%
        filter(label_raw != "Fit unavailable") %>%
        pull(subtitle_entry)

      subtitle_text <- if (length(subtitle_entries) > 0) {
        paste(subtitle_entries, collapse = " | ")
      } else {
        "Linear fits unavailable (>= 2 obs required per group)"
      }

      annotation_df <- annual_cluster_filtered %>%
        group_by(Cluster) %>%
        summarise(
          x = max(storage_sd, na.rm = TRUE),
          y = max(conc_sync_mean, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        left_join(select(cluster_stats, Cluster, annotation = label_raw), by = "Cluster") %>%
        filter(!is.na(annotation), annotation != "Fit unavailable")

      p <- annual_cluster_filtered %>%
        ggplot(aes(x = storage_sd, y = conc_sync_mean, color = Cluster)) +
        geom_point(alpha = 0.8, size = 2.8) +
        geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
        facet_wrap(~ Cluster) +
        scale_color_cluster() +
        labs(
          x = paste0("SD of ", primary_label, " across sites (annual)"),
          y = get_sync_label("conc_sync_allpairs"),
          caption = paste0(
            "n = ", nrow(annual_cluster_filtered),
            " | Synchrony metric = conc_sync_allpairs (Abbott)"
          )
        ) +
        theme_sync +
        theme(legend.position = "none") +
        geom_text(
          data = annotation_df,
          aes(x = x, y = y, label = annotation, color = Cluster),
          inherit.aes = FALSE,
          hjust = 1.03,
          vjust = 1.1,
          fontface = "italic",
          show.legend = FALSE
        )

      save_plot(p, "01c_annual_storage_divergence_conc_sync_by_cluster.png", width = 9, height = 7)

      message("  Annual cluster faceted plot saved (n groups = ", nrow(cluster_stats), ")")
    }
  }

  # =============================================================================
  # SECTION 1d: UNIVARIATE VS MULTIVARIATE SYNCHRONY RESPONSES TO STORAGE
  # =============================================================================

  message("\n=== UNIVARIATE VS MULTIVARIATE SYNCHRONY ===\n")

  # Calculate univariate synchrony (averaged across solutes) for each group
  univariate_by_group <- annual_joined %>%
    filter(!is.na(solute_type)) %>%
    group_by(Stream_Name, water_year, solute_type) %>%
    summarise(
      conc_sync_univariate = mean(conc_sync_allpairs, na.rm = TRUE),
      storage_val = first(.data[[metric]]),
      n_solutes = n(),
      .groups = "drop"
    ) %>%
    filter(!is.na(conc_sync_univariate), !is.na(storage_val))

  # Calculate storage divergence and mean synchrony by group
  storage_div_univariate <- univariate_by_group %>%
    group_by(water_year, solute_type) %>%
    summarise(
      storage_sd = sd(storage_val, na.rm = TRUE),
      sync_mean = mean(conc_sync_univariate, na.rm = TRUE),
      n_obs = n(),
      .groups = "drop"
    ) %>%
    filter(is.finite(storage_sd), is.finite(sync_mean), n_obs >= 2) %>%
    mutate(method = "Univariate\n(averaged)")

  # Now for multivariate: Calculate storage divergence with multivariate synchrony
  multivariate_by_group <- annual_joined %>%
    select(Stream_Name, water_year, all_of(metric)) %>%
    distinct() %>%
    left_join(sync_multivariate, by = c("Stream_Name", "water_year")) %>%
    filter(!is.na(.data[[metric]]))

  # Pivot multivariate synchrony to long format for easier processing
  storage_div_multivariate <- multivariate_by_group %>%
    group_by(water_year) %>%
    summarise(
      storage_sd = sd(.data[[metric]], na.rm = TRUE),
      sync_geogenic = first(concentration_mv_geogenic_sync_allpairs),
      sync_biogenic = first(concentration_mv_biogenic_sync_allpairs),
      sync_nutrient = first(concentration_mv_nutrient_sync_allpairs),
      n_obs = n(),
      .groups = "drop"
    ) %>%
    filter(n_obs >= 2, is.finite(storage_sd)) %>%
    tidyr::pivot_longer(
      cols = starts_with("sync_"),
      names_to = "solute_type",
      values_to = "sync_mean",
      names_prefix = "sync_"
    ) %>%
    mutate(
      solute_type = factor(
        tools::toTitleCase(solute_type),
        levels = c("Geogenic", "Biogenic", "Nutrient")
      ),
      method = "Multivariate\n(all together)"
    ) %>%
    filter(!is.na(sync_mean))

  # Combine univariate and multivariate
  storage_div_combined <- bind_rows(
    storage_div_univariate,
    storage_div_multivariate
  ) %>%
    mutate(method = factor(method, levels = c("Univariate\n(averaged)", "Multivariate\n(all together)")))

  # Filter to groups with sufficient data
  valid_groups_combined <- storage_div_combined %>%
    group_by(solute_type, method) %>%
    summarise(n_points = n(), .groups = "drop") %>%
    filter(n_points >= 3)

  if (nrow(valid_groups_combined) > 0) {
    storage_div_plot <- storage_div_combined %>%
      semi_join(valid_groups_combined, by = c("solute_type", "method"))

    # Compute stats for each solute_type × method combination
    stats_combined <- storage_div_plot %>%
      group_by(solute_type, method) %>%
      summarise(
        lm_stats = list(compute_lm_stats(sync_mean ~ storage_sd, cur_data())),
        .groups = "drop"
      ) %>%
      mutate(
        label_raw = vapply(
          lm_stats,
          function(stat) {
            label <- format_lm_label(stat)
            if (is.null(label)) "" else label
          },
          character(1)
        )
      ) %>%
      filter(label_raw != "")

    # Create faceted plot
    solute_levels <- levels(droplevels(storage_div_plot$solute_type))
    palette_vals <- solute_type_colors[solute_levels]

    # Annotation positions
    annotation_df <- storage_div_plot %>%
      group_by(solute_type, method) %>%
      summarise(
        x = max(storage_sd, na.rm = TRUE),
        y = max(sync_mean, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      left_join(select(stats_combined, solute_type, method, annotation = label_raw),
                by = c("solute_type", "method")) %>%
      filter(!is.na(annotation), annotation != "")

    p <- storage_div_plot %>%
      ggplot(aes(x = storage_sd, y = sync_mean, color = solute_type, linetype = method)) +
      geom_point(alpha = 0.7, size = 2.5) +
      geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 0.8) +
      facet_wrap(~ solute_type, scales = "free_y") +
      scale_color_manual(values = palette_vals) +
      scale_linetype_manual(values = c("Univariate\n(averaged)" = "dashed",
                                       "Multivariate\n(all together)" = "solid")) +
      labs(
        title = "Storage Divergence Effects: Univariate vs Multivariate Synchrony",
        x = paste0("SD of ", primary_label, " across sites (annual)"),
        y = "Abbott Synchrony (concentration)",
        color = "Solute Group",
        linetype = "Synchrony Method",
        caption = paste0(
          "Each point = one water year for that solute group. ",
          "Univariate = average of individual solute synchronies. ",
          "Multivariate = synchrony using all solutes together as a chemical profile."
        )
      ) +
      theme_sync +
      theme(legend.position = "right")

    # Add annotations if available
    if (nrow(annotation_df) > 0) {
      p <- p + geom_text(
        data = annotation_df,
        aes(x = x, y = y, label = annotation, color = solute_type),
        inherit.aes = FALSE,
        hjust = 1.05,
        vjust = 1.1,
        size = 2.8,
        fontface = "italic",
        show.legend = FALSE
      )
    }

    save_plot(p, "01d_storage_divergence_univariate_vs_multivariate.png", width = 11, height = 7)

    message("  Univariate vs multivariate comparison plot saved")

    # Print correlation between methods for each group
    for (stype in solute_levels) {
      uni_data <- storage_div_univariate %>%
        filter(solute_type == stype) %>%
        arrange(water_year)

      multi_data <- storage_div_multivariate %>%
        filter(solute_type == stype) %>%
        arrange(water_year)

      if (nrow(uni_data) > 0 && nrow(multi_data) > 0) {
        # Match by water_year
        merged <- inner_join(
          select(uni_data, water_year, storage_sd, sync_uni = sync_mean),
          select(multi_data, water_year, sync_multi = sync_mean),
          by = "water_year"
        )

        if (nrow(merged) >= 2) {
          corr <- cor(merged$sync_uni, merged$sync_multi, use = "complete.obs")
          message(sprintf("  %s: r(univariate, multivariate) = %.3f (n=%d water years)",
                         stype, corr, nrow(merged)))
        }
      }
    }
  } else {
    message("  Insufficient data for univariate vs multivariate comparison")
  }

} else if (length(primary_storage_metric) == 0) {
  warning("Primary storage metric not found in annual dataset; skipping annual divergence analysis")
}

# 1b. Seasonal: inter-site storage divergence vs synchrony

message("\nComputing seasonal storage divergence...")

if (length(primary_storage_metric) == 1 &&
    all(c("Stream_Name", "water_year", "hydrologic_season", primary_storage_metric) %in% names(seasonal)) &&
    "cq_slope" %in% names(seasonal)) {
  metric <- primary_storage_metric[[1]]
  
  seasonal_divergence <- seasonal %>%
    group_by(water_year, hydrologic_season) %>%
    summarise(
      storage_sd = sd(.data[[metric]], na.rm = TRUE),
      cq_sync_mean = mean(cq_slope, na.rm = TRUE),  # Mean within-season CQ
      .groups = "drop"
    ) %>%
    filter(!is.na(storage_sd))
  
  if (nrow(seasonal_divergence) > 4) {
    df_plot <- seasonal_divergence %>%
      filter(is.finite(storage_sd), is.finite(cq_sync_mean))
    
    if (nrow(df_plot) > 0) {
      season_labels <- split(df_plot, df_plot$hydrologic_season)
      season_labels <- lapply(season_labels, function(season_df) {
        if (nrow(season_df) < 2) {
          return(NA_character_)
        }
        stats <- compute_lm_stats(cq_sync_mean ~ storage_sd, season_df)
        if (is.null(stats) || !is.finite(stats$r2)) {
          return(NA_character_)
        }
        season_name <- as.character(season_df$hydrologic_season[1])
        sprintf("%s – R² = %.2f", season_name, stats$r2)
      })
      season_labels <- Filter(function(label) !is.null(label) && !is.na(label), season_labels)
      season_labels <- unlist(season_labels, use.names = FALSE)
      
      subtitle_text <- "Linear fits by hydrologic season (95% CI ribbons)"
      if (length(season_labels) > 0) {
        subtitle_text <- paste(subtitle_text, paste(season_labels, collapse = " | "), sep = " | ")
      }

      caption_text <- paste(
        "Season-year points:", nrow(df_plot),
        "| Synchrony metric = Seasonal mean CQ slope (cq_slope)",
        paste0("| Storage metric = SD of ", primary_label, " across sites")
      )
      
      p <- df_plot %>%
        ggplot(aes(x = storage_sd, y = cq_sync_mean, color = hydrologic_season)) +
        geom_point(alpha = 0.7, size = 3) +
        geom_smooth(method = "lm", se = TRUE, alpha = 0.15) +
        scale_color_season() +
        labs(
          x = paste0("SD of ", primary_label, " across sites (by season)"),
          y = "Mean CQ slope (within season)",
          color = "Season",
          title = "Storage Divergence vs CQ Slope by Season",
          subtitle = subtitle_text,
          caption = caption_text
        )
      
      save_plot(p, "03_seasonal_storage_divergence_cq.png", width = 9, height = 7)
      
      message("  Seasonal storage divergence plot created")
    }
  }
} else if (length(primary_storage_metric) == 0) {
  warning("Primary storage metric not found in seasonal dataset; skipping seasonal divergence analysis")
}

# =============================================================================
# SECTION 2: CQ REGIME & SYNCHRONY
# =============================================================================

message("\n=== CQ REGIME & SYNCHRONY ===\n")

if ("cq_slope" %in% names(annual_joined)) {
  
  # Classify each observation by CQ regime
  annual_reg <- annual_joined %>%
    mutate(
      cq_regime = case_when(
        cq_slope > 0.1 ~ "Mobilizing",
        cq_slope < -0.1 ~ "Diluting",
        is.finite(cq_slope) ~ "Chemostatic",
        TRUE ~ NA_character_
      ),
      cq_regime = factor(cq_regime, levels = c("Chemostatic", "Mobilizing", "Diluting"))
    )
  
  # Plot 2a: Concentration synchrony by regime
  if ("conc_sync_allpairs" %in% names(annual_reg)) {
    df_reg <- annual_reg %>% filter(!is.na(cq_regime), !is.na(conc_sync_allpairs))
    
    if (nrow(df_reg) > 0) {
      p <- df_reg %>%
        ggplot(aes(x = cq_regime, y = conc_sync_allpairs, fill = cq_regime)) +
        geom_boxplot(outlier.alpha = 0.4, alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
        scale_fill_manual(values = c("Chemostatic" = "grey60", "Mobilizing" = "green", "Diluting" = "blue")) +
        labs(
          x = "C–Q regime (annual)",
          y = get_sync_label("conc_sync_allpairs"),
          fill = "CQ Regime",
          title = "Concentration Synchrony by CQ Regime"
        )

      save_plot(p, "04_annual_sync_by_cq_regime_conc.png", width = 8, height = 6)
      
      message("  Conc sync by regime: n = ", nrow(df_reg))
    }
  }
  
  # Plot 2b: CQ-slope synchrony by regime
  if ("cqslope_sync_allpairs" %in% names(annual_reg)) {
    df_reg <- annual_reg %>% filter(!is.na(cq_regime), !is.na(cqslope_sync_allpairs))
    
    if (nrow(df_reg) > 0) {
      p <- df_reg %>%
        ggplot(aes(x = cq_regime, y = cqslope_sync_allpairs, fill = cq_regime)) +
        geom_boxplot(outlier.alpha = 0.4, alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
        scale_fill_manual(values = c("Chemostatic" = "grey60", "Mobilizing" = "green", "Diluting" = "blue")) +
        labs(
          x = "C–Q regime (annual)",
          y = get_sync_label("cqslope_sync_allpairs"),
          fill = "CQ Regime",
          title = "CQ Slope Synchrony by CQ Regime"
        )

      save_plot(p, "05_annual_sync_by_cq_regime_cqslope.png", width = 8, height = 6)
      
      message("  CQ-slope sync by regime: n = ", nrow(df_reg))
    }
  }
}

# =============================================================================
# SECTION 3: GAM - NONLINEAR STORAGE–SYNCHRONY RELATIONSHIP (Optional)
# =============================================================================

message("\n=== NONLINEAR RELATIONSHIP (GAM) ===\n")

if (quietly_loaded_mgcv && nrow(annual_divergence) >= 6) {
  
  message("Fitting GAM for nonlinear storage–synchrony relationship...")
  
  # Use annual data
  df_gam <- annual_divergence %>%
    filter(!is.na(storage_sd), !is.na(conc_sync_mean))
  
  if (nrow(df_gam) >= 6) {
    
    mod <- mgcv::gam(
      conc_sync_mean ~ s(storage_sd, k = 4),
      data = df_gam
    )
    
    # Prediction grid
    newd <- tibble(
      storage_sd = seq(min(df_gam$storage_sd), max(df_gam$storage_sd), length.out = 100)
    )
    
    pred <- predict(mod, newdata = newd, se.fit = TRUE, type = "response")
    newd$fit <- pred$fit
    newd$se <- pred$se.fit
    
    # Plot
    p <- ggplot() +
      geom_point(
        data = df_gam,
        aes(x = storage_sd, y = conc_sync_mean),
        alpha = 0.8, size = 3
      ) +
      geom_line(
        data = newd,
        aes(x = storage_sd, y = fit),
        color = "darkblue", linewidth = 1.2
      ) +
      geom_ribbon(
        data = newd,
        aes(x = storage_sd, ymin = fit - 2*se, ymax = fit + 2*se),
        fill = "lightblue", alpha = 0.3
      ) +
      labs(
        x = paste0("SD of ", primary_label, " across sites (annual)"),
        y = get_sync_label("conc_sync_allpairs"),
        title = "GAM: Nonlinear Storage-Synchrony Relationship"
      )

    save_plot(p, "06_gam_storage_synchrony.png", width = 8, height = 6)

    message("  GAM fitted successfully")
    print(summary(mod))
  }
}

message("\n=== SYNCHRONY DRIVERS COMPLETE ===\n")
