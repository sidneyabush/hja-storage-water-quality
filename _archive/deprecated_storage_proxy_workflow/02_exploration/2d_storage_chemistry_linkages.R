# =============================================================================
# 6a_storage_chemistry_linkages.R
# =============================================================================
# CORE ANALYSIS: How does storage similarity relate to chemistry similarity?
#
# RESEARCH QUESTION:
#   Do sites with similar dynamic storage ranges show similar CQ behavior?
#   Does this differ seasonally?
#
# APPROACH:
#   For each solute × water_year (annual) and season (seasonal):
#     - Compare all site pairs
#     - Calculate: pairwise storage distance (multi-metric)
#     - Calculate: pairwise CQ slope difference
#     - Examine relationship: storage distance vs CQ difference
#
# OUTPUT:
#   - Scatter plots: pairwise relationships
#   - Comparison: annual vs seasonal patterns
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})
if (!exists("%||%", inherits = TRUE)) {
  `%||%` <- function(x, y) if (!is.null(x)) x else y
}

# Source helpers from repo root
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
try(source(file.path(repo_dir, "00_helpers", "plot_prefs.R")), silent = TRUE)
try(source(file.path(repo_dir, "00_helpers", "workflow_config.R")), silent = TRUE)
try(source(file.path(repo_dir, "00_helpers", "plot_theme_set.R")), silent = TRUE)

# =============================================================================
# SETUP
# =============================================================================

base_dir    <- "/Users/sidneybush/Library/CloudStorage/Box-Box"
project_dir <- file.path(base_dir, "Sidney_Bush", "HJA_Water_Quality")
out_dir     <- file.path(project_dir, "outputs")

fig_dir     <- file.path(project_dir, "exploratory_plots", "02_exploration", "2d_storage_chemistry")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

supp_fig_dir <- file.path(fig_dir, "supplemental_storage")
dir.create(supp_fig_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, filename, width = 9, height = 7, supplemental = FALSE) {
  target_dir <- if (supplemental) supp_fig_dir else fig_dir
  ggplot2::ggsave(
    file.path(target_dir, filename),
    plot = p, width = width, height = height, dpi = 300
  )
}

# Build pairwise storage vs CQ difference table for a given metric
make_pairwise_table <- function(df, metric, group_vars = c("solute", "water_year")) {
  if (!is.data.frame(df) || !all(c("Stream_Name", metric, "cq_slope", group_vars) %in% names(df))) {
    return(tibble())
  }
  df_sub <- df %>%
    filter(!is.na(.data[[metric]]), !is.na(cq_slope))
  if (nrow(df_sub) < 2) return(tibble())
  df_sub %>%
    group_by(across(all_of(group_vars))) %>%
    group_modify(~ {
      if (nrow(.x) < 2) return(tibble())
      combs <- utils::combn(seq_len(nrow(.x)), 2)
      idx1 <- combs[1, ]
      idx2 <- combs[2, ]
      tibble(
        Stream1 = .x$Stream_Name[idx1],
        Stream2 = .x$Stream_Name[idx2],
        storage_metric = metric,
        storage_diff = abs(.x[[metric]][idx1] - .x[[metric]][idx2]),
        cq_diff = abs(.x$cq_slope[idx1] - .x$cq_slope[idx2])
      )
    }) %>%
    ungroup()
}

plot_pairwise_scatter <- function(tbl, metric, filename, title_text, supplemental = FALSE) {
  if (!is.data.frame(tbl) || nrow(tbl) == 0) return(invisible(NULL))
  metric_label <- get_label(metric)
  tbl <- flag_outlet_pairs(tbl, c("Stream1", "Stream2"))
  non_outlet <- tbl[is.na(tbl$is_outlet_pair) | !tbl$is_outlet_pair, , drop = FALSE]
  outlet_tbl <- tbl[tbl$is_outlet_pair %in% TRUE, , drop = FALSE]
  p <- ggplot(non_outlet, aes(x = storage_diff, y = cq_diff)) +
    geom_point(alpha = 0.45, size = 2, color = "#4B5563")
  if (nrow(outlet_tbl) > 0) {
    p <- p + geom_point(
      data = outlet_tbl,
      aes(x = storage_diff, y = cq_diff),
      inherit.aes = FALSE,
      alpha = 0.9,
      size = 3.5,
      shape = 17,
      color = "#111827"
    )
  }
  p <- p +
    geom_smooth(data = tbl, aes(x = storage_diff, y = cq_diff),
                inherit.aes = FALSE, method = "lm", se = TRUE,
                color = "#1D4ED8", fill = "#93C5FD") +
    labs(
      x = paste0("|Δ ", metric_label, "| between sites"),
      y = "|Δ CQ slope| between sites",
    ) +
    theme_linkage
  save_plot(p, filename, supplemental = supplemental)
  invisible(p)
}
theme_linkage <- theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

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
  )

annual <- readr::read_csv(
  file.path(out_dir, "HJA_clean_annual.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    Stream_Name = as.factor(Stream_Name),
    solute = as.factor(solute)
  )

# =============================================================================
# 1. ANNUAL: Pairwise storage similarity vs CQ similarity
# =============================================================================
# Unit: site × solute × water_year
# For each solute–year, compare all site pairs

message("\n=== ANNUAL: Pairwise analysis ===\n")

primary_storage_metric <- intersect(PRIMARY_STORAGE_METRIC, names(annual))
supp_storage_metrics <- setdiff(
  intersect(SUPPLEMENTAL_STORAGE_METRICS, names(annual)),
  primary_storage_metric
)

if (length(primary_storage_metric) == 1) {
  primary_pairs <- make_pairwise_table(annual, primary_storage_metric)
  if (nrow(primary_pairs) > 0) {
    message("Found ", nrow(primary_pairs), " site pairs (annual, primary metric)")
    plot_pairwise_scatter(
      primary_pairs, primary_storage_metric,
      filename = "01_annual_storage_vs_cq_primary.png",
      title_text = "ANNUAL: Storage similarity vs chemistry similarity (primary metric)"
    )
  } else {
    message("No complete pairs for primary storage metric in annual data")
  }
} else {
  warning("Primary storage metric not found in annual dataset")
}

if (length(supp_storage_metrics) > 0) {
  for (metric in supp_storage_metrics) {
    supp_pairs <- make_pairwise_table(annual, metric)
    if (nrow(supp_pairs) == 0) next
    message("Supplemental metric ", metric, ": ", nrow(supp_pairs), " pairs")
    plot_pairwise_scatter(
      supp_pairs, metric,
      filename = paste0("01_annual_storage_vs_cq_", metric, ".png"),
      title_text = paste("ANNUAL: Storage similarity vs chemistry similarity (", metric, ")"),
      supplemental = TRUE
    )
  }
}

# =============================================================================
# 2. SEASONAL: Pairwise storage similarity vs CQ similarity
# =============================================================================
# Unit: site × solute × water_year × hydrologic_season

message("\n=== SEASONAL: Pairwise analysis ===\n")

if (all(c("Stream_Name", "solute", "water_year", "hydrologic_season", "cq_slope") %in% names(seasonal))) {
  group_vars <- c("solute", "water_year", "hydrologic_season")
  if (length(primary_storage_metric) == 1) {
    seasonal_pairs <- make_pairwise_table(seasonal, primary_storage_metric, group_vars = group_vars)
    if (nrow(seasonal_pairs) > 0) {
      message("Found ", nrow(seasonal_pairs), " site pairs (seasonal, primary metric)")
      p <- seasonal_pairs %>%
        mutate(hydrologic_season = factor(hydrologic_season, levels = c("Wet", "Dry"))) %>%
        flag_outlet_pairs(c("Stream1", "Stream2"))
      non_outlet <- p[is.na(p$is_outlet_pair) | !p$is_outlet_pair, , drop = FALSE]
      outlet_tbl <- p[p$is_outlet_pair %in% TRUE, , drop = FALSE]
      g <- ggplot(non_outlet, aes(x = storage_diff, y = cq_diff, color = hydrologic_season)) +
        geom_point(alpha = 0.45, size = 2) +
        scale_color_manual(values = c("Wet" = "steelblue", "Dry" = "darkorange"))
      if (nrow(outlet_tbl) > 0) {
        g <- g + geom_point(
          data = outlet_tbl,
          aes(x = storage_diff, y = cq_diff),
          inherit.aes = FALSE,
          alpha = 0.95,
          size = 3,
          shape = 21,
          stroke = 1,
          color = "black",
          fill = "gold"
        )
      }
      g <- g +
        geom_smooth(data = p, aes(x = storage_diff, y = cq_diff, color = hydrologic_season),
                    inherit.aes = FALSE, method = "lm", se = TRUE, alpha = 0.15) +
        labs(
          x = paste0("|Δ ", get_label(primary_storage_metric), "| between sites"),
          y = "|Δ CQ slope| between sites",
          color = "Season",
    }
  }
  if (length(supp_storage_metrics) > 0) {
    for (metric in supp_storage_metrics) {
      seasonal_pairs_supp <- make_pairwise_table(seasonal, metric, group_vars = group_vars)
      if (nrow(seasonal_pairs_supp) == 0) next
      message("Seasonal supplemental metric ", metric, ": ", nrow(seasonal_pairs_supp), " pairs")
      plot_pairwise_scatter(
        seasonal_pairs_supp,
        metric,
        filename = paste0("02_seasonal_storage_vs_cq_", metric, ".png"),
        title_text = paste("SEASONAL: Storage similarity vs chemistry similarity (", metric, ")"),
        supplemental = TRUE
      )
    }
  }
}

message("\n=== STORAGE-CHEMISTRY LINKAGES COMPLETE ===\n")
