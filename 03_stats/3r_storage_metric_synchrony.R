#!/usr/bin/env Rscript
# =============================================================================
# 3r: Storage-metric controls on pairwise synchrony
# =============================================================================
# Tests whether pairs of catchments with more similar annual storage-paper
# metrics behave more similarly in the same water year. Concentration synchrony
# remains the primary response; C-Q agreement and cluster agreement are retained
# as compact sensitivity checks. Annual analyses are restricted to the overlap
# between the water-quality record and the finalized storage-paper metrics.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

rm(list = ls())

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_flag <- "--file="
  matches <- grep(file_flag, cmd_args)
  if (length(matches) > 0) {
    script_path <- sub(file_flag, "", cmd_args[matches[1]])
    return(dirname(normalizePath(script_path)))
  }
  for (i in rev(seq_along(sys.calls()))) {
    call_i <- sys.calls()[[i]]
    if (identical(call_i[[1]], as.name("source"))) {
      file_arg <- tryCatch(as.character(eval(call_i[[2]], envir = sys.frame(i))), error = function(...) NA_character_)
      if (is.character(file_arg) && length(file_arg) > 0 && file.exists(file_arg[1])) {
        return(dirname(normalizePath(file_arg[1])))
      }
    }
  }
  normalizePath(getwd())
}

find_repo_root <- function(start_dir) {
  current <- normalizePath(start_dir)
  repeat {
    if (dir.exists(file.path(current, "00_helpers")) || dir.exists(file.path(current, ".git"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Unable to locate project root from: ", start_dir)
    }
    current <- parent
  }
}

repo_dir <- find_repo_root(get_script_dir())
source(file.path(repo_dir, "00_helpers", "workflow_config.R"))
source(file.path(repo_dir, "00_helpers", "plot_prefs.R"))
theme_file <- file.path(repo_dir, "00_helpers", "plot_theme_set.R")
if (file.exists(theme_file)) source(theme_file)

paths <- get_project_paths()
out_dir <- paths$out_dir
res_dir <- file.path(out_dir, "03_stats", "storage_metric_synchrony")
fig_dir <- file.path(paths$fig_root, "03_stats", "storage_metric_synchrony")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

theme_storage_sync <- function(base_size = 11) {
  if (exists("theme_hja")) {
    theme_hja(base_size = base_size)
  } else {
    theme_bw(base_size = base_size) +
      theme(panel.grid = element_blank())
  }
}

pair_file <- file.path(out_dir, "HJA_pair_sync_metrics.csv")
cluster_file <- file.path(out_dir, "ClusterStreams_allSolutes_byWaterYear.csv")
storage_file <- file.path(out_dir, "HJA_storage_framework_annual.csv")

if (!file.exists(pair_file)) {
  stop("Missing pairwise synchrony file: ", pair_file, "\nRun 01_data_prep/1_MASTER_DATA_PREP.R first.")
}
if (!file.exists(cluster_file)) {
  stop("Missing annual cluster file: ", cluster_file, "\nRun 01_data_prep/1_MASTER_DATA_PREP.R first.")
}
if (!file.exists(storage_file)) {
  stop("Missing storage framework file: ", storage_file, "\nRun 01_data_prep/1j_import_storage_framework.R first.")
}

storage_metrics <- PAPER_FACING_STORAGE_METRICS
sync_responses <- c("Abbott_S", "prop_sync_wymore", "cluster_agreement")
sync_response_labels <- c(
  Abbott_S = "Concentration synchrony",
  prop_sync_wymore = "C-Q agreement",
  cluster_agreement = "Cluster pattern agreement"
)

pair_raw <- readr::read_csv(pair_file, show_col_types = FALSE) %>%
  mutate(
    Stream1 = standardize_wq_stream(Stream1),
    Stream2 = standardize_wq_stream(Stream2),
    hydrologic_season = ifelse(is.na(hydrologic_season), "annual", hydrologic_season)
  ) %>%
  filter(
    time_scale == "annual",
    hydrologic_season == "annual",
    water_year >= STORAGE_CHEMISTRY_YEAR_START,
    water_year <= STORAGE_CHEMISTRY_YEAR_END
  )

response_cols <- intersect(c("Abbott_S", "prop_sync_wymore"), names(pair_raw))
if (length(response_cols) == 0) {
  stop("No recognized synchrony response columns found in: ", pair_file)
}

pair_summary <- pair_raw %>%
  group_by(
    time_scale,
    hydrologic_season,
    water_year,
    Stream1,
    Stream2,
    is_outlet_pair
  ) %>%
  summarise(
    across(all_of(response_cols), ~ {
      out <- mean(.x, na.rm = TRUE)
      ifelse(is.nan(out), NA_real_, out)
    }),
    n_solutes = n_distinct(solute),
    n_sync_records = n(),
    .groups = "drop"
  )

cluster_annual <- readr::read_csv(cluster_file, show_col_types = FALSE) %>%
  mutate(Stream_Name = standardize_wq_stream(Stream_Name)) %>%
  transmute(
    water_year = as.integer(water_year),
    Stream_Name,
    solute = chemical,
    Cluster_climRef = as.character(Cluster_climRef)
  ) %>%
  filter(
    !is.na(water_year),
    !is.na(Stream_Name),
    !is.na(solute),
    !is.na(Cluster_climRef),
    water_year >= STORAGE_CHEMISTRY_YEAR_START,
    water_year <= STORAGE_CHEMISTRY_YEAR_END
  )

cluster_pair_summary <- cluster_annual %>%
  inner_join(
    cluster_annual,
    by = c("water_year", "solute"),
    suffix = c("_1", "_2"),
    relationship = "many-to-many"
  ) %>%
  filter(Stream_Name_1 < Stream_Name_2) %>%
  transmute(
    time_scale = "annual",
    hydrologic_season = "annual",
    water_year,
    Stream1 = Stream_Name_1,
    Stream2 = Stream_Name_2,
    solute,
    cluster_match = as.integer(Cluster_climRef_1 == Cluster_climRef_2)
  ) %>%
  group_by(time_scale, hydrologic_season, water_year, Stream1, Stream2) %>%
  summarise(
    cluster_agreement = mean(cluster_match, na.rm = TRUE),
    n_cluster_solutes = n_distinct(solute),
    n_cluster_records = n(),
    .groups = "drop"
  ) %>%
  mutate(is_outlet_pair = Stream1 == OUTLET_SITE | Stream2 == OUTLET_SITE)

pair_summary <- pair_summary %>%
  full_join(
    cluster_pair_summary,
    by = c("time_scale", "hydrologic_season", "water_year", "Stream1", "Stream2", "is_outlet_pair")
  )

response_cols <- intersect(sync_responses, names(pair_summary))

storage_annual <- readr::read_csv(storage_file, show_col_types = FALSE) %>%
  mutate(
    Stream_Name = standardize_wq_stream(Stream_Name),
    water_year = as.integer(water_year)
  ) %>%
  filter(
    !is.na(water_year),
    water_year >= STORAGE_CHEMISTRY_YEAR_START,
    water_year <= STORAGE_CHEMISTRY_YEAR_END
  ) %>%
  select(Stream_Name, water_year, any_of(storage_metrics))

metrics_available <- intersect(storage_metrics, names(storage_annual))
if (length(metrics_available) == 0) {
  stop("No paper-facing storage metrics found in: ", storage_file)
}

storage_1 <- storage_annual %>%
  rename_with(~ paste0(.x, "_1"), all_of(metrics_available))
storage_2 <- storage_annual %>%
  rename_with(~ paste0(.x, "_2"), all_of(metrics_available))

pair_storage <- pair_summary %>%
  left_join(storage_1, by = c("Stream1" = "Stream_Name", "water_year" = "water_year")) %>%
  left_join(storage_2, by = c("Stream2" = "Stream_Name", "water_year" = "water_year"))

for (metric in metrics_available) {
  x1 <- pair_storage[[paste0(metric, "_1")]]
  x2 <- pair_storage[[paste0(metric, "_2")]]
  pair_storage[[paste0(metric, "_abs_diff")]] <- abs(x1 - x2)
  metric_mean <- rowMeans(cbind(x1, x2), na.rm = TRUE)
  metric_mean[is.nan(metric_mean)] <- NA_real_
  pair_storage[[paste0(metric, "_pair_mean")]] <- metric_mean
}

readr::write_csv(
  pair_storage,
  file.path(res_dir, "pairwise_synchrony_with_storage_metrics.csv")
)

diff_cols <- paste0(metrics_available, "_abs_diff")

pair_storage_long <- pair_storage %>%
  select(
    time_scale,
    hydrologic_season,
    water_year,
    Stream1,
    Stream2,
    is_outlet_pair,
    n_solutes,
    n_sync_records,
    all_of(response_cols),
    all_of(diff_cols)
  ) %>%
  pivot_longer(
    cols = all_of(response_cols),
    names_to = "sync_response",
    values_to = "sync_value"
  ) %>%
  pivot_longer(
    cols = all_of(diff_cols),
    names_to = "storage_metric",
    values_to = "storage_abs_diff"
  ) %>%
  mutate(
    storage_metric = sub("_abs_diff$", "", storage_metric),
    storage_label = vapply(storage_metric, get_storage_label, character(1), short = TRUE),
    sync_response_label = dplyr::recode(sync_response, !!!sync_response_labels)
  )

readr::write_csv(
  pair_storage_long,
  file.path(res_dir, "pairwise_synchrony_storage_metric_long.csv")
)

safe_cor <- function(x, y) {
  idx <- is.finite(x) & is.finite(y)
  n <- sum(idx)
  if (n < 5 || sd(x[idx]) == 0 || sd(y[idx]) == 0) {
    return(tibble(n = n, r = NA_real_, p = NA_real_))
  }
  ct <- suppressWarnings(cor.test(x[idx], y[idx]))
  tibble(n = n, r = unname(ct$estimate), p = ct$p.value)
}

cor_tbl <- pair_storage_long %>%
  group_by(time_scale, hydrologic_season, sync_response, sync_response_label, storage_metric, storage_label) %>%
  summarise(safe_cor(storage_abs_diff, sync_value), .groups = "drop") %>%
  mutate(
    abs_r = abs(r),
    p_rank = rank(p, ties.method = "first", na.last = "keep")
  ) %>%
  arrange(time_scale, hydrologic_season, sync_response, desc(abs_r))

readr::write_csv(
  cor_tbl,
  file.path(res_dir, "pairwise_synchrony_storage_metric_correlations.csv")
)

top_hits <- cor_tbl %>%
  filter(n >= 8, is.finite(r)) %>%
  arrange(desc(abs_r)) %>%
  group_by(time_scale, hydrologic_season, sync_response) %>%
  slice_head(n = 5) %>%
  ungroup()

readr::write_csv(
  top_hits,
  file.path(res_dir, "pairwise_synchrony_storage_metric_top_hits.csv")
)

plot_tbl <- cor_tbl %>%
  filter(
    time_scale == "annual",
    hydrologic_season == "annual",
    is.finite(r)
  ) %>%
  mutate(
    storage_label = factor(
      storage_label,
      levels = vapply(storage_metrics, get_storage_label, character(1), short = TRUE)
    ),
    sync_response_label = factor(
      sync_response_label,
      levels = rev(unname(sync_response_labels))
    )
  )

if (nrow(plot_tbl) > 0) {
  p_heat <- ggplot(plot_tbl, aes(x = storage_label, y = sync_response_label, fill = r)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(n >= 8, sprintf("%.2f", r), "")), size = 3) +
    scale_fill_gradient2(
      low = diverging_low_color,
      mid = diverging_mid_color,
      high = diverging_high_color,
      midpoint = 0,
      limits = c(-1, 1),
      name = "r"
    ) +
    labs(
      x = "Absolute difference in storage metric",
      y = NULL
    ) +
    theme_storage_sync(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(color = "grey15"),
      legend.position = "bottom"
    )

  ggsave(
    file.path(fig_dir, "annual_pairwise_synchrony_storage_metric_heatmap.png"),
    p_heat,
    width = 9,
    height = 4,
    dpi = PLOT_DPI,
    bg = "white"
  )
}

message("Pairwise synchrony + storage metric outputs written to: ", res_dir)
message("Pairwise synchrony + storage metric figures written to: ", fig_dir)
