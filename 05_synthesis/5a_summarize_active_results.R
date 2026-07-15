#!/usr/bin/env Rscript
# =============================================================================
# 5a: Summarize active storage-chemistry workflow results
# =============================================================================
# Reads the active ordination and synchrony/agreement outputs and writes compact
# tables plus a short text summary for paper planning.
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

paths <- get_project_paths()
out_dir <- paths$out_dir
res_dir <- file.path(out_dir, "05_synthesis", "active_workflow")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

require_file <- function(path) {
  if (!file.exists(path)) {
    stop("Missing required active workflow output: ", path)
  }
  path
}

response_labels <- c(
  cq_slope = "C-Q slope chemistry profile",
  cq_CVc_CVq = "CVc/CVq chemistry profile",
  conc_sync_allpairs = "Concentration synchrony profile",
  wymore_crosssite_allpairs = "Wymore C-Q agreement profile",
  annual_stream_chemistry = "Annual stream chemistry",
  Abbott_S = "Concentration synchrony",
  prop_sync_wymore = "C-Q agreement",
  cluster_agreement = "Cluster pattern agreement"
)

label_response <- function(x) {
  out <- unname(response_labels[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

pca_dir <- file.path(out_dir, "04_PCA", "storage_metric_ordination")
annual_pca_dir <- file.path(out_dir, "04_PCA", "annual_chemistry_storage_ordination")
sync_dir <- file.path(out_dir, "03_stats", "storage_metric_synchrony")

pca_summary_file <- require_file(file.path(pca_dir, "storage_metric_ordination_summary.csv"))
annual_vector_file <- require_file(file.path(annual_pca_dir, "annual_stream_chemistry_storage_vectors.csv"))
annual_scores_file <- require_file(file.path(annual_pca_dir, "annual_stream_chemistry_pca_scores.csv"))
annual_variance_file <- require_file(file.path(annual_pca_dir, "annual_stream_chemistry_pca_variance.csv"))
sync_cor_file <- require_file(file.path(sync_dir, "pairwise_synchrony_storage_metric_correlations.csv"))
sync_pair_file <- require_file(file.path(sync_dir, "pairwise_synchrony_with_storage_metrics.csv"))
import_summary_file <- require_file(file.path(out_dir, "storage_framework_import_summary.csv"))
join_summary_file <- require_file(file.path(out_dir, "storage_framework_join_summary.csv"))

pca_summary <- readr::read_csv(pca_summary_file, show_col_types = FALSE) %>%
  mutate(
    response_label = label_response(response),
    top_storage_label = vapply(top_storage_metric, get_storage_label, character(1), short = TRUE),
    top_framework_axis_label = vapply(top_framework_axis, get_storage_label, character(1), short = TRUE)
  ) %>%
  select(
    response,
    response_label,
    solute_group,
    n_sites,
    n_solutes,
    PC1_variance,
    PC2_variance,
    top_storage_metric,
    top_storage_label,
    top_storage_metric_r2,
    top_framework_axis,
    top_framework_axis_label,
    top_framework_axis_r2
  )

readr::write_csv(
  pca_summary,
  file.path(res_dir, "ordination_summary.csv")
)

vector_files <- list.files(
  pca_dir,
  pattern = "_storage_vectors\\.csv$",
  full.names = TRUE
)

vector_tbl <- purrr::map_dfr(vector_files, readr::read_csv, show_col_types = FALSE) %>%
  mutate(response_label = label_response(response))

top_storage_vectors <- vector_tbl %>%
  filter(vector_type == "storage_metric", is.finite(vector_r2)) %>%
  arrange(response, desc(vector_r2)) %>%
  group_by(response, response_label, solute_group) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  transmute(
    response,
    response_label,
    solute_group,
    storage_metric = variable,
    storage_label = label,
    n,
    PC1_r,
    PC1_p,
    PC2_r,
    PC2_p,
    vector_r2
  )

readr::write_csv(
  top_storage_vectors,
  file.path(res_dir, "ordination_top_storage_vectors.csv")
)

top_framework_axes <- vector_tbl %>%
  filter(vector_type == "framework_axis", is.finite(vector_r2)) %>%
  arrange(response, desc(vector_r2)) %>%
  group_by(response, response_label, solute_group) %>%
  slice_head(n = 3) %>%
  ungroup() %>%
  transmute(
    response,
    response_label,
    solute_group,
    framework_axis = variable,
    framework_axis_label = label,
    n,
    PC1_r,
    PC1_p,
    PC2_r,
    PC2_p,
    vector_r2
  )

readr::write_csv(
  top_framework_axes,
  file.path(res_dir, "ordination_top_framework_axes_sensitivity.csv")
)

annual_scores <- readr::read_csv(annual_scores_file, show_col_types = FALSE)
annual_variance <- readr::read_csv(annual_variance_file, show_col_types = FALSE)
annual_vectors <- readr::read_csv(annual_vector_file, show_col_types = FALSE)

annual_top_storage_vectors <- annual_vectors %>%
  filter(vector_type == "storage_metric", is.finite(vector_r2)) %>%
  arrange(desc(vector_r2)) %>%
  slice_head(n = 5) %>%
  transmute(
    response = "annual_stream_chemistry",
    response_label = label_response(response),
    storage_metric = variable,
    storage_label = label,
    n,
    PC1_r,
    PC1_p,
    PC2_r,
    PC2_p,
    vector_r2
  )

readr::write_csv(
  annual_top_storage_vectors,
  file.path(res_dir, "annual_chemistry_top_storage_vectors.csv")
)

sync_cor <- readr::read_csv(sync_cor_file, show_col_types = FALSE) %>%
  mutate(response_label = label_response(sync_response))

top_sync_cor <- sync_cor %>%
  filter(n >= 8, is.finite(r)) %>%
  arrange(sync_response, desc(abs_r)) %>%
  group_by(sync_response, response_label) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  select(
    sync_response,
    response_label,
    storage_metric,
    storage_label,
    n,
    r,
    p,
    abs_r
  )

readr::write_csv(
  top_sync_cor,
  file.path(res_dir, "synchrony_top_storage_correlations.csv")
)

pair_storage <- readr::read_csv(sync_pair_file, show_col_types = FALSE)
import_summary <- readr::read_csv(import_summary_file, show_col_types = FALSE)
join_summary <- readr::read_csv(join_summary_file, show_col_types = FALSE)

coverage_summary <- tibble(
  item = c(
    "storage_chemistry_overlap_start",
    "storage_chemistry_overlap_end",
    "pair_year_rows",
    "site_pairs",
    "water_years",
    "annual_chemistry_site_year_rows",
    "annual_chemistry_sites",
    "annual_chemistry_water_years",
    "rows_with_concentration_synchrony",
    "rows_with_cq_agreement",
    "rows_with_cluster_agreement"
  ),
  value = c(
    STORAGE_CHEMISTRY_YEAR_START,
    STORAGE_CHEMISTRY_YEAR_END,
    nrow(pair_storage),
    n_distinct(paste(pair_storage$Stream1, pair_storage$Stream2, sep = "--")),
    n_distinct(pair_storage$water_year),
    nrow(annual_scores),
    n_distinct(annual_scores$Stream_Name),
    n_distinct(annual_scores$water_year),
    sum(is.finite(pair_storage$Abbott_S)),
    sum(is.finite(pair_storage$prop_sync_wymore)),
    sum(is.finite(pair_storage$cluster_agreement))
  )
)

readr::write_csv(
  coverage_summary,
  file.path(res_dir, "active_workflow_coverage_summary.csv")
)

combined_top_results <- bind_rows(
  pca_summary %>%
    transmute(
      result_family = "ordination",
      response = response,
      response_label = response_label,
      storage_metric = top_storage_metric,
      storage_label = top_storage_label,
      n = n_sites,
      statistic = "vector_r2",
      estimate = top_storage_metric_r2,
      p = NA_real_,
      note = paste0("PC1+PC2 variance = ", round(PC1_variance + PC2_variance, 1), "%")
    ),
  annual_top_storage_vectors %>%
    slice_head(n = 1) %>%
    transmute(
      result_family = "annual_ordination",
      response = response,
      response_label = response_label,
      storage_metric = storage_metric,
      storage_label = storage_label,
      n = n,
      statistic = "vector_r2",
      estimate = vector_r2,
      p = NA_real_,
      note = paste0(
        "Annual site-year chemistry; PC1+PC2 variance = ",
        round(sum(annual_variance$variance_explained[annual_variance$PC %in% c("PC1", "PC2")]), 1),
        "%"
      )
    ),
  top_sync_cor %>%
    group_by(sync_response, response_label) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    transmute(
      result_family = "pairwise_agreement",
      response = sync_response,
      response_label = response_label,
      storage_metric = storage_metric,
      storage_label = storage_label,
      n = n,
      statistic = "pearson_r",
      estimate = r,
      p = p,
      note = "Correlation between annual storage-metric distance and chemistry agreement"
    )
)

readr::write_csv(
  combined_top_results,
  file.path(res_dir, "active_workflow_top_results.csv")
)

top_result_lines <- combined_top_results %>%
  mutate(
    estimate_txt = ifelse(is.na(estimate), "NA", sprintf("%.2f", estimate)),
    line = paste0(response_label, ": ", storage_metric, " (", statistic, " = ", estimate_txt, ")")
  ) %>%
  pull(line)

summary_lines <- c(
  "Active HJA storage-chemistry workflow summary",
  paste0("Generated: ", Sys.Date()),
  "",
  "Method notes",
  "- Builds from the storage paper's multi-metric framework rather than treating storage as one high-to-low gradient.",
  paste0("- Annual storage-chemistry analyses use the overlapping water years ", STORAGE_CHEMISTRY_YEAR_START, "-", STORAGE_CHEMISTRY_YEAR_END, "."),
  "- PCA storage arrows are post hoc fitted correlations with chemistry ordination axes, not constrained ordination axes or eigenvectors.",
  "- Cluster agreement is treated as seasonal-pattern agreement, not the primary synchrony definition.",
  "",
  "Workflow coverage",
  paste0("- ", coverage_summary$item, ": ", coverage_summary$value),
  "",
  "Storage-framework import",
  paste0("- ", import_summary$item, ": ", import_summary$value),
  "",
  "Storage-framework joins",
  paste0(
    "- ", join_summary$output,
    ": ", join_summary$rows, " rows; ",
    join_summary$n_sites, " sites"
  ),
  "",
  "Top active results",
  paste0("- ", top_result_lines)
)

writeLines(summary_lines, file.path(res_dir, "active_workflow_results_summary.txt"))

message("Active workflow synthesis written to: ", res_dir)
