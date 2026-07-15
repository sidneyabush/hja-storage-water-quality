#!/usr/bin/env Rscript
# =============================================================================
# 6a: Candidate main-paper figures from active storage-chemistry outputs
# =============================================================================
# Builds a compact set of manuscript-candidate figures from the active
# storage-chemistry analysis outputs. Existing ordination plots are copied into
# the candidate figure folder and summary figures are regenerated from tables.
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
fig_root <- paths$fig_root
res_dir <- file.path(out_dir, "06_figures", "main_paper_candidates")
fig_dir <- file.path(fig_root, "06_figures", "main_paper_candidates")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

require_file <- function(path) {
  if (!file.exists(path)) {
    stop("Missing required figure input: ", path)
  }
  path
}

theme_candidate <- function(base_size = 11) {
  if (exists("theme_hja")) {
    theme_hja(base_size = base_size)
  } else {
    theme_bw(base_size = base_size) +
      theme(panel.grid = element_blank())
  }
}

analysis_family_labels <- c(
  site_profile_ordination = "Site-level chemistry profile",
  annual_ordination = "Annual chemistry",
  pairwise_agreement = "Pairwise synchrony/agreement"
)

response_order <- c(
  "cq_slope",
  "cq_CVc_CVq",
  "conc_sync_allpairs",
  "wymore_crosssite_allpairs",
  "annual_stream_chemistry",
  "Abbott_S",
  "prop_sync_wymore",
  "cluster_agreement"
)

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

response_matrix_file <- require_file(file.path(
  out_dir,
  "05_synthesis",
  "storage_chemistry_links",
  "storage_chemistry_response_matrix.csv"
))
top_by_response_file <- require_file(file.path(
  out_dir,
  "05_synthesis",
  "storage_chemistry_links",
  "storage_chemistry_top_by_response.csv"
))
pairwise_long_file <- require_file(file.path(
  out_dir,
  "03_stats",
  "storage_metric_synchrony",
  "pairwise_synchrony_storage_metric_long.csv"
))

response_matrix <- readr::read_csv(response_matrix_file, show_col_types = FALSE) %>%
  mutate(
    response = factor(response, levels = response_order),
    response_label = factor(response_label, levels = rev(unname(response_labels[response_order]))),
    analysis_family = factor(analysis_family, levels = names(analysis_family_labels)),
    analysis_family_label = factor(analysis_family_label, levels = unname(analysis_family_labels)),
    storage_label = factor(
      storage_label,
      levels = vapply(PAPER_FACING_STORAGE_METRICS, get_storage_label, character(1), short = TRUE)
    )
  )

p_matrix <- ggplot(response_matrix, aes(x = storage_label, y = response_label, fill = effect_strength)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(
    aes(label = ifelse(statistic == "pearson_r", sprintf("%+.2f", signed_effect), sprintf("%.2f", effect_strength))),
    size = 2.8,
    color = "grey15"
  ) +
  facet_grid(analysis_family_label ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient(
    low = "#F7F7F2",
    high = "#2F6B9A",
    limits = c(0, 1),
    name = "Strength"
  ) +
  labs(x = "Storage-paper metric", y = NULL) +
  theme_candidate(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y = element_text(angle = 0, hjust = 0),
    legend.position = "bottom",
    panel.spacing.y = unit(0.55, "lines")
  )

ggsave(
  file.path(fig_dir, "Fig1_storage_chemistry_response_matrix.png"),
  p_matrix,
  width = 9.8,
  height = 7.2,
  dpi = PLOT_DPI,
  bg = "white"
)

top_by_response <- readr::read_csv(top_by_response_file, show_col_types = FALSE) %>%
  filter(within_response_rank <= 3) %>%
  mutate(
    response = factor(response, levels = response_order),
    response_label = factor(response_label, levels = response_labels[response_order]),
    storage_metric = factor(storage_metric, levels = PAPER_FACING_STORAGE_METRICS),
    estimate_label = ifelse(
      statistic == "pearson_r",
      sprintf("r = %+.2f", signed_effect),
      sprintf("R2 = %.2f", effect_strength)
    )
  )

p_top <- ggplot(
  top_by_response,
  aes(x = effect_strength, y = fct_rev(storage_metric), color = storage_domain)
) +
  geom_segment(aes(x = 0, xend = effect_strength, yend = fct_rev(storage_metric)), linewidth = 0.6, alpha = 0.75) +
  geom_point(size = 2.4) +
  geom_text(aes(label = estimate_label), hjust = -0.12, size = 2.7, color = "grey20") +
  facet_wrap(~ response_label, ncol = 2, scales = "free_y") +
  scale_x_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.25)) +
  scale_color_manual(
    values = c(
      "Dynamic and extended-dynamic storage" = "#2F6B9A",
      "Flow-path partitioning" = "#4F7F52",
      "Mobile mixing and tracer storage" = "#8A5A83"
    ),
    name = NULL
  ) +
  labs(x = "Link strength", y = NULL) +
  theme_candidate(base_size = 10.5) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(hjust = 0),
    panel.spacing = unit(0.9, "lines")
  )

ggsave(
  file.path(fig_dir, "Fig2_top_storage_links_by_response.png"),
  p_top,
  width = 10.5,
  height = 8.0,
  dpi = PLOT_DPI,
  bg = "white"
)

pairwise_top <- readr::read_csv(top_by_response_file, show_col_types = FALSE) %>%
  filter(analysis_family == "pairwise_agreement", within_response_rank == 1) %>%
  select(sync_response = response, response_label, storage_metric, top_storage_label = storage_label)

pairwise_long <- readr::read_csv(pairwise_long_file, show_col_types = FALSE) %>%
  inner_join(pairwise_top, by = c("sync_response", "storage_metric")) %>%
  filter(is.finite(storage_abs_diff), is.finite(sync_value)) %>%
  mutate(
    response_metric_label = paste0(response_label, "\n", top_storage_label, " difference"),
    response_metric_label = factor(response_metric_label, levels = unique(response_metric_label))
  )

p_pair <- ggplot(pairwise_long, aes(x = storage_abs_diff, y = sync_value)) +
  geom_point(alpha = 0.18, size = 1.1, color = "#374151") +
  geom_smooth(method = "lm", se = TRUE, color = "#2F6B9A", fill = "#A9C3D7", linewidth = 0.8) +
  facet_wrap(~ response_metric_label, scales = "free", nrow = 1) +
  labs(
    x = "Absolute annual storage-metric difference",
    y = "Annual chemistry similarity/agreement"
  ) +
  theme_candidate(base_size = 11) +
  theme(
    strip.text = element_text(hjust = 0),
    panel.spacing = unit(1.0, "lines")
  )

ggsave(
  file.path(fig_dir, "Fig3_pairwise_storage_similarity_scatter.png"),
  p_pair,
  width = 11.5,
  height = 4.0,
  dpi = PLOT_DPI,
  bg = "white"
)

existing_figures <- tibble::tribble(
  ~candidate_file, ~source_file, ~figure_role,
  "Fig4_annual_stream_chemistry_storage_pca.png",
  file.path(fig_root, "04_PCA", "4j_annual_chemistry_storage_ordination", "annual_stream_chemistry_storage_pca.png"),
  "Annual site-year chemistry ordination with storage vectors",
  "FigS1_cq_slope_storage_pca.png",
  file.path(fig_root, "04_PCA", "4i_storage_metric_ordination", "cq_slope_all_storage_metrics_biplot.png"),
  "Site-level C-Q slope chemistry profile ordination",
  "FigS2_cvcq_storage_pca.png",
  file.path(fig_root, "04_PCA", "4i_storage_metric_ordination", "cq_CVc_CVq_all_storage_metrics_biplot.png"),
  "Site-level CVc/CVq chemistry profile ordination",
  "FigS3_concentration_synchrony_storage_pca.png",
  file.path(fig_root, "04_PCA", "4i_storage_metric_ordination", "conc_sync_allpairs_all_storage_metrics_biplot.png"),
  "Site-level concentration synchrony profile ordination",
  "FigS4_wymore_cq_agreement_storage_pca.png",
  file.path(fig_root, "04_PCA", "4i_storage_metric_ordination", "wymore_crosssite_allpairs_all_storage_metrics_biplot.png"),
  "Site-level Wymore C-Q agreement profile ordination",
  "FigS5_pairwise_synchrony_storage_heatmap.png",
  file.path(fig_root, "03_stats", "storage_metric_synchrony", "annual_pairwise_synchrony_storage_metric_heatmap.png"),
  "Pairwise annual synchrony/agreement storage-distance heatmap"
)

figure_manifest <- existing_figures %>%
  mutate(
    copied = file.exists(source_file),
    target_file = file.path(fig_dir, candidate_file)
  )

for (i in seq_len(nrow(figure_manifest))) {
  if (isTRUE(figure_manifest$copied[i])) {
    file.copy(figure_manifest$source_file[i], figure_manifest$target_file[i], overwrite = TRUE)
  }
}

generated_manifest <- tibble::tribble(
  ~candidate_file, ~source_file, ~figure_role, ~copied, ~target_file,
  "Fig1_storage_chemistry_response_matrix.png", response_matrix_file,
  "Cross-analysis storage-chemistry response matrix", FALSE,
  file.path(fig_dir, "Fig1_storage_chemistry_response_matrix.png"),
  "Fig2_top_storage_links_by_response.png", top_by_response_file,
  "Top storage metrics by chemistry response", FALSE,
  file.path(fig_dir, "Fig2_top_storage_links_by_response.png"),
  "Fig3_pairwise_storage_similarity_scatter.png", pairwise_long_file,
  "Pairwise annual storage-difference scatterplots", FALSE,
  file.path(fig_dir, "Fig3_pairwise_storage_similarity_scatter.png")
)

figure_manifest <- bind_rows(generated_manifest, figure_manifest) %>%
  mutate(exists = file.exists(target_file))

readr::write_csv(figure_manifest, file.path(res_dir, "candidate_main_figure_manifest.csv"))

message("Candidate main-paper figures written to: ", fig_dir)
message("Candidate main-paper figure manifest written to: ", res_dir)
