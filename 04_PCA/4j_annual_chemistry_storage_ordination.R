#!/usr/bin/env Rscript
# =============================================================================
# 4j: Annual stream-chemistry ordination with annual storage-paper metrics
# =============================================================================
# Rows are site-water-year observations; columns are annual mean log chemistry
# by solute. Storage-paper metrics are fitted as post hoc annual vectors.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
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
fig_dir <- file.path(paths$fig_root, "04_PCA", "4j_annual_chemistry_storage_ordination")
res_dir <- file.path(out_dir, "04_PCA", "annual_chemistry_storage_ordination")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

cq_file <- file.path(out_dir, "HJA_CQ_master.csv")
storage_file <- file.path(out_dir, "HJA_storage_framework_annual.csv")
if (!file.exists(cq_file)) stop("Missing chemistry master file: ", cq_file)
if (!file.exists(storage_file)) stop("Missing annual storage file: ", storage_file)

min_annual_solute_samples <- 3L

theme_annual_ord <- function(base_size = 12) {
  if (exists("theme_hja")) {
    theme_hja(base_size = base_size)
  } else {
    theme_bw(base_size = base_size) +
      theme(panel.grid = element_blank())
  }
}

impute_column_means <- function(mat) {
  for (j in seq_len(ncol(mat))) {
    idx <- !is.finite(mat[, j])
    if (any(idx)) mat[idx, j] <- mean(mat[, j], na.rm = TRUE)
  }
  mat
}

safe_cor_vec <- function(x, y) {
  idx <- is.finite(x) & is.finite(y)
  if (sum(idx) < 5 || sd(x[idx]) == 0 || sd(y[idx]) == 0) {
    return(c(r = NA_real_, p = NA_real_, n = sum(idx)))
  }
  ct <- suppressWarnings(cor.test(x[idx], y[idx]))
  c(r = unname(ct$estimate), p = ct$p.value, n = sum(idx))
}

annual_chem_long <- readr::read_csv(cq_file, show_col_types = FALSE) %>%
  mutate(
    Stream_Name = standardize_wq_stream(Stream_Name),
    water_year = if_else(Month >= 10L, Year + 1L, Year),
    variable = as.character(variable),
    value = as.numeric(value)
  ) %>%
  filter(
    variable %in% solute_order,
    is.finite(value),
    value > 0,
    water_year >= STORAGE_CHEMISTRY_YEAR_START,
    water_year <= STORAGE_CHEMISTRY_YEAR_END
  ) %>%
  group_by(Stream_Name, water_year, variable) %>%
  summarise(
    mean_conc = mean(value, na.rm = TRUE),
    mean_log10_conc = mean(log10(value), na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  )

readr::write_csv(
  annual_chem_long,
  file.path(res_dir, "annual_stream_chemistry_site_year_long.csv")
)

annual_wide <- annual_chem_long %>%
  filter(n_samples >= min_annual_solute_samples) %>%
  select(Stream_Name, water_year, variable, mean_log10_conc) %>%
  pivot_wider(names_from = variable, values_from = mean_log10_conc)

storage_annual <- readr::read_csv(storage_file, show_col_types = FALSE) %>%
  mutate(
    Stream_Name = standardize_wq_stream(Stream_Name),
    water_year = as.integer(water_year)
  ) %>%
  filter(
    water_year >= STORAGE_CHEMISTRY_YEAR_START,
    water_year <= STORAGE_CHEMISTRY_YEAR_END
  ) %>%
  select(Stream_Name, water_year, any_of(c(PAPER_FACING_STORAGE_METRICS, STORAGE_FRAMEWORK_AXIS_METRICS)))

annual_joined <- annual_wide %>%
  inner_join(storage_annual, by = c("Stream_Name", "water_year")) %>%
  mutate(row_id = paste(Stream_Name, water_year, sep = "_")) %>%
  arrange(Stream_Name, water_year)

readr::write_csv(
  annual_joined,
  file.path(res_dir, "annual_stream_chemistry_storage_joined.csv")
)

chem_cols <- intersect(solute_order, names(annual_joined))
if (length(chem_cols) < 2) {
  stop("Annual chemistry matrix has fewer than two solute columns after filtering.")
}

mat <- annual_joined %>%
  select(all_of(chem_cols)) %>%
  as.matrix()
rownames(mat) <- annual_joined$row_id

row_keep <- rowMeans(!is.finite(mat)) <= 0.5
mat <- mat[row_keep, , drop = FALSE]
annual_joined <- annual_joined[row_keep, , drop = FALSE]

col_var <- apply(mat, 2, function(x) var(x, na.rm = TRUE))
col_keep <- is.finite(col_var) & col_var > 0
mat <- mat[, col_keep, drop = FALSE]
chem_cols <- colnames(mat)

if (nrow(mat) < 5 || ncol(mat) < 2) {
  stop("Annual chemistry PCA matrix is too small after filtering.")
}

mat <- impute_column_means(mat)
pca <- prcomp(mat, center = TRUE, scale. = TRUE)
var_explained <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

scores <- as_tibble(pca$x[, 1:min(4, ncol(pca$x)), drop = FALSE], rownames = "row_id") %>%
  left_join(
    annual_joined %>%
      select(row_id, Stream_Name, water_year, any_of(c(PAPER_FACING_STORAGE_METRICS, STORAGE_FRAMEWORK_AXIS_METRICS))),
    by = "row_id"
  ) %>%
  mutate(Stream_Name = factor(Stream_Name, levels = site_order))

site_score_summary <- scores %>%
  group_by(Stream_Name) %>%
  summarise(
    site_label = standardize_storage_site(as.character(dplyr::first(Stream_Name))),
    PC1_mean = mean(PC1, na.rm = TRUE),
    PC1_sd = sd(PC1, na.rm = TRUE),
    PC2_mean = mean(PC2, na.rm = TRUE),
    PC2_sd = sd(PC2, na.rm = TRUE),
    n_years = n_distinct(water_year),
    .groups = "drop"
  ) %>%
  mutate(
    PC1_sd = replace_na(PC1_sd, 0),
    PC2_sd = replace_na(PC2_sd, 0),
    PC1_min = PC1_mean - PC1_sd,
    PC1_max = PC1_mean + PC1_sd,
    PC2_min = PC2_mean - PC2_sd,
    PC2_max = PC2_mean + PC2_sd,
    Stream_Name = factor(Stream_Name, levels = site_order)
  )

loadings <- as_tibble(pca$rotation[, 1:min(4, ncol(pca$rotation)), drop = FALSE], rownames = "solute")

vector_vars <- c(PAPER_FACING_STORAGE_METRICS, STORAGE_FRAMEWORK_AXIS_METRICS)
vector_vars <- vector_vars[vector_vars %in% names(scores)]
vector_tbl <- purrr::map_dfr(vector_vars, function(v) {
  pc1 <- safe_cor_vec(scores[[v]], scores$PC1)
  pc2 <- safe_cor_vec(scores[[v]], scores$PC2)
  tibble(
    variable = v,
    label = get_storage_label(v, short = TRUE),
    vector_type = ifelse(v %in% PAPER_FACING_STORAGE_METRICS, "storage_metric", "framework_axis"),
    PC1_r = pc1[["r"]],
    PC1_p = pc1[["p"]],
    PC2_r = pc2[["r"]],
    PC2_p = pc2[["p"]],
    n = min(pc1[["n"]], pc2[["n"]]),
    vector_r2 = PC1_r^2 + PC2_r^2
  )
}) %>%
  arrange(desc(vector_r2))

readr::write_csv(scores, file.path(res_dir, "annual_stream_chemistry_pca_scores.csv"))
readr::write_csv(site_score_summary, file.path(res_dir, "annual_stream_chemistry_pca_site_summary.csv"))
readr::write_csv(loadings, file.path(res_dir, "annual_stream_chemistry_pca_loadings.csv"))
readr::write_csv(vector_tbl, file.path(res_dir, "annual_stream_chemistry_storage_vectors.csv"))
readr::write_csv(
  tibble(
    PC = paste0("PC", seq_along(var_explained)),
    variance_explained = var_explained,
    cumulative_variance = cumsum(var_explained)
  ),
  file.path(res_dir, "annual_stream_chemistry_pca_variance.csv")
)

plot_extent <- max(
  abs(site_score_summary$PC1_min),
  abs(site_score_summary$PC1_max),
  abs(site_score_summary$PC2_min),
  abs(site_score_summary$PC2_max),
  na.rm = TRUE
)

loading_scale <- plot_extent /
  max(abs(loadings$PC1), abs(loadings$PC2), na.rm = TRUE) * 0.55
loadings_plot <- loadings %>%
  mutate(PC1_plot = PC1 * loading_scale, PC2_plot = PC2 * loading_scale)

vector_plot <- vector_tbl %>%
  filter(vector_type == "storage_metric", is.finite(PC1_r), is.finite(PC2_r))

if (nrow(vector_plot) > 0) {
  vector_scale <- plot_extent /
    max(sqrt(vector_plot$PC1_r^2 + vector_plot$PC2_r^2), na.rm = TRUE) * 0.85
  vector_plot <- vector_plot %>%
    mutate(PC1_plot = PC1_r * vector_scale, PC2_plot = PC2_r * vector_scale)
}

p <- ggplot() +
  geom_hline(yintercept = 0, color = "grey88", linewidth = 0.35) +
  geom_vline(xintercept = 0, color = "grey88", linewidth = 0.35) +
  geom_segment(
    data = site_score_summary,
    aes(x = PC1_min, xend = PC1_max, y = PC2_mean, yend = PC2_mean, color = Stream_Name),
    linewidth = 0.7,
    alpha = 0.65
  ) +
  geom_segment(
    data = site_score_summary,
    aes(x = PC1_mean, xend = PC1_mean, y = PC2_min, yend = PC2_max, color = Stream_Name),
    linewidth = 0.7,
    alpha = 0.65
  ) +
  geom_point(
    data = site_score_summary,
    aes(x = PC1_mean, y = PC2_mean, fill = Stream_Name, color = Stream_Name),
    shape = 21,
    size = 4.1,
    stroke = 0.5,
    alpha = 0.95
  ) +
  geom_text_repel(
    data = site_score_summary,
    aes(x = PC1_mean, y = PC2_mean, label = site_label, color = Stream_Name),
    size = 3.1,
    max.overlaps = 20,
    segment.alpha = 0.35
  ) +
  geom_segment(
    data = loadings_plot,
    aes(x = 0, y = 0, xend = PC1_plot, yend = PC2_plot),
    arrow = arrow(length = unit(0.016, "npc")),
    color = ordination_solute_vector_color,
    linewidth = 0.45,
    alpha = 0.55
  ) +
  geom_text_repel(
    data = loadings_plot,
    aes(x = PC1_plot, y = PC2_plot, label = solute),
    size = 3,
    color = ordination_solute_vector_color,
    max.overlaps = 30,
    segment.alpha = 0.15
  ) +
  geom_segment(
    data = vector_plot,
    aes(x = 0, y = 0, xend = PC1_plot, yend = PC2_plot),
    arrow = arrow(length = unit(0.024, "npc")),
    color = ordination_storage_vector_color,
    linewidth = 0.9
  ) +
  geom_text_repel(
    data = vector_plot,
    aes(x = PC1_plot, y = PC2_plot, label = label),
    size = 3,
    color = ordination_storage_label_color,
    fontface = "bold",
    max.overlaps = 30,
    segment.alpha = 0.25
  ) +
  labs(
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
    color = NULL
  ) +
  scale_color_site(name = NULL, drop = FALSE) +
  scale_fill_site(name = NULL, drop = FALSE) +
  theme_annual_ord(base_size = 12) +
  theme(
    legend.position = "none",
    plot.margin = margin(10, 14, 10, 10)
  )

ggsave(
  file.path(fig_dir, "annual_stream_chemistry_storage_pca.png"),
  p,
  width = 10,
  height = 8,
  dpi = PLOT_DPI,
  bg = "white"
)

message("Annual chemistry-storage ordination outputs written to: ", res_dir)
message("Annual chemistry-storage ordination figure written to: ", fig_dir)
