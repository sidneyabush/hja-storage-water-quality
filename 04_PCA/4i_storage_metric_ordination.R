# =============================================================================
# 4i: Storage-paper metric ordination of HJA stream chemistry
# =============================================================================
# Builds site x solute chemistry matrices and overlays finalized storage-paper
# metrics as fitted vectors. This keeps the chemistry PCA unconstrained
# while asking which storage dimensions align with the dominant chemistry axes.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
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
fig_dir <- file.path(paths$fig_root, "04_PCA", "4i_storage_metric_ordination")
res_dir <- file.path(out_dir, "04_PCA", "storage_metric_ordination")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

site_file <- file.path(out_dir, "HJA_master_site_means_storage_framework.csv")
if (!file.exists(site_file)) {
  stop(
    "Missing storage-framework joined site file: ", site_file,
    "\nRun 01_data_prep/1j_import_storage_framework.R first."
  )
}

site_means <- readr::read_csv(site_file, show_col_types = FALSE) %>%
  mutate(
    site = if ("site" %in% names(.)) standardize_storage_site(site) else standardize_storage_site(Stream_Name),
    Stream_Name = if ("Stream_Name" %in% names(.)) Stream_Name else standardize_wq_stream(site)
  ) %>%
  add_solute_type("solute", three_way = TRUE)

if (!"cq_slope" %in% names(site_means) && "cq_slope_windowed" %in% names(site_means)) {
  site_means <- site_means %>% mutate(cq_slope = cq_slope_windowed)
}
if (!"cq_CVc_CVq" %in% names(site_means) && "cq_CVc_CVq_windowed" %in% names(site_means)) {
  site_means <- site_means %>% mutate(cq_CVc_CVq = cq_CVc_CVq_windowed)
}

empty_annual_profiles <- tibble(
  site = character(),
  Stream_Name = character(),
  solute = character(),
  water_year = integer()
)

annual_cq_file <- file.path(out_dir, "HJA_master_annual_storage_framework.csv")
annual_sync_file <- file.path(out_dir, "HJA_composite_synchrony_annual.csv")
annual_wymore_file <- file.path(out_dir, "HJA_wymore_crosssite_sync.csv")

annual_cq <- if (file.exists(annual_cq_file)) {
  readr::read_csv(annual_cq_file, show_col_types = FALSE) %>%
    transmute(
      site = standardize_storage_site(Stream_Name),
      Stream_Name = standardize_wq_stream(Stream_Name),
      solute = as.character(solute),
      water_year = as.integer(water_year),
      cq_slope = cq_slope_windowed,
      cq_CVc_CVq = cq_CVc_CVq_windowed
    ) %>%
    filter(
      water_year >= STORAGE_CHEMISTRY_YEAR_START,
      water_year <= STORAGE_CHEMISTRY_YEAR_END
    )
} else {
  empty_annual_profiles
}

annual_sync <- if (file.exists(annual_sync_file)) {
  readr::read_csv(annual_sync_file, show_col_types = FALSE) %>%
    transmute(
      site = standardize_storage_site(Stream_Name),
      Stream_Name = standardize_wq_stream(Stream_Name),
      solute = as.character(solute),
      water_year = as.integer(water_year),
      conc_sync_allpairs = conc_sync_allpairs
    ) %>%
    filter(
      water_year >= STORAGE_CHEMISTRY_YEAR_START,
      water_year <= STORAGE_CHEMISTRY_YEAR_END
    )
} else {
  empty_annual_profiles
}

annual_wymore <- if (file.exists(annual_wymore_file)) {
  readr::read_csv(annual_wymore_file, show_col_types = FALSE) %>%
    mutate(
      Stream1 = standardize_wq_stream(Stream1),
      Stream2 = standardize_wq_stream(Stream2),
      water_year = as.integer(water_year),
      solute = as.character(solute)
    ) %>%
    filter(
      !is.na(sync),
      water_year >= STORAGE_CHEMISTRY_YEAR_START,
      water_year <= STORAGE_CHEMISTRY_YEAR_END
    ) %>%
    select(solute, water_year, Stream1, Stream2, sync) %>%
    pivot_longer(
      cols = c(Stream1, Stream2),
      names_to = "pair_position",
      values_to = "Stream_Name"
    ) %>%
    group_by(Stream_Name, solute, water_year) %>%
    summarise(
      wymore_crosssite_allpairs = mean(sync, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      site = standardize_storage_site(Stream_Name),
      Stream_Name = standardize_wq_stream(Stream_Name)
    ) %>%
    select(site, Stream_Name, solute, water_year, wymore_crosssite_allpairs)
} else {
  empty_annual_profiles
}

annual_profile_list <- list(annual_cq, annual_sync, annual_wymore)
annual_profile_list <- annual_profile_list[vapply(annual_profile_list, nrow, integer(1)) > 0]
annual_profiles <- if (length(annual_profile_list) > 0) {
  purrr::reduce(
    annual_profile_list,
    full_join,
    by = c("site", "Stream_Name", "solute", "water_year")
  ) %>%
    add_solute_type("solute", three_way = TRUE)
} else {
  empty_annual_profiles %>%
    add_solute_type("solute", three_way = TRUE)
}

theme_storage_ord <- function(base_size = 12) {
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
  if (sum(idx) < 3 || sd(x[idx]) == 0 || sd(y[idx]) == 0) {
    return(c(r = NA_real_, p = NA_real_, n = sum(idx)))
  }
  ct <- suppressWarnings(cor.test(x[idx], y[idx]))
  c(r = unname(ct$estimate), p = ct$p.value, n = sum(idx))
}

project_annual_profile_scores <- function(annual_df, response_col, solute_group, pca, pca_solutes) {
  if (!response_col %in% names(annual_df) || nrow(annual_df) == 0) {
    return(tibble())
  }

  df_group <- annual_df %>%
    filter(
      is.finite(.data[[response_col]]),
      !is.na(solute),
      !is.na(site),
      solute %in% pca_solutes
    )

  if (!identical(solute_group, "All")) {
    df_group <- df_group %>% filter(as.character(solute_type) == solute_group)
  }

  if (n_distinct(df_group$site) < 3 || nrow(df_group) == 0) {
    return(tibble())
  }

  wide <- df_group %>%
    group_by(site, Stream_Name, water_year, solute) %>%
    summarise(value = mean(.data[[response_col]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = solute, values_from = value)

  missing_solutes <- setdiff(pca_solutes, names(wide))
  for (solute_i in missing_solutes) {
    wide[[solute_i]] <- NA_real_
  }

  annual_mat <- wide %>%
    select(all_of(pca_solutes)) %>%
    as.matrix()

  row_keep <- rowMeans(!is.finite(annual_mat)) <= 0.5
  annual_mat <- annual_mat[row_keep, , drop = FALSE]
  wide <- wide[row_keep, , drop = FALSE]

  if (nrow(annual_mat) == 0) {
    return(tibble())
  }

  for (j in seq_along(pca_solutes)) {
    idx <- !is.finite(annual_mat[, j])
    if (any(idx)) annual_mat[idx, j] <- pca$center[[pca_solutes[j]]]
  }

  annual_scaled <- sweep(annual_mat, 2, pca$center[pca_solutes], "-")
  annual_scaled <- sweep(annual_scaled, 2, pca$scale[pca_solutes], "/")
  projected <- annual_scaled %*% pca$rotation[pca_solutes, 1:min(4, ncol(pca$rotation)), drop = FALSE]

  as_tibble(projected) %>%
    mutate(
      site = wide$site,
      Stream_Name = wide$Stream_Name,
      water_year = wide$water_year,
      response = response_col,
      solute_group = solute_group
    ) %>%
    relocate(site, Stream_Name, water_year, response, solute_group)
}

run_response_pca <- function(df, response_col, response_label, solute_group = "All") {
  if (!response_col %in% names(df)) {
    message("Skipping missing response: ", response_col)
    return(NULL)
  }

  df_group <- df %>%
    filter(is.finite(.data[[response_col]]), !is.na(solute), !is.na(site))

  if (!identical(solute_group, "All")) {
    df_group <- df_group %>% filter(as.character(solute_type) == solute_group)
  }

  if (n_distinct(df_group$site) < 3 || n_distinct(df_group$solute) < 2) {
    message("Skipping ", response_col, " / ", solute_group, ": insufficient sites or solutes.")
    return(NULL)
  }

  wide <- df_group %>%
    group_by(site, solute) %>%
    summarise(value = mean(.data[[response_col]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = solute, values_from = value)

  site_names <- wide$site
  mat <- wide %>% select(-site) %>% as.matrix()
  rownames(mat) <- site_names

  row_keep <- rowMeans(!is.finite(mat)) <= 0.5
  mat <- mat[row_keep, , drop = FALSE]
  site_names <- site_names[row_keep]

  col_var <- apply(mat, 2, function(x) var(x, na.rm = TRUE))
  col_keep <- is.finite(col_var) & col_var > 0
  mat <- mat[, col_keep, drop = FALSE]

  if (nrow(mat) < 3 || ncol(mat) < 2) {
    message("Skipping ", response_col, " / ", solute_group, ": matrix too small after filtering.")
    return(NULL)
  }

  mat <- impute_column_means(mat)
  pca <- prcomp(mat, center = TRUE, scale. = TRUE)
  var_explained <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

  site_storage <- df %>%
    select(
      site,
      Stream_Name,
      any_of(c(
        STORAGE_FRAMEWORK_AXIS_METRICS,
        PAPER_FACING_STORAGE_METRICS,
        "geology_class",
        "geomorphology_class"
      ))
    ) %>%
    distinct(site, .keep_all = TRUE)

  scores <- as_tibble(pca$x[, 1:min(4, ncol(pca$x)), drop = FALSE], rownames = "site") %>%
    left_join(site_storage, by = "site") %>%
    mutate(
      Stream_Name = ifelse(is.na(Stream_Name), standardize_wq_stream(site), Stream_Name),
      response = response_col,
      solute_group = solute_group
    )

  annual_projected_scores <- project_annual_profile_scores(
    annual_profiles,
    response_col,
    solute_group,
    pca,
    colnames(mat)
  )

  score_errorbars <- if (nrow(annual_projected_scores) > 0) {
    annual_projected_scores %>%
      group_by(site) %>%
      summarise(
        PC1_sd = sd(PC1, na.rm = TRUE),
        PC2_sd = sd(PC2, na.rm = TRUE),
        n_years = n_distinct(water_year),
        .groups = "drop"
      )
  } else {
    tibble(site = scores$site, PC1_sd = 0, PC2_sd = 0, n_years = NA_integer_)
  }

  scores <- scores %>%
    left_join(score_errorbars, by = "site") %>%
    mutate(
      PC1_sd = replace_na(PC1_sd, 0),
      PC2_sd = replace_na(PC2_sd, 0),
      PC1_min = PC1 - PC1_sd,
      PC1_max = PC1 + PC1_sd,
      PC2_min = PC2 - PC2_sd,
      PC2_max = PC2 + PC2_sd
    )

  loadings <- as_tibble(pca$rotation[, 1:min(4, ncol(pca$rotation)), drop = FALSE], rownames = "solute") %>%
    mutate(response = response_col, solute_group = solute_group)

  # Paper-facing ordination arrows use the named storage-paper metrics.
  # Figure 7 framework axes are still written to the vector tables for
  # sensitivity checks, but they are not the default biplot arrows.
  vector_vars <- c(PAPER_FACING_STORAGE_METRICS, STORAGE_FRAMEWORK_AXIS_METRICS)
  vector_vars <- vector_vars[vector_vars %in% names(scores)]
  vector_tbl <- purrr::map_dfr(vector_vars, function(v) {
    pc1 <- safe_cor_vec(scores[[v]], scores$PC1)
    pc2 <- safe_cor_vec(scores[[v]], scores$PC2)
    tibble(
      response = response_col,
      response_label = response_label,
      solute_group = solute_group,
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

  metric_vector_tbl <- vector_tbl %>%
    filter(vector_type == "storage_metric")

  framework_axis_tbl <- vector_tbl %>%
    filter(vector_type == "framework_axis")

  prefix <- paste(response_col, tolower(gsub("[^A-Za-z0-9]+", "_", solute_group)), sep = "_")
  readr::write_csv(scores, file.path(res_dir, paste0(prefix, "_scores.csv")))
  readr::write_csv(annual_projected_scores, file.path(res_dir, paste0(prefix, "_annual_projected_scores.csv")))
  readr::write_csv(loadings, file.path(res_dir, paste0(prefix, "_loadings.csv")))
  readr::write_csv(vector_tbl, file.path(res_dir, paste0(prefix, "_storage_vectors.csv")))
  readr::write_csv(
    tibble(
      response = response_col,
      solute_group = solute_group,
      PC = paste0("PC", seq_along(var_explained)),
      variance_explained = var_explained,
      cumulative_variance = cumsum(var_explained)
    ),
    file.path(res_dir, paste0(prefix, "_variance.csv"))
  )

  plot_extent <- max(
    abs(scores$PC1_min),
    abs(scores$PC1_max),
    abs(scores$PC2_min),
    abs(scores$PC2_max),
    na.rm = TRUE
  )

  loading_scale <- plot_extent /
    max(abs(loadings$PC1), abs(loadings$PC2), na.rm = TRUE) * 0.65
  loadings_plot <- loadings %>%
    mutate(PC1_plot = PC1 * loading_scale, PC2_plot = PC2 * loading_scale)

  vector_plot <- metric_vector_tbl %>%
    filter(is.finite(PC1_r), is.finite(PC2_r))

  if (nrow(vector_plot) > 0) {
    vector_scale <- plot_extent /
      max(sqrt(vector_plot$PC1_r^2 + vector_plot$PC2_r^2), na.rm = TRUE) * 0.9
    vector_plot <- vector_plot %>%
      mutate(PC1_plot = PC1_r * vector_scale, PC2_plot = PC2_r * vector_scale)
  }

  scores <- scores %>%
    mutate(Stream_Name = factor(Stream_Name, levels = site_order))

  p <- ggplot() +
    geom_hline(yintercept = 0, color = "grey88", linewidth = 0.35) +
    geom_vline(xintercept = 0, color = "grey88", linewidth = 0.35) +
    geom_segment(
      data = scores,
      aes(x = PC1_min, xend = PC1_max, y = PC2, yend = PC2, color = Stream_Name),
      linewidth = 0.7,
      alpha = 0.65
    ) +
    geom_segment(
      data = scores,
      aes(x = PC1, xend = PC1, y = PC2_min, yend = PC2_max, color = Stream_Name),
      linewidth = 0.7,
      alpha = 0.65
    ) +
    geom_point(
      data = scores,
      aes(x = PC1, y = PC2, fill = Stream_Name, color = Stream_Name),
      shape = 21,
      size = 4.1,
      stroke = 0.5,
      alpha = 0.95
    ) +
    geom_text_repel(
      data = scores,
      aes(x = PC1, y = PC2, label = site, color = Stream_Name),
      size = 3,
      max.overlaps = 20,
      segment.alpha = 0.35
    ) +
    geom_segment(
      data = loadings_plot,
      aes(x = 0, y = 0, xend = PC1_plot, yend = PC2_plot),
      arrow = arrow(length = unit(0.018, "npc")),
      color = ordination_solute_vector_color,
      linewidth = 0.55,
      alpha = 0.7
    ) +
    geom_text_repel(
      data = loadings_plot,
      aes(x = PC1_plot, y = PC2_plot, label = solute),
      size = 3.2,
      color = ordination_solute_vector_color,
      max.overlaps = 30,
      segment.alpha = 0.2
    ) +
    geom_segment(
      data = vector_plot,
      aes(x = 0, y = 0, xend = PC1_plot, yend = PC2_plot),
      arrow = arrow(length = unit(0.026, "npc")),
      color = ordination_storage_vector_color,
      linewidth = 1
    ) +
    geom_text_repel(
      data = vector_plot,
      aes(x = PC1_plot, y = PC2_plot, label = label),
      size = 3.2,
      color = ordination_storage_label_color,
      fontface = "bold",
      max.overlaps = 30,
      segment.alpha = 0.25
    ) +
    labs(
      x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
      y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
      fill = NULL
    ) +
    scale_color_site(name = NULL, drop = FALSE) +
    scale_fill_site(name = NULL, drop = FALSE) +
    theme_storage_ord(base_size = 12) +
    theme(
      legend.position = "none",
      plot.margin = margin(10, 14, 10, 10)
    )

  ggsave(file.path(fig_dir, paste0(prefix, "_storage_metrics_biplot.png")), p, width = 10, height = 8, dpi = PLOT_DPI, bg = "white")

  list(
    response = response_col,
    solute_group = solute_group,
    n_sites = nrow(scores),
    n_solutes = nrow(loadings),
    PC1_variance = var_explained[1],
    PC2_variance = var_explained[2],
    top_storage_metric = if (nrow(metric_vector_tbl) > 0) metric_vector_tbl$variable[1] else NA_character_,
    top_storage_metric_r2 = if (nrow(metric_vector_tbl) > 0) metric_vector_tbl$vector_r2[1] else NA_real_,
    top_framework_axis = if (nrow(framework_axis_tbl) > 0) framework_axis_tbl$variable[1] else NA_character_,
    top_framework_axis_r2 = if (nrow(framework_axis_tbl) > 0) framework_axis_tbl$vector_r2[1] else NA_real_
  )
}

response_specs <- tribble(
  ~response_col, ~response_label,
  "cq_slope", "C-Q slope chemistry profile",
  "cq_CVc_CVq", "CVc/CVq chemistry profile",
  "conc_sync_allpairs", "Concentration synchrony profile",
  "wymore_crosssite_allpairs", "Wymore C-Q agreement profile"
)

# Primary ordinations use all solutes together. Splitting by solute group makes
# the PCA variants harder to interpret and is not needed for the main paper.
solute_groups_to_run <- "All"

summary_rows <- purrr::pmap_dfr(response_specs, function(response_col, response_label) {
  purrr::map_dfr(solute_groups_to_run, function(group_i) {
    out <- run_response_pca(site_means, response_col, response_label, group_i)
    if (is.null(out)) return(NULL)
    as_tibble(out)
  })
})

readr::write_csv(summary_rows, file.path(res_dir, "storage_metric_ordination_summary.csv"))

message("Storage-metric ordination outputs written to: ", res_dir)
message("Storage-metric ordination figures written to: ", fig_dir)
