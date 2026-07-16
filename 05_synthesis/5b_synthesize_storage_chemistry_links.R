# =============================================================================
# 5b: Synthesize storage-chemistry links across active analyses
# =============================================================================
# Combines the site-level chemistry-profile ordinations, annual chemistry
# ordination, and pairwise synchrony analyses into a compact paper summary of
# which finalized storage-paper metrics organize water-quality behavior.
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
res_dir <- file.path(out_dir, "05_synthesis", "storage_chemistry_links")
fig_dir <- file.path(paths$fig_root, "05_synthesis", "storage_chemistry_links")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

require_file <- function(path) {
  if (!file.exists(path)) {
    stop("Missing required analysis output: ", path)
  }
  path
}

theme_storage_links <- function(base_size = 11) {
  if (exists("theme_hja")) {
    theme_hja(base_size = base_size)
  } else {
    ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(panel.grid = ggplot2::element_blank())
  }
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
response_label_order <- unname(response_labels[response_order])

analysis_family_labels <- c(
  site_profile_ordination = "Site-level chemistry profile",
  annual_ordination = "Annual chemistry",
  pairwise_agreement = "Pairwise synchrony/agreement"
)

label_response <- function(x) {
  out <- unname(response_labels[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

metric_domain <- function(metric) {
  dplyr::case_when(
    metric %in% c("RBI", "RCS", "FDC", "SD", "WB", "dynamic_storage_strength_z") ~
      "Dynamic and extended-dynamic storage",
    metric %in% c("BF", "flow_path_partitioning_z") ~
      "Flow-path partitioning",
    metric %in% c("DR", "Fyw", "MTT", "mobile_mixing_no_bf_z", "mobile_mixing_with_bf_z") ~
      "Mobile mixing and tracer storage",
    metric %in% c("unified_state_index") ~
      "Integrated storage state",
    metric %in% c("geology_pc1", "geology_pc2") ~
      "Physical template",
    TRUE ~ "Other"
  )
}

effect_note <- function(analysis_family, signed_effect) {
  dplyr::case_when(
    analysis_family != "pairwise_agreement" ~
      "PCA axis sign is arbitrary; interpret vector strength, not sign.",
    is.na(signed_effect) ~
      "Pairwise correlation sign unavailable.",
    signed_effect < 0 ~
      "Larger storage contrast is associated with lower chemistry similarity.",
    signed_effect > 0 ~
      "Larger storage contrast is associated with higher chemistry similarity.",
    TRUE ~
      "No directional pairwise association."
  )
}

pca_dir <- file.path(out_dir, "04_PCA", "storage_metric_ordination")
annual_pca_dir <- file.path(out_dir, "04_PCA", "annual_chemistry_storage_ordination")
sync_dir <- file.path(out_dir, "03_stats", "storage_metric_synchrony")

site_vector_files <- list.files(
  pca_dir,
  pattern = "_storage_vectors\\.csv$",
  full.names = TRUE
)
if (length(site_vector_files) == 0) {
  stop("No site-level storage-vector outputs found in: ", pca_dir)
}

site_vectors <- purrr::map_dfr(site_vector_files, readr::read_csv, show_col_types = FALSE) %>%
  filter(vector_type == "storage_metric") %>%
  transmute(
    analysis_family = "site_profile_ordination",
    analysis_family_label = unname(analysis_family_labels["site_profile_ordination"]),
    response,
    response_label,
    solute_group,
    storage_metric = variable,
    storage_label = label,
    storage_domain = metric_domain(storage_metric),
    n,
    statistic = "vector_r2",
    effect_strength = vector_r2,
    signed_effect = NA_real_,
    PC1_r,
    PC1_p,
    PC2_r,
    PC2_p,
    p = NA_real_
  )

annual_vectors <- readr::read_csv(
  require_file(file.path(annual_pca_dir, "annual_stream_chemistry_storage_vectors.csv")),
  show_col_types = FALSE
) %>%
  filter(vector_type == "storage_metric") %>%
  transmute(
    analysis_family = "annual_ordination",
    analysis_family_label = unname(analysis_family_labels["annual_ordination"]),
    response = "annual_stream_chemistry",
    response_label = label_response(response),
    solute_group = "All",
    storage_metric = variable,
    storage_label = label,
    storage_domain = metric_domain(storage_metric),
    n,
    statistic = "vector_r2",
    effect_strength = vector_r2,
    signed_effect = NA_real_,
    PC1_r,
    PC1_p,
    PC2_r,
    PC2_p,
    p = NA_real_
  )

pairwise_links <- readr::read_csv(
  require_file(file.path(sync_dir, "pairwise_synchrony_storage_metric_correlations.csv")),
  show_col_types = FALSE
) %>%
  filter(storage_metric %in% PAPER_FACING_STORAGE_METRICS) %>%
  transmute(
    analysis_family = "pairwise_agreement",
    analysis_family_label = unname(analysis_family_labels["pairwise_agreement"]),
    response = sync_response,
    response_label = sync_response_label,
    solute_group = "All",
    storage_metric,
    storage_label,
    storage_domain = metric_domain(storage_metric),
    n,
    statistic = "pearson_r",
    effect_strength = abs_r,
    signed_effect = r,
    PC1_r = NA_real_,
    PC1_p = NA_real_,
    PC2_r = NA_real_,
    PC2_p = NA_real_,
    p
  )

response_matrix <- bind_rows(site_vectors, annual_vectors, pairwise_links) %>%
  filter(is.finite(effect_strength)) %>%
  mutate(
    storage_metric = factor(storage_metric, levels = PAPER_FACING_STORAGE_METRICS),
    storage_label = factor(
      storage_label,
      levels = vapply(PAPER_FACING_STORAGE_METRICS, get_storage_label, character(1), short = TRUE)
    ),
    analysis_family = factor(analysis_family, levels = names(analysis_family_labels)),
    analysis_family_label = factor(analysis_family_label, levels = unname(analysis_family_labels)),
    response = factor(response, levels = response_order),
    response_label = factor(response_label, levels = rev(response_label_order)),
    effect_note = effect_note(as.character(analysis_family), signed_effect)
  ) %>%
  group_by(analysis_family, response) %>%
  arrange(desc(effect_strength), .by_group = TRUE) %>%
  mutate(within_response_rank = row_number()) %>%
  ungroup()

readr::write_csv(
  response_matrix,
  file.path(res_dir, "storage_chemistry_response_matrix.csv")
)

top_by_response <- response_matrix %>%
  arrange(analysis_family, response, within_response_rank) %>%
  group_by(analysis_family, analysis_family_label, response, response_label) %>%
  slice_head(n = 5) %>%
  ungroup()

readr::write_csv(
  top_by_response,
  file.path(res_dir, "storage_chemistry_top_by_response.csv")
)

top_one_by_response <- top_by_response %>%
  group_by(analysis_family, analysis_family_label, response, response_label) %>%
  slice_head(n = 1) %>%
  ungroup()

metric_consensus <- response_matrix %>%
  group_by(storage_metric, storage_label, storage_domain) %>%
  summarise(
    n_response_families = n_distinct(analysis_family),
    n_responses = n_distinct(response),
    mean_effect_strength = mean(effect_strength, na.rm = TRUE),
    median_effect_strength = median(effect_strength, na.rm = TRUE),
    max_effect_strength = max(effect_strength, na.rm = TRUE),
    n_top3_responses = sum(within_response_rank <= 3, na.rm = TRUE),
    n_top1_responses = sum(within_response_rank == 1, na.rm = TRUE),
    strongest_response = as.character(response_label[which.max(effect_strength)]),
    strongest_analysis_family = as.character(analysis_family_label[which.max(effect_strength)]),
    pairwise_negative_links = sum(analysis_family == "pairwise_agreement" & signed_effect < 0, na.rm = TRUE),
    pairwise_positive_links = sum(analysis_family == "pairwise_agreement" & signed_effect > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_top1_responses), desc(n_top3_responses), desc(mean_effect_strength))

readr::write_csv(
  metric_consensus,
  file.path(res_dir, "storage_chemistry_metric_consensus.csv")
)

domain_summary <- response_matrix %>%
  group_by(storage_domain) %>%
  summarise(
    n_metrics = n_distinct(storage_metric),
    n_responses = n_distinct(response),
    mean_effect_strength = mean(effect_strength, na.rm = TRUE),
    max_effect_strength = max(effect_strength, na.rm = TRUE),
    n_top3_responses = sum(within_response_rank <= 3, na.rm = TRUE),
    top_metric = as.character(storage_metric[which.max(effect_strength)]),
    top_response = as.character(response_label[which.max(effect_strength)]),
    .groups = "drop"
  ) %>%
  arrange(desc(n_top3_responses), desc(mean_effect_strength))

readr::write_csv(
  domain_summary,
  file.path(res_dir, "storage_chemistry_domain_summary.csv")
)

plot_tbl <- response_matrix %>%
  mutate(
    tile_label = ifelse(
      analysis_family == "pairwise_agreement",
      sprintf("%+.2f", signed_effect),
      sprintf("%.2f", effect_strength)
    )
  )

p_heat <- ggplot(plot_tbl, aes(x = storage_label, y = response_label, fill = effect_strength)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = tile_label), size = 2.9, color = "grey15") +
  facet_grid(analysis_family_label ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient(
    low = "#F7F7F2",
    high = "#2F6B9A",
    limits = c(0, 1),
    name = "Strength"
  ) +
  labs(
    x = "Storage-paper metric",
    y = NULL
  ) +
  theme_storage_links(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y = element_text(angle = 0, hjust = 0),
    legend.position = "bottom",
    panel.spacing.y = unit(0.5, "lines")
  )

ggsave(
  file.path(fig_dir, "storage_chemistry_response_matrix.png"),
  p_heat,
  width = 9.8,
  height = 7.2,
  dpi = PLOT_DPI,
  bg = "white"
)

top_lines <- top_one_by_response %>%
  mutate(
    estimate_txt = ifelse(
      statistic == "pearson_r",
      sprintf("r = %.2f", signed_effect),
      sprintf("R2 = %.2f", effect_strength)
    ),
    line = paste0(
      "- ", response_label, ": ", storage_metric, " (", storage_domain, "; ", estimate_txt, ")"
    )
  ) %>%
  pull(line)

consensus_lines <- metric_consensus %>%
  slice_head(n = 5) %>%
  mutate(
    line = paste0(
      "- ", storage_metric,
      ": top-ranked for ", n_top1_responses,
      " response(s), top-three for ", n_top3_responses,
      "; strongest in ", strongest_response
    )
  ) %>%
  pull(line)

domain_lines <- domain_summary %>%
  mutate(
    line = paste0(
      "- ", storage_domain,
      ": ", n_top3_responses,
      " top-three response link(s); strongest metric = ", top_metric,
      " for ", top_response
    )
  ) %>%
  pull(line)

summary_lines <- c(
  "Storage-chemistry link synthesis",
  paste0("Generated: ", Sys.Date()),
  "",
  "Scope",
  paste0("- Uses the storage-chemistry overlap years ", STORAGE_CHEMISTRY_YEAR_START, "-", STORAGE_CHEMISTRY_YEAR_END, "."),
  "- Combines unconstrained chemistry ordination vectors, annual chemistry PCA vectors, and pairwise storage-distance correlations.",
  "- Ordination values are vector R2 against PC1 and PC2; pairwise values are Pearson r between storage-metric distance and chemistry similarity/agreement.",
  "",
  "Top storage link by response",
  top_lines,
  "",
  "Consensus across storage metrics",
  consensus_lines,
  "",
  "Storage-domain summary",
  domain_lines,
  "",
  "Notes for reading these results",
  "- PCA axis signs are arbitrary, so ordination rows should be interpreted by vector strength and biplot direction, not by sign alone.",
  "- Negative pairwise r values support the storage-similarity hypothesis: catchments with more similar storage metrics have more similar chemistry behavior.",
  "- Positive pairwise r values flag responses where storage contrast does not translate into chemical dissimilarity in the expected direction."
)

writeLines(summary_lines, file.path(res_dir, "storage_chemistry_links_summary.txt"))

format_effect <- function(statistic, effect_strength, signed_effect) {
  ifelse(
    statistic == "pearson_r",
    sprintf("r = %+.2f", signed_effect),
    sprintf("link strength = %.2f", effect_strength)
  )
}

outline_top_lines <- top_one_by_response %>%
  mutate(
    estimate_txt = format_effect(statistic, effect_strength, signed_effect),
    line = paste0(
      "- ", response_label, ": strongest link with ", storage_label,
      " (", storage_domain, "; ", estimate_txt, ")."
    )
  ) %>%
  pull(line)

outline_consensus_lines <- metric_consensus %>%
  slice_head(n = 5) %>%
  mutate(
    line = paste0(
      "- ", storage_label, ": top-three storage link for ", n_top3_responses,
      " chemistry response(s); strongest for ", strongest_response, "."
    )
  ) %>%
  pull(line)

outline_domain_lines <- domain_summary %>%
  mutate(
    line = paste0(
      "- ", storage_domain, ": ", n_top3_responses,
      " top-three chemistry link(s); strongest current link is ",
      top_metric, " with ", top_response, "."
    )
  ) %>%
  pull(line)

outline_lines <- c(
  "# Storage-Water Quality Results Outline",
  "",
  paste0("Generated: ", Sys.Date()),
  "",
  "## Core Result",
  "",
  "- Finalized storage-paper metrics organize HJA stream-chemistry differences across site-level profiles, annual chemistry structure, and pairwise chemistry similarity.",
  "- Dynamic/extended-dynamic storage metrics provide the broadest current chemistry links; mobile-mixing and flow-path metrics add specific support for selected response families.",
  "- Pairwise storage-distance results are more variable than the ordination results, so they currently read best as supporting evidence rather than the main result.",
  "",
  "## Figure Shortlist",
  "",
  "- Fig1_storage_chemistry_response_matrix.png: main overview figure for storage-chemistry link strength across response families.",
  "- Fig2_top_storage_links_by_response.png: main narrowing figure for the strongest storage links by response.",
  "- Fig3_pairwise_storage_similarity_scatter.png: likely supplement unless the pairwise argument becomes central.",
  "- Fig4_annual_stream_chemistry_storage_pca.png: useful annual-chemistry support figure.",
  "",
  "## Strongest Link By Chemistry Response",
  "",
  outline_top_lines,
  "",
  "## Most Consistent Storage Metrics",
  "",
  outline_consensus_lines,
  "",
  "## Storage Domains",
  "",
  outline_domain_lines,
  "",
  "## Writing Next",
  "",
  "- Draft first results paragraph around Fig1: storage metrics align with multiple chemistry-response families, but the strength and metric identity differ by response type.",
  "- Draft second results paragraph around Fig2: RBI, FDC, SD, BF, MTT, and tracer-storage metrics separate different kinds of chemistry behavior.",
  "- Keep pairwise synchrony as supporting evidence until Fig3 is cleaned and the direction of each pairwise link is checked against the response matrix."
)

writeLines(outline_lines, file.path(res_dir, "storage_water_quality_results_outline.md"))

message("Storage-chemistry synthesis tables written to: ", res_dir)
message("Storage-chemistry synthesis figure written to: ", fig_dir)
