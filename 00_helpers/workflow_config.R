# =============================================================================
# workflow_config.R -- Centralized analysis configuration for HJA Stream project
# =============================================================================
# This script centralizes shared constants so every analysis script references
# the same storage metrics, outlet definitions, and helper utilities.
# Source this near the top of any script that needs consistent settings.
# =============================================================================

# =============================================================================
# TEMPORAL RANGE FOR PRIMARY ANALYSES
# =============================================================================
# All clustering, synchrony baselines, and site-mean calculations use this range
# Rationale: 1979 is first year with ≥50 site-solute combinations
#            2020 matches the storage-paper pre-fire endpoint; 2021-2024
#            chemistry can be retained for future post-fire sensitivity checks
ANALYSIS_YEAR_START <- 1979
ANALYSIS_YEAR_END   <- 2020

# Main annual storage-chemistry comparisons use the overlap between the long
# water-quality record and the finalized storage-paper annual metrics.
STORAGE_PAPER_YEAR_START <- 1997
STORAGE_PAPER_YEAR_END   <- 2020
STORAGE_CHEMISTRY_YEAR_START <- max(ANALYSIS_YEAR_START, STORAGE_PAPER_YEAR_START)
STORAGE_CHEMISTRY_YEAR_END   <- min(ANALYSIS_YEAR_END, STORAGE_PAPER_YEAR_END)

# =============================================================================
# STORAGE METRICS FROM THE STORAGE PAPER
# =============================================================================
# The water-quality paper now uses the finalized multi-metric storage framework
# from the storage paper rather than a single Q_dS_range_mm storage proxy.
PRIMARY_STORAGE_METRIC <- "SD"
PRIMARY_FLASHINESS_METRIC <- "RBI"
PRIMARY_MOBILE_STORAGE_METRIC <- "mobile_mixing_no_bf_z"
PRIMARY_STORAGE_FRAMEWORK_AXIS <- "unified_state_index"

STORAGE_PAPER_WORKFLOW_ROOT <- Sys.getenv(
  "HJA_STORAGE_FINAL_WORKFLOW_ROOT",
  unset = "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/final_workflow"
)

STORAGE_FRAMEWORK_SITE_FILE <- file.path(
  STORAGE_PAPER_WORKFLOW_ROOT,
  "outputs",
  "models",
  "unified_framework",
  "unified_framework_site_axes.csv"
)

STORAGE_PAPER_MASTER_SITE_FILE <- file.path(
  STORAGE_PAPER_WORKFLOW_ROOT,
  "outputs",
  "master",
  "master_site.csv"
)

STORAGE_PAPER_MASTER_ANNUAL_FILE <- file.path(
  STORAGE_PAPER_WORKFLOW_ROOT,
  "outputs",
  "master",
  "master_annual.csv"
)

STORAGE_FRAMEWORK_AXIS_METRICS <- c(
  "dynamic_storage_strength_z",
  "mobile_mixing_no_bf_z",
  "mobile_mixing_with_bf_z",
  "flow_path_partitioning_z",
  "unified_state_index",
  "geology_pc1",
  "geology_pc2"
)

STORAGE_PAPER_RAW_METRICS <- c(
  "RBI",
  "RCS",
  "FDC",
  "SD",
  "WB",
  "BF",
  "DR",
  "Fyw",
  "MTT"
)

PAPER_FACING_STORAGE_METRICS <- STORAGE_PAPER_RAW_METRICS

# Tier 1 predictors available for all 8 headwater sites
TIER1_PREDICTORS <- c(

  "dynamic_storage_strength_z", # Dynamic and extended-dynamic storage axis
  "mobile_mixing_no_bf_z",      # Tracer-only mobile-mixing axis
  "flow_path_partitioning_z",   # Ca-derived baseflow/flow-path partitioning
  "unified_state_index",        # Dynamic + mobile framework index
  "RBI",                        # Flashiness index
  "Elevation_mean_m",   # Topographic position
  "Slope_mean",         # Terrain steepness
  "Area_km2",           # Catchment size
  "Pyro_per",           # Pyroclastic geology %
  "Lava1_per",          # Basaltic lava %
  "Lava2_per",          # Secondary lava %
  "Ash_Per"             # Volcanic ash %
)

# =============================================================================
# SUPPLEMENTAL METRICS (Tier 2 - isotope subset, ~5-6 sites)
# =============================================================================
# These require isotope data not available at all sites
TIER2_PREDICTORS <- c(
  "MTT",                # Mean transit time
  "Fyw",                # Young water fraction
  "DR"                  # Damping ratio (isotope signal attenuation)
)

# Storage framework axes for sensitivity plots beyond the raw storage metrics
# used in the paper.
SUPPLEMENTAL_STORAGE_METRICS <- STORAGE_FRAMEWORK_AXIS_METRICS

# Raw storage-paper metrics for main paper analyses
SUPPLEMENTAL_HYDRO_METRICS <- PAPER_FACING_STORAGE_METRICS

# DEPRECATED: Do not use in new analyses
# - Q_dS_range_mm: pre-storage-paper storage proxy used in the exploratory workflow
# - WB_dS_range_mm: older water-quality water-balance proxy
# - DS_drawdown_*: Annual drawdown metrics (temporal mismatch with CQ data)
# - Abbott CQ-slope sync: Anti-correlated with other sync metrics

# Outlet definition (for highlighting hub-and-spoke behavior with GSLOOK)
OUTLET_SITE <- "GSLOOK"
OUTLET_SHAPE_TRIANGLE <- 17

# Convenience helpers ---------------------------------------------------------

#' Tag outlet streams in a site-level data frame
#' @param df data frame containing Stream_Name (or custom column)
#' @param col string column name for stream labels
#' @return df with logical column `is_outlet`
flag_outlet_stream <- function(df, col = "Stream_Name") {
  if (!is.data.frame(df) || !col %in% names(df)) return(df)
  df$is_outlet <- df[[col]] == OUTLET_SITE
  df$outlet_marker <- ifelse(df$is_outlet, "GSLOOK", "Other sites")
  df$outlet_marker <- factor(df$outlet_marker, levels = c("Other sites", "GSLOOK"))
  df
}

#' Tag outlet pairs in pairwise tables
#' Looks for columns named `site1`/`site2`, `Stream_Name.x`/`.y`, or similar.
#' @param df data frame
#' @param cols character vector of length 2 naming site columns
#' @return df with `is_outlet_pair` logical (TRUE if either member is outlet)
flag_outlet_pairs <- function(df, cols = NULL) {
  if (!is.data.frame(df)) return(df)
  guess_cols <- cols
  if (is.null(guess_cols)) {
    site_column_options <- list(
      c("Stream_Name.x", "Stream_Name.y"),
      c("site1", "site2"),
      c("site_a", "site_b"),
      c("Stream_Name_1", "Stream_Name_2")
    )
    for (site_column_option in site_column_options) {
      if (all(site_column_option %in% names(df))) {
        guess_cols <- site_column_option
        break
      }
    }
  }
  if (is.null(guess_cols)) return(df)
  df$is_outlet_pair <- df[[guess_cols[1]]] == OUTLET_SITE |
    df[[guess_cols[2]]] == OUTLET_SITE
  df$outlet_pair_marker <- ifelse(df$is_outlet_pair, "Contains GSLOOK", "All other pairs")
  df$outlet_pair_marker <- factor(df$outlet_pair_marker, levels = c("All other pairs", "Contains GSLOOK"))
  df
}

#' Filter to sites with complete data for supplemental metrics
#' @param df data frame
#' @param metrics vector of metric names (defaults to Tier 2 predictors)
#' @return list with `$full` complete cases and `$summary` availability table
split_by_metric_availability <- function(df, metrics = TIER2_PREDICTORS) {
  if (!is.data.frame(df)) {
    return(list(full = data.frame(), summary = data.frame()))
  }
  metrics <- metrics[metrics %in% names(df)]
  if (length(metrics) == 0) {
    return(list(full = data.frame(), summary = data.frame()))
  }
  complete_idx <- stats::complete.cases(df[, metrics, drop = FALSE])
  full <- df[complete_idx, , drop = FALSE]
  summary <- data.frame(
    metric = metrics,
    n_available = vapply(metrics, function(m) sum(!is.na(df[[m]])), numeric(1)),
    n_total = nrow(df),
    stringsAsFactors = FALSE
  )
  list(full = full, summary = summary)
}

HJA_VERBOSE <- tolower(Sys.getenv("HJA_VERBOSE", "false")) %in% c("1", "true", "yes", "y")
if (HJA_VERBOSE) {
  message("[config] Storage metrics used in paper: ", paste(PAPER_FACING_STORAGE_METRICS, collapse = ", "))
  message("[config] Analysis years: ", ANALYSIS_YEAR_START, "-", ANALYSIS_YEAR_END)
  message("[config] Storage-chemistry overlap: ", STORAGE_CHEMISTRY_YEAR_START, "-", STORAGE_CHEMISTRY_YEAR_END)
}

# -----------------------------------------------------------------------------
# Convenience helpers used across scripts
# -----------------------------------------------------------------------------

legend_bottom <- function() ggplot2::theme(legend.position = "bottom")
legend_right  <- function() ggplot2::theme(legend.position = "right")
guide_none    <- function() ggplot2::guides(color = "none", fill = "none", shape = "none", linetype = "none")

get_project_paths <- function() {
  base_dir <- Sys.getenv(
    "HJA_BOX_ROOT",
    unset = "/Users/sidneybush/Library/CloudStorage/Box-Box"
  )
  project_dir <- Sys.getenv(
    "HJA_WQ_PROJECT_DIR",
    unset = file.path(base_dir, "Sidney_Bush", "HJA_Water_Quality")
  )
  list(
    base_dir = base_dir,
    project_dir = project_dir,
    data_dir = file.path(project_dir, "data"),
    raw_dir = file.path(project_dir, "raw_data"),
    out_dir = file.path(project_dir, "outputs"),
    fig_root = file.path(project_dir, "exploratory_plots"),
    storage_final_workflow_root = STORAGE_PAPER_WORKFLOW_ROOT,
    storage_framework_site_file = STORAGE_FRAMEWORK_SITE_FILE,
    storage_paper_master_site_file = STORAGE_PAPER_MASTER_SITE_FILE,
    storage_paper_master_annual_file = STORAGE_PAPER_MASTER_ANNUAL_FILE
  )
}

standardize_storage_site <- function(site) {
  site_chr <- as.character(site)
  site_chr <- dplyr::case_when(
    site_chr %in% c("GSLOOK", "LOOK", "Lookout", "Look") ~ "Look",
    site_chr %in% c("GSMACK", "GSWSMC", "MACK", "Mack") ~ "Mack",
    grepl("^GSWS", site_chr) ~ sub("^GS", "", site_chr),
    TRUE ~ site_chr
  )
  site_chr <- ifelse(grepl("^WS[0-9]$", site_chr), sub("^WS", "WS0", site_chr), site_chr)
  site_chr
}

standardize_wq_stream <- function(site) {
  site_chr <- standardize_storage_site(site)
  dplyr::case_when(
    site_chr == "Look" ~ "GSLOOK",
    site_chr == "Mack" ~ "GSMACK",
    grepl("^WS", site_chr) ~ paste0("GS", site_chr),
    TRUE ~ site_chr
  )
}

get_primary_storage <- function(df, metric = PRIMARY_STORAGE_METRIC) {
  if (!metric %in% names(df)) return(rep(NA_real_, nrow(df)))
  df[[metric]]
}

get_flashiness <- function(df) {
  if ("RBI" %in% names(df)) return(df$RBI)
  rep(NA_real_, nrow(df))
}
