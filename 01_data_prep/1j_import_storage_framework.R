#!/usr/bin/env Rscript
# =============================================================================
# Step 1j: Import finalized storage-paper framework
# =============================================================================
# Uses the storage paper final workflow as the canonical source for HJA storage
# metrics, Figure 7 storage axes, and catchment characteristics. This step
# leaves legacy water-quality outputs untouched and writes storage-framework
# joined versions for the next-paper analyses.
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

paths <- get_project_paths()
out_dir <- paths$out_dir
data_dir <- paths$data_dir
storage_copy_dir <- file.path(data_dir, "storage_paper_framework")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(storage_copy_dir, recursive = TRUE, showWarnings = FALSE)

message("=== STEP 1j: IMPORT STORAGE-PAPER FRAMEWORK ===")
message("Water-quality project: ", paths$project_dir)
message("Storage-paper workflow: ", paths$storage_final_workflow_root)

required_files <- c(
  paths$storage_framework_site_file,
  paths$storage_paper_master_site_file,
  paths$storage_paper_master_annual_file
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required storage-paper files:\n", paste(missing_files, collapse = "\n"))
}

copy_storage_input <- function(from, name) {
  to <- file.path(storage_copy_dir, name)
  ok <- file.copy(from, to, overwrite = TRUE)
  if (!ok) stop("Could not copy storage-paper file to: ", to)
  to
}

copied_files <- c(
  unified_framework_site_axes = copy_storage_input(
    paths$storage_framework_site_file,
    "storage_paper_unified_framework_site_axes.csv"
  ),
  master_site = copy_storage_input(
    paths$storage_paper_master_site_file,
    "storage_paper_master_site.csv"
  ),
  master_annual = copy_storage_input(
    paths$storage_paper_master_annual_file,
    "storage_paper_master_annual.csv"
  )
)

storage_site_raw <- readr::read_csv(paths$storage_paper_master_site_file, show_col_types = FALSE) %>%
  mutate(
    site = standardize_storage_site(site),
    Stream_Name = standardize_wq_stream(site)
  )

storage_axes <- readr::read_csv(paths$storage_framework_site_file, show_col_types = FALSE) %>%
  mutate(
    site = standardize_storage_site(site),
    Stream_Name = standardize_wq_stream(site)
  )

storage_site_metrics <- storage_site_raw %>%
  transmute(
    site,
    Stream_Name,
    n_years_storage_paper = n_years,
    storage_paper_year_mean = year_mean,
    RBI = RBI_mean,
    RCS = RCS_mean,
    FDC = FDC_mean,
    SD = SD_mean,
    WB = WB_mean,
    BF = BF_mean,
    DR = DR,
    Fyw = Fyw,
    MTT = MTT,
    Q5norm = Q5norm_mean,
    CV_Q5norm = CV_Q5norm_mean,
    T_7DMax = T_7DMax_mean,
    Q_7Q5 = Q_7Q5_mean,
    Pws = Pws_mean,
    precip_nov_may_mm = precip_nov_may_mm_mean,
    Area_km2,
    Elevation_min_m,
    Elevation_mean_m,
    Elevation_max_m,
    Slope_mean,
    Slope_Std,
    Aspect_Mean_deg,
    Harvest,
    Age,
    Landslide_Young,
    Landslide_Mod,
    Landslide_Old,
    Landslide_Total,
    Lava1_per,
    Lava2_per,
    Ash_Per,
    Pyro_per
  )

axis_cols <- c(
  "site",
  "Stream_Name",
  "dynamic_storage_strength",
  "mobile_mixing",
  "mobile_mixing_with_bf",
  "mobile_mixing_no_bf",
  "BF_mean",
  "n_dynamic_components",
  "n_mobile_components_with_bf",
  "n_mobile_components_no_bf",
  STORAGE_FRAMEWORK_AXIS_METRICS,
  "geology_class",
  "geomorphology_class",
  "Landslide_Total",
  "Lava1_per",
  "Lava2_per",
  "Ash_Per",
  "Pyro_per"
)
axis_cols <- intersect(axis_cols, names(storage_axes))

storage_site_framework <- storage_site_metrics %>%
  left_join(
    storage_axes %>% select(all_of(axis_cols), -any_of("Stream_Name")),
    by = "site",
    suffix = c("", "_axis")
  ) %>%
  mutate(Stream_Name = standardize_wq_stream(site)) %>%
  relocate(Stream_Name, .after = site)

readr::write_csv(storage_site_framework, file.path(out_dir, "HJA_storage_framework_site.csv"))

storage_annual_framework <- readr::read_csv(paths$storage_paper_master_annual_file, show_col_types = FALSE) %>%
  mutate(
    site = standardize_storage_site(site),
    Stream_Name = standardize_wq_stream(site),
    water_year = as.integer(year)
  ) %>%
  select(
    site,
    Stream_Name,
    water_year,
    storage_paper_year = year,
    any_of(c("RCS", "RBI", "SD", "FDC", "Q99", "Q50", "Q01", "Q5norm", "CV_Q5norm", "BF", "WB", "T_7DMax", "Q_7Q5", "Pws", "precip_nov_may_mm"))
  ) %>%
  left_join(
    storage_site_framework %>%
      select(
        site,
        any_of(c(
          "dynamic_storage_strength",
          "mobile_mixing",
          "mobile_mixing_with_bf",
          "mobile_mixing_no_bf",
          "n_dynamic_components",
          "n_mobile_components_with_bf",
          "n_mobile_components_no_bf",
          STORAGE_FRAMEWORK_AXIS_METRICS,
          "geology_class",
          "geomorphology_class",
          "Area_km2",
          "Elevation_min_m",
          "Elevation_mean_m",
          "Elevation_max_m",
          "Slope_mean",
          "Slope_Std",
          "Aspect_Mean_deg",
          "Harvest",
          "Age",
          "Landslide_Young",
          "Landslide_Mod",
          "Landslide_Old",
          "Landslide_Total",
          "Lava1_per",
          "Lava2_per",
          "Ash_Per",
          "Pyro_per",
          "DR",
          "Fyw",
          "MTT"
        ))
      ),
    by = "site"
  )

readr::write_csv(storage_annual_framework, file.path(out_dir, "HJA_storage_framework_annual.csv"))

rename_conflicts_for_join <- function(df, incoming, by_cols) {
  conflicts <- intersect(names(df), setdiff(names(incoming), by_cols))
  if (length(conflicts) == 0) return(df)
  rename_map <- setNames(conflicts, paste0("legacy_", conflicts))
  rename(df, !!!rename_map)
}

add_site_key <- function(df) {
  if ("Stream_Name" %in% names(df)) {
    df %>% mutate(site = standardize_storage_site(Stream_Name))
  } else if ("site" %in% names(df)) {
    df %>% mutate(site = standardize_storage_site(site))
  } else {
    stop("Input file does not contain Stream_Name or site.")
  }
}

join_storage_framework <- function(input_name, output_name) {
  input_file <- file.path(out_dir, input_name)
  if (!file.exists(input_file)) {
    message("Skipping missing output: ", input_name)
    return(NULL)
  }

  df <- readr::read_csv(input_file, show_col_types = FALSE) %>%
    add_site_key()

  has_year <- "water_year" %in% names(df)
  if (has_year) {
    incoming <- storage_annual_framework %>% select(-Stream_Name)
    by_cols <- c("site", "water_year")
  } else {
    incoming <- storage_site_framework %>% select(-Stream_Name)
    by_cols <- "site"
  }

  df_join <- df %>%
    rename_conflicts_for_join(incoming, by_cols) %>%
    left_join(incoming, by = by_cols)

  out_file <- file.path(out_dir, output_name)
  readr::write_csv(df_join, out_file)

  tibble(
    input = input_name,
    output = output_name,
    rows = nrow(df_join),
    cols = ncol(df_join),
    n_sites = n_distinct(df_join$site),
    n_with_dynamic_axis = sum(is.finite(df_join$dynamic_storage_strength_z)),
    n_with_mobile_axis = sum(is.finite(df_join$mobile_mixing_no_bf_z)),
    n_with_unified_index = sum(is.finite(df_join$unified_state_index))
  )
}

join_targets <- tribble(
  ~input, ~output,
  "HJA_master_site_means.csv", "HJA_master_site_means_storage_framework.csv",
  "HJA_clean_site_means.csv", "HJA_clean_site_means_storage_framework.csv",
  "HJA_master_annual.csv", "HJA_master_annual_storage_framework.csv",
  "HJA_clean_annual.csv", "HJA_clean_annual_storage_framework.csv",
  "HJA_master_seasonal.csv", "HJA_master_seasonal_storage_framework.csv",
  "HJA_clean_seasonal.csv", "HJA_clean_seasonal_storage_framework.csv",
  "HJA_master_rolling_windows.csv", "HJA_master_rolling_windows_storage_framework.csv",
  "HJA_clean_windows.csv", "HJA_clean_windows_storage_framework.csv"
)

join_summary <- purrr::map2_dfr(join_targets$input, join_targets$output, join_storage_framework)

storage_summary <- tibble(
  item = c(
    "storage_framework_site_rows",
    "storage_framework_annual_rows",
    "copied_storage_files",
    "joined_output_files"
  ),
  value = c(
    as.character(nrow(storage_site_framework)),
    as.character(nrow(storage_annual_framework)),
    as.character(length(copied_files)),
    as.character(nrow(join_summary))
  )
)

readr::write_csv(join_summary, file.path(out_dir, "storage_framework_join_summary.csv"))
readr::write_csv(storage_summary, file.path(out_dir, "storage_framework_import_summary.csv"))

message("Storage framework site rows: ", nrow(storage_site_framework))
message("Storage framework annual rows: ", nrow(storage_annual_framework))
message("Joined output files: ", nrow(join_summary))
message("Copied canonical storage inputs to: ", storage_copy_dir)
message("=== STEP 1j COMPLETE ===")
