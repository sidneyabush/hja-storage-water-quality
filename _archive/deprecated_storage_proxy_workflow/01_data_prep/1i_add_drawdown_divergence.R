# =============================================================================
# STEP 01i: ADD DS_DRAWDOWN DIVERGENCE TO SYNCHRONY DATA
# =============================================================================
# Purpose: Calculate inter-site DS_drawdown divergence for pair-level analysis
#
# This script:
# 1. Loads annual DS_drawdown data
# 2. Calculates pair-level divergence metrics (mean, SD, range)
# 3. Merges with existing synchrony data
# 4. Updates annual synchrony files with drawdown divergence
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

rm(list = ls())

# Source helpers
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
try(source(file.path(repo_dir, "00_helpers", "plot_prefs.R")), silent = TRUE)
try(source(file.path(repo_dir, "00_helpers", "workflow_config.R")), silent = TRUE)

# Paths
paths <- get_project_paths()
out_dir <- paths$out_dir
data_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality/data"

invisible(NULL)

# =============================================================================
# LOAD DATA
# =============================================================================
invisible(NULL)

# DS drawdown (annual by site)
ds_draw <- read_csv(file.path(data_dir, "DS_drawdown_annual.csv"),
                    show_col_types = FALSE) %>%
  rename(Stream_Name = SITECODE, water_year = waterYear, DS_drawdown = DS_sum)

invisible(NULL)

# Annual synchrony data (try multiple possible filenames)
sync_files <- c("HJA_annual_synchrony.csv", "HJA_composite_synchrony_annual.csv",
                "HJA_Abbott_synchrony_windows.csv")
sync_file <- NULL
for (f in sync_files) {
  if (file.exists(file.path(out_dir, f))) {
    sync_file <- f
    break
  }
}

if (is.null(sync_file)) {
  stop("Could not find annual synchrony file. Checked: ", paste(sync_files, collapse = ", "))
}

sync_annual <- read_csv(file.path(out_dir, sync_file),
                        show_col_types = FALSE)

invisible(NULL)

# =============================================================================
# CALCULATE PAIR-LEVEL DS_DRAWDOWN DIVERGENCE
# =============================================================================
invisible(NULL)

# Create all site pairs with drawdown data
ds_pairs <- ds_draw %>%
  inner_join(ds_draw, by = "water_year", suffix = c("1", "2")) %>%
  filter(Stream_Name1 < Stream_Name2) %>%  # Unique pairs only
  mutate(
    DS_drawdown_diff = abs(DS_drawdown1 - DS_drawdown2),
    DS_drawdown_mean_pair = (DS_drawdown1 + DS_drawdown2) / 2
  ) %>%
  select(Stream_Name1, Stream_Name2, water_year,
         DS_drawdown1, DS_drawdown2, DS_drawdown_diff, DS_drawdown_mean_pair)

invisible(NULL)

# Aggregate by pair across all years
ds_pairs_summary <- ds_pairs %>%
  group_by(Stream_Name1, Stream_Name2) %>%
  summarise(
    DS_drawdown_diff_mean = mean(DS_drawdown_diff, na.rm = TRUE),
    DS_drawdown_diff_sd = sd(DS_drawdown_diff, na.rm = TRUE),
    DS_drawdown_diff_min = min(DS_drawdown_diff, na.rm = TRUE),
    DS_drawdown_diff_max = max(DS_drawdown_diff, na.rm = TRUE),
    DS_drawdown_corr = cor(DS_drawdown1, DS_drawdown2, use = "complete.obs"),
    n_years_drawdown = n(),
    .groups = "drop"
  )

invisible(NULL)

# =============================================================================
# MERGE WITH ANNUAL SYNCHRONY
# =============================================================================
invisible(NULL)

# Backup original
backup_file <- file.path(out_dir, paste0(tools::file_path_sans_ext(sync_file), "_BACKUP_before_drawdown.csv"))
if (!file.exists(backup_file)) {
  write_csv(sync_annual, backup_file)
  invisible(NULL)
}

# Merge annual pair-level drawdown
sync_annual_updated <- sync_annual %>%
  left_join(ds_pairs, by = c("Stream1" = "Stream_Name1",
                              "Stream2" = "Stream_Name2",
                              "water_year"))

# Also merge long-term divergence summary
sync_annual_updated <- sync_annual_updated %>%
  left_join(ds_pairs_summary, by = c("Stream1" = "Stream_Name1",
                                      "Stream2" = "Stream_Name2"))

invisible(NULL)
invisible(
    sum(!is.na(sync_annual_updated$DS_drawdown_diff)),
    "of", nrow(sync_annual_updated), "\n\n")

# Save
write_csv(sync_annual_updated, file.path(out_dir, sync_file))
