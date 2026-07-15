# =============================================================================
# Step 01: Master Data Prep Runner
# =============================================================================
# Order:
#   1a - CQ Data Wrangle
#   1b - Define Hydro Seasons
#   1c - CQ Rolling Analysis (wet75/dry150)
#   1d - Define Chemistry Clusters for annual pattern-agreement check
#   1e - Synchrony (Abbott + Wymore)
#   1f - Composite Synchrony
#   1g - Combine Master + Clean Tables (chemistry outputs)
#   1h - Import finalized storage-paper framework and join storage axes
#
# Cluster visualization/testing branches are archived; the active workflow only
# carries clusters forward as annual pattern-agreement annotations.
#
# Outputs written to Box: /Users/sidneybush/Library/CloudStorage/Box-Box/
#                         Sidney_Bush/HJA_Water_Quality/outputs/
# =============================================================================

suppressPackageStartupMessages({
 library(tidyverse)
})

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

script_dir <- get_script_dir()
find_repo_root <- function(start_dir) {
	current <- normalizePath(start_dir)
	sentinel <- ".git"
	repeat {
		helper_dir <- file.path(current, "00_helpers")
		git_dir    <- file.path(current, sentinel)
		if (dir.exists(helper_dir) || dir.exists(git_dir)) {
			return(current)
		}
		parent <- dirname(current)
		if (identical(parent, current)) {
			stop("Unable to locate project root from: ", start_dir)
		}
		current <- parent
	}
}

repo_dir <- find_repo_root(script_dir)

source(file.path(repo_dir, "00_helpers", "workflow_config.R"))
source(file.path(repo_dir, "00_helpers", "plot_prefs.R"))
theme_file <- file.path(repo_dir, "00_helpers", "plot_theme_set.R")
if (file.exists(theme_file)) source(theme_file)

paths <- get_project_paths()
out_dir <- paths$out_dir

source_repo <- function(...) {
	script_path <- file.path(repo_dir, ...)
	if (!file.exists(script_path)) {
		stop("Script not found: ", script_path)
	}
	env <- new.env(parent = globalenv())
	sys.source(script_path, envir = env)
}

message("=== STEP 01: DATA PREPARATION ===\n")

# 1a: Core CQ data wrangle
source_repo("01_data_prep", "1a_CQ_Data_Wrangle.R")

# 1b: Define hydrologic seasons
source_repo("01_data_prep", "1b_define_hydro_seasons.R")

# 1c: CQ Rolling Analysis (asymmetric windows: wet=75d, dry=150d)
source_repo("01_data_prep", "1c_CQ_Rolling_Analysis_windows_wet75_dry150.R")

# 1d: Define chemistry clusters
cluster_output <- file.path(out_dir, "ClusterStreams_allSolutes_modalClusters.csv")
if (requireNamespace("dtw", quietly = TRUE) && requireNamespace("dtwclust", quietly = TRUE)) {
	source_repo("01_data_prep", "1e_clusters_define.R")
} else if (file.exists(cluster_output)) {
	message("Skipping 1e_clusters_define.R because packages 'dtw' and/or 'dtwclust' are not installed; using existing cluster output: ", cluster_output)
} else {
	stop("Packages 'dtw' and 'dtwclust' are required to build chemistry clusters, and no existing cluster output was found: ", cluster_output)
}

# 1e: Abbott + Wymore synchrony calculations
source_repo("01_data_prep", "1f_sync_abbott_and_wymore.R")

# 1f: Composite synchrony aggregation
source_repo("01_data_prep", "1g_sync_composite.R")

# 1g: Combine chemistry outputs into master + clean tables
source_repo("01_data_prep", "1h_combine_master_clean.R")

# 1h: Import storage-paper final workflow metrics and Figure 7 framework axes
source_repo("01_data_prep", "1j_import_storage_framework.R")
