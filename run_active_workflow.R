# =============================================================================
# Run active HJA storage-chemistry analysis
# =============================================================================
# Runs the current paper analysis in order. Set HJA_SKIP_PREP=true to rerun only
# ordination, synchrony, and synthesis from existing prepared outputs.
#
# Optional: set HJA_RUN_STEPS to a comma-separated subset of:
# prep,pca,annual_chemistry,synchrony,summary,links,figures
# =============================================================================

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

repo_dir <- get_script_dir()
rscript <- file.path(R.home("bin"), "Rscript")

is_true <- function(x) {
  tolower(x) %in% c("1", "true", "yes", "y")
}

steps <- list(
  list(
    id = "prep",
    label = "Prepare chemistry tables and import storage-paper metrics",
    script = file.path(repo_dir, "01_data_prep", "1_MASTER_DATA_PREP.R")
  ),
  list(
    id = "pca",
    label = "Run all-solute chemistry ordinations with storage vectors",
    script = file.path(repo_dir, "04_PCA", "4i_storage_metric_ordination.R")
  ),
  list(
    id = "annual_chemistry",
    label = "Run annual stream-chemistry ordination with annual storage metrics",
    script = file.path(repo_dir, "04_PCA", "4j_annual_chemistry_storage_ordination.R")
  ),
  list(
    id = "synchrony",
    label = "Run annual pairwise synchrony and cluster-agreement analysis",
    script = file.path(repo_dir, "03_stats", "3r_storage_metric_synchrony.R")
  ),
  list(
    id = "summary",
    label = "Summarize current results",
    script = file.path(repo_dir, "05_synthesis", "5a_summarize_active_results.R")
  ),
  list(
    id = "links",
    label = "Synthesize storage-chemistry links across active analyses",
    script = file.path(repo_dir, "05_synthesis", "5b_synthesize_storage_chemistry_links.R")
  ),
  list(
    id = "figures",
    label = "Build prelim main-paper figures",
    script = file.path(repo_dir, "06_figures", "6a_prelim_main_figures.R")
  )
)

skip_prep <- is_true(Sys.getenv("HJA_SKIP_PREP", "false"))
requested_steps <- Sys.getenv("HJA_RUN_STEPS", "")
if (nzchar(requested_steps)) {
  requested_steps <- trimws(strsplit(requested_steps, ",", fixed = TRUE)[[1]])
  steps <- Filter(function(step) step$id %in% requested_steps, steps)
}
if (skip_prep) {
  steps <- Filter(function(step) step$id != "prep", steps)
}

if (length(steps) == 0) {
  stop("No analysis steps selected.")
}

run_step <- function(step) {
  if (!file.exists(step$script)) {
    stop("Missing analysis script: ", step$script)
  }

  start_time <- Sys.time()
  message("\n=== ", step$label, " ===")
  status <- system2(rscript, args = normalizePath(step$script))
  elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 2)

  if (!identical(status, 0L)) {
    stop("Analysis step failed: ", step$id, " (", step$script, ")")
  }
  message("Completed ", step$id, " in ", elapsed, " minutes.")
}

message("Running active HJA storage-chemistry analysis from: ", repo_dir)
invisible(lapply(steps, run_step))
message("\nActive analysis complete.")
