# =============================================================================
# STEP 2a: CLUSTER CREATION + CALENDAR YEAR ASSIGNMENTS
# =============================================================================
# - Create climatological reference clusters
# - Assign each calendar-year series to nearest climatology centroid
# - Compute cluster stability + modal clusters
# - Write all cluster CSVs used by downstream scripts
# =============================================================================

suppressPackageStartupMessages({
  library(dtw)
  library(dtwclust)
  library(tidyverse)
  library(reshape2)
})

rm(list = ls())

# Source shared workflow settings and plot preferences
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
source(file.path(repo_dir, "00_helpers", "workflow_config.R"))
source(file.path(repo_dir, "00_helpers", "plot_prefs.R"))

# =============================================================================
# TEMPORAL RANGE FOR ANALYSIS (from workflow_config.R)
# =============================================================================
message("Clustering temporal range: ", ANALYSIS_YEAR_START, "-", ANALYSIS_YEAR_END)

# =============================================================================
# Solute groups defined in plot_prefs.R (GEO_SOLUTES, BIO_SOLUTES)
# =============================================================================
geo_solutes <- GEO_SOLUTES
bio_solutes <- BIO_SOLUTES

# =============================================================================
# PATHS: data directory + outputs 
# =============================================================================
data_dir   <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality/data"
output_dir <- file.path(dirname(data_dir), "outputs")

# Save original working directory and restore at end
original_wd <- getwd()
on.exit(setwd(original_wd), add = TRUE)
setwd(data_dir)

# =============================================================================
# 1. CLIMATOLOGICAL REFERENCE WORKFLOW (NO WATER YEAR)
# =============================================================================

# read climatological monthly means
chem_clim <- readr::read_csv(
  file.path(output_dir, "HJA_CF002_monthly_means.csv"),
  show_col_types = FALSE
)
colnames(chem_clim)[3:14] <- as.character(seq(1, 12, 1))

# pivot longer, normalize within variables and sites
month_conc_norm_clim <- chem_clim %>%
  pivot_longer(cols = 3:14, values_to = "concentration", names_to = "month") %>%
  dplyr::group_by(variable, Stream_Name) %>%
  dplyr::mutate(norm_conc = scale(concentration)) %>%
  dplyr::ungroup()

month_conc_norm_clim$month <- as.numeric(month_conc_norm_clim$month)

# cast to wide matrix: rows = months, cols = site–solute
month_cast_clim <- dcast(
  month_conc_norm_clim,
  formula = month ~ Stream_Name + variable,
  value.var = "norm_conc"
)

month_norm_clim   <- month_cast_clim
month_norm_t_clim <- as.data.frame(t(month_norm_clim[, 2:ncol(month_norm_clim)]))

# =============================================================================
# (Optional) Cluster Validity Indices for k = 2–6
# =============================================================================
index      <- c("Sil","D","COP","DB","DBstar","CH","SF")
cvi_ideal  <- c("max","max","min","min","min","max","max")
cvi_index  <- data.frame(index, cvi_ideal)

clust.dba_cv_clim <- tsclust(
  month_norm_t_clim,
  type        = "partitional",
  centroid    = "dba",
  distance    = "dtw",
  window.size = 1L,
  k           = 2L:6L,
  seed        = 8
)

cvi_df_clim            <- lapply(clust.dba_cv_clim, cvi)
cluster_stats_clim      <- do.call(rbind, cvi_df_clim)
cluster_stats_melt_clim <- melt(cluster_stats_clim)

cluster_stats_melt_clim$Var1 <- cluster_stats_melt_clim$Var1 + 1
colnames(cluster_stats_melt_clim)[2] <- "index"

cluster_stats_melt_clim <- merge(cluster_stats_melt_clim, cvi_index, by = "index")
colnames(cluster_stats_melt_clim)[2] <- "number_of_clusters"
cluster_stats_melt_clim$cvi_goal <- paste0(
  cluster_stats_melt_clim$index,
  "-",
  cluster_stats_melt_clim$cvi_ideal
)

# =============================================================================
# FINAL CLIMATOLOGICAL CLUSTERING WITH k = 4
# =============================================================================
k_clim <- 4L

clust.dba_clim <- tsclust(
  month_norm_t_clim,
  type        = "partitional",
  centroid    = "dba",
  distance    = "dtw",
  window.size = 1L,  # Sakoe-Chiba band to match reference methodology
  k           = k_clim,
  seed        = 8
)

# pull out individual lines
mydata_clim <- clust.dba_clim@datalist

# convert to df
data_df_clim <- t(do.call(cbind, mydata_clim))

# pull out cluster information and add site and chemical back in
month_clusters_clim <- as.data.frame(cbind(data_df_clim, clust.dba_clim@cluster))
colnames(month_clusters_clim) <- c(paste(seq(1:12)), "Cluster_climRef")
month_clusters_clim <- tibble::rownames_to_column(month_clusters_clim, "Site")

month_clusters_clim <- month_clusters_clim %>% 
  tidyr::extract(
    Site,
    into  = c("Stream_Name", "chemical"),
    regex = "(.*)_([^_]+)$"
  )

# ADD water_year COLUMN (NA) SO IT MATCHES WY TABLES
month_clusters_clim <- month_clusters_clim %>%
  dplyr::mutate(water_year = NA_integer_) %>%
  dplyr::select(
    water_year,
    Stream_Name,
    chemical,
    dplyr::all_of(as.character(1:12)),
    Cluster_climRef
  )

# =============================================================================
# RENUMBER CLIMATOLOGICAL CLUSTERS ONCE: Cluster_climRef = 1 IS LARGEST
# =============================================================================
cluster_sizes_clim <- month_clusters_clim %>%
  dplyr::group_by(Cluster_climRef) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(n)) %>%
  dplyr::mutate(new_cluster = dplyr::row_number())

message("\nCluster size reordering (old -> new):\n")
print(cluster_sizes_clim)

# map old -> new labels in climatology table
month_clusters_clim <- month_clusters_clim %>%
  dplyr::left_join(
    cluster_sizes_clim %>% dplyr::select(Cluster_climRef, new_cluster),
    by = "Cluster_climRef"
  ) %>%
  dplyr::mutate(Cluster_climRef = new_cluster) %>%
  dplyr::select(-new_cluster)

# reorder centroids to match new cluster labels
# OLD cluster IDs are in cluster_sizes_clim$Cluster_climRef (before renumbering)
# We need centroids in NEW order: position i should have old cluster cluster_sizes_clim$Cluster_climRef[i]
centroids_clim_list_raw <- lapply(clust.dba_clim@centroids, function(x) as.numeric(x))
old_cluster_order <- cluster_sizes_clim$Cluster_climRef  # e.g., c(2, 1, 4, 3) means new cluster 1 = old 2
centroids_clim_list <- centroids_clim_list_raw[old_cluster_order]
k_clim <- length(centroids_clim_list)

# =============================================================================
# WRITE CLIMATOLOGICAL WORKFLOW CSV (WITH water_year COLUMN)
# =============================================================================
readr::write_csv(
  month_clusters_clim,
  file.path(output_dir, "ClusterStreams_allSolutes.csv")
)

# =============================================================================
# 2. CALENDAR-YEAR WORKFLOW (ASSIGN EACH CY SERIES TO CLIMATOLOGICAL REGIMES)
# =============================================================================
# IMPORTANT: We use CALENDAR YEAR, not water year, for cluster assignment.
# Each cluster should be assigned to the actual calendar year the data came from.
# FILTERED TO: ANALYSIS_YEAR_START - ANALYSIS_YEAR_END

chem_cy <- readr::read_csv(
  file.path(output_dir, "HJA_CF002_monthly_means_byCY.csv"),
  show_col_types = FALSE
)

# APPLY TEMPORAL FILTER FROM WORKFLOW CONFIG
chem_cy <- chem_cy %>%
  dplyr::filter(calendar_year >= ANALYSIS_YEAR_START & calendar_year <= ANALYSIS_YEAR_END)

# Identify month columns
month_cols_cy <- setdiff(names(chem_cy), c("Stream_Name", "variable", "calendar_year"))

# LONG → STANDARDIZE → WIDE MATRIX FOR CLUSTERING
# NOTE: Z-score normalize by variable & stream ONLY (across all calendar years)
# This matches the reference methodology: group_by(stream) %>% mutate(scale(...))
# Each year gets normalized using the mean/SD calculated across ALL years
month_conc_norm_cy <- chem_cy %>%
  tidyr::pivot_longer(
    cols      = tidyselect::all_of(month_cols_cy),
    names_to  = "month",
    values_to = "concentration"
  ) %>%
  dplyr::filter(!is.na(concentration)) %>%
  dplyr::group_by(variable, Stream_Name) %>%
  dplyr::mutate(
    norm_conc = {
      mu  <- mean(concentration, na.rm = TRUE)
      sdv <- stats::sd(concentration, na.rm = TRUE)
      if (is.na(sdv) || sdv == 0) 0 else (concentration - mu) / sdv
    }
  ) %>%
  dplyr::ungroup()

month_conc_norm_cy$month <- as.numeric(month_conc_norm_cy$month)

month_cast_cy <- reshape2::dcast(
  month_conc_norm_cy,
  formula = month ~ calendar_year + Stream_Name + variable,
  value.var = "norm_conc"
)

month_norm_cy   <- month_cast_cy
month_norm_t_cy <- as.data.frame(t(month_norm_cy[, 2:ncol(month_norm_cy)]))

# Drop any series with missing values
valid_rows_cy <- stats::complete.cases(month_norm_t_cy)
if (!all(valid_rows_cy)) {
}
month_norm_t_cy <- month_norm_t_cy[valid_rows_cy, , drop = FALSE]

# =============================================================================
# 2A. ASSIGN EACH CY SERIES TO NEAREST CLIMATOLOGICAL REGIME
# =============================================================================
# NOTE: Using Sakoe-Chiba band with window.size=1L to match reference methodology
# Reference: dtwDist(..., window.type="sakoechiba", window.size=1L)
dtw_to_clim_centroids <- function(series_vec, centroids_list) {
  v <- as.numeric(series_vec)
  sapply(centroids_list, function(cntr) {
    dtw(v, cntr,
        distance.only = TRUE,
        window.type = "sakoechiba",
        window.size = 1L)$distance
  })
}

# distance matrix: rows = CY series, cols = climatology clusters 1..k_clim
dist_mat_clim <- t(apply(
  month_norm_t_cy,
  1,
  dtw_to_clim_centroids,
  centroids_list = centroids_clim_list
))

Cluster_climRef_idx <- max.col(-dist_mat_clim)      # index of min distance
dist_climRef        <- apply(dist_mat_clim, 1, min)

# =============================================================================
# 2B. BUILD CY TABLE USING CLIMATOLOGY-REFERENCE CLUSTERS (FIXED IDS)
# =============================================================================
# NOTE: Clusters are now assigned to CALENDAR YEAR, not water year

data_df_cy <- month_norm_t_cy
colnames(data_df_cy) <- as.character(seq_len(ncol(data_df_cy)))
month_cols_clust_cy <- colnames(data_df_cy)

month_clusters_cy <- as.data.frame(data_df_cy) %>%
  tibble::rownames_to_column("Site") %>%
  tidyr::extract(
    Site,
    into  = c("calendar_year", "Stream_Name", "chemical"),
    regex = "^(.+?)_(.+?)_([^_]+)$"
  ) %>%
  dplyr::mutate(
    calendar_year   = as.integer(calendar_year),
    Cluster_climRef = Cluster_climRef_idx,
    dist_climRef    = dist_climRef
  ) %>%
  dplyr::select(
    calendar_year,
    Stream_Name,
    chemical,
    dplyr::all_of(month_cols_clust_cy),
    Cluster_climRef,
    dist_climRef
  )

# =============================================================================
# 3. CLUSTER STABILITY + MODAL CLUSTERS BY CALENDAR YEAR
#    (USING FIXED CLIMATOLOGY-REFERENCE CLUSTERS)
# =============================================================================

cluster_changes_cy <- month_clusters_cy %>%
  dplyr::arrange(Stream_Name, chemical, calendar_year) %>%
  dplyr::group_by(Stream_Name, chemical) %>%
  dplyr::mutate(
    Cluster_climRef_prev = dplyr::lag(Cluster_climRef),
    changed              = Cluster_climRef != Cluster_climRef_prev
  ) %>%
  dplyr::ungroup()

cluster_change_summary_cy <- cluster_changes_cy %>%
  dplyr::group_by(Stream_Name, chemical) %>%
  dplyr::summarise(
    n_years     = dplyr::n(),
    n_changes   = sum(changed, na.rm = TRUE),
    prop_change = ifelse(n_years > 1, n_changes / (n_years - 1), NA_real_),
    .groups     = "drop"
  ) %>%
  dplyr::mutate(stability = 1 - prop_change)

print(cluster_change_summary_cy)

# Save cluster stability metrics to CSV (used by characterization / by-site scripts)
readr::write_csv(
  cluster_change_summary_cy,
  file.path(output_dir, "ClusterStreams_stability_metrics.csv")
)

# Helper for modal cluster
mode_val <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# MODAL CLIMATOLOGY-REFERENCE CLUSTER PER STREAM × SOLUTE (across calendar years)
cluster_modal <- month_clusters_cy %>%
  dplyr::group_by(Stream_Name, chemical) %>%
  dplyr::summarise(
    Cluster_mode = mode_val(Cluster_climRef),
    n_years      = dplyr::n(),
    .groups      = "drop"
  )

# =============================================================================
# SECOND REORDERING: ENSURE CLUSTER 1 = MOST COMMON MODAL CLUSTER
# BUT PRESERVE ALL 4 CLUSTERS FROM CLIMATOLOGY
# =============================================================================

# First, get all clusters from climatology (should be 1-4)
all_clim_clusters <- sort(unique(month_clusters_clim$Cluster_climRef))

# Then get modal cluster frequencies
modal_cluster_sizes <- cluster_modal %>%
  dplyr::group_by(Cluster_mode) %>%
  dplyr::summarise(n_modal = dplyr::n(), .groups = "drop") %>%
  # Add any missing clusters from climatology with n_modal = 0
  dplyr::full_join(
    tibble(Cluster_mode = all_clim_clusters),
    by = "Cluster_mode"
  ) %>%
  dplyr::mutate(n_modal = dplyr::coalesce(n_modal, 0L)) %>%
  dplyr::arrange(dplyr::desc(n_modal)) %>%
  dplyr::mutate(final_cluster = dplyr::row_number())


# Apply final reordering to ALL cluster columns
# Remap: cluster_modal$Cluster_mode
cluster_modal <- cluster_modal %>%
  dplyr::left_join(
    modal_cluster_sizes %>% dplyr::select(Cluster_mode, final_cluster),
    by = "Cluster_mode"
  ) %>%
  dplyr::mutate(Cluster_mode = final_cluster) %>%
  dplyr::select(-final_cluster)

# Remap: month_clusters_clim$Cluster_climRef
month_clusters_clim <- month_clusters_clim %>%
  dplyr::left_join(
    modal_cluster_sizes %>% dplyr::select(Cluster_mode, final_cluster),
    by = c("Cluster_climRef" = "Cluster_mode")
  ) %>%
  dplyr::mutate(Cluster_climRef = final_cluster) %>%
  dplyr::select(-final_cluster)

# Remap: month_clusters_cy$Cluster_climRef
month_clusters_cy <- month_clusters_cy %>%
  dplyr::left_join(
    modal_cluster_sizes %>% dplyr::select(Cluster_mode, final_cluster),
    by = c("Cluster_climRef" = "Cluster_mode")
  ) %>%
  dplyr::mutate(Cluster_climRef = final_cluster) %>%
  dplyr::select(-final_cluster)

# Reorder centroids one more time to match final cluster IDs
old_to_final_map <- modal_cluster_sizes$Cluster_mode  # e.g., c(2, 1, 4, 3)
centroids_clim_list <- centroids_clim_list[old_to_final_map]

# attach modal cluster back to CY table (for convenience)
month_clusters_cy <- month_clusters_cy %>%
  dplyr::left_join(cluster_modal, by = c("Stream_Name", "chemical"))

# =============================================================================
# 3A. MAJORITY CLUSTERS BY SITE / STREAM / GEO vs BIO (USING MODAL CLUSTERS)
# =============================================================================

cluster_modal_site_solute <- cluster_modal %>%
  dplyr::mutate(
    solute_type = dplyr::case_when(
      chemical %in% geo_solutes ~ "Geogenic",
      chemical %in% bio_solutes ~ "Biogenic",
      TRUE                      ~ "Other"
    ),
    solute_type = factor(solute_type,
                         levels = c("Geogenic","Biogenic","Other"))
  )

# Majority cluster per STREAM (across all solutes)
majority_cluster_byStream <- cluster_modal_site_solute %>%
  dplyr::group_by(Stream_Name, Cluster_mode) %>%
  dplyr::summarise(n_solutes = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(Stream_Name) %>%
  dplyr::slice_max(n_solutes, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Stream_Name)


# Majority cluster per STREAM, split by solute type (geo vs bio)
majority_cluster_byStream_soluteType <- cluster_modal_site_solute %>%
  dplyr::filter(solute_type != "Other") %>%
  dplyr::group_by(Stream_Name, solute_type, Cluster_mode) %>%
  dplyr::summarise(n_solutes = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(Stream_Name, solute_type) %>%
  dplyr::slice_max(n_solutes, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(solute_type, Stream_Name)

# Full distribution of modal clusters per stream (for stacked bars)
site_cluster_distribution <- cluster_modal_site_solute %>%
  dplyr::group_by(Stream_Name, Cluster_mode) %>%
  dplyr::summarise(n_solutes = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(Stream_Name) %>%
  dplyr::mutate(
    total_solutes = sum(n_solutes),
    prop_solutes  = n_solutes / total_solutes
  ) %>%
  dplyr::ungroup()

# =============================================================================
# 4. WRITE CALENDAR-YEAR CSV OUTPUTS (USING FIXED CLIMATOLOGY-REFERENCE CLUSTERS)
# =============================================================================
# NOTE: Output file is now named "byCalendarYear" to be explicit about what it contains.
# The file "byWaterYear" is kept for backward compatibility but is actually CY data.

# CY table with normalized months, cluster IDs, distance, and modal cluster per site–solute
readr::write_csv(
  month_clusters_cy,
  file.path(output_dir, "ClusterStreams_allSolutes_byCalendarYear.csv")
)

# ALSO write with old name for backward compatibility (scripts reference byWaterYear)
# TODO: Update all downstream scripts to use byCalendarYear, then remove this
readr::write_csv(
  month_clusters_cy %>%
    dplyr::rename(water_year = calendar_year),  # Rename for compatibility
  file.path(output_dir, "ClusterStreams_allSolutes_byWaterYear.csv")
)


# Modal cluster per stream–solute
readr::write_csv(
  cluster_modal,
  file.path(output_dir, "ClusterStreams_allSolutes_modalClusters.csv")
)

# site–solute modal with solute_type tag
readr::write_csv(
  cluster_modal_site_solute,
  file.path(output_dir, "ClusterStreams_modal_byStreamSolute.csv")
)

# majority modal cluster per stream (all solutes)
readr::write_csv(
  majority_cluster_byStream,
  file.path(output_dir, "ClusterStreams_majorityCluster_byStream.csv")
)

# majority modal cluster per stream, split by solute type (geo vs bio)
readr::write_csv(
  majority_cluster_byStream_soluteType,
  file.path(output_dir, "ClusterStreams_majorityCluster_byStream_geoBio.csv")
)

# full modal cluster distribution per stream
readr::write_csv(
  site_cluster_distribution,
  file.path(output_dir, "ClusterStreams_clusterDistribution_byStream.csv")
)

# Restore working directory
setwd(original_wd)

# =============================================================================
# END 
# =============================================================================
