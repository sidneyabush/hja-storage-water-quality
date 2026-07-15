# =============================================================================
# 1c CLUSTER TESTS: TIMING vs MECHANISM (non-plotting)
# =============================================================================
# Purpose: Test whether DTW clusters capture timing of flowpath activation
# or geochemical mechanism (deep vs shallow flowpaths)
#
# Approach:
# 1. Z-normalized DTW (emphasizes timing/shape)
# 2. Raw-scaled DTW (retains magnitude/composition)
# 3. Feature-based clustering (derives mean, CV, IQR, skewness, kurtosis)
# 4. Compute Adjusted Rand Index (ARI) for stability
# 5. Indicator species analysis (which solutes define each cluster?)
# 6. Integrate with hydrologic state and synchrony metrics
#
# Outputs:
#   - outputs/clusters/ari_by_year.csv (stability: ARI between global/per-year)
#   - outputs/clusters/indicator_global.csv (which solutes define clusters)
#   - outputs/clusters/cluster_site_year.csv (cluster membership by site/year)
#   - outputs/clusters/clustering_summary.txt (summary statistics)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(dtwclust)
  library(dtw)
  library(mclust)      # for adjustedRandIndex
  library(indicspecies)
  library(moments)     # for skewness/kurtosis
})

rm(list = ls())

# Paths (use repo-local outputs from prep)
repo_root <- getwd()
prep_dir <- file.path(repo_root, "outputs", "prep")
output_dir <- file.path(repo_root, "outputs")
cluster_dir <- file.path(output_dir, "clusters")
dir.create(prep_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cluster_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(42)

message("=== CLUSTER ANALYSIS: TESTING TIMING VS MECHANISM ===\n\n")

# =============================================================================
# 1. LOAD CQ DAILY DATA & HYDRO SEASONS
# =============================================================================
message("Step 1: Loading CQ timeseries...\n")

cq_path <- file.path(prep_dir, "HJA_CQ_daily_for_analysis.csv")
if (!file.exists(cq_path)) {
  stop(paste0("Missing prep output: ", cq_path,
              "\nRun 01_data_prep steps (1a–1d) to generate it."))
}

cq_daily <- readr::read_csv(
  cq_path,
  show_col_types = FALSE
) %>%
  mutate(date = as.Date(date)) %>%
  arrange(Stream_Name, date)

seasons_path <- file.path(prep_dir, "HJA_daily_Q_with_seasons.csv")
if (!file.exists(seasons_path)) {
  stop(paste0("Missing prep output: ", seasons_path,
              "\nRun 01_data_prep step 1b to generate it."))
}

seasons <- readr::read_csv(
  seasons_path,
  show_col_types = FALSE
) %>%
  mutate(date = as.Date(date)) %>%
  select(date, hydrologic_season)

cq_daily <- cq_daily %>%
  left_join(seasons, by = "date")

solutes <- cq_daily %>%
  select(-c(Stream_Name, date, Qcms, water_year, month, hydrologic_season)) %>%
  names()

sites <- unique(cq_daily$Stream_Name)
years <- sort(unique(cq_daily$water_year))

message("  Sites:", paste(sites, collapse = ", "), "\n")
message("  Solutes:", paste(solutes, collapse = ", "), "\n")
message("  Years:", min(years), "to", max(years), "\n")
message("  Records:", nrow(cq_daily), "\n\n")

# =============================================================================
# 2. BUILD TIMESERIES & FEATURE MATRICES FOR CLUSTERING
# =============================================================================
message("Step 2: Preparing timeseries for clustering...\n")

ts_global_znorm <- list()
ts_global_raw <- list()
features_global <- tibble()

for (sol in solutes) {
  for (site in sites) {
    
    ts <- cq_daily %>%
      filter(Stream_Name == site, !is.na(.data[[sol]])) %>%
      arrange(date) %>%
      pull(!!sym(sol))
    
    if (length(ts) > 50) {
      name <- paste0(site, "_", sol)
      
      # Z-normalized (for DTW emphasizing timing)
      ts_znorm <- as.numeric(scale(ts))
      ts_global_znorm[[name]] <- ts_znorm
      
      # Raw (for DTW retaining magnitude)
      ts_global_raw[[name]] <- ts
      
      # Features (for amplitude/composition-based clustering)
      features_global <- features_global %>%
        bind_rows(tibble(
          name = name,
          site = site,
          solute = sol,
          mean_conc = mean(ts, na.rm = TRUE),
          sd_conc = sd(ts, na.rm = TRUE),
          cv_conc = sd(ts, na.rm = TRUE) / pmax(mean(ts, na.rm = TRUE), 0.001),
          median_conc = median(ts, na.rm = TRUE),
          p25 = quantile(ts, 0.25, na.rm = TRUE),
          p75 = quantile(ts, 0.75, na.rm = TRUE),
          iqr = IQR(ts, na.rm = TRUE),
          min_conc = min(ts, na.rm = TRUE),
          max_conc = max(ts, na.rm = TRUE),
          range_conc = max(ts, na.rm = TRUE) - min(ts, na.rm = TRUE)
        ))
    }
  }
}

message("  Created", length(ts_global_znorm), "global timeseries\n")
message("  Feature matrix:", nrow(features_global), "rows\n\n")

# =============================================================================
# 3. GLOBAL DTW CLUSTERING: Z-NORMALIZED (TIMING-SENSITIVE)
# =============================================================================
message("Step 3a: Z-normalized DTW clustering (emphasizes timing)...\n")

ts_znorm_mat <- do.call(rbind, ts_global_znorm)

# Choose k by silhouette
silh_znorm <- numeric(5)
for (k_test in 2:6) {
  clust_test <- tsclust(ts_znorm_mat, type = "partitional",
                        k = k_test,
                        distance = "dtw_basic",
                        centroid = "partition",
                        seed = 42,
                        trace = FALSE)
  # Silhouette score
  sil <- cluster::silhouette(clust_test@cluster, 
                             proxy::dist(ts_znorm_mat, method = "dtw_basic"))
  silh_znorm[k_test - 1] <- mean(sil[, "sil_width"])
}

k_opt_znorm <- which.max(silh_znorm) + 1
message("  Silhouette (k=2:6):", round(silh_znorm, 3), "\n")
message("  Optimal k:", k_opt_znorm, "\n")

clust_znorm <- tsclust(ts_znorm_mat, type = "partitional",
                       k = k_opt_znorm,
                       distance = "dtw_basic",
                       centroid = "partition",
                       seed = 42,
                       trace = FALSE)
labels_znorm <- setNames(clust_znorm@cluster, rownames(ts_znorm_mat))

message("  Cluster sizes:", table(labels_znorm), "\n\n")

# =============================================================================
# 3b. GLOBAL DTW CLUSTERING: RAW-SCALED (MAGNITUDE-SENSITIVE)
# =============================================================================
message("Step 3b: Raw-scaled DTW clustering (retains magnitude/mechanism)...\n")

ts_raw_mat <- do.call(rbind, ts_global_raw)
# Min-max scale to [0, 1] for comparability across solutes
ts_raw_scaled <- t(apply(ts_raw_mat, 1, function(x) {
  (x - min(x, na.rm = TRUE)) / pmax(max(x, na.rm = TRUE) - min(x, na.rm = TRUE), 0.001)
}))

silh_raw <- numeric(5)
for (k_test in 2:6) {
  clust_test <- tsclust(ts_raw_scaled, type = "partitional",
                        k = k_test,
                        distance = "dtw_basic",
                        centroid = "partition",
                        seed = 42,
                        trace = FALSE)
  sil <- cluster::silhouette(clust_test@cluster,
                             proxy::dist(ts_raw_scaled, method = "dtw_basic"))
  silh_raw[k_test - 1] <- mean(sil[, "sil_width"])
}

k_opt_raw <- which.max(silh_raw) + 1
message("  Optimal k:", k_opt_raw, "\n")

clust_raw <- tsclust(ts_raw_scaled, type = "partitional",
                     k = k_opt_raw,
                     distance = "dtw_basic",
                     centroid = "partition",
                     seed = 42,
                     trace = FALSE)
labels_raw <- setNames(clust_raw@cluster, rownames(ts_raw_mat))

message("  Cluster sizes:", table(labels_raw), "\n\n")

# =============================================================================
# 3c. FEATURE-BASED CLUSTERING (EUCLIDEAN ON DERIVED FEATURES)
# =============================================================================
message("Step 3c: Feature-based clustering (composition/amplitude)...\n")

feat_mat <- features_global %>%
  select(mean_conc, cv_conc, median_conc, iqr, range_conc) %>%
  as.matrix()
rownames(feat_mat) <- features_global$name

feat_scaled <- scale(feat_mat)
hc_feat <- hclust(dist(feat_scaled), method = "ward.D2")
labels_feat <- setNames(cutree(hc_feat, k = k_opt_znorm), 
                        rownames(feat_mat))

message("  Cluster sizes:", table(labels_feat), "\n")
message("  (Using k =", k_opt_znorm, "to match z-norm DTW)\n\n")

# =============================================================================
# 4. PER-YEAR CLUSTERING & ARI (STABILITY METRIC)
# =============================================================================
message("Step 4: Computing per-year clustering and ARI...\n")

ari_results <- tibble()
labels_per_year <- list()

for (yr in years) {
  
  ts_year_znorm <- list()
  
  for (sol in solutes) {
    for (site in sites) {
      
      ts <- cq_daily %>%
        filter(Stream_Name == site, water_year == yr, !is.na(.data[[sol]])) %>%
        arrange(date) %>%
        pull(!!sym(sol))
      
      if (length(ts) > 30) {
        name <- paste0(site, "_", sol)
        ts_year_znorm[[name]] <- as.numeric(scale(ts))
      }
    }
  }
  
  if (length(ts_year_znorm) > 1) {
    
    ts_year_mat <- do.call(rbind, ts_year_znorm)
    
    clust_year <- tsclust(ts_year_mat, type = "partitional",
                          k = k_opt_znorm,
                          distance = "dtw_basic",
                          centroid = "partition",
                          seed = 42,
                          trace = FALSE)
    
    labels_year <- setNames(clust_year@cluster, rownames(ts_year_mat))
    labels_per_year[[as.character(yr)]] <- labels_year
    
    # Compute ARI: global vs per-year for matching series
    matching <- intersect(names(labels_znorm), names(labels_year))
    
    if (length(matching) > 2) {
      ari <- mclust::adjustedRandIndex(labels_znorm[matching], 
                                       labels_year[matching])
    } else {
      ari <- NA_real_
    }
    
    ari_results <- ari_results %>%
      bind_rows(tibble(
        water_year = yr,
        n_series = length(ts_year_znorm),
        n_matching = length(matching),
        ari = ari
      ))
  }
}

write_csv(ari_results, file.path(cluster_dir, "ari_by_year.csv"))

message("  ARI summary:\n")
message("    Mean ARI:", round(mean(ari_results$ari, na.rm = TRUE), 3), "\n")
message("    SD ARI:", round(sd(ari_results$ari, na.rm = TRUE), 3), "\n")
message("    Range:", round(min(ari_results$ari, na.rm = TRUE), 3), "to",
    round(max(ari_results$ari, na.rm = TRUE), 3), "\n\n")

# =============================================================================
# 5. INDICATOR SPECIES ANALYSIS
# =============================================================================
message("Step 5: Identifying indicator solutes for each cluster...\n")

# Build solute x cluster contingency table
solute_clust_tab <- tibble(
  name = names(labels_znorm),
  cluster = labels_znorm
) %>%
  separate(name, into = c("site", "solute"), sep = "_", extra = "merge") %>%
  group_by(solute, cluster) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = cluster, values_from = count, values_fill = 0) %>%
  column_to_rownames("solute") %>%
  as.matrix()

# Indicator analysis: which solutes significantly associated with each cluster?
ind_res <- indicspecies::multipatt(solute_clust_tab,
                                   rep(seq_len(ncol(solute_clust_tab)), 
                                       nrow(solute_clust_tab)),
                                   control = how(nperm = 999),
                                   func = "r")

# Extract summary
ind_summary <- data.frame(
  solute = rownames(solute_clust_tab),
  cluster = apply(ind_res$sign, 1, function(x) which(x == 1)),
  stat = ind_res$A[,1],
  p_value = ind_res$B
) %>%
  as_tibble() %>%
  arrange(p_value)

write_csv(ind_summary, file.path(cluster_dir, "indicator_solutes.csv"))

message("  Top indicator solutes:\n")
print(head(ind_summary, 8))
message("\n")

# =============================================================================
# 6. SUMMARIZE CLUSTER MEMBERSHIP BY SITE & YEAR
# =============================================================================
message("Step 6: Computing modal cluster by site and water year...\n")

cluster_site_year <- tibble()

for (site in sites) {
  for (yr in years) {
    
    site_yr_series <- features_global %>%
      filter(site == !!site) %>%
      pull(name)
    
    if (length(site_yr_series) > 0) {
      
      clust_znorm_modal <- as.numeric(names(sort(table(labels_znorm[site_yr_series]), 
                                                 decreasing = TRUE)[1])[1])
      clust_raw_modal <- as.numeric(names(sort(table(labels_raw[site_yr_series]), 
                                               decreasing = TRUE)[1])[1])
      clust_feat_modal <- as.numeric(names(sort(table(labels_feat[site_yr_series]), 
                                                decreasing = TRUE)[1])[1])
      
      cluster_site_year <- cluster_site_year %>%
        bind_rows(tibble(
          site = site,
          water_year = yr,
          cluster_znorm_global = clust_znorm_modal,
          cluster_raw_global = clust_raw_modal,
          cluster_feat_global = clust_feat_modal,
          n_solutes = length(site_yr_series)
        ))
    }
  }
}

write_csv(cluster_site_year, file.path(cluster_dir, "cluster_site_year.csv"))

message("  Saved cluster assignments for", nrow(cluster_site_year), "site-year combinations\n\n")

