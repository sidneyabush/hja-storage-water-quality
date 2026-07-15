# =============================================================================
# STEP 2h: CLUSTER-MOBILIZATION INTEGRATION
# =============================================================================
# Purpose: Integrate DTW clusters (concentration synchrony) with Wymore
#          synchrony (mobilization synchrony) to test whether:
#
#   1. Solutes in the same DTW cluster have higher Wymore synchrony
#   2. Storage divergence modulates this relationship
#   3. Concentration synchrony (clusters) predicts mobilization synchrony (Wymore)
#
# Key Distinction:
#   - DTW clusters: Long-term seasonal concentration patterns (concentration synchrony)
#   - Wymore synchrony: Real-time CQ slope correlation (mobilization synchrony)
#
# Hypothesis:
#   Sites in the same cluster should have higher baseline Wymore synchrony,
#   BUT storage divergence can still produce asynchronous mobilization even
#   within clusters.
#
# Research Question Addressed:
#   RQ3: Do solute export regimes signal synchronized vs. independent behavior?
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
})

rm(list = ls())

# Source helpers
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
try(source(file.path(repo_dir, "00_helpers", "plot_prefs.R")), silent = TRUE)
try(source(file.path(repo_dir, "00_helpers", "workflow_config.R")), silent = TRUE)
try(source(file.path(repo_dir, "00_helpers", "plot_theme_set.R")), silent = TRUE)

# Paths
paths <- get_project_paths()
out_dir <- paths$out_dir
fig_dir <- file.path(paths$fig_root, "02_exploration", "2h_cluster_mobilization")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

message("\n╔════════════════════════════════════════════════════════════════╗\n")
message("║  CLUSTER-MOBILIZATION INTEGRATION                             ║\n")
message("╚════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# LOAD DATA
# =============================================================================
message("=== LOADING DATA ===\n")

# DTW Clusters (concentration synchrony)
clusters_modal <- read_csv(file.path(out_dir, "ClusterStreams_allSolutes_modalClusters.csv"),
                           show_col_types = FALSE) %>%
  mutate(
    Stream_Name = factor(Stream_Name, levels = site_order),
    solute = factor(chemical, levels = solute_order),
    Cluster_mode = factor(Cluster_mode, levels = cluster_levels)
  ) %>%
  select(Stream_Name, solute, Cluster_mode) %>%
  rename(Cluster = Cluster_mode)

# Wymore cross-site synchrony (mobilization synchrony)
wymore_crosssite <- read_csv(file.path(out_dir, "HJA_wymore_crosssite_sync.csv"),
                             show_col_types = FALSE) %>%
  mutate(
    Stream1 = factor(Stream1, levels = site_order),
    Stream2 = factor(Stream2, levels = site_order),
    solute = factor(solute, levels = solute_order)
  )

# Abbott synchrony (for comparison)
abbott_sync <- read_csv(file.path(out_dir, "HJA_Abbott_synchrony_windows.csv"),
                        show_col_types = FALSE) %>%
  filter(time_scale == "annual") %>%
  mutate(
    Stream1 = factor(Stream1, levels = site_order),
    Stream2 = factor(Stream2, levels = site_order),
    solute = factor(solute, levels = solute_order)
  )

# Composite synchrony (long-term)
composite_sync <- read_csv(file.path(out_dir, "HJA_composite_synchrony.csv"),
                           show_col_types = FALSE) %>%
  mutate(
    Stream_Name = factor(Stream_Name, levels = site_order),
    solute = factor(solute, levels = solute_order)
  )

# Storage data (for modulation analysis)
site_means <- read_csv(file.path(out_dir, "HJA_master_site_means.csv"),
                       show_col_types = FALSE) %>%
  mutate(
    Stream_Name = factor(Stream_Name, levels = site_order),
    solute = factor(solute, levels = solute_order)
  )

message("  Clusters:", nrow(clusters_modal), "site-solute combinations\n")
message("  Wymore crosssite:", nrow(wymore_crosssite), "site-pair-window observations\n")
message("  Abbott synchrony:", nrow(abbott_sync), "site-pair-year observations\n")

# Identify storage metric
storage_metric <- intersect(PRIMARY_STORAGE_METRIC, names(site_means))
if (length(storage_metric) == 0) {
  storage_metric <- intersect(c("Q_dS_range_mm", "WB_dS_range_mm"), names(site_means))
}
storage_metric <- storage_metric[[1]]
message("  Using storage metric:", storage_metric, "\n")

# =============================================================================
# ANALYSIS 1: SAME-CLUSTER PAIRS HAVE HIGHER WYMORE SYNCHRONY?
# =============================================================================
message("\n=== ANALYSIS 1: Within-Cluster vs. Between-Cluster Wymore Synchrony ===\n")

# Add cluster assignments to site pairs
wymore_with_clusters <- wymore_crosssite %>%
  left_join(clusters_modal %>% rename(Cluster1 = Cluster),
            by = c("Stream1" = "Stream_Name", "solute")) %>%
  left_join(clusters_modal %>% rename(Cluster2 = Cluster),
            by = c("Stream2" = "Stream_Name", "solute")) %>%
  mutate(
    same_cluster = !is.na(Cluster1) & !is.na(Cluster2) & Cluster1 == Cluster2,
    cluster_pair = if_else(same_cluster,
                           paste0("Within-", Cluster1),
                           paste0(Cluster1, "-", Cluster2))
  )

# Summary statistics
cluster_sync_summary <- wymore_with_clusters %>%
  filter(!is.na(sync)) %>%
  group_by(same_cluster) %>%
  summarise(
    n = n(),
    mean_sync = mean(sync, na.rm = TRUE),
    sd_sync = sd(sync, na.rm = TRUE),
    se_sync = sd_sync / sqrt(n),
    .groups = "drop"
  )

message("\nWymore synchrony by cluster pairing:\n")
print(cluster_sync_summary)

# Statistical test
if (nrow(cluster_sync_summary) == 2) {
  within_sync <- wymore_with_clusters %>% filter(same_cluster == TRUE, !is.na(sync)) %>% mutate(sync_num = as.numeric(sync)) %>% pull(sync_num)
  between_sync <- wymore_with_clusters %>% filter(same_cluster == FALSE, !is.na(sync)) %>% mutate(sync_num = as.numeric(sync)) %>% pull(sync_num)

  if (length(within_sync) > 0 && length(between_sync) > 0) {
    test_result <- wilcox.test(within_sync, between_sync)
    message("\nWilcoxon test (within vs. between cluster):\n")
    message("  W =", test_result$statistic, ", p =", format.pval(test_result$p.value, digits = 4), "\n")

    if (test_result$p.value < 0.05) {
      message("  → Significant difference! Same-cluster pairs are",
          if_else(mean(within_sync, na.rm = TRUE) > mean(between_sync, na.rm = TRUE),
                  "MORE", "LESS"),
          "synchronous.\n")
    } else {
      message("  → No significant difference.\n")
    }
  }
}

# Visualization 1: Boxplot
p1 <- wymore_with_clusters %>%
  filter(!is.na(sync), !is.na(same_cluster)) %>%
  mutate(cluster_type = if_else(same_cluster, "Within-cluster", "Between-cluster")) %>%
  ggplot(aes(x = cluster_type, y = sync, fill = cluster_type)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
  geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4, fill = "red") +
  scale_fill_manual(values = c("Within-cluster" = cluster_colors[["1"]],
                                "Between-cluster" = "grey60")) +
  labs(
    x = "Cluster Pairing Type",
    y = "Wymore Synchrony (proportion synchronous)",

ggsave(file.path(fig_dir, "01_within_vs_between_cluster_wymore.png"), p1,
       width = 10, height = 7, dpi = 300)

# Visualization 2: By-cluster breakdown
p2 <- wymore_with_clusters %>%
  filter(!is.na(sync), !is.na(Cluster1), !is.na(Cluster2)) %>%
  ggplot(aes(x = Cluster1, y = sync, fill = same_cluster)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
  facet_wrap(~Cluster2, labeller = labeller(Cluster2 = \(x) paste("Pair with Cluster", x))) +
  scale_fill_manual(
    values = c("TRUE" = cluster_colors[["1"]], "FALSE" = "grey60"),
    name = "Same Cluster",
    labels = c("FALSE" = "No", "TRUE" = "Yes")
  ) +
  labs(
    x = "Cluster 1",
    y = "Wymore Synchrony",

ggsave(file.path(fig_dir, "01b_wymore_by_cluster_pairing.png"), p2,
       width = 12, height = 9, dpi = 300)

message("\n✓ Analysis 1 complete: Within-cluster synchrony\n")

# =============================================================================
# ANALYSIS 2: STORAGE MODULATES CLUSTER-WYMORE RELATIONSHIP?
# =============================================================================
message("\n=== ANALYSIS 2: Storage Divergence Modulates Cluster-Synchrony Link ===\n")

# Compute storage divergence for each site pair
storage_for_pairs <- site_means %>%
  select(Stream_Name, solute, all_of(storage_metric))

wymore_with_storage <- wymore_with_clusters %>%
  left_join(storage_for_pairs %>% rename(storage1 = all_of(storage_metric)),
            by = c("Stream1" = "Stream_Name", "solute")) %>%
  left_join(storage_for_pairs %>% rename(storage2 = all_of(storage_metric)),
            by = c("Stream2" = "Stream_Name", "solute")) %>%
  mutate(
    storage_diff = abs(storage1 - storage2),
    storage_mean = (storage1 + storage2) / 2
  )

# Aggregate by window to get mean synchrony per window
wymore_window_agg <- wymore_with_storage %>%
  group_by(solute, window_center, water_year, hydrologic_season, same_cluster) %>%
  summarise(
    mean_sync = mean(sync, na.rm = TRUE),
    mean_storage_diff = mean(storage_diff, na.rm = TRUE),
    n_pairs = n(),
    .groups = "drop"
  ) %>%
  filter(n_pairs >= 3, is.finite(mean_sync), is.finite(mean_storage_diff))

message("  Aggregated to", nrow(wymore_window_agg), "windows\n")

# Linear model: sync ~ storage_diff * same_cluster
if (nrow(wymore_window_agg) >= 20) {

  lm_interaction <- lm(mean_sync ~ mean_storage_diff * same_cluster,
                       data = wymore_window_agg)

  message("\nLinear model: Wymore sync ~ storage_diff * same_cluster\n")
  print(summary(lm_interaction))

  # Interaction test
  interaction_p <- anova(lm_interaction)["mean_storage_diff:same_clusterTRUE", "Pr(>F)"]

  if (!is.na(interaction_p) && interaction_p < 0.05) {
    message("\n  → Significant interaction (p =", round(interaction_p, 4), ")\n")
    message("     Storage divergence affects within-cluster and between-cluster pairs differently!\n")
  } else {
    message("\n  → No significant interaction\n")
  }

  # Visualization: Storage divergence vs synchrony, colored by cluster pairing
  p3 <- ggplot(wymore_window_agg, aes(x = mean_storage_diff, y = mean_sync,
                                      color = same_cluster)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1.2) +
    scale_color_manual(
      values = c("TRUE" = cluster_colors[["1"]], "FALSE" = "grey60"),
      name = "Same Cluster",
      labels = c("FALSE" = "Between-cluster", "TRUE" = "Within-cluster")
    ) +
    labs(
      x = paste0("Mean |Δ ", get_storage_label(storage_metric, short = TRUE), "| between site pairs"),
      y = "Mean Wymore Synchrony (per window)",
                      "Ribbons = 95% CI for linear fit")
    ) +
    theme_hja() +
    legend_bottom()

  ggsave(file.path(fig_dir, "02_storage_modulation_cluster_sync.png"), p3,
         width = 12, height = 8, dpi = 300)
}

message("\n✓ Analysis 2 complete: Storage modulation\n")

# =============================================================================
# ANALYSIS 3: CLUSTER MEMBERSHIP PREDICTS LONG-TERM WYMORE SYNCHRONY?
# =============================================================================
message("\n=== ANALYSIS 3: Long-term Cluster Prediction of Wymore Synchrony ===\n")

# Aggregate Wymore synchrony to site-solute level (long-term)
wymore_longterm <- wymore_with_clusters %>%
  group_by(Stream1, Stream2, solute, same_cluster, Cluster1, Cluster2) %>%
  summarise(
    mean_sync = mean(sync, na.rm = TRUE),
    n_windows = n(),
    .groups = "drop"
  ) %>%
  filter(n_windows >= 10)  # At least 10 windows

message("  Long-term pairs:", nrow(wymore_longterm), "\n")

# Join with Abbott synchrony for comparison
if ("cqslope_sync_allpairs" %in% names(composite_sync)) {

  # Reshape composite sync to pairwise format (approximation)
  abbott_by_cluster <- composite_sync %>%
    left_join(clusters_modal, by = c("Stream_Name", "solute")) %>%
    select(Stream_Name, solute, Cluster, cqslope_sync_allpairs)

  # Compare within-cluster vs between-cluster Abbott synchrony
  cluster_abbott_comparison <- abbott_by_cluster %>%
    group_by(Cluster, solute) %>%
    summarise(
      mean_abbott_sync = mean(cqslope_sync_allpairs, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )

  message("\nAbbott CQ-slope synchrony by cluster:\n")
  print(cluster_abbott_comparison)
}

# Visualization 3: Cluster × Wymore synchrony matrix
cluster_wymore_matrix <- wymore_longterm %>%
  group_by(Cluster1, Cluster2) %>%
  summarise(
    mean_sync = mean(mean_sync, na.rm = TRUE),
    n_pairs = n(),
    .groups = "drop"
  )

p4 <- ggplot(cluster_wymore_matrix, aes(x = Cluster1, y = Cluster2, fill = mean_sync)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = round(mean_sync, 2)), color = "white", size = 5, fontface = "bold") +
  scale_fill_viridis_c(option = "plasma", name = "Mean\nWymore\nSync") +
  labs(
    x = "Cluster 1",
    y = "Cluster 2",

ggsave(file.path(fig_dir, "03_cluster_wymore_matrix.png"), p4,
       width = 10, height = 9, dpi = 300)

message("\n✓ Analysis 3 complete: Long-term cluster prediction\n")

# =============================================================================
# ANALYSIS 4: CONCENTRATION SYNCHRONY (Abbott) vs MOBILIZATION SYNCHRONY (Wymore)
# =============================================================================
message("\n=== ANALYSIS 4: Abbott vs. Wymore Synchrony ===\n")

# Compare Abbott concentration sync with Wymore CQ sync
if ("conc_sync_allpairs" %in% names(composite_sync) &&
    "wymore_crosssite_allpairs" %in% names(composite_sync)) {

  sync_comparison <- composite_sync %>%
    left_join(clusters_modal, by = c("Stream_Name", "solute")) %>%
    filter(!is.na(conc_sync_allpairs), !is.na(wymore_crosssite_allpairs)) %>%
    mutate(
      sync_diff = conc_sync_allpairs - wymore_crosssite_allpairs,
      sync_agreement = case_when(
        conc_sync_allpairs > 0.5 & wymore_crosssite_allpairs > 0.5 ~ "Both High",
        conc_sync_allpairs < 0.5 & wymore_crosssite_allpairs < 0.5 ~ "Both Low",
        TRUE ~ "Discordant"
      )
    )

  message("\nSynchrony agreement:\n")
  print(table(sync_comparison$sync_agreement, sync_comparison$Cluster))

  # Correlation
  cor_abbott_wymore <- cor.test(sync_comparison$conc_sync_allpairs,
                                 sync_comparison$wymore_crosssite_allpairs)
  message("\nCorrelation (Abbott conc sync vs Wymore crosssite):\n")
  message("  r =", round(cor_abbott_wymore$estimate, 3),
      ", p =", format.pval(cor_abbott_wymore$p.value, digits = 4), "\n")

  # Visualization
  p5 <- ggplot(sync_comparison, aes(x = conc_sync_allpairs, y = wymore_crosssite_allpairs,
                                    color = Cluster)) +
    geom_point(alpha = 0.7, size = 3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
    scale_color_cluster() +
    labs(
      x = "Abbott Concentration Synchrony (all-pairs)",
      y = "Wymore Cross-site Synchrony (all-pairs)",
                       ", p = ", format.pval(cor_abbott_wymore$p.value, digits = 3)),

  ggsave(file.path(fig_dir, "04_abbott_vs_wymore_synchrony.png"), p5,
         width = 11, height = 9, dpi = 300)

  message("\n✓ Analysis 4 complete: Abbott vs Wymore comparison\n")
}

# =============================================================================
# SAVE RESULTS
# =============================================================================
message("\n=== SAVING RESULTS ===\n")

# Save within-cluster vs between-cluster summary
write_csv(cluster_sync_summary,
          file.path(out_dir, "02_exploration/cluster_wymore_summary.csv"))

# Save cluster-wymore matrix
write_csv(cluster_wymore_matrix,
          file.path(out_dir, "02_exploration/cluster_wymore_matrix.csv"))

# Save synchrony comparison
if (exists("sync_comparison")) {
  write_csv(sync_comparison,
            file.path(out_dir, "02_exploration/abbott_wymore_comparison.csv"))
}

# =============================================================================
# SUMMARY
# =============================================================================
message("\n╔════════════════════════════════════════════════════════════════╗\n")
message("║  ANALYSIS COMPLETE                                            ║\n")
message("╚════════════════════════════════════════════════════════════════╝\n\n")

message("KEY FINDINGS:\n\n")

message("1. Within-cluster pairs:\n")
if (exists("cluster_sync_summary") && nrow(cluster_sync_summary) > 0) {
  within_mean <- cluster_sync_summary %>% filter(same_cluster == TRUE) %>% pull(mean_sync)
  between_mean <- cluster_sync_summary %>% filter(same_cluster == FALSE) %>% pull(mean_sync)
  if (length(within_mean) > 0 && length(between_mean) > 0) {
    message("   Within-cluster mean sync:", round(within_mean, 3), "\n")
    message("   Between-cluster mean sync:", round(between_mean, 3), "\n")
    message("   Difference:", round(within_mean - between_mean, 3), "\n")
  }
}

message("\n2. Storage modulation:\n")
if (exists("interaction_p") && !is.na(interaction_p)) {
  if (interaction_p < 0.05) {
    message("   Significant interaction (p =", round(interaction_p, 4), ")\n")
    message("   → Storage divergence affects cluster pairs differently\n")
  } else {
    message("   No significant interaction\n")
  }
}

message("\n3. Abbott vs Wymore:\n")
if (exists("cor_abbott_wymore")) {
  message("   Correlation: r =", round(cor_abbott_wymore$estimate, 3),
      ", p =", format.pval(cor_abbott_wymore$p.value, digits = 4), "\n")
  if (cor_abbott_wymore$p.value < 0.05) {
    if (cor_abbott_wymore$estimate > 0) {
      message("   → Positive correlation: Sites with high concentration synchrony also show high mobilization synchrony\n")
    } else {
      message("   → Negative correlation: Concentration and mobilization synchrony are decoupled\n")
    }
  } else {
    message("   → No significant correlation: Concentration and mobilization synchrony are independent\n")
  }
}

message("\nOutputs saved to:\n")
message("  Figures:", fig_dir, "\n")
message("  Tables:", file.path(out_dir, "02_exploration/"), "\n\n")
