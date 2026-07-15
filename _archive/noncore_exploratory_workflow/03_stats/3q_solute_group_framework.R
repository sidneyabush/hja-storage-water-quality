# =============================================================================
# 3q_solute_group_framework.R
# Comprehensive analysis framework stratified by solute groups
# 
# Key Questions:
# 1. Does cluster type predict synchrony?
# 2. Do catchment properties → sync differ by solute group?
# 3. Do catchment properties → cluster differ by solute group?
# =============================================================================

library(tidyverse)
library(lme4)
library(MuMIn)
library(nnet)
library(randomForest)
library(car)

# Paths
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
box_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality"
out_dir <- file.path(box_dir, "outputs")
fig_dir <- file.path(box_dir, "figures")  # Main figures folder
meeting_fig_dir <- file.path(box_dir, "updates_12042025")  # Meeting-specific copy

# Create directories if needed
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(meeting_fig_dir, showWarnings = FALSE, recursive = TRUE)

source(file.path(repo_dir, "00_helpers", "workflow_config.R"))
source(file.path(repo_dir, "00_helpers", "plot_prefs.R"))
try(source(file.path(repo_dir, "00_helpers", "plot_theme_set.R")), silent = TRUE)

# =============================================================================
# DEFINE SOLUTE GROUPS
# =============================================================================

solute_groups <- tribble(
  ~solute, ~group,
  "Ca", "Geogenic",
  "Mg", "Geogenic", 
  "Na", "Geogenic",
  "K", "Geogenic",
  "DSi", "Geogenic",
  "Cl", "Geogenic",
  "SO4", "Geogenic",
  "NO3", "Nutrients",
  "NH3", "Nutrients",
  "PO4", "Nutrients",
  "DOC", "Organic"
)

message("=============================================================================\n")
message("SOLUTE GROUP FRAMEWORK ANALYSIS\n")
message("=============================================================================\n\n")

message("Solute Groups:\n")
solute_groups %>% 
  group_by(group) %>% 
  summarise(solutes = paste(solute, collapse = ", "), n = n()) %>%
  print()

# =============================================================================
# LOAD DATA
# =============================================================================

# Synchrony data - this is pairwise with Stream1/Stream2
sync_raw <- read_csv(file.path(out_dir, "HJA_Abbott_synchrony_windows.csv"), show_col_types = FALSE)

# Modal clusters - use labels from plot_prefs.R
modal_clusters <- read_csv(file.path(out_dir, "ClusterStreams_allSolutes_modalClusters.csv"), show_col_types = FALSE) %>%
  rename(solute = chemical, modal_cluster = Cluster_mode) %>%
  mutate(modal_cluster = cluster_labels[as.character(modal_cluster)])

# Site characteristics
site_chars <- read_csv(file.path(out_dir, "Catchment_site_characteristics.csv"), show_col_types = FALSE)

# Add solute groups to modal clusters
modal_clusters <- modal_clusters %>%
  left_join(solute_groups, by = "solute")

message("\nData loaded:\n")
message("  Synchrony records:", nrow(sync_raw), "\n")
message("  Modal cluster records:", nrow(modal_clusters), "\n")
message("  Sites with characteristics:", nrow(site_chars), "\n")

# =============================================================================
# PREPARE SYNC DATA
# =============================================================================

# Filter to outlet pairs only for outlet sync analysis
# Note: Stream1 is GSLOOK (outlet), Stream2 is the tributary
outlet_sync <- sync_raw %>%
  filter(is_outlet_pair == TRUE) %>%
  left_join(solute_groups, by = "solute") %>%
  rename(Stream_Name = Stream2)  # The NON-outlet stream (tributary)

# Get pairwise sync (non-outlet)
pairwise_sync <- sync_raw %>%
  filter(is_outlet_pair == FALSE) %>%
  left_join(solute_groups, by = "solute")

message("\nOutlet sync records:", nrow(outlet_sync), "\n")
message("Pairwise sync records:", nrow(pairwise_sync), "\n")

# =============================================================================
# KEY FINDING: DISTRIBUTION OF r VALUES BY SOLUTE GROUP
# =============================================================================
# Mean r can hide cancellation between positive and negative correlations!
# |r| (absolute value) shows coupling STRENGTH regardless of direction

message("\n=============================================================================\n")
message("KEY FINDING: DISTRIBUTION OF PEARSON r BY SOLUTE GROUP\n")
message("=============================================================================\n")

# Calculate distribution stats
r_distribution <- outlet_sync %>%
  filter(!is.na(pearson_r), !is.na(group)) %>%
  group_by(group) %>%
  summarise(
    n = n(),
    pct_negative = round(mean(pearson_r < 0) * 100, 1),
    mean_r = round(mean(pearson_r), 3),
    mean_abs_r = round(mean(abs(pearson_r)), 3),
    sd_r = round(sd(pearson_r), 3),
    cancellation = round(mean_abs_r - abs(mean_r), 3),
    .groups = "drop"
  )

message("\nDistribution of Pearson r (outlet synchrony) by Solute Group:\n")
print(r_distribution)

message("\n*** CRITICAL: Nutrients show major cancellation! ***\n")
message("  - 37% negative r values for nutrients\n")
message("  - Mean r = 0.21 (looks 'desynchronized')\n") 
message("  - Mean |r| = 0.62 (actually moderately coupled!)\n")
message("  - The low mean isn't weak coupling - it's MIXED DIRECTIONALITY\n")

# Create violin plot
group_colors <- c("Geogenic" = "#2166AC", "Nutrients" = "#B2182B", "Organic" = "#4DAF4A")

p_violin <- outlet_sync %>%
  filter(!is.na(pearson_r), !is.na(group)) %>%
  mutate(group = factor(group, levels = c("Geogenic", "Nutrients", "Organic"))) %>%
  ggplot(aes(x = group, y = pearson_r, fill = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.5) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
  scale_fill_manual(values = group_colors) +
  labs(
    title = "Distribution of Outlet Synchrony (Pearson r) by Solute Group",
    subtitle = "Nutrients show bimodal distribution: some positive, some NEGATIVE correlations",
    x = "Solute Group",
    y = "Pearson r (outlet synchrony)",
    caption = "Diamond = mean. Nutrients: mean r = 0.21 but mean |r| = 0.62"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "r_distribution_violin_by_group.png"), p_violin,
       width = 9, height = 7, dpi = 300)
ggsave(file.path(meeting_fig_dir, "r_distribution_violin_by_group.png"), p_violin,
       width = 9, height = 7, dpi = 300)
message("\nSaved: r_distribution_violin_by_group.png (to figures/ and updates/)\n")

# Create histogram faceted by group
p_hist <- outlet_sync %>%
  filter(!is.na(pearson_r), !is.na(group)) %>%
  mutate(group = factor(group, levels = c("Geogenic", "Nutrients", "Organic"))) %>%
  ggplot(aes(x = pearson_r, fill = group)) +
  geom_histogram(bins = 30, alpha = 0.8, color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray20", linewidth = 0.8) +
  facet_wrap(~group, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = group_colors) +
  labs(
    title = "Histogram of Outlet Synchrony (Pearson r) by Solute Group",
    subtitle = "Nutrients have substantial mass below zero (negative correlations)",
    x = "Pearson r (outlet synchrony)",
    y = "Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 12))

ggsave(file.path(fig_dir, "r_distribution_histogram_by_group.png"), p_hist,
       width = 8, height = 9, dpi = 300)
ggsave(file.path(meeting_fig_dir, "r_distribution_histogram_by_group.png"), p_hist,
       width = 8, height = 9, dpi = 300)
message("Saved: r_distribution_histogram_by_group.png (to figures/ and updates/)\n")

# Create comparison: mean r vs mean |r|
comparison_df <- r_distribution %>%
  select(group, mean_r, mean_abs_r) %>%
  pivot_longer(cols = c(mean_r, mean_abs_r), 
               names_to = "metric", values_to = "value") %>%
  mutate(
    metric = case_when(
      metric == "mean_r" ~ "Mean r\n(can cancel)",
      metric == "mean_abs_r" ~ "Mean |r|\n(coupling strength)"
    ),
    group = factor(group, levels = c("Geogenic", "Nutrients", "Organic"))
  )

p_comparison <- ggplot(comparison_df, aes(x = group, y = value, fill = metric)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray40") +
  scale_fill_manual(values = c("Mean r\n(can cancel)" = "#7570B3", 
                               "Mean |r|\n(coupling strength)" = "#E7298A"),
                    name = "Metric") +
  labs(
    title = "Mean r vs Mean |r| by Solute Group",
    subtitle = "Nutrients: mean r (0.21) drastically underestimates coupling strength (|r| = 0.62)",
    x = "Solute Group",
    y = "Value"
  ) +
  theme_minimal(base_size = 14) +
  coord_cartesian(ylim = c(0, 1))

ggsave(file.path(fig_dir, "mean_r_vs_abs_r_comparison.png"), p_comparison,
       width = 9, height = 6, dpi = 300)
ggsave(file.path(meeting_fig_dir, "mean_r_vs_abs_r_comparison.png"), p_comparison,
       width = 9, height = 6, dpi = 300)
message("Saved: mean_r_vs_abs_r_comparison.png (to figures/ and updates/)\n")

# =============================================================================
# PART 1: DOES CLUSTER TYPE PREDICT SYNCHRONY?
# =============================================================================

message("\n=============================================================================\n")
message("PART 1: CLUSTER TYPE → SYNCHRONY\n")
message("=============================================================================\n")

# Join cluster info to outlet sync data  
outlet_with_cluster <- outlet_sync %>%
  left_join(modal_clusters %>% select(Stream_Name, solute, modal_cluster), 
            by = c("Stream_Name", "solute"))

message("\nOutlet sync with clusters:\n")
message("  Records with cluster:", sum(!is.na(outlet_with_cluster$modal_cluster)), "\n")
message("  Records without cluster:", sum(is.na(outlet_with_cluster$modal_cluster)), "\n")

# Filter to records with clusters
sync_cluster <- outlet_with_cluster %>%
  filter(!is.na(modal_cluster), !is.na(pearson_r))

# 1a. Simple correlation: Cluster → Outlet Sync
message("\n--- 1a. Cluster Type → Outlet Synchrony ---\n")

sync_by_cluster <- sync_cluster %>%
  group_by(modal_cluster) %>%
  summarise(
    n = n(),
    mean_sync = mean(pearson_r, na.rm = TRUE),
    sd_sync = sd(pearson_r, na.rm = TRUE),
    .groups = "drop"
  )

message("\nOutlet Sync by Cluster Type:\n")
print(sync_by_cluster)

# ANOVA
anova_cluster <- aov(pearson_r ~ modal_cluster, data = sync_cluster)
message("\nANOVA - Cluster → Outlet Sync:\n")
print(summary(anova_cluster))

# Effect size (eta-squared)
ss <- summary(anova_cluster)[[1]]
eta_sq <- ss["Sum Sq"][1,1] / sum(ss["Sum Sq"])
message("\nEta-squared (effect size):", round(eta_sq, 3), "\n")

# 1b. Mixed model: Cluster → Sync with random effects
message("\n--- 1b. Mixed Model: Cluster → Outlet Sync ---\n")

# Check if we have enough levels for random effects
n_streams <- n_distinct(sync_cluster$Stream_Name)
n_solutes <- n_distinct(sync_cluster$solute)
message("Unique streams:", n_streams, ", Unique solutes:", n_solutes, "\n")

if (n_streams > 1 && n_solutes > 1) {
  m_cluster_sync <- lmer(pearson_r ~ modal_cluster + (1|Stream_Name) + (1|solute), 
                         data = sync_cluster)
  message("\nMixed Model Summary:\n")
  print(summary(m_cluster_sync))
  
  r2_cluster <- r.squaredGLMM(m_cluster_sync)
  message("\nR² (marginal):", round(r2_cluster[1], 3), "\n")
  message("R² (conditional):", round(r2_cluster[2], 3), "\n")
} else if (n_streams > 1) {
  m_cluster_sync <- lmer(pearson_r ~ modal_cluster + (1|Stream_Name), 
                         data = sync_cluster)
  r2_cluster <- r.squaredGLMM(m_cluster_sync)
  message("R² (marginal):", round(r2_cluster[1], 3), "\n")
  message("R² (conditional):", round(r2_cluster[2], 3), "\n")
} else {
  # Fall back to simple model
  m_cluster_sync <- lm(pearson_r ~ modal_cluster, data = sync_cluster)
  r2_cluster <- c(summary(m_cluster_sync)$r.squared, summary(m_cluster_sync)$r.squared)
  message("R² (simple):", round(r2_cluster[1], 3), "\n")
}

# 1c. Does cluster predict sync WITHIN solute groups?
message("\n--- 1c. Cluster → Sync by Solute Group ---\n")

cluster_sync_by_group <- sync_cluster %>%
  filter(!is.na(group)) %>%
  group_by(group) %>%
  summarise(
    n = n(),
    n_clusters = n_distinct(modal_cluster),
    .groups = "drop"
  )
print(cluster_sync_by_group)

results_cluster_sync <- list()

for (grp in unique(na.omit(sync_cluster$group))) {
  message("\n===", grp, "===\n")
  dat <- sync_cluster %>% filter(group == grp)
  
  # Summary by cluster within group
  grp_summary <- dat %>%
    group_by(modal_cluster) %>%
    summarise(n = n(), mean_sync = mean(pearson_r, na.rm = TRUE), .groups = "drop")
  print(grp_summary)
  
  if (n_distinct(dat$modal_cluster) > 1 && nrow(dat) > 30) {
    m <- try(lmer(pearson_r ~ modal_cluster + (1|Stream_Name), data = dat), silent = TRUE)
    if (!inherits(m, "try-error")) {
      r2 <- r.squaredGLMM(m)
      message("R²m:", round(r2[1], 3), "| R²c:", round(r2[2], 3), "\n")
      
      # ANOVA within group
      aov_grp <- aov(pearson_r ~ modal_cluster, data = dat)
      p_val <- summary(aov_grp)[[1]]["Pr(>F)"][1,1]
      message("ANOVA p-value:", format.pval(p_val, digits = 3), "\n")
      
      results_cluster_sync[[grp]] <- list(R2m = r2[1], R2c = r2[2], p = p_val)
    }
  } else {
    message("Insufficient clusters or data in this group\n")
  }
}

# =============================================================================
# PART 2: CATCHMENT → OUTLET SYNC BY SOLUTE GROUP
# =============================================================================

message("\n=============================================================================\n")
message("PART 2: CATCHMENT PROPERTIES → OUTLET SYNC BY SOLUTE GROUP\n")
message("=============================================================================\n")

# Prepare data
sync_catch <- outlet_sync %>%
  filter(!is.na(pearson_r), !is.na(group)) %>%
  left_join(site_chars, by = "Stream_Name") %>%
  filter(!is.na(Elevation_mean_m))

# Scale predictors (use only columns that exist)
sync_catch <- sync_catch %>%
  mutate(
    Elevation_sc = scale(Elevation_mean_m)[,1],
    Slope_sc = scale(Slope_mean)[,1],
    Area_sc = scale(Area_km2)[,1],
    DR_sc = scale(DR_Overall)[,1],
    MTT_sc = scale(MTT_overall)[,1],
    Fyw_sc = scale(Fyw_overall)[,1]
  )

message("\nData for catchment analysis: n =", nrow(sync_catch), "\n")

message("\nFull model (all solutes):\n")
m_full <- lmer(pearson_r ~ Elevation_sc + Slope_sc + Area_sc + 
               DR_sc + MTT_sc + Fyw_sc + (1|Stream_Name) + (1|solute),
               data = sync_catch, na.action = na.omit)
r2_full <- r.squaredGLMM(m_full)
message("R²m:", round(r2_full[1], 3), "| R²c:", round(r2_full[2], 3), "\n")

message("\n--- Models by Solute Group ---\n")

results_outlet <- list()

for (grp in c("Geogenic", "Nutrients", "Organic")) {
  message("\n===", grp, "===\n")
  dat <- sync_catch %>% filter(group == grp)
  message("n =", nrow(dat), "\n")
  
  if (nrow(dat) > 20) {
    # Model
    m <- try(lmer(pearson_r ~ Elevation_sc + Slope_sc + DR_sc + MTT_sc + 
                  (1|Stream_Name), data = dat, na.action = na.omit), silent = TRUE)
    
    if (!inherits(m, "try-error")) {
      r2 <- r.squaredGLMM(m)
      message("R²m:", round(r2[1], 3), "| R²c:", round(r2[2], 3), "\n")
      
      # Fixed effects
      fe <- fixef(m)
      message("\nFixed Effects:\n")
      print(round(fe, 3))
      
      results_outlet[[grp]] <- list(
        R2m = r2[1],
        R2c = r2[2],
        coefs = fe
      )
    }
  }
}

# =============================================================================
# PART 3: CATCHMENT → CLUSTER BY SOLUTE GROUP
# =============================================================================

message("\n=============================================================================\n")
message("PART 3: CATCHMENT PROPERTIES → CLUSTER BY SOLUTE GROUP\n")
message("=============================================================================\n")

cluster_catch <- modal_clusters %>%
  filter(!is.na(group)) %>%
  left_join(site_chars, by = "Stream_Name") %>%
  filter(!is.na(Elevation_mean_m), !is.na(modal_cluster))

results_cluster <- list()

for (grp in c("Geogenic", "Nutrients", "Organic")) {
  message("\n===", grp, "===\n")
  dat <- cluster_catch %>% filter(group == grp)
  message("n =", nrow(dat), "\n")
  message("Clusters:", paste(unique(dat$modal_cluster), collapse = ", "), "\n")
  
  if (nrow(dat) > 15 && n_distinct(dat$modal_cluster) > 1) {
    # Multinomial logistic
    m <- try(multinom(modal_cluster ~ Elevation_mean_m + Slope_mean + 
                      DR_Overall + MTT_overall, 
                      data = dat, trace = FALSE), silent = TRUE)
    
    if (!inherits(m, "try-error")) {
      null_m <- multinom(modal_cluster ~ 1, data = dat, trace = FALSE)
      mcfadden <- 1 - logLik(m)/logLik(null_m)
      message("McFadden R²:", round(mcfadden, 3), "\n")
      
      # Random Forest
      rf_dat <- dat %>%
        select(modal_cluster, Elevation_mean_m, Slope_mean, Area_km2,
               DR_Overall, MTT_overall, Fyw_overall) %>%
        drop_na()
      
      if (nrow(rf_dat) > 10) {
        rf <- randomForest(as.factor(modal_cluster) ~ ., data = rf_dat, 
                          ntree = 500, importance = TRUE)
        oob_acc <- 1 - rf$err.rate[500, "OOB"]
        message("RF OOB Accuracy:", round(oob_acc * 100, 1), "%\n")
        
        results_cluster[[grp]] <- list(
          McFadden = as.numeric(mcfadden),
          RF_Accuracy = oob_acc
        )
      }
    }
  }
}

# =============================================================================
# PART 4: VISUALIZATION
# =============================================================================

message("\n=============================================================================\n")
message("PART 4: GENERATING FIGURES\n")
message("=============================================================================\n")

# 4a. Synchrony barplot by solute group
sync_summary <- outlet_sync %>%
  filter(!is.na(group), !is.na(pearson_r)) %>%
  group_by(group) %>%
  summarise(
    mean_sync = mean(pearson_r, na.rm = TRUE),
    se_sync = sd(pearson_r, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(group = factor(group, levels = c("Geogenic", "Nutrients", "Organic")))

p_sync_group <- ggplot(sync_summary, aes(x = group, y = mean_sync, fill = group)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = mean_sync - se_sync, ymax = mean_sync + se_sync),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = c("Geogenic" = "#2166AC", 
                               "Nutrients" = "#B2182B",
                               "Organic" = "#4DAF4A")) +
  labs(
    title = "Outlet Synchrony by Solute Group",
    subtitle = "Nutrients show near-zero synchrony with outlet",
    x = "Solute Group",
    y = "Mean Pearson r with Outlet",
    caption = paste0("Geogenic: n=", sync_summary$n[sync_summary$group=="Geogenic"],
                    ", Nutrients: n=", sync_summary$n[sync_summary$group=="Nutrients"],
                    ", Organic: n=", sync_summary$n[sync_summary$group=="Organic"])
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(-0.1, 1))

ggsave(file.path(fig_dir, "sync_by_solute_group.png"), p_sync_group,
       width = 8, height = 6, dpi = 300)
ggsave(file.path(meeting_fig_dir, "sync_by_solute_group.png"), p_sync_group,
       width = 8, height = 6, dpi = 300)
message("Saved: sync_by_solute_group.png (to figures/ and updates/)\n")

# 4b. Cluster distribution by solute group
cluster_dist <- modal_clusters %>%
  filter(!is.na(group), !is.na(modal_cluster)) %>%
  count(group, modal_cluster) %>%
  group_by(group) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("Geogenic", "Nutrients", "Organic")))

p_cluster_group <- ggplot(cluster_dist, aes(x = group, y = pct, fill = modal_cluster)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = cluster_colors_labeled, name = "Cluster") +
  labs(
    title = "Seasonal Pattern Cluster Distribution by Solute Group",
    subtitle = "Clusters based on DTW of monthly concentration z-scores",
    x = "Solute Group",
    y = "Percentage"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "cluster_by_solute_group.png"), p_cluster_group,
       width = 9, height = 6, dpi = 300)
ggsave(file.path(meeting_fig_dir, "cluster_by_solute_group.png"), p_cluster_group,
       width = 9, height = 6, dpi = 300)
message("Saved: cluster_by_solute_group.png (to figures/ and updates/)\n")

# 4b2. FLIPPED: Solute group distribution by cluster (clusters on x-axis)
cluster_dist_flipped <- modal_clusters %>%
  filter(!is.na(group), !is.na(modal_cluster)) %>%
  count(modal_cluster, group) %>%
  group_by(modal_cluster) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("Geogenic", "Nutrients", "Organic")))

p_group_by_cluster <- ggplot(cluster_dist_flipped, aes(x = modal_cluster, y = pct, fill = group)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = group_colors, name = "Solute Group") +
  labs(
    title = "Solute Group Composition by Seasonal Pattern Cluster",
    subtitle = "Clusters based on DTW of monthly concentration z-scores",
    x = "Cluster",
    y = "Percentage"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 15, hjust = 1))

ggsave(file.path(fig_dir, "group_by_cluster.png"), p_group_by_cluster,
       width = 10, height = 6, dpi = 300)
ggsave(file.path(meeting_fig_dir, "group_by_cluster.png"), p_group_by_cluster,
       width = 10, height = 6, dpi = 300)
message("Saved: group_by_cluster.png (FLIPPED version - clusters on x-axis)\n")

# 4c. Sync by Cluster Type
if (nrow(sync_cluster) > 0) {
  sync_cluster_summary <- sync_cluster %>%
    filter(!is.na(modal_cluster)) %>%
    group_by(modal_cluster) %>%
    summarise(
      mean_sync = mean(pearson_r, na.rm = TRUE),
      se_sync = sd(pearson_r, na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = "drop"
    )
  
  p_sync_cluster <- ggplot(sync_cluster_summary, aes(x = modal_cluster, y = mean_sync, 
                                                     fill = modal_cluster)) +
    geom_col(width = 0.7) +
    geom_errorbar(aes(ymin = mean_sync - se_sync, ymax = mean_sync + se_sync),
                  width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_fill_manual(values = cluster_colors_labeled) +
    labs(
      title = "Outlet Synchrony by Seasonal Pattern Cluster",
      subtitle = paste0("ANOVA p < ", format.pval(summary(anova_cluster)[[1]]["Pr(>F)"][1,1], digits = 2),
                       ", η² = ", round(eta_sq, 2)),
      x = "Modal Cluster",
      y = "Mean Pearson r with Outlet"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none") +
    coord_cartesian(ylim = c(-0.1, 1))
  
  ggsave(file.path(fig_dir, "sync_by_cluster_type.png"), p_sync_cluster,
         width = 8, height = 6, dpi = 300)
  ggsave(file.path(meeting_fig_dir, "sync_by_cluster_type.png"), p_sync_cluster,
         width = 8, height = 6, dpi = 300)
  message("Saved: sync_by_cluster_type.png (to figures/ and updates/)\n")
}

# 4d. Sync × Cluster × Solute Group
if (nrow(sync_cluster) > 0) {
  sync_cluster_group <- sync_cluster %>%
    filter(!is.na(group), !is.na(modal_cluster)) %>%
    group_by(group, modal_cluster) %>%
    summarise(
      mean_sync = mean(pearson_r, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(group = factor(group, levels = c("Geogenic", "Nutrients", "Organic")))
  
  p_sync_cluster_group <- ggplot(sync_cluster_group, 
                                 aes(x = modal_cluster, y = mean_sync, fill = group)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = c("Geogenic" = "#2166AC", 
                                 "Nutrients" = "#B2182B",
                                 "Organic" = "#4DAF4A"),
                      name = "Solute Group") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      title = "Outlet Synchrony by Seasonal Cluster and Solute Group",
      subtitle = "Nutrients show mixed synchrony (positive + negative r) regardless of cluster",
      x = "Modal Cluster",
      y = "Mean Pearson r with Outlet"
    ) +
    theme_minimal(base_size = 14) +
    coord_cartesian(ylim = c(-0.2, 1))
  
  ggsave(file.path(fig_dir, "sync_by_cluster_and_group.png"), p_sync_cluster_group,
         width = 10, height = 6, dpi = 300)
  ggsave(file.path(meeting_fig_dir, "sync_by_cluster_and_group.png"), p_sync_cluster_group,
         width = 10, height = 6, dpi = 300)
  message("Saved: sync_by_cluster_and_group.png (to figures/ and updates/)\n")
}

# =============================================================================
# SUMMARY TABLE
# =============================================================================

message("\n=============================================================================\n")
message("SUMMARY OF RESULTS\n")
message("=============================================================================\n")

message("\n--- Cluster → Sync ---\n")
message("Overall: R²m =", round(r2_cluster[1], 3), "\n")
message("ANOVA p-value:", format.pval(summary(anova_cluster)[[1]]["Pr(>F)"][1,1], digits = 2), "\n")
message("Effect size (η²):", round(eta_sq, 3), "\n")

message("\n--- Cluster → Sync by Group ---\n")
for (grp in names(results_cluster_sync)) {
  message(grp, ": R²m =", round(results_cluster_sync[[grp]]$R2m, 3), 
      ", p =", format.pval(results_cluster_sync[[grp]]$p, digits = 2), "\n")
}

message("\n--- Catchment → Sync by Group ---\n")
for (grp in names(results_outlet)) {
  message(grp, ": R²m =", round(results_outlet[[grp]]$R2m, 3), "\n")
}

message("\n--- Catchment → Cluster by Group ---\n")
for (grp in names(results_cluster)) {
  message(grp, ": McFadden R² =", round(results_cluster[[grp]]$McFadden, 3),
      ", RF Accuracy =", round(results_cluster[[grp]]$RF_Accuracy * 100, 1), "%\n")
}

message("\n=============================================================================\n")
message("ANALYSIS COMPLETE\n")
message("=============================================================================\n")
