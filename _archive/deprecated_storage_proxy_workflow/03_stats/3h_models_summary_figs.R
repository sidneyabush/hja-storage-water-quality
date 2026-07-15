suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
source(file.path(repo_dir, "00_helpers", "plot_prefs.R"))

base_dir    <- "/Users/sidneybush/Library/CloudStorage/Box-Box"
project_dir <- file.path(base_dir, "Sidney_Bush", "HJA_Water_Quality")
out_stats   <- file.path(project_dir, "outputs", "03_stats")
fig_root    <- file.path(project_dir, "exploratory_plots", "03_stats", "3h_mixed_effects")
dir.create(fig_root, showWarnings = FALSE, recursive = TRUE)

# Load model summaries
summary_file <- file.path(out_stats, "model_summary_all.csv")
model1_file  <- file.path(out_stats, "model1_conc_sync_comparison.csv")
model2_file  <- file.path(out_stats, "model2_cqslope_sync_comparison.csv")

if (!file.exists(summary_file)) stop("model_summary_all.csv not found in outputs/03_stats")

summary_all <- readr::read_csv(summary_file, show_col_types = FALSE)
model1_cmp  <- if (file.exists(model1_file)) readr::read_csv(model1_file, show_col_types = FALSE) else NULL
model2_cmp  <- if (file.exists(model2_file)) readr::read_csv(model2_file, show_col_types = FALSE) else NULL

# Identify top models per response by lowest AIC (or highest R2_marginal when available)
top_by_response <- summary_all %>%
  filter(!is.na(Response)) %>%
  group_by(Response) %>%
  arrange(AIC, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

readr::write_csv(top_by_response, file.path(out_stats, "model_top_by_response.csv"))

# Figure A: Bar showing AIC and R2m for top models per response
figA <- top_by_response %>%
  mutate(Model = forcats::fct_inorder(Model)) %>%
  ggplot(aes(x = Model, y = -AIC, fill = Response)) +
  geom_col() +
  geom_text(aes(label = paste0("R2m=", sprintf("%.2f", R2_marginal))),
            vjust = -0.5, size = 3) +
  labs(title = "Top Mixed-Effects Models per Response",
       subtitle = "Bars show -AIC (higher is better); labels show marginal R2",
       x = "Model", y = "-AIC") +
  theme_hja() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggplot2::ggsave(file.path(fig_root, "figA_top_models_per_response.png"), figA, width = 10, height = 6, dpi = 300)

# Figure B: Comparison tables saved as PNG via simple ggplot text grob (optional)
tbl_text <- top_by_response %>%
  mutate(line = sprintf("%s: %s (AIC=%.0f, R2m=%.2f)", Response, Model, AIC, R2_marginal)) %>%
  mutate(idx = row_number())

figB <- ggplot(tbl_text, aes(y = rev(idx), x = 0, label = line)) +
  geom_text(hjust = 0, family = "mono") +
  labs(title = "Top Models Summary", x = NULL, y = NULL) +
  theme_void() +
  xlim(0, 1)

ggplot2::ggsave(file.path(fig_root, "figB_top_models_summary_text.png"), figB, width = 8, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# CQ slope vs Cluster by Solute Group (Geogenic / Biogenic / Nutrient)
# -----------------------------------------------------------------------------

# Load site means or seasonal data with cq_slope and Cluster
seasonal_file <- file.path(project_dir, "outputs", "HJA_master_seasonal.csv")
site_means_file <- file.path(project_dir, "outputs", "HJA_clean_site_means.csv")

df <- if (file.exists(site_means_file)) {
  readr::read_csv(site_means_file, show_col_types = FALSE)
} else if (file.exists(seasonal_file)) {
  readr::read_csv(seasonal_file, show_col_types = FALSE)
} else {
  stop("No source data found for CQ slope vs Cluster plots")
}

# Prefer Cluster_mode_wy if Cluster not present
cluster_col <- if ("Cluster" %in% names(df)) "Cluster" else if ("Cluster_mode_wy" %in% names(df)) "Cluster_mode_wy" else NA
if (is.na(cluster_col)) stop("Required cluster column missing (Cluster or Cluster_mode_wy)")
if (!all(c("cq_slope", "solute") %in% names(df))) stop("Required columns missing: cq_slope or solute")

# Add solute group labels from plot_prefs helpers
df <- df %>% mutate(solute_group = categorize_solute(solute))

# Compute site-mean per cluster/solute
df_grouped <- df %>%
  rename(Cluster = !!cluster_col) %>%
  group_by(Cluster, solute_group, solute) %>%
  summarise(cq_slope_mean = mean(cq_slope, na.rm = TRUE), .groups = "drop")

# Order clusters and groups
df_grouped$Cluster <- factor(df_grouped$Cluster, levels = c("1", "2", "3", "4"))
df_grouped$solute_group <- factor(df_grouped$solute_group, levels = c("Geogenic", "Biogenic", "Nutrient"))

figC <- ggplot(df_grouped, aes(x = solute, y = cq_slope_mean, fill = solute_group)) +
  geom_col() +
  facet_grid(solute_group ~ Cluster, labeller = labeller(
    Cluster = c("1" = "1-Baseflow Enriched", "2" = "2-Chemostatic",
                "3" = "3-Spring/Early Summer Enriched", "4" = "4-Winter Flushing")
  )) +
  scale_fill_manual(values = solute_type_colors, name = "Solute Group") +
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey60") +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey60") +
  labs(title = "CQ Slope by Cluster Across Solute Groups",
       x = "Solute", y = "Mean CQ Slope (site/season means)") +
  theme_hja() +
  theme(strip.text = element_text(face = "bold"))

ggplot2::ggsave(file.path(fig_root, "figC_cqslope_by_cluster_and_group.png"), figC, width = 12, height = 8, dpi = 300)

## -----------------------------------------------------------------------------
## Concentration synchrony vs storage by solute group
## -----------------------------------------------------------------------------

site_means <- readr::read_csv(site_means_file, show_col_types = FALSE)
if (all(c("conc_sync_outlet", "Q_dS_range_mm", "solute") %in% names(site_means))) {
  site_means <- site_means %>% mutate(solute_group = categorize_solute(solute))
  # Compute group-wise r and p for labels
  stats_by_group <- site_means %>%
    group_by(solute_group) %>%
    summarise(
      r = cor(conc_sync_outlet, Q_dS_range_mm, use = "pairwise.complete.obs"),
      n = sum(!is.na(conc_sync_outlet) & !is.na(Q_dS_range_mm)),
      .groups = "drop"
    )
  label_df <- stats_by_group %>% mutate(label = sprintf("r=%.2f (n=%d)", r, n))

  figD <- ggplot(site_means, aes(x = Q_dS_range_mm, y = conc_sync_outlet, color = solute_group)) +
    geom_point(alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    facet_wrap(~ solute_group, ncol = 3) +
    scale_color_manual(values = solute_type_colors, name = "Solute Group") +
    labs(title = "Concentration Synchrony vs Storage by Solute Group",
         x = get_label("Q_dS_range_mm"),
         y = "Outlet Concentration Synchrony") +
    theme_hja()

  ggplot2::ggsave(file.path(fig_root, "figD_conc_sync_vs_storage_by_group.png"), figD, width = 12, height = 6, dpi = 300)
}

## -----------------------------------------------------------------------------
## Concentration synchrony vs storage by solute group AND cluster
## -----------------------------------------------------------------------------

if (all(c("conc_sync_outlet", "Q_dS_range_mm", "solute") %in% names(site_means))) {
  cluster_col2 <- if ("Cluster" %in% names(site_means)) "Cluster" else if ("Cluster_mode_wy" %in% names(site_means)) "Cluster_mode_wy" else NA
  if (!is.na(cluster_col2)) {
    site_means2 <- site_means %>%
      mutate(solute_group = categorize_solute(solute)) %>%
      rename(Cluster = !!cluster_col2) %>%
      mutate(Cluster = factor(Cluster, levels = c("1", "2", "3", "4")))

    figE <- ggplot(site_means2, aes(x = Q_dS_range_mm, y = conc_sync_outlet, color = solute_group)) +
      geom_point(alpha = 0.8, size = 2) +
      geom_smooth(method = "lm", se = TRUE, color = "black") +
      facet_grid(solute_group ~ Cluster, labeller = labeller(
        Cluster = c("1" = "1-Baseflow Enriched", "2" = "2-Chemostatic",
                    "3" = "3-Spring/Early Summer Enriched", "4" = "4-Winter Flushing")
      )) +
      scale_color_manual(values = solute_type_colors, name = "Solute Group") +
      labs(title = "Concentration Synchrony vs Storage by Group and Cluster",
           x = get_label("Q_dS_range_mm"),
           y = "Outlet Concentration Synchrony") +
      theme_hja()

    ggplot2::ggsave(file.path(fig_root, "figE_conc_sync_vs_storage_by_group_and_cluster.png"), figE, width = 12, height = 8, dpi = 300)
  }
}

message("Figures written to:", fig_root, "\n")
