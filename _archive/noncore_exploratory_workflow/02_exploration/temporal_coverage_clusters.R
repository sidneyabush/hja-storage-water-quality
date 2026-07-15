# =============================================================================
# Temporal Coverage and Cluster Distribution Figure
# =============================================================================
# Creates multi-panel figure showing:
# A) Temporal data coverage (1969-2017)
# B) Cluster distribution by decade
# C) Era comparison (pre/post 2000)
# D) Site × Year data availability

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(lubridate)
})

# Source helpers
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
source(file.path(repo_dir, "00_helpers", "plot_prefs.R"))

# Paths
base_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality"
output_dir <- file.path(base_dir, "outputs")
plot_dir <- file.path(base_dir, "exploratory_plots", "01_data_prep", "conc_sample_bias")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
cq_master <- read_csv(file.path(output_dir, "HJA_CQ_master.csv"), show_col_types = FALSE)
cluster_data <- read_csv(file.path(output_dir, "ClusterStreams_allSolutes_byWaterYear.csv"),
                        show_col_types = FALSE)

# =============================================================================
# PANEL A: Temporal Data Coverage (site-solute combinations per year)
# =============================================================================

# Count site-solute combinations per calendar year
coverage_by_year <- cq_master %>%
  mutate(calendar_year = Year) %>%
  group_by(calendar_year, Stream_Name, variable) %>%
  summarise(n_samples = n(), .groups = "drop") %>%
  group_by(calendar_year) %>%
  summarise(n_combinations = n(), .groups = "drop")

# Define coverage tiers based on number of site-solute combinations
# Max possible = 9 sites × 11 solutes = 99 combinations
coverage_by_year <- coverage_by_year %>%
  mutate(
    coverage_tier = case_when(
      n_combinations >= 90 ~ "Full (90+)",
      n_combinations >= 50 ~ "Good (50-89)",
      n_combinations >= 35 ~ "Moderate (35-49)",
      TRUE ~ "Sparse (<35)"
    ),
    coverage_tier = factor(coverage_tier,
                          levels = c("Full (90+)", "Good (50-89)",
                                   "Moderate (35-49)", "Sparse (<35)"))
  )

# Panel A plot
p_coverage <- ggplot(coverage_by_year, aes(x = calendar_year, y = n_combinations, fill = coverage_tier)) +
  geom_col() +
  geom_hline(yintercept = c(50, 90), linetype = "dashed", color = "gray40") +
  scale_fill_manual(
    values = c("Full (90+)" = "#08519c",
               "Good (50-89)" = "#6baed6",
               "Moderate (35-49)" = "#c6dbef",
               "Sparse (<35)" = "#fcbba1"),
    name = "Coverage Tier"
  ) +
  scale_x_continuous(breaks = seq(1970, 2020, 10)) +
  scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 105)) +
  labs(
    x = "Calendar Year",
    y = "Site-Solute Combinations"
  ) +
  theme_hja() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

# =============================================================================
# PANEL B: Cluster Distribution by Decade
# =============================================================================

cluster_by_decade <- cluster_data %>%
  mutate(
    decade = paste0(floor(water_year / 10) * 10, "s"),
    Cluster = Cluster_mode  # Rename for consistency
  ) %>%
  count(decade, Cluster) %>%
  group_by(decade) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    Cluster = factor(Cluster, levels = cluster_levels),
    decade = factor(decade, levels = c("1960s", "1970s", "1980s", "1990s", "2000s", "2010s"))
  )

p_decade <- ggplot(cluster_by_decade, aes(x = decade, y = prop, fill = Cluster)) +
  geom_col() +
  scale_fill_manual(values = cluster_colors, name = "Cluster") +
  scale_y_continuous(labels = scales::percent_format(), breaks = seq(0, 1, 0.25)) +
  labs(
    x = "Decade",
    y = "Proportion"
  ) +
  theme_hja() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

# =============================================================================
# PANEL C: Era Comparison (Pre/Post 2000)
# =============================================================================

cluster_by_era <- cluster_data %>%
  mutate(
    era = if_else(water_year < 2000, "Pre-2000\n(1969-1999)", "Post-2000\n(2000-2017)"),
    Cluster = Cluster_mode  # Rename for consistency
  ) %>%
  count(era, Cluster) %>%
  group_by(era) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    Cluster = factor(Cluster, levels = cluster_levels),
    era = factor(era, levels = c("Post-2000\n(2000-2017)", "Pre-2000\n(1969-1999)"))
  )

# Chi-square test
era_matrix <- cluster_data %>%
  mutate(
    era = if_else(water_year < 2000, "Pre-2000", "Post-2000"),
    Cluster = Cluster_mode  # Rename for consistency
  ) %>%
  count(era, Cluster) %>%
  pivot_wider(names_from = Cluster, values_from = n, values_fill = 0) %>%
  select(-era) %>%
  as.matrix()

chi_test <- chisq.test(era_matrix)

p_era <- ggplot(cluster_by_era, aes(x = era, y = prop, fill = Cluster)) +
  geom_col() +
  scale_fill_manual(values = cluster_colors, name = "Cluster") +
  scale_y_continuous(labels = scales::percent_format(), breaks = seq(0, 1, 0.2)) +
  labs(
    x = NULL,
    y = "Proportion"
  ) +
  theme_hja() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

# =============================================================================
# PANEL D: Site × Year Data Availability (heatmap)
# =============================================================================

site_year_availability <- cq_master %>%
  mutate(calendar_year = Year) %>%
  group_by(Stream_Name, calendar_year) %>%
  summarise(n_solutes = n_distinct(variable), .groups = "drop") %>%
  mutate(
    Stream_Name = factor(Stream_Name, levels = site_order)
  )

p_heatmap <- ggplot(site_year_availability, aes(x = calendar_year, y = Stream_Name, fill = n_solutes)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(
    option = "viridis",
    name = "N Solutes",
    breaks = seq(0, 11, 1),
    limits = c(0, 11)
  ) +
  scale_x_continuous(breaks = seq(1970, 2020, 10)) +
  labs(
    x = "Calendar Year",
    y = "Site"
  ) +
  theme_hja() +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 10)
  )

# =============================================================================
# Combine all panels
# =============================================================================

layout <- "
AABB
CCDD
"

p_combined <- p_coverage + p_decade + p_era + p_heatmap +
  plot_layout(design = layout, guides = "collect") &
  theme(legend.position = "bottom")

# Assign final plot
p_final <- p_combined

# Save
ggsave(
  file.path(plot_dir, "temporal_coverage_clusters.png"),
  plot = p_final,
  width = 14,
  height = 10,
  dpi = 300,
  bg = "white"
)

message("\n✓ Temporal coverage and cluster distribution figure saved")
