# =============================================================================
# STEP 02b: CLUSTER VISUALIZATIONS
# =============================================================================
# Analysis of concentration climatology clusters:
#   - Section 1: Cluster Patterns (climatology, water-year, modal)
#   - Section 2: Cluster Stability Analysis (inter-annual persistence)
#   - Section 3: Cluster Controls (catchment, hydrology, storage)
#   - Section 4: Cluster Transitions (what drives cluster shifts?)
#   - Section 5: Cluster × CQ × Synchrony Relationships
#   - Section 6: Cluster Patterns by Site
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(forcats)
  library(patchwork)
  library(tidyr)
  library(grid)
  library(corrplot)
})

rm(list = ls())

# Source helpers from repo root
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
source(file.path(repo_dir, "00_helpers", "plot_prefs.R"))

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
base_dir    <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality"
data_path   <- file.path(base_dir, "data")
output_dir  <- file.path(base_dir, "outputs")

# Create plot directories organized by clustering type
plot_base     <- file.path(base_dir, "exploratory_plots", "02_exploration", "2b_clusters")

# CLIMATOLOGICAL clustering (average across all years)
plot_climatological <- file.path(plot_base, "climatological")

# ANNUAL clustering (by calendar year) - all plots in one folder
plot_annual <- file.path(plot_base, "annual")

# Legacy variables for compatibility (all point to annual folder)
plot_stability <- plot_annual
plot_controls <- plot_annual
plot_transitions <- plot_annual
plot_relationships <- plot_annual
plot_bysite   <- plot_annual

dir.create(plot_climatological, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_annual, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Shared theme and palettes
# =============================================================================
month_labels <- c("J","F","M","A","M","J","J","A","S","O","N","D")

# Solute palettes (standardized lists from helpers)
geo_solutes <- if (exists("GEO_SOLUTES")) GEO_SOLUTES else c("Ca","Mg","Na","K","DSi","Cl","SO4")
bio_solutes <- if (exists("BIO_SOLUTES")) BIO_SOLUTES else c("DOC","NH3","NO3","PO4")

geo_palette <- c(Ca = "#F9D5A7", Mg = "#F8A978", Na = "#F4976C",
                 K = "#D7867E", DSi = "#A75D5D", Cl = "#85586F", SO4 = "#E3B778")
bio_palette <- c(DOC = "#A7D0CD", NH3 = "#74B49B", NO3 = "#508CA4", PO4 = "#87A8A4")

# Cluster colors pulled from the updated reference palette (helper fallback)
if (!exists("cluster_colors")) {
  cluster_colors <- c(
    "1" = "#CFA980",
    "2" = "#98B89F",
    "3" = "#5E8AA1",
    "4" = "#526B8E"
  )
}

if (!exists("site_order")) {
  site_order <- c("GSWS09", "GSWS10", "GSWS01", "GSLOOK", "GSWS02", "GSWS06", "GSWS07", "GSWS08", "GSMACK")
}
if (!exists("solute_order")) {
  solute_order <- c("Ca", "Mg", "Na", "K", "DSi", "Cl", "SO4", "DOC", "NH3", "NO3", "PO4")
}
solute_colors <- c(geo_palette, bio_palette)

if (!exists("cluster_levels")) {
  cluster_levels <- c("1", "2", "3", "4")
}

ensure_cluster_levels <- function(df, cols) {
  cols <- cols[cols %in% names(df)]
  if (length(cols) == 0) {
    return(df)
  }
  df[cols] <- lapply(df[cols], function(x) factor(as.character(x), levels = cluster_levels))
  df
}

# =============================================================================
# LOAD ALL DATA
# =============================================================================
message("\n=== LOADING DATA ===\n")

# Climatology clusters
raw <- readr::read_csv(file.path(output_dir, "ClusterStreams_allSolutes.csv"), show_col_types = FALSE) %>%
  select(-water_year) %>%
  rename(Jan = `1`, Feb = `2`, Mar = `3`, Apr = `4`, May = `5`, Jun = `6`,
         Jul = `7`, Aug = `8`, Sep = `9`, Oct = `10`, Nov = `11`, Dec = `12`,
         Cluster = Cluster_climRef) %>%
  { df <- ensure_cluster_levels(., "Cluster"); apply_factor_orders(df) }

# Water-year clusters
cluster_wy <- readr::read_csv(file.path(output_dir, "ClusterStreams_allSolutes_byWaterYear.csv"), show_col_types = FALSE) %>%
  rename(Cluster = Cluster_climRef) %>%
  { df <- ensure_cluster_levels(., "Cluster"); apply_factor_orders(df) }

# Modal clusters
cluster_modal <- readr::read_csv(file.path(output_dir, "ClusterStreams_allSolutes_modalClusters.csv"), show_col_types = FALSE) %>%
  { df <- ensure_cluster_levels(., "Cluster_mode"); apply_factor_orders(df) }

# Stability metrics
cluster_stability <- readr::read_csv(file.path(output_dir, "ClusterStreams_stability_metrics.csv"), show_col_types = FALSE) %>%
  { df <- ensure_cluster_levels(., "Cluster"); apply_factor_orders(df) }

# Catchment characteristics
catchment <- readr::read_csv(file.path(output_dir, "Catchment_site_characteristics.csv"), show_col_types = FALSE)

# Site-level data (with hydro metrics)
site_data <- tryCatch({
  readr::read_csv(file.path(output_dir, "HJA_exploratory_site.csv"), show_col_types = FALSE)
}, error = function(e) NULL)

# Season boundaries
season_bounds <- tryCatch({
  readr::read_csv(file.path(output_dir, "season_boundaries.csv"), show_col_types = FALSE)
}, error = function(e) NULL)

# CQ rolling data
cq_data <- tryCatch({
  readr::read_csv(file.path(output_dir, "CQ_rolling_window_results.csv"), show_col_types = FALSE)
}, error = function(e) NULL)

# Mega dataset (hydro + CQ + static)
mega <- tryCatch({
  readr::read_csv(file.path(output_dir, "HJA_mega_90d_windows_CQ_hydro_static.csv"), show_col_types = FALSE)
}, error = function(e) NULL)

message("  Data loaded successfully\n")

# =============================================================================
# SECTION 1: CLUSTER PATTERNS (CLIMATOLOGY)
# =============================================================================
message("\n=== SECTION 1: CLUSTER PATTERNS ===\n")

# Convert to long format
long <- raw %>%
  pivot_longer(Jan:Dec, names_to = "month", values_to = "z") %>%
  mutate(month = factor(month, levels = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")),
         month_num = as.numeric(month))

present_solutes <- unique(long$chemical)
all_solutes <- c(geo_solutes[geo_solutes %in% present_solutes], bio_solutes[bio_solutes %in% present_solutes])
long <- long %>% mutate(chemical = factor(chemical, levels = all_solutes))

# Mean by cluster (overall)
mean_by_cluster <- long %>%
  group_by(Cluster, month_num) %>%
  summarise(mean_z = mean(z, na.rm = TRUE), .groups = "drop")

p_cluster_pattern <- ggplot(mean_by_cluster, aes(x = month_num, y = mean_z, color = Cluster, fill = Cluster)) +
  geom_ribbon(aes(ymin = mean_z - 0.1, ymax = mean_z + 0.1), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = 1:12, labels = month_labels) +
  scale_color_cluster() +
  scale_fill_cluster() +
  labs(x = "Month", y = "Mean Normalized Concentration") +
  theme_clean()

ggsave(file.path(plot_climatological, "climatological_cluster_patterns.png"), p_cluster_pattern, width = 10, height = 7, dpi = 300)

# Mean by solute × cluster
mean_solute_cluster <- long %>%
  group_by(chemical, Cluster, month_num) %>%
  summarise(mean_z = mean(z, na.rm = TRUE), .groups = "drop")

p_solute_cluster <- ggplot(mean_solute_cluster, aes(x = month_num, y = mean_z, color = chemical)) +
  geom_line(linewidth = 1, alpha = 0.8) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:12, labels = month_labels) +
  scale_color_solute() +
  facet_wrap(~ Cluster, ncol = 2, labeller = labeller(Cluster = function(x) paste("Cluster", x))) +
  labs(x = "Month", y = "Normalized Concentration", color = "Solute") +
  theme_clean()

ggsave(file.path(plot_climatological, "climatological_patterns_by_solute.png"), p_solute_cluster, width = 12, height = 10, dpi = 300)

# Spaghetti plot: all individual site-solute traces colored by cluster
p_spaghetti <- ggplot(
  long,
  aes(
    x = month_num,
    y = z,
    group = interaction(Stream_Name, chemical),
    color = Cluster
  )
) +
  geom_line(alpha = 0.45, linewidth = 1.2) +
  scale_x_continuous(breaks = 1:12, labels = month_labels) +
  scale_color_cluster() +
  facet_wrap(~ chemical, ncol = 3) +
  labs(
    x = "Month",
    y = "Normalized Concentration",
    color = "Cluster"
  )

ggsave(file.path(plot_climatological, "climatological_spaghetti_all_traces.png"), p_spaghetti, width = 14, height = 12, dpi = 300)

# Climatological patterns by site
p_site_patterns <- ggplot(long, aes(x = month_num, y = z, color = chemical, group = chemical)) +
  geom_line(linewidth = 1.2, alpha = 0.7) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_x_continuous(breaks = 1:12, labels = month_labels) +
  scale_color_solute() +
  facet_wrap(~ Stream_Name, ncol = 3) +
  labs(
    x = "Month",
    y = "Normalized Concentration",
    color = "Solute"
  )

ggsave(file.path(plot_climatological, "climatological_patterns_by_site.png"), p_site_patterns, width = 14, height = 11, dpi = 300)

# Cluster distribution by site (climatological)
cluster_by_site_clim <- raw %>%
  count(Stream_Name, Cluster) %>%
  mutate(
    Stream_Name = factor(Stream_Name, levels = site_order),
    Cluster = factor(Cluster, levels = cluster_levels)
  ) %>%
  tidyr::complete(
    Stream_Name,
    Cluster = factor(cluster_levels, levels = cluster_levels),
    fill = list(n = 0)
  ) %>%
  arrange(Stream_Name, desc(Cluster)) %>%
  group_by(Stream_Name) %>%
  mutate(
    total = sum(n),
    prop = ifelse(total > 0, n / total, 0),
    # Calculate label position (cumulative proportion - half of current segment)
    label_y = cumsum(prop) - 0.5 * prop,
    # Format percentage label
    label = ifelse(prop > 0.05, paste0(round(prop * 100), "%"), "")
  ) %>%
  ungroup()

p_site_clusters <- ggplot(cluster_by_site_clim, aes(x = Stream_Name, y = prop, fill = Cluster)) +
  geom_col(position = "stack", alpha = 0.8) +
  geom_text(aes(y = label_y, label = label), size = 3.5, color = "white") +
  scale_fill_cluster() +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "Site",
    y = "Proportion",
    fill = "Cluster"
  )

ggsave(file.path(plot_climatological, "climatological_cluster_distribution_by_site.png"), p_site_clusters, width = 10, height = 7, dpi = 300)

# Solute × Site matrix (tile plot) - climatological
p_matrix_clim <- ggplot(raw, aes(x = chemical, y = Stream_Name, fill = Cluster)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Cluster), color = "white", size = 4) +
  scale_fill_cluster() +
  labs(
    x = "Solute",
    y = "Site",
    fill = "Cluster",
  ) +
  theme_clean() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plot_climatological, "climatological_solute_by_site_cluster_matrix.png"), p_matrix_clim, width = 10, height = 8, dpi = 300)

# Combined cluster characterization: temporal patterns + solute composition
# Create panel for each cluster
cluster_panels_clim <- list()

for (cl in levels(long$Cluster)) {
  cl_data <- long %>% filter(Cluster == cl)
  cl_pct <- round(100 * nrow(cl_data %>% distinct(Stream_Name, chemical)) /
                        nrow(raw %>% distinct(Stream_Name, chemical)), 1)

  # Left panel: temporal pattern (spaghetti plot colored by solute)
  p_left <- ggplot(cl_data, aes(x = month_num, y = z,
                                group = interaction(Stream_Name, chemical),
                                color = chemical)) +
    geom_line(alpha = 0.45, linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
    scale_x_continuous(breaks = 1:12, labels = month_labels) +
    scale_color_solute() +
    labs(
      x = "Month",
      y = "Normalized Concentration",
      color = "Solute"
    ) +
    theme_clean(base_size = 11) +
    theme(legend.position = "none")

  # Right panel: solute composition (count bar chart)
  solute_counts <- cl_data %>%
    distinct(Stream_Name, chemical) %>%
    count(chemical) %>%
    complete(chemical = factor(solute_order, levels = solute_order), fill = list(n = 0)) %>%
    mutate(chemical = factor(chemical, levels = solute_order))

  p_right <- ggplot(solute_counts, aes(x = chemical, y = n, fill = chemical)) +
    geom_col(alpha = 0.8) +
    scale_fill_solute() +
    labs(x = "Solute", y = "Count") +
    theme_clean(base_size = 11) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
    )

  # Combine left and right with adjusted widths (spaghetti gets more space)
  cluster_panels_clim[[cl]] <- p_left + p_right + plot_layout(widths = c(2, 1))
}

# Combine all clusters vertically with shared legend
combined_clim <- wrap_plots(cluster_panels_clim, ncol = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(
    theme = theme(
    )
  )

ggsave(file.path(plot_climatological, "climatological_cluster_characterization_combined.png"),
       combined_clim, width = 12, height = 16, dpi = 300)

# =============================================================================
# SECTION 2: CLUSTER STABILITY ANALYSIS
# =============================================================================

# Stability by solute
# Define shapes for each site (9 sites need 9 distinct shapes)
site_shapes <- c(
  "GSWS09" = 21,  # circle
  "GSWS10" = 22,  # square
  "GSWS01" = 23,  # diamond
  "GSLOOK" = 24,  # triangle up
  "GSWS02" = 25,  # triangle down
  "GSWS06" = 21,  # circle (repeated)
  "GSWS07" = 22,  # square (repeated)
  "GSWS08" = 23,  # diamond (repeated)
  "GSMACK" = 24   # triangle up (repeated)
)

p_stab_solute <- cluster_stability %>%
  filter(!is.na(Stream_Name)) %>%
  mutate(
    Stream_Name = forcats::fct_relevel(Stream_Name, site_order),
    chemical = factor(chemical, levels = solute_order)
  ) %>%
  ggplot(aes(x = chemical, y = stability)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, fill = "gray90", color = "gray40", linewidth = 0.5) +
  geom_jitter(aes(fill = Stream_Name, shape = Stream_Name), width = 0.18, alpha = 0.85, size = 4, stroke = 0.3, color = "black") +
  scale_fill_site() +
  scale_shape_manual(values = site_shapes, name = "Site") +
  labs(x = "Solute", y = "Cluster stability (1 = never changes)", fill = "Site") +
  theme_clean() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(clip = "off") +
  guides(
    fill = guide_legend(override.aes = list(size = 3)),
    shape = guide_legend(override.aes = list(size = 3))
  )

ggsave(file.path(plot_stability, "annual_stability_by_solute.png"), p_stab_solute, width = 11, height = 8, dpi = 300)

# Stability by site
# Define shapes for each solute (11 solutes need 11 distinct shapes)
solute_shapes <- c(
  "Ca" = 16,   # filled circle
  "Mg" = 15,   # filled square
  "Na" = 17,   # filled triangle
  "K" = 18,    # filled diamond
  "DSi" = 19,  # larger filled circle
  "Cl" = 1,    # open circle
  "SO4" = 0,   # open square
  "DOC" = 2,   # open triangle
  "NH3" = 5,   # open diamond
  "NO3" = 6,   # inverted triangle
  "PO4" = 8    # star
)

p_stab_site <- ggplot(cluster_stability, aes(x = Stream_Name, y = stability)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, fill = "gray90", color = "gray40", linewidth = 0.5) +
  geom_jitter(aes(color = chemical, shape = chemical), width = 0.2, alpha = 0.7, size = 4) +
  scale_color_solute() +
  scale_shape_manual(values = solute_shapes, name = "Solute") +
  labs(x = "Site", y = "Cluster stability (1 = never changes)") +
  theme_clean() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plot_stability, "annual_stability_by_site.png"), p_stab_site, width = 12, height = 7, dpi = 300)

# Stability by site-solute (faceted by site)
# Horizontal bars for easier cross-site comparison
# Reference lines: 0.5 = stable half the time, 0.75 = stable 3/4 of the time
p_stab_site_solute <- cluster_stability %>%
  filter(!is.na(Stream_Name)) %>%
  mutate(
    Stream_Name = factor(Stream_Name, levels = site_order),
    chemical = factor(chemical, levels = rev(solute_order))  # Reverse so Ca is at top
  ) %>%
  ggplot(aes(x = stability, y = chemical, fill = chemical)) +
  geom_col(alpha = 0.8) +
  geom_vline(xintercept = c(0.25, 0.5, 0.75), linetype = "dashed", color = "gray50", linewidth = 0.3) +
  scale_fill_solute() +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  facet_wrap(~ Stream_Name, ncol = 3) +
  labs(
    y = "Solute",
    x = "Cluster stability (1 = never changes)",
  ) +
  theme_clean() +
  theme(
    axis.text.y = element_text(size = 9),
    strip.text = element_text(),
    legend.position = "none"
  )

ggsave(file.path(plot_stability, "annual_stability_by_site_solute_faceted.png"), p_stab_site_solute, width = 20, height = 12, dpi = 300)

# Transition matrix
cluster_transitions <- cluster_wy %>%
  arrange(Stream_Name, chemical, water_year) %>%
  group_by(Stream_Name, chemical) %>%
  mutate(Cluster_prev = lag(Cluster)) %>%
  filter(!is.na(Cluster_prev)) %>%
  ungroup() %>%
  mutate(
    Cluster = factor(as.character(Cluster), levels = cluster_levels),
    Cluster_prev = factor(as.character(Cluster_prev), levels = cluster_levels)
  )

transition_matrix <- cluster_transitions %>%
  count(Cluster_prev, Cluster, .drop = FALSE) %>%
  group_by(Cluster_prev) %>%
  mutate(
    total = sum(n),
    proportion = ifelse(total > 0, n / total, NA_real_)
  ) %>%
  ungroup()

stay_rates <- transition_matrix %>%
  filter(Cluster_prev == Cluster) %>%
  transmute(
    Cluster_prev,
    stay_prop = proportion,
    stay_label = paste0(Cluster_prev, ": ", scales::percent(stay_prop, accuracy = 1))
  )

transition_delta <- crossing(
  Cluster_prev = factor(cluster_levels, levels = cluster_levels),
  Cluster = factor(cluster_levels, levels = cluster_levels)
) %>%
  filter(Cluster_prev != Cluster) %>%
  left_join(
    cluster_transitions %>%
      filter(Cluster != Cluster_prev) %>%
      count(Cluster_prev, Cluster, name = "n_transitions"),
    by = c("Cluster_prev", "Cluster")
  ) %>%
  left_join(
    cluster_transitions %>%
      filter(Cluster != Cluster_prev) %>%
      count(Cluster_prev, name = "total_transitions"),
    by = "Cluster_prev"
  ) %>%
  mutate(
    n_transitions = replace_na(n_transitions, 0L),
    total_transitions = replace_na(total_transitions, 0L),
    share = dplyr::case_when(
      total_transitions == 0L ~ NA_real_,
      TRUE ~ n_transitions / total_transitions
    )
  )

delta_subtitle <- if (nrow(stay_rates) > 0) {
  paste0(
    "Stay-within cluster rates — ",
    paste(stay_rates$stay_label, collapse = "; ")
  )
} else {
  "No stay-within cluster transitions recorded"
}

p_trans_delta <- ggplot(transition_delta, aes(x = Cluster, y = Cluster_prev, fill = share)) +
  geom_tile(color = "white", linewidth = 0.8, na.rm = FALSE) +
  geom_text(
    aes(label = ifelse(is.na(share), "—", scales::percent(share, accuracy = 1))),
    color = "black",
    size = 4
  ) +
  scale_fill_gradient(
    low = "#F2F2F7",
    high = "#49006A",
    na.value = "white",
    labels = scales::percent_format(accuracy = 1),
    name = "% of transitions"
  ) +
  labs(
    x = "Cluster (year t)",
    y = "Cluster (year t-1)"
  )

ggsave(
  file.path(plot_stability, "annual_cluster_transition_matrix.png"),
  p_trans_delta,
  width = 8.5,
  height = 7,
  dpi = 300
)

# Cluster timeline - ensure all 4 clusters are represented
cluster_timeline <- cluster_wy %>%
  filter(!is.na(Cluster)) %>%
  mutate(Cluster = factor(as.character(Cluster), levels = cluster_levels)) %>%
  count(water_year, Cluster, .drop = FALSE) %>%
  group_by(water_year) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  replace_na(list(prop = 0))

p_timeline <- ggplot(cluster_timeline, aes(x = water_year, y = prop, fill = Cluster)) +
  geom_area(alpha = 0.8, color = "white", linewidth = 0.3) +
  scale_fill_cluster() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  labs(x = "Calendar year", y = "Proportion of site-solute combos", fill = "Cluster") +
  theme_clean()

ggsave(file.path(plot_stability, "annual_cluster_timeline.png"), p_timeline, width = 12, height = 6, dpi = 300)

cluster_timeline_site <- cluster_wy %>%
  filter(!is.na(Cluster)) %>%
  mutate(Cluster = factor(as.character(Cluster), levels = cluster_levels)) %>%
  count(Stream_Name, water_year, Cluster, .drop = FALSE) %>%
  group_by(Stream_Name, water_year) %>%
  tidyr::complete(
    Cluster = factor(cluster_levels, levels = cluster_levels),
    fill = list(n = 0)
  ) %>%
  mutate(
    total = sum(n),
    prop = ifelse(total > 0, n / total, 0)
  ) %>%
  ungroup() %>%
  select(-total) %>%
  apply_factor_orders()

p_timeline_site <- ggplot(cluster_timeline_site, aes(x = water_year, y = prop, fill = Cluster)) +
  geom_col(position = "stack", color = "white", linewidth = 0.2) +
  scale_fill_cluster() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  facet_wrap(~ Stream_Name, ncol = 3) +
  labs(
    x = "Calendar year",
    y = "Proportion of site-solute combos",
    fill = "Cluster"
  ) +
  theme_clean()

ggsave(file.path(plot_stability, "annual_cluster_timeline_by_site.png"), p_timeline_site, width = 12, height = 10, dpi = 300)

cluster_heatmap_input <- cluster_wy %>%
  filter(!is.na(Cluster)) %>%
  mutate(
    Cluster = factor(as.character(Cluster), levels = cluster_levels),
    water_year = as.integer(water_year)
  ) %>%
  apply_factor_orders()

if (nrow(cluster_heatmap_input) > 0) {
  year_range <- range(cluster_heatmap_input$water_year, na.rm = TRUE)
  cluster_heatmap <- cluster_heatmap_input %>%
    group_by(Stream_Name, chemical, water_year) %>%
    summarise(Cluster = dplyr::first(Cluster), .groups = "drop") %>%
    tidyr::complete(
      Stream_Name = site_order,
      chemical = solute_order,
      water_year = seq.int(year_range[1], year_range[2]),
      fill = list(Cluster = NA)
    ) %>%
    mutate(
      Stream_Name = factor(Stream_Name, levels = site_order),
      chemical = factor(chemical, levels = solute_order)
    )

  p_heatmap_site <- ggplot(cluster_heatmap, aes(x = water_year, y = chemical, fill = Cluster)) +
    geom_tile(color = "white", linewidth = 0.15, na.rm = FALSE) +
    scale_fill_manual(
      values = cluster_colors,
      na.value = "#FFFFFF",
      na.translate = FALSE,
      name = "Cluster"
    ) +
    scale_x_continuous(expand = c(0, 0), breaks = scales::pretty_breaks(n = 8)) +
    scale_y_discrete(limits = rev) +
    labs(
      x = "Calendar year",
      y = "Solute",
      fill = "Cluster",
    ) +
    facet_wrap(~ Stream_Name, ncol = 3) +
    theme_clean() +
    theme(
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(),
      legend.position = "bottom"
    )

  ggsave(
    file.path(plot_stability, "Cluster_calendar_year_heatmap_by_site.png"),
    p_heatmap_site,
    width = 16,
    height = 14,
    dpi = 300
  )
} else {
  warning("No cluster-year assignments available to plot calendar-year heatmap")
}


# =============================================================================
# SECTION 3: CLUSTER CONTROLS (Catchment + Hydrology)
# =============================================================================

# Join catchment with modal clusters
catchment_cluster <- cluster_modal %>%
  left_join(catchment, by = "Stream_Name") %>%
  rename(Cluster = Cluster_mode) %>%
  filter(!is.na(Cluster)) %>%
  ensure_cluster_levels("Cluster") %>%
  apply_factor_orders()

# Add RBI from site_data if available
if (!is.null(site_data) && "RBI" %in% names(site_data)) {
  rbi_data <- site_data %>%
    select(Stream_Name, solute, RBI) %>%
    rename(chemical = solute)

  catchment_cluster <- catchment_cluster %>%
    left_join(rbi_data, by = c("Stream_Name", "chemical"))
}

# Key catchment variables
catchment_vars <- c("Area_km2", "Elevation_mean_m", "Slope_mean", "Harvest", "Age",
                    "Lava1_per", "Ash_Per", "DR_Overall", "MTT_overall", "Fyw_overall")
catchment_labels <- c(
  "Area_km2" = "Basin area (km²)",
  "Elevation_mean_m" = "Mean elevation (m)",
  "Slope_mean" = "Mean slope (degrees)",
  "Aspec_Mean_deg" = "Mean aspect (degrees)",
  "Harvest" = "Harvest (%)",
  "Age" = "Stand age (years)",
  "Landslide_Total" = "Landslide density",
  "Lava1_per" = "Lava geology 1 (%)",
  "Lava2_per" = "Lava geology 2 (%)",
  "Ash_Per" = "Ash geology (%)",
  "Pyro_per" = "Pyroclastic geology (%)",
  "DR_Overall" = "Damping ratio",
  "MTT_overall" = "Mean transit time (years)",
  "Fyw_overall" = "Young water fraction",
  "RBI" = "Richards-Baker Flashiness Index"
)

# Filter to available variables
catchment_vars <- catchment_vars[catchment_vars %in% names(catchment_cluster)]

# Individual plots removed - only keeping full combined grid below

# FULL Combined controls panel (STATIC catchment characteristics only)
all_vars <- c("Area_km2", "Elevation_mean_m", "Slope_mean", "Aspec_Mean_deg", "Harvest", "Age",
              "Landslide_Total", "Lava1_per", "Lava2_per", "Ash_Per", "Pyro_per")
all_vars <- all_vars[all_vars %in% names(catchment_cluster)]

if (length(all_vars) > 0) {
  controls_all_long <- catchment_cluster %>%
    select(Stream_Name, chemical, Cluster, all_of(all_vars)) %>%
    pivot_longer(cols = all_of(all_vars), names_to = "metric", values_to = "value") %>%
    filter(!is.na(value)) %>%
    mutate(metric = factor(metric, levels = all_vars, labels = catchment_labels[all_vars]))
  
  # Calculate number of rows needed (3 columns)
  n_vars <- length(all_vars)
  n_rows <- ceiling(n_vars / 3)
  fig_height <- max(10, n_rows * 3)
  
  p_controls_all <- controls_all_long %>%
    filter(!is.na(Stream_Name)) %>%
    mutate(Stream_Name = forcats::fct_relevel(Stream_Name, intersect(site_order, unique(Stream_Name)))) %>%
    ggplot(aes(x = Cluster, y = value)) +
    geom_boxplot(alpha = 0.3, outlier.shape = NA, fill = "gray90", color = "gray40", linewidth = 0.5) +
    geom_jitter(aes(fill = Stream_Name, shape = Stream_Name), width = 0.15, alpha = 0.85, size = 2.5, stroke = 0.3, color = "black") +
    scale_fill_site(name = "Site",
            breaks = intersect(site_order, unique(controls_all_long$Stream_Name)),
            guide = "legend") +
    scale_shape_manual(values = site_shapes,
               name = "Site",
               breaks = intersect(site_order, unique(controls_all_long$Stream_Name)),
               guide = "legend") +
    facet_wrap(~ metric, scales = "free_y", ncol = 3) +
    labs(
      x = "Cluster",
      y = "Value"
    ) +
    theme_clean() +
    coord_cartesian(clip = "off") +
    theme(
      legend.position = c(0.75, 0.08),
      legend.box = "horizontal",
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11)
    ) +
    guides(
      color = guide_legend(override.aes = list(fill = cluster_colors, color = NA, shape = 22, size = 3), order = 2),
      fill = guide_legend(order = 1),
      shape = guide_legend(order = 1)
    )
  
  ggsave(file.path(plot_controls, "annual_controls_combined_ALL.png"), p_controls_all,
         width = 16, height = 14, dpi = 300)
}

# =============================================================================
# SECTION 4: CLUSTER TRANSITIONS (What drives cluster shifts?)
# =============================================================================

# Identify transition years
transitions <- cluster_wy %>%
  arrange(Stream_Name, chemical, water_year) %>%
  group_by(Stream_Name, chemical) %>%
  mutate(
    Cluster_prev = lag(Cluster),
    transition = Cluster != Cluster_prev & !is.na(Cluster_prev),
    transition_type = ifelse(transition, paste(Cluster_prev, "→", Cluster), "No change")
  ) %>%
  ungroup()

# Count transitions by year
transitions_by_year <- transitions %>%
  group_by(water_year) %>%
  summarise(
    n_transitions = sum(transition, na.rm = TRUE),
    n_total = n(),
    pct_transition = n_transitions / n_total * 100,
    .groups = "drop"
  )

p_trans_time <- ggplot(transitions_by_year, aes(x = water_year, y = pct_transition)) +
  geom_line(linewidth = 1, color = "#2166AC") +
  geom_point(size = 2, color = "#2166AC") +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.2, color = "#B2182B") +
  labs(x = "Water year", y = "% of site-solute combos that changed cluster") +
  theme_clean()

ggsave(file.path(plot_transitions, "annual_transition_rate_over_time.png"), p_trans_time, width = 12, height = 6, dpi = 300)

# If we have season boundaries, check if wet season timing relates to transitions
if (!is.null(season_bounds) && "wet_start_doy" %in% names(season_bounds)) {
  trans_season <- transitions_by_year %>%
    left_join(season_bounds %>% select(water_year, wet_start_doy, wet_duration), by = "water_year")
  
  if (nrow(trans_season) > 10 && sum(!is.na(trans_season$wet_start_doy)) > 10) {
    p_trans_wet <- ggplot(trans_season, aes(x = wet_start_doy, y = pct_transition)) +
      geom_point(size = 3, alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE, color = "#2166AC") +
      labs(x = "Wet season start (DOY)", y = "% transitions") +
      theme_clean()
    
    ggsave(file.path(plot_transitions, "annual_transitions_vs_wet_start.png"), p_trans_wet, width = 8, height = 6, dpi = 300)
  }
}

# Transitions by solute type
transitions_solute <- transitions %>%
  filter(transition) %>%
  mutate(solute_type = ifelse(chemical %in% geo_solutes, "Geogenic", "Biogenic")) %>%
  count(solute_type, transition_type) %>%
  group_by(solute_type) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

if (nrow(transitions_solute) > 0) {
  available_solute_type_colors <- solute_type_colors[names(solute_type_colors) %in% unique(transitions_solute$solute_type)]
  p_trans_solute <- ggplot(transitions_solute, aes(x = transition_type, y = pct, fill = solute_type)) +
    geom_col(position = "dodge", alpha = 0.8) +
    scale_fill_manual(values = available_solute_type_colors) +
    labs(x = "Transition type", y = "% of transitions", fill = "Solute type") +
    theme_clean() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(plot_transitions, "annual_transition_types_by_solute.png"), p_trans_solute, width = 12, height = 7, dpi = 300)
}


# =============================================================================
# SECTION 5: CLUSTER × CQ × SYNCHRONY RELATIONSHIPS
# =============================================================================

if (!is.null(mega)) {
  # Join mega data with cluster info
  # First check if mega already has Cluster column
  mega_for_join <- mega
  if ("Cluster" %in% names(mega_for_join)) {
    mega_for_join <- mega_for_join %>% select(-Cluster)
  }
  
  mega_cluster <- mega_for_join %>%
    left_join(cluster_modal %>% select(Stream_Name, chemical, Cluster_mode),
              by = c("Stream_Name", "solute" = "chemical")) %>%
    rename(Cluster = Cluster_mode) %>%
    { df <- ensure_cluster_levels(., "Cluster"); apply_factor_orders(df) } %>%
    filter(!is.na(Cluster))
  
  # CQ slope by cluster
  if ("cq_slope.x" %in% names(mega_cluster)) {
    p_cq_cluster <- ggplot(mega_cluster, aes(x = Cluster, y = cq_slope.x, fill = Cluster)) +
      geom_boxplot(alpha = 0.3, outlier.alpha = 0.1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      scale_fill_cluster() +
      labs(x = "Cluster", y = "CQ slope") +
      theme_clean() +
      theme(legend.position = "none")
    
    ggsave(file.path(plot_relationships, "annual_CQ_slope_by_cluster.png"), p_cq_cluster, width = 8, height = 6, dpi = 300)
    
    # CQ slope by cluster and season
    if ("hydrologic_season" %in% names(mega_cluster)) {
      p_cq_season <- ggplot(mega_cluster, aes(x = Cluster, y = cq_slope.x, fill = Cluster)) +
        geom_boxplot(alpha = 0.3, outlier.alpha = 0.1) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        scale_fill_cluster() +
        facet_wrap(~ hydrologic_season) +
        labs(x = "Cluster", y = "CQ slope") +
        theme_clean() +
        theme(legend.position = "none")
      
      ggsave(file.path(plot_relationships, "annual_CQ_slope_by_cluster_season.png"), p_cq_season, width = 12, height = 6, dpi = 300)
    }
  }
  
  # Synchrony by cluster
  if ("cq_sync" %in% names(mega_cluster)) {
    # Note: cq_sync originates from the Wymore quadrant flag (windows landing in Q1 or Q3)
    sync_cluster <- mega_cluster %>%
      group_by(Cluster) %>%
      summarise(
        pct_sync = mean(cq_sync == "sync", na.rm = TRUE) * 100,
        n = n(),
        .groups = "drop"
      )
    
    p_sync_cluster <- ggplot(sync_cluster, aes(x = Cluster, y = pct_sync, fill = Cluster)) +
      geom_col(alpha = 0.8) +
      geom_text(aes(label = sprintf("%.1f%%", pct_sync)), vjust = -0.5, size = 4) +
      scale_fill_cluster() +
      labs(
        x = "Cluster",
        y = "% windows synchronous (Wymore quadrant)",
      ) +
      ylim(0, 100) +
      theme_clean() +
      theme(legend.position = "none")
    
    ggsave(file.path(plot_relationships, "annual_synchrony_by_cluster.png"), p_sync_cluster, width = 8, height = 6, dpi = 300)
  }
  
  # Storage by cluster (Q-dS vs WB-dS) using patchwork
  storage_plots <- list()
  storage_labels <- character()
  
  if ("Q_dS_range_mm" %in% names(mega_cluster)) {
    q_ds_data <- mega_cluster %>% filter(!is.na(Q_dS_range_mm), Q_dS_range_mm > 0)
    if (nrow(q_ds_data) > 0) {
      p_q_ds <- ggplot(q_ds_data, aes(x = Cluster, y = Q_dS_range_mm, fill = Cluster)) +
      geom_boxplot(alpha = 0.3, outlier.alpha = 0.1) +
      scale_fill_cluster() +
      scale_y_log10() +
      labs(x = "Cluster", y = "Storage range (mm)") +
      theme_clean() +
      theme(legend.position = "none")
      storage_plots <- c(storage_plots, list(p_q_ds))
      storage_labels <- c(storage_labels, "Q–dS metric")
    }
  }
  
  if ("WB_dS_range_mm" %in% names(mega_cluster)) {
    wb_ds_data <- mega_cluster %>% filter(!is.na(WB_dS_range_mm), WB_dS_range_mm > 0)
    if (nrow(wb_ds_data) > 0) {
      p_wb_ds <- ggplot(wb_ds_data, aes(x = Cluster, y = WB_dS_range_mm, fill = Cluster)) +
      geom_boxplot(alpha = 0.3, outlier.alpha = 0.1) +
      scale_fill_cluster() +
      scale_y_log10() +
      labs(x = "Cluster", y = "Storage range (mm)") +
      theme_clean() +
      theme(legend.position = "none")
      storage_plots <- c(storage_plots, list(p_wb_ds))
      storage_labels <- c(storage_labels, "WB–dS metric")
    }
  }
  
  if (length(storage_plots) > 0) {
    panel_cols <- ifelse(length(storage_plots) > 1, 2, 1)
    storage_subtitle <- if (length(storage_plots) == 2) {
      "Left: Q–dS metric, Right: WB–dS metric"
    } else {
      paste("Only", storage_labels[1], "available")
    }
    combined_storage <- patchwork::wrap_plots(plotlist = storage_plots, ncol = panel_cols) +
      patchwork::plot_annotation(
      )
    plot_width <- ifelse(length(storage_plots) > 1, 12, 8)
    ggsave(
      file.path(plot_relationships, "Storage_by_cluster_Q_vs_WB.png"),
      combined_storage,
      width = plot_width,
      height = 6,
      dpi = 300
    )
  }
  
  # COMPREHENSIVE CORRELATION MATRIX
  message("  Creating correlation matrices...\n")
  
  # Select key variables for correlation
  corr_vars <- c("cq_slope.x", "Q_dS_range_mm", "RBI", "DR_Overall", "MTT_final", "FYw_final")
  corr_vars <- corr_vars[corr_vars %in% names(mega_cluster)]
  
  if (length(corr_vars) >= 3) {
    # Overall correlation
    corr_data <- mega_cluster %>%
      select(all_of(corr_vars)) %>%
      drop_na()
    
    if (nrow(corr_data) > 50) {
      corr_mat <- cor(corr_data, use = "pairwise.complete.obs")
      
      png(file.path(plot_relationships, "Correlation_matrix_overall.png"), width = 10, height = 10, units = "in", res = 300)
      corrplot(corr_mat, method = "color", type = "upper",
               addCoef.col = "black", number.cex = 0.8,
               tl.col = "black", tl.srt = 45)
      dev.off()
      
      # By cluster
      for (cl in levels(mega_cluster$Cluster)) {
        cl_data <- mega_cluster %>%
          filter(Cluster == cl) %>%
          select(all_of(corr_vars)) %>%
          drop_na()
        
        if (nrow(cl_data) > 30) {
          cl_mat <- cor(cl_data, use = "pairwise.complete.obs")
          
          png(file.path(plot_relationships, paste0("Correlation_matrix_cluster", cl, ".png")),
              width = 10, height = 10, units = "in", res = 300)
          corrplot(cl_mat, method = "color", type = "upper",
                   addCoef.col = "black", number.cex = 0.7,
                   tl.col = "black", tl.srt = 45)
          dev.off()
        }
      }
      
      # By season
      if ("hydrologic_season" %in% names(mega_cluster)) {
        for (season in unique(mega_cluster$hydrologic_season)) {
          s_data <- mega_cluster %>%
            filter(hydrologic_season == season) %>%
            select(all_of(corr_vars)) %>%
            drop_na()
          
          if (nrow(s_data) > 30) {
            s_mat <- cor(s_data, use = "pairwise.complete.obs")
            
            png(file.path(plot_relationships, paste0("Correlation_matrix_", season, ".png")),
                width = 10, height = 10, units = "in", res = 300)
            corrplot(s_mat, method = "color", type = "upper",
                     addCoef.col = "black", number.cex = 0.7,
                     tl.col = "black", tl.srt = 45)
            dev.off()
          }
        }
      }
    }
  }
}


# =============================================================================
# SECTION 6: CLUSTER PATTERNS BY SITE
# =============================================================================

# Cluster distribution by site
cluster_by_site <- cluster_modal %>%
  rename(Cluster = Cluster_mode) %>%
  count(Stream_Name, Cluster) %>%
  mutate(Cluster = factor(as.character(Cluster), levels = cluster_levels)) %>%
  tidyr::complete(
    Stream_Name,
    Cluster = factor(cluster_levels, levels = cluster_levels),
    fill = list(n = 0)
  ) %>%
  apply_factor_orders() %>%
  arrange(Stream_Name, desc(Cluster)) %>%
  group_by(Stream_Name) %>%
  mutate(
    total = sum(n),
    prop = ifelse(total > 0, n / total, 0),
    # Format percentage label
    label = ifelse(prop > 0.05, paste0(round(prop * 100), "%"), "")
  ) %>%
  ungroup() %>%
  select(-total)

# Proportional bar
p_site2 <- ggplot(cluster_by_site, aes(x = Stream_Name, y = prop, fill = Cluster)) +
  geom_col(position = "stack", alpha = 0.8) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3.5, color = "white") +
  scale_fill_cluster() +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Site", y = "Proportion", fill = "Cluster") +
  theme_clean() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plot_bysite, "annual_cluster_proportions_by_site.png"), p_site2, width = 10, height = 7, dpi = 300)

# Solute × Site matrix (tile plot)
p_matrix <- ggplot(cluster_modal, aes(x = chemical, y = Stream_Name, fill = Cluster_mode)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Cluster_mode), color = "white", size = 4) +
  scale_fill_cluster() +
  labs(x = "Solute", y = "Site", fill = "Cluster") +
  theme_clean() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plot_bysite, "annual_solute_by_site_cluster_matrix.png"), p_matrix, width = 10, height = 8, dpi = 300)

# Cluster by site HEATMAP (count-based)
p_heatmap <- ggplot(cluster_by_site, aes(x = Cluster, y = Stream_Name, fill = n)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient(
    low = "#F5F4FB",
    high = "#2D1E70",
    name = "Modal count",
    breaks = scales::pretty_breaks(n = 5)
  ) +
  geom_text(
    aes(label = n, colour = n == 0),
    size = 4
  ) +
  scale_colour_manual(values = c(`TRUE` = "grey40", `FALSE` = "white"), guide = "none") +
  labs(
    x = "Cluster",
    y = "Site",
    fill = "Modal count"
  )

ggsave(file.path(plot_bysite, "annual_cluster_by_site_heatmap.png"), p_heatmap, width = 10, height = 6, dpi = 300)

# Combined cluster characterization for ANNUAL (using modal clusters)
# Need to get monthly data for modal cluster assignments
modal_long <- raw %>%
  select(-Cluster) %>%  # Remove climatological cluster to avoid duplicate
  inner_join(cluster_modal %>% select(Stream_Name, chemical, Cluster_mode),
             by = c("Stream_Name", "chemical")) %>%
  pivot_longer(Jan:Dec, names_to = "month", values_to = "z") %>%
  mutate(
    month = factor(month, levels = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")),
    month_num = as.numeric(month),
    chemical = factor(chemical, levels = solute_order)
  ) %>%
  rename(Cluster = Cluster_mode)

cluster_panels_annual <- list()

for (cl in levels(modal_long$Cluster)) {
  cl_data <- modal_long %>% filter(Cluster == cl)
  cl_pct <- round(100 * nrow(cl_data %>% distinct(Stream_Name, chemical)) /
                        nrow(cluster_modal %>% distinct(Stream_Name, chemical)), 1)

  # Left panel: temporal pattern (spaghetti plot colored by solute)
  p_left <- ggplot(cl_data, aes(x = month_num, y = z,
                                group = interaction(Stream_Name, chemical),
                                color = chemical)) +
    geom_line(alpha = 0.45, linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
    scale_x_continuous(breaks = 1:12, labels = month_labels) +
    scale_color_solute() +
    labs(
      x = "Month",
      y = "Normalized Concentration",
      color = "Solute"
    ) +
    theme_clean(base_size = 11) +
    theme(legend.position = "none")

  # Right panel: solute composition (count bar chart)
  solute_counts <- cl_data %>%
    distinct(Stream_Name, chemical) %>%
    count(chemical) %>%
    complete(chemical = factor(solute_order, levels = solute_order), fill = list(n = 0)) %>%
    mutate(chemical = factor(chemical, levels = solute_order))

  p_right <- ggplot(solute_counts, aes(x = chemical, y = n, fill = chemical)) +
    geom_col(alpha = 0.8) +
    scale_fill_solute() +
    labs(x = "Solute", y = "Count") +
    theme_clean(base_size = 11) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
    )

  # Combine left and right with adjusted widths (spaghetti gets more space)
  cluster_panels_annual[[cl]] <- p_left + p_right + plot_layout(widths = c(2, 1))
}

# Combine all clusters vertically with shared legend
combined_annual <- wrap_plots(cluster_panels_annual, ncol = 1) +
  plot_layout(guides = "collect")

ggsave(file.path(plot_annual, "annual_cluster_characterization_combined.png"),
       combined_annual, width = 12, height = 16, dpi = 300)

# =============================================================================
# CLEANUP AND SUMMARY
# =============================================================================

# Delete old folders if they exist
old_folders <- c("heatmaps", "characterization", "hydrology")
for (folder in old_folders) {
  old_dir <- file.path(plot_base, folder)
  if (dir.exists(old_dir)) {
    unlink(old_dir, recursive = TRUE)
    message(paste0("  [✓] Removed old ", folder, " folder\n"))
  }
}
