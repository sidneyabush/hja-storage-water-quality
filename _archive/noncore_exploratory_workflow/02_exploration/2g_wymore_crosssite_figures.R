# =============================================================================
# Wymore Cross-Site Synchrony Quadrant Plots
# 
# Visualize cross-site CQ slope synchrony using quadrant classification
# Q1: Both sites mobilizing (+ slopes) = sync
# Q2: Site 1 mobilizing, Site 2 diluting = async
# Q3: Both sites diluting (- slopes) = sync
# Q4: Site 1 diluting, Site 2 mobilizing = async
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

rm(list = ls())

if (!exists("%||%", inherits = TRUE)) {
  `%||%` <- function(x, y) if (!is.null(x)) x else y
}

# Source helpers from repo root
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
try(source(file.path(repo_dir, "00_helpers", "plot_prefs.R")), silent = TRUE)
try(source(file.path(repo_dir, "00_helpers", "workflow_config.R")), silent = TRUE)
try(source(file.path(repo_dir, "00_helpers", "plot_theme_set.R")), silent = TRUE)

# Paths
base_dir    <- "/Users/sidneybush/Library/CloudStorage/Box-Box"
project_dir <- file.path(base_dir, "Sidney_Bush", "HJA_Water_Quality")
out_dir     <- file.path(project_dir, "outputs")
fig_dir     <- file.path(project_dir, "exploratory_plots", "02_exploration", "2g_wymore")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Load Wymore cross-site sync data
wymore <- readr::read_csv(file.path(out_dir, "HJA_wymore_crosssite_sync.csv"),
                          show_col_types = FALSE) %>%
  filter(!is.na(quadrant))

message("Loaded ", nrow(wymore), " site-pair × window observations")
message("Quadrant distribution:")
print(table(wymore$quadrant))

# Order factors for consistent faceting
pair_levels <- combn(site_order, 2, function(x) paste(x[1], "vs", x[2]))
pair_levels <- unname(pair_levels)

cluster_wy <- readr::read_csv(
  file.path(out_dir, "ClusterStreams_allSolutes_byWaterYear.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    Cluster_climRef = factor(as.character(Cluster_climRef), levels = cluster_levels)
  ) %>%
  select(Stream_Name, chemical, water_year, Cluster_climRef)

wymore <- wymore %>%
  filter(Stream1 %in% site_order, Stream2 %in% site_order) %>%
  mutate(
    solute   = factor(solute, levels = solute_order),
    Stream1  = factor(Stream1, levels = site_order),
    Stream2  = factor(Stream2, levels = site_order),
    pair_label = factor(paste(Stream1, "vs", Stream2), levels = pair_levels)
  ) %>%
  left_join(
    cluster_wy %>% rename(cluster1 = Cluster_climRef),
    by = c("Stream1" = "Stream_Name", "solute" = "chemical", "water_year")
  ) %>%
  left_join(
    cluster_wy %>% rename(cluster2 = Cluster_climRef),
    by = c("Stream2" = "Stream_Name", "solute" = "chemical", "water_year")
  ) %>%
  drop_na(pair_label)

wymore <- add_solute_type(wymore, solute_col = "solute", three_way = TRUE)

# Ensure solute is a factor before droplevels
if (!is.factor(wymore$solute)) {
  wymore$solute <- factor(wymore$solute, levels = solute_order)
}
solutes_for_plots <- levels(droplevels(wymore$solute))
if (length(solutes_for_plots) == 0) {
  stop("No overlapping solutes found for Wymore plots.")
}

# Theme
theme_wymore <- theme_bw(base_size = 16) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray95", colour = NA),
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 13)
  )

# Colors
quadrant_colors <- c(
  "Q1" = "#2E86AB",  # Blue - sync (both mobilizing)
  "Q2" = "#F4A261",  # Orange - async
  "Q3" = "#2A9D8F",  # Teal - sync (both diluting)
  "Q4" = "#E76F51"   # Red-orange - async
)

facet_pair_plot <- function(solute_name) {
  solute_df <- wymore %>% filter(solute == solute_name) %>%
    mutate(pair_label = forcats::fct_drop(pair_label))
  if (nrow(solute_df) == 0) return(NULL)
  slope_range <- max(abs(c(solute_df$slope1, solute_df$slope2)), na.rm = TRUE)
  if (!is.finite(slope_range) || slope_range == 0) slope_range <- 1
  lim <- slope_range * 1.1 * c(-1, 1)
  inset_span <- diff(lim) * 0.25
  inset_xmin <- lim[2] - inset_span
  inset_ymin <- lim[2] - inset_span
  quadrant_levels <- c("Q1", "Q2", "Q3", "Q4")

  inset_bg <- solute_df %>%
    distinct(pair_label) %>%
    mutate(
      xmin = inset_xmin,
      xmax = lim[2],
      ymin = inset_ymin,
      ymax = lim[2]
    )

  counts <- solute_df %>%
    count(pair_label, quadrant, name = "n") %>%
    complete(pair_label, quadrant = factor(quadrant_levels, levels = quadrant_levels), fill = list(n = 0)) %>%
    group_by(pair_label) %>%
    mutate(
      n_total = sum(n),
      n_max = max(n),
      n_max = ifelse(is.finite(n_max) & n_max > 0, n_max, 1),
      idx = as.integer(factor(quadrant, levels = quadrant_levels)),
      bar_width = inset_span / length(quadrant_levels),
      xleft = inset_xmin + (idx - 1) * bar_width + bar_width * 0.12,
      xright = inset_xmin + idx * bar_width - bar_width * 0.12,
      ybottom = inset_ymin,
      ytop = inset_ymin + (n / n_max) * inset_span,
      label_x = (xleft + xright) / 2,
      label_y = ytop + inset_span * 0.05,
      label_y = pmin(label_y, lim[2])
    ) %>%
    ungroup()

  ggplot(
    solute_df,
    aes(
      x = slope1,
      y = slope2,
      color = quadrant,
      shape = if_else(is_outlet_pair, "Outlet pairs (GSLOOK)", "All other pairs")
    )
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray70", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray70", linewidth = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "gray55") +
    geom_rect(
      data = inset_bg,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "white",
      alpha = 0.85,
      color = NA
    ) +
    geom_rect(
      data = counts,
      aes(xmin = xleft, xmax = xright, ymin = ybottom, ymax = ytop, fill = quadrant),
      inherit.aes = FALSE,
      color = NA,
      alpha = 0.9
    ) +
    geom_text(
      data = counts,
      aes(x = label_x, y = label_y, label = n),
      inherit.aes = FALSE,
      size = 3.5,
      fontface = "bold",
      color = "#1F2933"
    ) +
    geom_point(alpha = 0.4, size = 1.4) +
    scale_color_manual(
      values = quadrant_colors,
      breaks = quadrant_levels,
      labels = c("Q1: +/+", "Q2: +/-", "Q3: -/-", "Q4: -/+"),
      name = "Quadrant"
    ) +
    scale_fill_manual(
      values = quadrant_colors,
      breaks = quadrant_levels,
      guide = "none"
    ) +
    scale_shape_manual(
      values = c("All other pairs" = 16, "Outlet pairs (GSLOOK)" = 17),
      name = "Pair type"
    ) +
    coord_equal(xlim = lim, ylim = lim, clip = "off") +
    facet_wrap(~ pair_label, ncol = 4, drop = TRUE) +
    labs(
      x = "Site 1 CQ slope",
      y = "Site 2 CQ slope",
      strip.text = element_text(size = 12),
      plot.margin = margin(10, 20, 10, 10)
    )
}

for (solute_name in solutes_for_plots) {
  pair_plot <- facet_pair_plot(solute_name)
  if (!is.null(pair_plot)) {
    ggsave(
      file.path(fig_dir, paste0("wymore_quadrant_pairs_", solute_name, ".png")),
      plot = pair_plot,
      width = 16,
      height = 18,
      dpi = 300
    )
  }
}

# --- Plot 2: Synchrony fraction by solute and season ---
sync_by_solute <- wymore %>%
  group_by(solute, hydrologic_season) %>%
  summarise(
    n = n(),
    sync_fraction = mean(sync, na.rm = TRUE),
    .groups = "drop"
  )

p_solute <- sync_by_solute %>%
  mutate(solute = forcats::fct_drop(solute)) %>%
  ggplot(aes(x = solute, y = sync_fraction, fill = hydrologic_season)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.75) +
  geom_text(
    aes(label = scales::percent(sync_fraction, accuracy = 1)),
    position = position_dodge(width = 0.85),
    vjust = -0.6,
    size = 4.5,
    fontface = "bold"
  ) +
  geom_text(
    aes(label = paste0("n = ", scales::comma(n))),
    position = position_dodge(width = 0.85),
    vjust = -2.1,
    size = 3.8
  ) +
  scale_fill_manual(values = season_colors) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1.05),
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    x = "Solute",
    y = "Cross-Site Synchrony Fraction",
    fill = "Season",

# --- Plot 3: Synchrony fraction by cluster and season ---
cluster_long <- wymore %>%
  select(hydrologic_season, sync, cluster1, cluster2) %>%
  pivot_longer(c(cluster1, cluster2), values_to = "cluster", names_to = "cluster_role") %>%
  filter(!is.na(cluster)) %>%
  mutate(cluster = factor(cluster, levels = cluster_levels))

sync_by_cluster <- cluster_long %>%
  group_by(cluster, hydrologic_season) %>%
  summarise(n = n(), sync_fraction = mean(sync, na.rm = TRUE), .groups = "drop")

p_cluster <- ggplot(sync_by_cluster, aes(x = cluster, y = sync_fraction, fill = hydrologic_season)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_text(
    aes(label = scales::percent(sync_fraction, accuracy = 1)),
    position = position_dodge(width = 0.75),
    vjust = -0.7,
    size = 4.5,
    fontface = "bold"
  ) +
  geom_text(
    aes(label = paste0("n = ", scales::comma(n))),
    position = position_dodge(width = 0.75),
    vjust = -2.0,
    size = 3.6
  ) +
  scale_fill_manual(values = season_colors) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1.05),
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    x = "Cluster (calendar-year assignment)",
    y = "Cross-Site Synchrony Fraction",
    fill = "Season",

p_type <- ggplot(sync_by_type, aes(x = solute_type, y = sync_fraction, fill = hydrologic_season)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_text(
    aes(label = scales::percent(sync_fraction, accuracy = 1)),
    position = position_dodge(width = 0.75),
    vjust = -0.7,
    size = 4.5,
    fontface = "bold"
  ) +
  geom_text(
    aes(label = paste0("n = ", scales::comma(n))),
    position = position_dodge(width = 0.75),
    vjust = -2.0,
    size = 3.6
  ) +
  scale_fill_manual(values = season_colors) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1.05),
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    x = "Solute group",
    y = "Cross-Site Synchrony Fraction",
    fill = "Season",

p_quadrant <- ggplot(quadrant_by_season, aes(x = hydrologic_season, y = prop, fill = quadrant)) +
  geom_col(position = "fill") +
  geom_text(
    aes(label = scales::percent(prop, accuracy = 1)),
    position = position_fill(vjust = 0.5),
    color = "white",
    fontface = "bold",
    size = 4
  ) +
  scale_fill_manual(
    values = quadrant_colors,
    labels = c(
      "Q1" = "Q1: Both mobilizing",
      "Q2" = "Q2: Mixed",
      "Q3" = "Q3: Both diluting",
      "Q4" = "Q4: Mixed"
    )
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Season",
    y = "Proportion of site pairs",
    fill = "Quadrant",
    n = n(),
    .groups = "drop"
  )

p_outlet <- sync_outlet %>%
  mutate(solute = factor(solute, levels = solute_order)) %>%
  ggplot(aes(x = solute, y = sync_fraction, fill = pair_type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(
    aes(label = scales::percent(sync_fraction, accuracy = 1)),
    position = position_dodge(width = 0.8),
    vjust = -0.7,
    size = 4.2,
    fontface = "bold"
  ) +
  geom_text(
    aes(label = paste0("n = ", scales::comma(n))),
    position = position_dodge(width = 0.8),
    vjust = -2.0,
    size = 3.4
  ) +
  scale_fill_manual(values = c("Outlet pairs\n(includes GSLOOK)" = "#2E86AB",
                               "All other pairs" = "#E07A5F"), name = "Pair type") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1.05),
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    x = "Solute",
    y = "Cross-site synchrony fraction",

summary_layout <- (p_quadrant | p_solute) / (p_cluster | p_type) /
  p_outlet +
  plot_annotation(

ggsave(
  file.path(fig_dir, "wymore_summary_panels.png"),
  plot = summary_layout,
  width = 16,
  height = 15,
  dpi = 300
)
ggsave(
  file.path(fig_dir, "wymore_sync_by_solute_season.png"),
  plot = p_solute,
  width = 12,
  height = 7,
  dpi = 300
)
ggsave(
  file.path(fig_dir, "wymore_sync_by_cluster_season.png"),
  plot = p_cluster,
  width = 10,
  height = 7,
  dpi = 300
)
ggsave(
  file.path(fig_dir, "wymore_sync_by_solute_type_season.png"),
  plot = p_type,
  width = 10,
  height = 7,
  dpi = 300
)
ggsave(
  file.path(fig_dir, "wymore_outlet_comparison.png"),
  plot = p_outlet,
  width = 12,
  height = 7,
  dpi = 300
)

message("\nSaved Wymore cross-site synchrony figures to: ", fig_dir)

# Print summary
message("\n=== Summary Statistics ===")
message("Overall cross-site synchrony: ", 
        round(mean(wymore$sync, na.rm = TRUE) * 100, 1), "%")
message("Wet season sync: ", 
        round(mean(wymore$sync[wymore$hydrologic_season == "Wet"], na.rm = TRUE) * 100, 1), "%")
message("Dry season sync: ",
        round(mean(wymore$sync[wymore$hydrologic_season == "Dry"], na.rm = TRUE) * 100, 1), "%")
message("Outlet pairs sync: ",
        round(mean(wymore$sync[wymore$is_outlet_pair], na.rm = TRUE) * 100, 1), "%")
