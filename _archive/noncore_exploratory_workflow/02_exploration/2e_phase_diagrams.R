# =============================================================================
# 6b_phase_diagrams.R
# =============================================================================
# PHASE SPACE ANALYSIS: Where do sites/seasons operate in storage–chemistry space?
#
# RESEARCH QUESTION:
#   Do sites move through different phases of the storage–CQ relationship?
#   How does this vary seasonally vs annually, and across sites?
#
# APPROACH:
#   Create phase diagrams plotting:
#     X-axis: storage-paper metrics or storage-metric ordination scores
#     Y-axis: CQ slope
#   These reveal whether sites are mobilizing, diluting, or chemostatic
#   and how storage varies across those regimes.
#
# SCALES:
#   1. Seasonal phase diagram (main): shows within-year variability
#   2. Site-mean phase diagram (supplemental): long-term patterns
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
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

# =============================================================================
# SETUP
# =============================================================================

base_dir    <- "/Users/sidneybush/Library/CloudStorage/Box-Box"
project_dir <- file.path(base_dir, "Sidney_Bush", "HJA_Water_Quality")
out_dir     <- file.path(project_dir, "outputs")

fig_dir     <- file.path(project_dir, "exploratory_plots", "02_exploration", "2e_phase_diagrams")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, filename, width = 10, height = 7) {
  ggplot2::ggsave(
    file.path(fig_dir, filename),
    plot = p, width = width, height = height, dpi = 300
  )
}

theme_phase <- theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# =============================================================================
# LOAD DATA
# =============================================================================

message("Loading data...")

seasonal <- readr::read_csv(
  file.path(out_dir, "HJA_clean_seasonal.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    hydrologic_season = factor(hydrologic_season, levels = c("Wet", "Dry")),
    Stream_Name = as.factor(Stream_Name),
    solute = as.factor(solute)
  ) %>%
  apply_factor_orders()
if (exists("ALL_SOLUTES")) {
  seasonal <- seasonal %>% mutate(solute = forcats::fct_relevel(solute, ALL_SOLUTES))
}
if (exists("site_order")) {
  seasonal <- seasonal %>% mutate(Stream_Name = forcats::fct_relevel(Stream_Name, site_order))
}
if (!"solute_group" %in% names(seasonal) && exists("get_solute_group")) {
  seasonal <- seasonal %>% mutate(solute_group = get_solute_group(as.character(solute)))
}

site_means <- readr::read_csv(
  file.path(out_dir, "HJA_clean_site_means.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    Stream_Name = as.factor(Stream_Name),
    solute = as.factor(solute)
  ) %>%
  apply_factor_orders()
if (exists("ALL_SOLUTES")) {
  site_means <- site_means %>% mutate(solute = forcats::fct_relevel(solute, ALL_SOLUTES))
}
if (exists("site_order")) {
  site_means <- site_means %>% mutate(Stream_Name = forcats::fct_relevel(Stream_Name, site_order))
}
if (!"solute_group" %in% names(site_means) && exists("get_solute_group")) {
  site_means <- site_means %>% mutate(solute_group = get_solute_group(as.character(solute)))
}

primary_storage_metric <- intersect(PRIMARY_STORAGE_METRIC, names(seasonal))
label_metric <- if (length(primary_storage_metric) == 1) primary_storage_metric else PRIMARY_STORAGE_METRIC
primary_label <- get_label(label_metric)

legend_shapes <- c(
  "Wet" = 19,
  "Dry" = 17,
  "Outlet (GSLOOK)" = OUTLET_SHAPE_TRIANGLE
)

# =============================================================================
# 1. SEASONAL PHASE DIAGRAM (main)
# =============================================================================
# Visualize how 90-day windows vary in storage–CQ space

message("\n=== SEASONAL PHASE DIAGRAM ===\n")

if (length(primary_storage_metric) == 1 &&
    all(c(primary_storage_metric, "cq_slope") %in% names(seasonal))) {
  metric <- primary_storage_metric[[1]]
  df_phase <- seasonal %>%
    filter(!is.na(.data[[metric]]), !is.na(cq_slope))
  df_phase <- flag_outlet_stream(df_phase)
  
  if (nrow(df_phase) > 0) {
    non_outlet <- df_phase %>% filter(!is_outlet)
    outlet_pts <- df_phase %>% filter(is_outlet)
    
    p_seasonal <- ggplot(non_outlet, aes(
      x = .data[[metric]],
      y = cq_slope,
      color = Stream_Name,
      shape = hydrologic_season
    )) +
      # Reference lines for CQ regimes
      geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey50", alpha = 0.5, linewidth = 0.5) +
      geom_hline(yintercept = -0.1, linetype = "dashed", color = "grey50", alpha = 0.5, linewidth = 0.5) +
      geom_hline(yintercept = 0, linetype = "solid", color = "grey70", alpha = 0.3, linewidth = 0.3) +
      # Zones
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.1, ymax = Inf,
               fill = "green", alpha = 0.05) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = -0.1, ymax = 0.1,
               fill = "grey", alpha = 0.05) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -0.1,
               fill = "blue", alpha = 0.05) +
      # Regime labels
      annotate("text", x = Inf, y = 0.3, label = "Mobilizing",
               hjust = 1.1, vjust = 1, color = "darkgreen", size = 3, alpha = 0.6) +
      annotate("text", x = Inf, y = 0, label = "Chemostatic",
               hjust = 1.1, vjust = 1, color = "grey40", size = 3, alpha = 0.6) +
      annotate("text", x = Inf, y = -0.3, label = "Diluting",
               hjust = 1.1, vjust = 1, color = "darkblue", size = 3, alpha = 0.6) +
      # Data
      geom_point(alpha = 0.75, size = 2.2) +
      scale_color_site(name = "Site") +
      scale_shape_manual(values = legend_shapes, name = "Season / Outlet") +
      guides(
        color = guide_legend(order = 1, override.aes = list(size = 3, alpha = 1)),
        shape = guide_legend(
          order = 2,
          override.aes = list(
            color = c("grey35", "grey35", "#111827"),
            size = c(2.4, 2.4, 2.8)
          )
        )
      ) +
      labs(
        x = paste0("Seasonal ", primary_label),
        y = "Seasonal C–Q slope",
          aes(x = .data[[metric]], y = cq_slope, shape = "Outlet (GSLOOK)"),
          inherit.aes = FALSE,
          size = 2.6,
          color = "#111827",
          alpha = 0.95,
          stroke = 0.6,
          show.legend = TRUE
        )
    }

    save_plot(p_seasonal, "01_seasonal_phase_diagram.png", width = 11, height = 8)

    if ("solute" %in% names(df_phase)) {
      save_plot(p_seasonal + facet_wrap(~ solute, ncol = 3),
                "01_seasonal_phase_diagram_by_solute.png", width = 12, height = 8)
      save_plot(p_seasonal + facet_wrap(~ Stream_Name, ncol = 3),
                "01_seasonal_phase_diagram_by_site.png", width = 12, height = 9)
    }

    message("Created seasonal phase diagrams with ", nrow(df_phase), " observations")
  }
} else {
  warning("Primary storage metric not found in seasonal dataset")
}

# =============================================================================
# 2. SITE-MEAN PHASE DIAGRAM (supplemental)
# =============================================================================
# Longer-term view: where do sites sit in the phase space on average?

message("\n=== SITE-MEAN PHASE DIAGRAM ===\n")

if (length(primary_storage_metric) == 1 &&
    all(c(primary_storage_metric, "cq_slope") %in% names(site_means))) {
  metric <- primary_storage_metric[[1]]
  df_phase_site <- site_means %>%
    filter(!is.na(.data[[metric]]), !is.na(cq_slope))
  df_phase_site <- flag_outlet_stream(df_phase_site)
  
  if (nrow(df_phase_site) > 0) {
    non_outlet <- df_phase_site %>% filter(!is_outlet)
    outlet_pts <- df_phase_site %>% filter(is_outlet)
    
    p_site_mean <- non_outlet %>%
      ggplot(aes(x = .data[[metric]], y = cq_slope, color = Stream_Name)) +
      # Reference lines
      geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey50", alpha = 0.5, linewidth = 0.5) +
      geom_hline(yintercept = -0.1, linetype = "dashed", color = "grey50", alpha = 0.5, linewidth = 0.5) +
      geom_hline(yintercept = 0, linetype = "solid", color = "grey70", alpha = 0.3, linewidth = 0.3) +
      # Zones
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.1, ymax = Inf,
               fill = "green", alpha = 0.05) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = -0.1, ymax = 0.1,
               fill = "grey", alpha = 0.05) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -0.1,
               fill = "blue", alpha = 0.05) +
      # Data
      geom_point(alpha = 0.85, size = 3, show.legend = TRUE) +
      scale_color_site(name = "Site") +
      scale_shape_manual(
        values = c("Outlet (GSLOOK)" = OUTLET_SHAPE_TRIANGLE),
        name = NULL
      ) +
      guides(
        color = guide_legend(order = 1, override.aes = list(size = 3.2, alpha = 1)),
        shape = guide_legend(order = 2, override.aes = list(color = "black", size = 3.2))
      ) +
      labs(
        x = paste0("Long-term ", primary_label),
        y = "Site-mean C–Q slope",
        aes(x = .data[[metric]], y = cq_slope, shape = "Outlet (GSLOOK)"),
        inherit.aes = FALSE,
        size = 3.2,
        color = "black",
        alpha = 0.95,
        stroke = 0.7,
        show.legend = TRUE
      )
    }

    save_plot(p_site_mean, "02_site_mean_phase_diagram.png", width = 11, height = 8)

    save_plot(
      p_site_mean + facet_wrap(~ solute, ncol = 3),
      "02_site_mean_phase_diagram_by_solute.png",
      width = 12,
      height = 8
    )

    save_plot(
      p_site_mean + facet_wrap(~ Stream_Name, ncol = 3),
      "02_site_mean_phase_diagram_by_site.png",
      width = 12,
      height = 9
    )

    message("Created site-mean phase diagrams with ", nrow(df_phase_site), " observations")
  }
} else {
  warning("Primary storage metric not found in site-mean dataset")
}

message("\n=== PHASE DIAGRAMS COMPLETE ===\n")
