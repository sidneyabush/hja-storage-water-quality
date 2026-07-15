# =============================================================================
# STEP 03f: STORAGE METHOD COMPARISON & SYNCHRONY ANALYSIS
# =============================================================================
# Goal: Compare two storage estimation methods and their relationships with
#       CQ behavior and synchrony metrics
#
# Storage methods:
#   - Q_dS_range_mm: Discharge-based storage (Q - dS)
#   - WB_dS_range_mm: Water balance storage (P - ET - Q) - requires precip data
#
# Synchrony metrics compared:
#   - Abbott absolute concentration synchrony (conc_sync_allpairs, conc_sync_outlet)
#   - Abbott absolute CQ-slope synchrony (cqslope_sync_allpairs, cqslope_sync_outlet)
#   - Wymore cross-site consistency (wymore_cvcq_consistency, wymore_crosssite)
#
# Also incorporates annual storage metrics with baseflow from:
#   - HJA_StorageMetrics_Annual.csv (mean_bf, rbfi, recession_curve_slope)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(corrplot)
})
try(source(file.path("/Users/sidneybush/Documents/GitHub/hja-water-quality", "00_helpers", "plot_theme_set.R")), silent = TRUE)

rm(list = ls())

# Source helpers
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
try(source(file.path(repo_dir, "00_helpers", "plot_prefs.R")), silent = TRUE)
try(source(file.path(repo_dir, "00_helpers", "workflow_config.R")), silent = TRUE)

# Paths
base_dir   <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality"
data_dir   <- file.path(base_dir, "data")
output_dir <- file.path(base_dir, "outputs")
plot_dir   <- file.path(base_dir, "exploratory_plots", "03_stats", "3f_storage_sync_comparison")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "03_stats"), showWarnings = FALSE, recursive = TRUE)

# Theme
theme_clean <- function(base_size = BASE_SIZE) {
  theme_hja(base_size = base_size) +
    theme(
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "white", color = NA),
      legend.title = element_text(face = "bold")
    )
}

message("\n╔════════════════════════════════════════════════════════════════╗\n")
message("║  STORAGE METHOD & SYNCHRONY COMPARISON                        ║\n")
message("╚════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# LOAD DATA
# =============================================================================
message("=== LOADING DATA ===\n")

# Main window-level data
mega <- read_csv(file.path(output_dir, "HJA_mega_90d_windows_CQ_hydro_static.csv"), 
                 show_col_types = FALSE)

# Site-level synchrony data
site_data <- read_csv(file.path(output_dir, "HJA_exploratory_site.csv"), 
                      show_col_types = FALSE)

# Annual storage metrics with baseflow
annual_storage <- tryCatch({
  read_csv(file.path(data_dir, "HJA_StorageMetrics_Annual.csv"), show_col_types = FALSE) %>%
    rename(Stream_Name = site, water_year = year)
}, error = function(e) {
  message("  Note: Annual storage metrics file not found\n")
  NULL
})

message("  Mega data:", nrow(mega), "rows\n")
message("  Site data:", nrow(site_data), "rows\n")
if (!is.null(annual_storage)) message("  Annual storage:", nrow(annual_storage), "rows\n")

# Check available columns
storage_cols <- c(PRIMARY_STORAGE_METRIC, "Q_dS_range_mm", "WB_dS_range_mm", "Q_dS_net_mm", "WB_dS_net_mm")
storage_cols <- unique(storage_cols)
storage_avail <- storage_cols[storage_cols %in% names(mega)]
message("\n  Available storage columns:", paste(storage_avail, collapse = ", "), "\n")

sync_cols <- c("conc_sync_allpairs", "cqslope_sync_allpairs", "conc_sync_outlet", 
               "cqslope_sync_outlet", "wymore_cvcq_consistency", 
               "wymore_crosssite_allpairs", "wymore_crosssite_outlet")
sync_avail <- sync_cols[sync_cols %in% names(site_data)]
message("  Available sync columns:", paste(sync_avail, collapse = ", "), "\n")

# =============================================================================
# SECTION 1: COMPARE Q-dS vs WB-dS STORAGE METHODS
# =============================================================================
message("\n=== SECTION 1: Q-dS vs WB-dS COMPARISON ===\n")

if ("Q_dS_range_mm" %in% names(mega) && "WB_dS_range_mm" %in% names(mega)) {
  
  wb_coverage <- sum(!is.na(mega$WB_dS_range_mm)) / nrow(mega) * 100
  q_coverage <- sum(!is.na(mega$Q_dS_range_mm)) / nrow(mega) * 100
  message("  Q_dS coverage:", round(q_coverage, 1), "%\n")
  message("  WB_dS coverage:", round(wb_coverage, 1), "%\n")

  both_valid <- mega %>%
    filter(!is.na(Q_dS_range_mm), !is.na(WB_dS_range_mm)) %>%
    filter(Q_dS_range_mm > 0, WB_dS_range_mm > 0) %>%
    mutate(Stream_Name = factor(Stream_Name, levels = site_order))

  if (nrow(both_valid) > 100) {
    cor_methods <- cor.test(both_valid$Q_dS_range_mm, both_valid$WB_dS_range_mm)
    message("\n  Correlation between methods:\n")
    message("    r =", round(cor_methods$estimate, 3), ", p =", format.pval(cor_methods$p.value), "\n")
    message("    n =", nrow(both_valid), "windows with both values\n")

    scatter_sample <- both_valid %>% sample_frac(min(1, 12000 / nrow(both_valid)))
    p_methods <- ggplot(scatter_sample, aes(x = Q_dS_range_mm, y = WB_dS_range_mm)) +
      geom_point(aes(color = Stream_Name), alpha = 0.35, size = 1) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
      geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7) +
      scale_color_site() +
      scale_x_log10(labels = scales::label_number()) +
      scale_y_log10(labels = scales::label_number()) +
      labs(
        x = get_storage_label("Q_dS_range_mm"),
        y = get_storage_label("WB_dS_range_mm"),
        title = paste0("Storage Method Comparison (r = ", round(cor_methods$estimate, 2), ")"),
        subtitle = "Points = 90-day windows; color-coded by site",
        caption = "Dashed line shows 1:1 agreement; solid line is a least-squares fit"
      ) +
      theme_clean() +
      legend_bottom()

    save_plot(p_methods, "01_storage_method_comparison.png", plot_dir, width = 11, height = 9)

    if ("hydrologic_season" %in% names(both_valid)) {
      p_methods_season <- both_valid %>%
        filter(!is.na(hydrologic_season)) %>%
        mutate(hydrologic_season = dplyr::case_when(
          hydrologic_season == "wet" ~ "Wet",
          hydrologic_season == "dry" ~ "Dry",
          TRUE ~ stringr::str_to_title(hydrologic_season)
        )) %>%
        ggplot(aes(x = Q_dS_range_mm, y = WB_dS_range_mm)) +
        geom_point(aes(color = hydrologic_season), alpha = 0.4, size = 1) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
        geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.6) +
        scale_color_manual(values = c("Wet" = "#2166AC", "Dry" = "#B2182B"), name = "Season") +
        scale_x_log10(labels = scales::label_number()) +
        scale_y_log10(labels = scales::label_number()) +
        facet_wrap(~hydrologic_season) +
        labs(
          x = get_storage_label("Q_dS_range_mm"),
          y = get_storage_label("WB_dS_range_mm"),
          title = "Storage Method Comparison by Season"
        ) +
        theme_clean()

      save_plot(p_methods_season, "02_storage_method_by_season.png", plot_dir, width = 12, height = 6)
    }
  }
}

# =============================================================================
# SECTION 2: CQ SLOPE vs BOTH STORAGE METHODS
# =============================================================================
message("\n=== SECTION 2: CQ SLOPE vs STORAGE METHODS ===\n")

if ("cq_slope.x" %in% names(mega)) {
  
  # Correlations with each storage method
  storage_cq_cors <- tibble(
    Storage_Method = character(),
    Correlation = numeric(),
    p_value = numeric(),
    n = integer()
  )
  
  for (s_col in storage_avail) {
    valid_data <- mega %>% filter(!is.na(cq_slope.x) & !is.na(.data[[s_col]]))
    if (nrow(valid_data) > 100) {
      ct <- cor.test(valid_data$cq_slope.x, valid_data[[s_col]])
      storage_cq_cors <- bind_rows(storage_cq_cors,
                                    tibble(Storage_Method = s_col,
                                           Correlation = ct$estimate,
                                           p_value = ct$p.value,
                                           n = nrow(valid_data)))
    }
  }
  
  message("\nCQ slope correlations with storage methods:\n")
  print(storage_cq_cors)
  write_csv(storage_cq_cors, file.path(output_dir, "03_stats/storage_cq_correlations.csv"))
  
  # Side-by-side plots
  plot_list <- list()
  storage_subset <- storage_avail[seq_len(min(2, length(storage_avail)))]
  for (s_col in storage_subset) {
    if (!s_col %in% names(mega)) next
    plot_data <- mega %>%
      filter(!is.na(cq_slope.x), !is.na(.data[[s_col]]), .data[[s_col]] > 0) %>%
      mutate(Stream_Name = factor(Stream_Name, levels = site_order))
    if (nrow(plot_data) == 0) next
    plot_sample <- plot_data %>% sample_frac(min(1, 15000 / nrow(plot_data)))
    plot_list[[s_col]] <- ggplot(plot_sample, aes(x = .data[[s_col]], y = cq_slope.x)) +
      geom_point(aes(color = solute), alpha = 0.25, size = 0.9) +
      geom_smooth(method = "lm", se = FALSE, color = "grey35", linewidth = 0.6) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.5) +
      scale_x_log10(labels = scales::label_number()) +
      scale_color_solute() +
      facet_wrap(~Stream_Name, ncol = 3) +
      labs(
        x = get_storage_label(s_col),
        y = "CQ Slope",
        title = get_storage_label(s_col, short = TRUE),
        subtitle = "Faceted by site; points are 90-day windows"
      ) +
      theme_clean(base_size = BASE_SIZE - 1) +
      legend_bottom()
  }

  if (length(plot_list) > 0) {
    ncol <- min(2, length(plot_list))
    p_combined <- wrap_plots(plot_list, ncol = ncol) +
      plot_annotation(title = "CQ Slope vs Storage: Method Comparison")
    save_plot(p_combined, "03_cq_vs_storage_methods.png", plot_dir, width = 14, height = 7)
  }
}

# =============================================================================
# SECTION 3: ABBOTT vs WYMORE SYNCHRONY COMPARISON
# =============================================================================
message("\n=== SECTION 3: ABBOTT vs WYMORE SYNCHRONY ===\n")

# Abbott metrics: conc_sync, cqslope_sync (correlation-based spatial coherence)
# Wymore metrics: wymore_cvcq_consistency, wymore_crosssite (CVc/CVq ratio similarity)

if (length(sync_avail) >= 2) {
  
  # Correlation matrix of synchrony metrics
  sync_data <- site_data %>%
    select(Stream_Name, solute, all_of(sync_avail)) %>%
    drop_na()
  
  if (nrow(sync_data) > 20) {
    sync_corr <- cor(sync_data %>% select(all_of(sync_avail)), use = "pairwise.complete.obs")
    
    # Higher resolution PNG for crisp figures (300 DPI)
    png(file.path(plot_dir, "04_sync_method_correlations.png"), 
        width = 10, height = 9, units = "in", res = 300)
    corrplot(sync_corr, method = "color", type = "upper",
             addCoef.col = "black", number.cex = 0.8,
             tl.col = "black", tl.srt = 45,
             title = "Synchrony Metric Correlations\n(Abbott vs Wymore)", mar = c(0,0,3,0))
    dev.off()
    
    message("\nSynchrony metric correlations saved\n")
  }
  
  # Compare concentration sync vs CQ-slope sync (Abbott)
  if ("conc_sync_allpairs" %in% names(site_data) && "cqslope_sync_allpairs" %in% names(site_data)) {
    abbott_sync <- site_data %>%
      filter(!is.na(conc_sync_allpairs), !is.na(cqslope_sync_allpairs)) %>%
      flag_outlet_stream() %>%
      mutate(Stream_Name = factor(Stream_Name, levels = site_order))

    if (nrow(abbott_sync) > 0) {
      p_abbott <- ggplot(abbott_sync, aes(x = conc_sync_allpairs, y = cqslope_sync_allpairs)) +
        geom_point(aes(color = solute, shape = outlet_marker), size = 2.5, alpha = 0.85) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
        scale_color_solute() +
        scale_shape_manual(values = c(16, OUTLET_SHAPE_TRIANGLE), name = "Contains GSLOOK") +
        facet_wrap(~Stream_Name, ncol = 3) +
        labs(
          x = get_sync_label("conc_sync_allpairs"),
          y = get_sync_label("cqslope_sync_allpairs"),
          title = "Abbott synchrony comparison",
          subtitle = "Faceted by site; colors denote solute"
        ) +
        theme_clean() +
        legend_bottom()

      save_plot(p_abbott, "05_abbott_conc_vs_cq_sync.png", plot_dir, width = 12, height = 9)
    }
  }
}

# =============================================================================
# SECTION 4: STORAGE vs SYNCHRONY RELATIONSHIPS
# =============================================================================
message("\n=== SECTION 4: STORAGE vs SYNCHRONY ===\n")

# Aggregate window-level storage to match site-level sync
window_storage <- mega %>%
  group_by(Stream_Name, solute) %>%
  summarise(
    across(all_of(storage_avail), ~mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
    mean_cq = mean(cq_slope.x, na.rm = TRUE),
    pct_sync = mean(cq_sync == "sync", na.rm = TRUE) * 100,
    n_windows = n(),
    .groups = "drop"
  )

# Join with site-level sync metrics
storage_sync <- window_storage %>%
  left_join(site_data %>% select(Stream_Name, solute, all_of(sync_avail)),
            by = c("Stream_Name", "solute")) %>%
  flag_outlet_stream() %>%
  mutate(Stream_Name = factor(Stream_Name, levels = site_order))

storage_mean_cols <- paste0("mean_", storage_avail)
storage_mean_cols <- storage_mean_cols[storage_mean_cols %in% names(storage_sync)]

message("\nStorage vs Synchrony correlations:\n")
sync_storage_cors <- tibble()

for (sync_col in sync_avail) {
  for (storage_col in storage_mean_cols) {
    valid <- storage_sync %>%
      filter(!is.na(.data[[sync_col]]), !is.na(.data[[storage_col]]))
    if (nrow(valid) > 10) {
      ct <- cor.test(valid[[sync_col]], valid[[storage_col]])
      sync_storage_cors <- bind_rows(
        sync_storage_cors,
        tibble(
          Sync_Metric = sync_col,
          Storage_Method = storage_col,
          r = ct$estimate,
          p = ct$p.value,
          n = nrow(valid)
        )
      )
    }
  }
}

if (nrow(sync_storage_cors) > 0) {
  sync_storage_cors <- sync_storage_cors %>%
    mutate(significant = p < 0.05, abs_r = abs(r)) %>%
    arrange(desc(abs_r))

  print(sync_storage_cors, n = 20)
  write_csv(sync_storage_cors, file.path(output_dir, "03_stats/sync_storage_correlations.csv"))

  top_rel <- sync_storage_cors %>% slice(1)
  storage_metric_top <- stringr::str_remove(top_rel$Storage_Method, "^mean_")

  top_data <- storage_sync %>%
    filter(
      !is.na(.data[[top_rel$Sync_Metric]]),
      !is.na(.data[[top_rel$Storage_Method]])
    )

  if (nrow(top_data) > 0) {
    p_top <- ggplot(top_data, aes(x = .data[[top_rel$Storage_Method]], y = .data[[top_rel$Sync_Metric]])) +
      geom_point(aes(color = solute, shape = outlet_marker), size = 2.5, alpha = 0.85) +
      geom_smooth(method = "lm", se = FALSE, color = "grey35", linewidth = 0.6) +
      facet_wrap(~Stream_Name, ncol = 3) +
      scale_color_solute() +
      scale_shape_manual(values = c(16, OUTLET_SHAPE_TRIANGLE), name = "Contains GSLOOK") +
      labs(
        x = get_storage_label(storage_metric_top),
        y = get_sync_label(top_rel$Sync_Metric),
        title = "Storage vs synchrony by site",
        subtitle = paste0("Strongest relationship: r = ", round(top_rel$r, 2), 
                          ", n = ", top_rel$n),
        caption = "Points represent solute-level means; lines show site-level linear fits"
      ) +
      theme_clean() +
      legend_bottom()

    save_plot(p_top, "06_strongest_storage_sync.png", plot_dir, width = 12, height = 8)
  }
}

# =============================================================================
# SECTION 5: ANNUAL BASEFLOW METRICS
# =============================================================================
message("\n=== SECTION 5: ANNUAL BASEFLOW & STORAGE METRICS ===\n")

if (!is.null(annual_storage)) {
  message("\nAnnual storage columns:", paste(names(annual_storage), collapse = ", "), "\n")
  
  # Classify water years
  annual_storage <- annual_storage %>%
    group_by(Stream_Name) %>%
    mutate(
      Q5_percentile = percent_rank(Q5norm),
      wy_type = case_when(
        Q5_percentile < 0.2 ~ "Very Dry",
        Q5_percentile < 0.4 ~ "Dry",
        Q5_percentile < 0.6 ~ "Normal",
        Q5_percentile < 0.8 ~ "Wet",
        TRUE ~ "Very Wet"
      ),
      wy_type = factor(wy_type, levels = c("Very Dry", "Dry", "Normal", "Wet", "Very Wet"))
    ) %>%
    ungroup()
  
  message("\nWater year type distribution:\n")
  print(table(annual_storage$wy_type))
  
  # Save classified annual data
  write_csv(annual_storage, file.path(output_dir, "03_stats/annual_storage_classified.csv"))
  
  # Join with mega data
  mega_annual <- mega %>%
    left_join(annual_storage %>% select(Stream_Name, water_year, mean_bf, rbfi, 
                                         recession_curve_slope, Q5norm, wy_type),
              by = c("Stream_Name", "water_year"))
  
  # Baseflow vs CQ slope
  if ("mean_bf" %in% names(mega_annual) && sum(!is.na(mega_annual$mean_bf)) > 100) {
    bf_cq_cor <- cor.test(mega_annual$mean_bf, mega_annual$cq_slope.x, use = "complete.obs")
    message("\nBaseflow vs CQ slope: r =", round(bf_cq_cor$estimate, 3),
        ", p =", format.pval(bf_cq_cor$p.value), "\n")

    bf_plot_data <- mega_annual %>%
      filter(!is.na(mean_bf), !is.na(cq_slope.x)) %>%
      mutate(Stream_Name = factor(Stream_Name, levels = site_order))

    p_bf <- ggplot(bf_plot_data, aes(x = mean_bf, y = cq_slope.x)) +
      geom_point(aes(color = Stream_Name), alpha = 0.35, size = 1.2) +
      geom_smooth(method = "lm", se = FALSE, color = "grey35", linewidth = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
      scale_color_site(name = "Site") +
      labs(
        x = "Mean Baseflow",
        y = "CQ slope",
        title = paste0("Baseflow vs CQ slope (r = ", round(bf_cq_cor$estimate, 2), ")"),
        subtitle = "Points represent 90-day windows with paired baseflow estimates"
      ) +
      theme_clean() +
      legend_bottom()

    save_plot(p_bf, "07_baseflow_vs_cq.png", plot_dir, width = 11, height = 8)
  }
  
  # CQ by water year type
  if ("wy_type" %in% names(mega_annual)) {
    cq_by_wytype <- mega_annual %>%
      filter(!is.na(wy_type) & !is.na(cq_slope.x)) %>%
      group_by(wy_type) %>%
      summarise(
        mean_cq = mean(cq_slope.x, na.rm = TRUE),
        sd_cq = sd(cq_slope.x, na.rm = TRUE),
        n = n(),
        .groups = "drop"
      )
    
    message("\nCQ slope by water year type:\n")
    print(cq_by_wytype)
    
    wytype_plot_data <- mega_annual %>%
      filter(!is.na(wy_type), !is.na(cq_slope.x)) %>%
      mutate(Stream_Name = factor(Stream_Name, levels = site_order))

    p_wytype <- ggplot(wytype_plot_data, aes(x = wy_type, y = cq_slope.x)) +
      geom_jitter(width = 0.15, height = 0, alpha = 0.45, color = "#526B8E", size = 1.2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
      facet_wrap(~Stream_Name, ncol = 3) +
      labs(
        x = "Water year type",
        y = "CQ slope",
        title = "CQ slope by water year type",
        subtitle = "Per-year values; panels separate sites"
      ) +
      theme_clean() +
      theme(axis.text.x = element_text(angle = 12, hjust = 1))

    save_plot(p_wytype, "08_cq_by_wy_type_5cat.png", plot_dir, width = 12, height = 8)
  }
}

# =============================================================================
# SECTION 6: SUMMARY TABLE
# =============================================================================
message("\n=== GENERATING SUMMARY ===\n")

summary_table <- tibble(
  Metric = c("Windows with Q_dS", "Windows with WB_dS", "Windows with both",
             "Storage method correlation", "Best CQ-storage correlation"),
  Value = c(
    as.character(sum(!is.na(mega$Q_dS_range_mm))),
    as.character(sum(!is.na(mega$WB_dS_range_mm))),
    as.character(sum(!is.na(mega$Q_dS_range_mm) & !is.na(mega$WB_dS_range_mm))),
    if (exists("cor_methods")) as.character(round(cor_methods$estimate, 3)) else "NA",
    if (nrow(storage_cq_cors) > 0) paste0(storage_cq_cors$Storage_Method[1], ": r=", round(storage_cq_cors$Correlation[1], 3)) else "NA"
  )
)

write_csv(summary_table, file.path(output_dir, "03_stats/storage_sync_summary.csv"))

# =============================================================================
# SUMMARY
# =============================================================================
message("\n╔════════════════════════════════════════════════════════════════╗\n")
message("║  STORAGE & SYNCHRONY COMPARISON COMPLETE                      ║\n")
message("╚════════════════════════════════════════════════════════════════╝\n\n")

message("Key findings:\n")
if (exists("cor_methods")) {
  message("  - Q-dS vs WB-dS correlation: r =", round(cor_methods$estimate, 2), "\n")
}
if (nrow(storage_cq_cors) > 0) {
  message("  - Best CQ-storage relationship:", storage_cq_cors$Storage_Method[1], 
      "(r =", round(storage_cq_cors$Correlation[1], 2), ")\n")
}
message("\nOutputs saved to:\n")
message("  Plots:", plot_dir, "\n")
message("  Tables:", file.path(output_dir, "03_stats"), "\n\n")
