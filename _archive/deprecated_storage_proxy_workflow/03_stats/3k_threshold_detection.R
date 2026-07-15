# =============================================================================
# STEP 03k: THRESHOLD DETECTION ANALYSIS
# =============================================================================
# Purpose: Detect thresholds, breakpoints, and nonlinear responses in:
#   - Storage-synchrony relationships
#   - Storage-CQ slope relationships
#   - Critical storage divergence values where synchrony collapses
#
# Methods:
#   1. Segmented regression (breakpoint detection)
#   2. Threshold GAM (smooth transitions)
#   3. CART/Decision trees (natural splits)
#   4. Piecewise linear regression
#
# Research Questions:
#   - Is there a critical storage divergence threshold?
#   - Do CQ regime transitions occur at specific storage values?
#   - Are relationships linear or threshold-driven?
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(segmented)  # Segmented regression
  library(mgcv)       # GAMs
  library(rpart)      # Decision trees
  library(rpart.plot) # Tree plotting
  library(patchwork)
})

rm(list = ls())

# Source helpers
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
try(source(file.path(repo_dir, "00_helpers", "plot_prefs.R")), silent = TRUE)
try(source(file.path(repo_dir, "00_helpers", "workflow_config.R")), silent = TRUE)

# Paths
paths <- get_project_paths()
out_dir <- paths$out_dir
fig_dir <- file.path(paths$fig_root, "03_stats", "3k_thresholds")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "03_stats"), showWarnings = FALSE, recursive = TRUE)

message("\n╔════════════════════════════════════════════════════════════════╗\n")
message("║  THRESHOLD DETECTION: Breakpoints in Storage-Synchrony       ║\n")
message("╚════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# LOAD DATA
# =============================================================================
message("=== LOADING DATA ===\n")

# Annual data with synchrony
annual <- read_csv(file.path(out_dir, "HJA_master_annual.csv"), show_col_types = FALSE)
sync_annual <- read_csv(file.path(out_dir, "HJA_composite_synchrony_annual.csv"), show_col_types = FALSE)
clusters_modal <- read_csv(file.path(out_dir, "ClusterStreams_allSolutes_modalClusters.csv"), show_col_types = FALSE)

# Join
annual_full <- annual %>%
  left_join(sync_annual, by = c("Stream_Name", "solute", "water_year")) %>%
  left_join(
    clusters_modal %>% dplyr::select(Stream_Name, chemical, Cluster_mode),
    by = c("Stream_Name", "solute" = "chemical")
  ) %>%
  mutate(
    Stream_Name = as.factor(Stream_Name),
    solute = as.factor(solute),
    Cluster = factor(Cluster_mode, levels = cluster_levels)
  )

# Storage metric
storage_metric <- intersect(PRIMARY_STORAGE_METRIC, names(annual_full))
if (length(storage_metric) == 0) {
  storage_metric <- intersect(c("Q_dS_range_mm", "WB_dS_range_mm"), names(annual_full))
}
storage_metric <- storage_metric[[1]]
message("  Storage metric:", storage_metric, "\n")

# Compute storage divergence
annual_divergence <- annual_full %>%
  group_by(water_year, solute) %>%
  summarise(
    storage_sd = sd(.data[[storage_metric]], na.rm = TRUE),
    storage_mean = mean(.data[[storage_metric]], na.rm = TRUE),
    n_sites = sum(!is.na(.data[[storage_metric]])),
    .groups = "drop"
  ) %>%
  filter(n_sites >= 3)

# Join back
model_data <- annual_full %>%
  left_join(annual_divergence, by = c("water_year", "solute")) %>%
  filter(!is.na(storage_sd))

message("  Model data:", nrow(model_data), "observations\n")

# =============================================================================
# METHOD 1: SEGMENTED REGRESSION (Breakpoint Detection)
# =============================================================================
message("\n=== METHOD 1: SEGMENTED REGRESSION ===\n")

response_vars <- c("conc_sync_allpairs", "cqslope_sync_allpairs")
response_vars <- intersect(response_vars, names(model_data))

threshold_results <- list()

for (resp in response_vars) {

  message("\n--- Response:", resp, "---\n")

  # Prepare data
  df_seg <- model_data %>%
    filter(!is.na(.data[[resp]]), !is.na(storage_sd), is.finite(storage_sd)) %>%
    dplyr::select(storage_sd, y = all_of(resp)) %>%
    arrange(storage_sd)

  if (nrow(df_seg) < 20) {
    message("  Insufficient data (n =", nrow(df_seg), "); skipping\n")
    next
  }

  message("  N =", nrow(df_seg), "observations\n")

  # Fit linear model first
  lm_base <- lm(y ~ storage_sd, data = df_seg)

  # Fit segmented regression (automatic breakpoint detection)
  seg_fit <- tryCatch({
    segmented(lm_base, seg.Z = ~storage_sd, psi = NA, control = seg.control(display = FALSE))
  }, error = function(e) {
    message("  Segmented regression failed; trying with initial guess\n")
    # Try with initial guess at median
    psi_guess <- median(df_seg$storage_sd, na.rm = TRUE)
    tryCatch({
      segmented(lm_base, seg.Z = ~storage_sd, psi = psi_guess, control = seg.control(display = FALSE))
    }, error = function(e) NULL)
  })

  if (is.null(seg_fit)) {
    message("  Could not fit segmented model; using linear model only\n")
    threshold_results[[resp]] <- tibble(
      response = resp,
      method = "segmented",
      breakpoint = NA_real_,
      breakpoint_se = NA_real_,
      slope_below = coef(lm_base)["storage_sd"],
      slope_above = NA_real_,
      davies_test_p = NA_real_
    )
    next
  }

  # Extract breakpoint
  bp <- summary(seg_fit)$psi
  breakpoint <- bp[1, "Est."]
  breakpoint_se <- bp[1, "St.Err"]

  # Extract slopes
  slopes <- slope(seg_fit, APC = FALSE)
  slope_below <- slopes$storage_sd[1, 1]
  slope_above <- slopes$storage_sd[2, 1]

  # Davies test for breakpoint existence
  davies_p <- tryCatch({
    davies.test(lm_base, seg.Z = ~storage_sd)$p.value
  }, error = function(e) NA_real_)

  message("  Breakpoint:", round(breakpoint, 3), "±", round(breakpoint_se, 3), "\n")
  message("  Slope below breakpoint:", round(slope_below, 4), "\n")
  message("  Slope above breakpoint:", round(slope_above, 4), "\n")
  message("  Davies test p-value:", round(davies_p, 4), "\n")

  # Store results
  threshold_results[[resp]] <- tibble(
    response = resp,
    method = "segmented",
    breakpoint = breakpoint,
    breakpoint_se = breakpoint_se,
    slope_below = slope_below,
    slope_above = slope_above,
    davies_test_p = davies_p
  )

  # Visualize
  df_seg$fitted_seg <- fitted(seg_fit)

  p_seg <- ggplot(df_seg, aes(x = storage_sd, y = y)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_line(aes(y = fitted_seg), color = "red", linewidth = 1.2) +
    geom_vline(xintercept = breakpoint, linetype = "dashed", color = "darkred", linewidth = 1) +
    annotate("text", x = breakpoint, y = max(df_seg$y, na.rm = TRUE),
             label = paste0("Breakpoint = ", round(breakpoint, 2)),
             hjust = -0.1, vjust = 1, color = "darkred", fontface = "bold") +
    labs(
      x = paste0("Storage Divergence (SD of ", get_storage_label(storage_metric, short = TRUE), ")"),
      y = get_sync_label(resp),
      title = paste("Segmented Regression:", get_sync_label(resp)),
      subtitle = paste0("Breakpoint: ", round(breakpoint, 2), " ± ", round(breakpoint_se, 2),
                       " | Davies test p = ", format.pval(davies_p, digits = 3)),
      caption = paste0("Slope below: ", round(slope_below, 3), " | Slope above: ", round(slope_above, 3))
    ) +
    theme_hja()

  ggsave(file.path(fig_dir, paste0("segmented_", resp, ".png")), p_seg,
         width = 11, height = 8, dpi = 300)

  message("  ✓ Segmented plot saved\n")
}

# =============================================================================
# METHOD 2: THRESHOLD GAM (Smooth Nonlinearity)
# =============================================================================
message("\n=== METHOD 2: THRESHOLD GAM ===\n")

for (resp in response_vars) {

  message("\n--- Response:", resp, "---\n")

  # Prepare data
  df_gam <- model_data %>%
    filter(!is.na(.data[[resp]]), !is.na(storage_sd), is.finite(storage_sd)) %>%
    dplyr::select(storage_sd, y = all_of(resp))

  if (nrow(df_gam) < 20) {
    message("  Insufficient data; skipping\n")
    next
  }

  message("  N =", nrow(df_gam), "observations\n")

  # Fit GAM with penalized spline
  gam_fit <- gam(y ~ s(storage_sd, k = 5), data = df_gam, method = "REML")

  # Summary
  message("\n  GAM summary:\n")
  print(summary(gam_fit))

  # Deviance explained
  dev_exp <- summary(gam_fit)$dev.expl * 100
  message("\n  Deviance explained:", round(dev_exp, 1), "%\n")

  # Prediction grid
  pred_grid <- tibble(storage_sd = seq(min(df_gam$storage_sd), max(df_gam$storage_sd), length.out = 200))
  pred_fit <- predict(gam_fit, newdata = pred_grid, se.fit = TRUE, type = "response")
  pred_grid$fit <- pred_fit$fit
  pred_grid$se <- pred_fit$se.fit
  pred_grid$lwr <- pred_grid$fit - 1.96 * pred_grid$se
  pred_grid$upr <- pred_grid$fit + 1.96 * pred_grid$se

  # Detect potential threshold (where slope changes most)
  # Compute first derivative
  deriv <- diff(pred_grid$fit) / diff(pred_grid$storage_sd)
  max_change_idx <- which.max(abs(deriv))
  threshold_gam <- pred_grid$storage_sd[max_change_idx]

  message("  Potential threshold (max slope change):", round(threshold_gam, 3), "\n")

  # Visualize
  p_gam <- ggplot(df_gam, aes(x = storage_sd, y = y)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_ribbon(data = pred_grid, aes(y = fit, ymin = lwr, ymax = upr),
                inherit.aes = FALSE, fill = "blue", alpha = 0.2) +
    geom_line(data = pred_grid, aes(y = fit), color = "blue", linewidth = 1.2) +
    geom_vline(xintercept = threshold_gam, linetype = "dashed", color = "darkblue") +
    annotate("text", x = threshold_gam, y = max(df_gam$y, na.rm = TRUE),
             label = paste0("Max slope change = ", round(threshold_gam, 2)),
             hjust = -0.1, vjust = 1, color = "darkblue", fontface = "bold") +
    labs(
      x = paste0("Storage Divergence (SD of ", get_storage_label(storage_metric, short = TRUE), ")"),
      y = get_sync_label(resp),
      title = paste("Threshold GAM:", get_sync_label(resp)),
      subtitle = paste0("Deviance explained: ", round(dev_exp, 1), "%"),
      caption = "Shaded area = 95% CI | Dashed line = point of maximum slope change"
    ) +
    theme_hja()

  ggsave(file.path(fig_dir, paste0("gam_", resp, ".png")), p_gam,
         width = 11, height = 8, dpi = 300)

  message("  ✓ GAM plot saved\n")
}

# =============================================================================
# METHOD 3: CART (Decision Tree)
# =============================================================================
message("\n=== METHOD 3: DECISION TREE (CART) ===\n")

for (resp in response_vars) {

  message("\n--- Response:", resp, "---\n")

  # Prepare data with additional predictors
  df_tree <- model_data %>%
    filter(!is.na(.data[[resp]]), !is.na(storage_sd), !is.na(Cluster)) %>%
    dplyr::select(y = all_of(resp), storage_sd, Cluster, storage_mean) %>%
    drop_na()

  if (nrow(df_tree) < 20) {
    message("  Insufficient data; skipping\n")
    next
  }

  message("  N =", nrow(df_tree), "observations\n")

  # Fit regression tree
  tree_fit <- rpart(y ~ storage_sd + Cluster + storage_mean,
                    data = df_tree,
                    method = "anova",
                    control = rpart.control(minsplit = 10, cp = 0.01))

  # Print tree
  message("\n  Decision tree:\n")
  print(tree_fit)

  # Variable importance
  var_imp <- tree_fit$variable.importance
  if (!is.null(var_imp)) {
    message("\n  Variable importance:\n")
    print(var_imp)
  }

  # Plot tree
  png(file.path(fig_dir, paste0("tree_", resp, ".png")),
      width = 12, height = 9, units = "in", res = 300)
  rpart.plot(tree_fit,
             main = paste("Decision Tree:", get_sync_label(resp)),
             box.palette = "auto",
             shadow.col = "gray",
             nn = TRUE)
  dev.off()

  message("  ✓ Tree plot saved\n")

  # Extract thresholds (splits on storage_sd)
  splits <- tree_fit$splits
  if (!is.null(splits) && nrow(splits) > 0) {
    storage_splits <- splits[grepl("storage", rownames(splits)), , drop = FALSE]
    if (nrow(storage_splits) > 0) {
      message("\n  Storage splits detected:\n")
      print(storage_splits)
    }
  }
}

# =============================================================================
# METHOD 4: PIECEWISE LINEAR (Two Regimes)
# =============================================================================
message("\n=== METHOD 4: PIECEWISE LINEAR REGRESSION ===\n")

for (resp in response_vars) {

  message("\n--- Response:", resp, "---\n")

  df_pw <- model_data %>%
    filter(!is.na(.data[[resp]]), !is.na(storage_sd)) %>%
    dplyr::select(storage_sd, y = all_of(resp)) %>%
    arrange(storage_sd)

  if (nrow(df_pw) < 20) {
    message("  Insufficient data; skipping\n")
    next
  }

  # Test a grid of potential breakpoints
  grid_breaks <- quantile(df_pw$storage_sd, probs = seq(0.2, 0.8, by = 0.1), na.rm = TRUE)
  aic_vals <- numeric(length(grid_breaks))

  for (i in seq_along(grid_breaks)) {
    bp <- grid_breaks[i]
    df_pw$regime <- if_else(df_pw$storage_sd <= bp, "low", "high")

    # Fit piecewise model
    pw_model <- lm(y ~ storage_sd * regime, data = df_pw)
    aic_vals[i] <- AIC(pw_model)
  }

  # Best breakpoint
  best_idx <- which.min(aic_vals)
  best_bp <- grid_breaks[best_idx]
  best_aic <- aic_vals[best_idx]

  message("  Best breakpoint (AIC):", round(best_bp, 3), "(AIC =", round(best_aic, 1), ")\n")

  # Fit final model
  df_pw$regime <- if_else(df_pw$storage_sd <= best_bp, "low", "high")
  pw_final <- lm(y ~ storage_sd * regime, data = df_pw)

  message("\n  Piecewise model:\n")
  print(summary(pw_final))

  # Plot
  df_pw$fitted_pw <- fitted(pw_final)

  p_pw <- ggplot(df_pw, aes(x = storage_sd, y = y, color = regime)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_line(aes(y = fitted_pw, group = regime), linewidth = 1.2) +
    geom_vline(xintercept = best_bp, linetype = "dashed", color = "black") +
    scale_color_manual(values = c("low" = cluster_colors[["1"]], "high" = cluster_colors[["3"]]),
                       name = "Regime") +
    annotate("text", x = best_bp, y = max(df_pw$y, na.rm = TRUE),
             label = paste0("Breakpoint = ", round(best_bp, 2)),
             hjust = -0.1, vjust = 1, fontface = "bold") +
    labs(
      x = paste0("Storage Divergence (SD of ", get_storage_label(storage_metric, short = TRUE), ")"),
      y = get_sync_label(resp),
      title = paste("Piecewise Linear:", get_sync_label(resp)),
      subtitle = paste0("Best breakpoint: ", round(best_bp, 2), " (AIC = ", round(best_aic, 1), ")"),
      caption = "Two regimes: low storage divergence vs. high storage divergence"
    ) +
    theme_hja() +
    legend_bottom()

  ggsave(file.path(fig_dir, paste0("piecewise_", resp, ".png")), p_pw,
         width = 11, height = 8, dpi = 300)

  message("  ✓ Piecewise plot saved\n")
}

# =============================================================================
# SAVE THRESHOLD RESULTS
# =============================================================================
message("\n=== SAVING RESULTS ===\n")

if (length(threshold_results) > 0) {
  threshold_table <- bind_rows(threshold_results)
  write_csv(threshold_table, file.path(out_dir, "03_stats/threshold_detection_results.csv"))
  message("  Threshold results saved\n")
}

# =============================================================================
# SUMMARY
# =============================================================================
message("\n╔════════════════════════════════════════════════════════════════╗\n")
message("║  ANALYSIS COMPLETE                                            ║\n")
message("╚════════════════════════════════════════════════════════════════╝\n\n")

message("INTERPRETATION:\n\n")

message("Segmented regression:\n")
message("  → Detects sharp breakpoints in the relationship\n")
message("  → Davies test: p < 0.05 means breakpoint is significant\n")
message("  → Slopes before and after breakpoint can differ\n\n")

message("Threshold GAM:\n")
message("  → Detects smooth, gradual transitions\n")
message("  → More flexible than linear models\n")
message("  → Point of max slope change = potential threshold\n\n")

message("Decision tree:\n")
message("  → Finds natural splits in data\n")
message("  → Includes interactions with other variables (Cluster, mean storage)\n")
message("  → Easy interpretation (if-then rules)\n\n")

message("Piecewise linear:\n")
message("  → Two linear regimes with different slopes\n")
message("  → AIC-based selection of breakpoint\n\n")

if (exists("threshold_table")) {
  message("DETECTED THRESHOLDS:\n\n")
  print(threshold_table, width = Inf)
}

message("\nOutputs saved to:\n")
message("  Tables:", file.path(out_dir, "03_stats/"), "\n")
message("  Figures:", fig_dir, "\n\n")

message("Key files:\n")
message("  - threshold_detection_results.csv\n")
message("  - segmented_*.png (breakpoint models)\n")
message("  - gam_*.png (smooth nonlinearity)\n")
message("  - tree_*.png (decision trees)\n")
message("  - piecewise_*.png (two-regime models)\n\n")
