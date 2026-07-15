# =============================================================================
# 3a_site_significance_tests.R
# =============================================================================
# Statistical significance tests comparing metrics across sites
# Tests: CQ behavior composition, synchrony patterns, cluster membership
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# Source helpers and shared configuration
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
try(source(file.path(repo_dir, "00_helpers", "plot_prefs.R")), silent = TRUE)
try(source(file.path(repo_dir, "00_helpers", "workflow_config.R")), silent = TRUE)

paths    <- get_project_paths()
out_dir  <- paths$out_dir
fig_dir  <- file.path(paths$fig_root, "03_stats", "3a_site_significance")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
annual <- readr::read_csv(file.path(out_dir, "HJA_clean_annual.csv"), show_col_types = FALSE)
seasonal <- readr::read_csv(file.path(out_dir, "HJA_clean_seasonal.csv"), show_col_types = FALSE)
sync <- suppressWarnings(readr::read_csv(file.path(out_dir, "HJA_composite_synchrony.csv"), show_col_types = FALSE))

annual <- annual %>%
  mutate(
    Stream_Name = factor(Stream_Name, levels = site_order),
    solute = factor(solute, levels = solute_order)
  )

seasonal <- seasonal %>%
  mutate(
    Stream_Name = factor(Stream_Name, levels = site_order),
    solute = factor(solute, levels = solute_order)
  )

sync <- sync %>%
  mutate(Stream_Name = factor(Stream_Name, levels = site_order))

# Helpers
save_csv <- function(df, name) {
  readr::write_csv(df, file.path(out_dir, "03_stats", paste0(name, ".csv")))
}

# Ensure output dir exists
dir.create(file.path(out_dir, "03_stats"), showWarnings = FALSE, recursive = TRUE)

# 1) CQ behavior composition differences across sites per water year (per solute)
# Use proportions of behavior (prop_sync, prop_enrich, prop_dilute, prop_chemostat)
behav_tests <- annual %>%
  pivot_longer(cols = c(prop_sync, prop_enrich, prop_dilute, prop_chemostat), names_to = "behavior", values_to = "prop") %>%
  filter(!is.na(prop)) %>%
  group_by(solute, water_year, behavior) %>%
  group_modify(function(df, key) {
    if (n_distinct(df$Stream_Name) < 2) return(tibble(test = NA_character_, stat = NA_real_, p = NA_real_))
    fit <- try(aov(prop ~ Stream_Name, data = df), silent = TRUE)
    if (inherits(fit, "try-error")) {
      kw <- suppressWarnings(kruskal.test(prop ~ Stream_Name, data = df))
      tibble(test = "Kruskal", stat = unname(kw$statistic), p = kw$p.value)
    } else {
      sm <- summary(fit)[[1]]
      tibble(test = "ANOVA", stat = unname(sm$`F value`[1]), p = sm$`Pr(>F)`[1])
    }
  }) %>%
  ungroup()
save_csv(behav_tests, "wy_behavior_composition_site_differences")

# 2) CQ slope mean differences across sites per water year (per solute)
cq_tests <- annual %>%
  filter(!is.na(cq_slope)) %>%
  group_by(solute, water_year) %>%
  group_modify(function(df, key) {
    if (n_distinct(df$Stream_Name) < 2) return(tibble(test = NA_character_, stat = NA_real_, p = NA_real_))
    fit <- try(aov(cq_slope ~ Stream_Name, data = df), silent = TRUE)
    if (inherits(fit, "try-error")) {
      kw <- suppressWarnings(kruskal.test(cq_slope ~ Stream_Name, data = df))
      tibble(test = "Kruskal", stat = unname(kw$statistic), p = kw$p.value)
    } else {
      sm <- summary(fit)[[1]]
      tibble(test = "ANOVA", stat = unname(sm$`F value`[1]), p = sm$`Pr(>F)`[1])
    }
  }) %>%
  ungroup()
save_csv(cq_tests, "wy_cqslope_site_differences")

# 3) Synchrony differences across water years (Abbott vs Wymore)
# Expect sync to have columns: method (e.g., "Abbott_conc", "Wymore_cqslope"), water_year, value, optional site/pair
sync_long <- sync %>%
  pivot_longer(cols = matches("sync"), names_to = "metric", values_to = "value")

sync_tests <- sync_long %>%
  filter(!is.na(value), !is.na(Stream_Name)) %>%
  group_by(metric) %>%
  group_modify(function(df, key) {
    if (n_distinct(df$Stream_Name) < 2) {
      return(tibble(
        n_years = if ("n_years_sync" %in% names(df)) mean(df$n_years_sync, na.rm = TRUE) else NA_real_,
        test = NA_character_,
        stat = NA_real_,
        p = NA_real_
      ))
    }
    kt <- tryCatch(suppressWarnings(kruskal.test(value ~ Stream_Name, data = df)), error = function(e) NULL)
    tibble(
      n_years = if ("n_years_sync" %in% names(df)) mean(df$n_years_sync, na.rm = TRUE) else NA_real_,
      test = if (!is.null(kt)) "Kruskal" else NA_character_,
      stat = if (!is.null(kt)) unname(kt$statistic) else NA_real_,
      p = if (!is.null(kt)) kt$p.value else NA_real_
    )
  }) %>%
  ungroup()
save_csv(sync_tests, "wy_synchrony_differences")

# Minimal figure: volcano-style p-values
# Standardize columns before binding to avoid NA propagation
pvals <- bind_rows(
  behav_tests %>% 
    filter(!is.na(p)) %>%
    mutate(domain = "Behavior", metric = behavior) %>% 
    select(solute, water_year, metric, domain, test, stat, p),
  cq_tests %>% 
    filter(!is.na(p)) %>%
    mutate(domain = "CQ slope", metric = "cq_slope") %>% 
    select(solute, water_year, metric, domain, test, stat, p),
  sync_tests %>% 
    filter(!is.na(p)) %>%
    mutate(domain = "Synchrony", solute = NA_character_, water_year = NA_integer_) %>% 
    select(solute, water_year, metric, domain, test, stat, p)
) %>%
  mutate(
    p = pmax(p, .Machine$double.xmin),
    neg_log_p = -log10(p),
    label = case_when(
      !is.na(solute) & !is.na(water_year) ~ paste(domain, solute, water_year, metric, sep = " | "),
      !is.na(solute) ~ paste(domain, solute, metric, sep = " | "),
      TRUE ~ paste(domain, metric, sep = " | ")
    ),
    domain = factor(domain, levels = c("Behavior", "CQ slope", "Synchrony"))
  )

if (nrow(pvals) > 0) {
  domain_colors <- c(
    "Behavior" = cluster_colors[["4"]],
    "CQ slope" = cluster_colors[["3"]],
    "Synchrony" = cluster_colors[["2"]]
  )

  pvals <- pvals %>%
    mutate(label = forcats::fct_reorder(label, neg_log_p, .fun = max, .desc = TRUE))

  p <- ggplot(pvals, aes(x = label, y = neg_log_p, fill = domain)) +
    geom_col(width = 0.78) +
    coord_flip() +
    scale_fill_manual(values = domain_colors, name = "Test domain") +
    labs(
      title = "Per–water-year significance tests",
      subtitle = "-log10(p) grouped by metric domain",
      y = "-log10(p)",
      x = "Test group"
    ) +
    theme_hja() +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = BASE_SIZE - 2)
    ) +
    legend_bottom()

  save_plot(p, "wy_significance_overview.png", fig_dir, width = 11, height = 7)
}

message("Significance tests complete: CSVs in outputs/03_stats, figure in exploratory_plots/03_stats/3a_site_significance/")
