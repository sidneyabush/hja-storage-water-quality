# =============================================================================
# STEP 05: FORMAL MIXED-EFFECTS MODELS (moved from 03_stats)
# =============================================================================

suppressPackageStartupMessages({
	library(tidyverse)
	library(lme4)
	library(lmerTest)
	library(MuMIn)
	library(broom.mixed)
	library(patchwork)
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
fig_dir <- file.path(paths$fig_root, "05_mixed_effects")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "05_mixed_effects"), showWarnings = FALSE, recursive = TRUE)

# Load data
annual <- readr::read_csv(file.path(out_dir, "HJA_master_annual.csv"), show_col_types = FALSE)
site_means <- readr::read_csv(file.path(out_dir, "HJA_master_site_means.csv"), show_col_types = FALSE)
sync_annual <- readr::read_csv(file.path(out_dir, "HJA_composite_synchrony_annual.csv"), show_col_types = FALSE)
clusters_modal <- readr::read_csv(file.path(out_dir, "ClusterStreams_allSolutes_modalClusters.csv"), show_col_types = FALSE)

# Join data
annual_full <- annual %>%
	left_join(sync_annual, by = c("Stream_Name", "solute", "water_year")) %>%
	left_join(
		clusters_modal %>% dplyr::select(Stream_Name, chemical, Cluster_mode),
		by = c("Stream_Name", "solute" = "chemical")
	) %>%
	mutate(
		Stream_Name = as.factor(Stream_Name),
		solute = as.factor(solute),
		solute_type = factor(categorize_solute(solute), levels = c("Geogenic", "Biogenic", "Nutrient")),
		Cluster = factor(Cluster_mode, levels = cluster_levels)
	) %>%
	apply_factor_orders()

# Storage metric
storage_metric <- intersect(PRIMARY_STORAGE_METRIC, names(annual_full))
if (length(storage_metric) == 0) {
	storage_metric <- intersect(c("Q_dS_range_mm", "WB_dS_range_mm"), names(annual_full))
	if (length(storage_metric) == 0) stop("No storage metric found!")
}
storage_metric <- storage_metric[[1]]

# Storage divergence
annual_divergence <- annual_full %>%
	group_by(water_year, solute) %>%
	summarise(
		storage_sd = sd(.data[[storage_metric]], na.rm = TRUE),
		storage_mean = mean(.data[[storage_metric]], na.rm = TRUE),
		n_sites = sum(!is.na(.data[[storage_metric]])),
		.groups = "drop"
	) %>%
	filter(n_sites >= 3, is.finite(storage_sd))

model_data <- annual_full %>%
	left_join(annual_divergence, by = c("water_year", "solute")) %>%
	filter(!is.na(storage_sd), !is.na(Cluster)) %>%
	mutate(
		storage_sd_c = as.numeric(scale(storage_sd, center = TRUE, scale = FALSE)),
		storage_mean_c = as.numeric(scale(storage_mean, center = TRUE, scale = FALSE))
	)

# Model 1: Concentration synchrony
if ("conc_sync_allpairs" %in% names(model_data)) {
	df_conc <- model_data %>% filter(!is.na(conc_sync_allpairs), !is.na(storage_sd_c)) %>% droplevels()
	m1a_null <- lmer(conc_sync_allpairs ~ 1 + (1|water_year) + (1|solute) + (1|Stream_Name), data = df_conc, REML = TRUE)
	m1b_storage <- lmer(conc_sync_allpairs ~ storage_sd_c + (1|water_year) + (1|solute) + (1|Stream_Name), data = df_conc, REML = TRUE)
	m1c_cluster <- lmer(conc_sync_allpairs ~ storage_sd_c + Cluster + (1|water_year) + (1|solute) + (1|Stream_Name), data = df_conc, REML = TRUE)
	m1d_interaction <- lmer(conc_sync_allpairs ~ storage_sd_c * Cluster + (1|water_year) + (1|solute) + (1|Stream_Name), data = df_conc, REML = TRUE)
	anova_conc <- anova(m1a_null, m1b_storage, m1c_cluster, m1d_interaction)
	best_conc <- list(m1a_null, m1b_storage, m1c_cluster, m1d_interaction)[[which.min(anova_conc$AIC)]]
	fixed_conc <- broom.mixed::tidy(best_conc, effects = "fixed") %>% mutate(across(where(is.numeric), ~round(., 4)))
	random_conc <- broom.mixed::tidy(best_conc, effects = "ran_pars") %>% mutate(across(where(is.numeric), ~round(., 4)))
	r2_conc <- MuMIn::r.squaredGLMM(best_conc)
	readr::write_csv(fixed_conc, file.path(out_dir, "05_mixed_effects/model1_conc_sync_fixed_effects.csv"))
	readr::write_csv(random_conc, file.path(out_dir, "05_mixed_effects/model1_conc_sync_random_effects.csv"))
	readr::write_csv(as_tibble(anova_conc, rownames = "model"), file.path(out_dir, "05_mixed_effects/model1_conc_sync_comparison.csv"))
}

# Model 2: CQ-slope synchrony
if ("cqslope_sync_allpairs" %in% names(model_data)) {
	df_cqslope <- model_data %>% filter(!is.na(cqslope_sync_allpairs), !is.na(storage_sd_c)) %>% droplevels()
	m2a_null <- lmer(cqslope_sync_allpairs ~ 1 + (1|water_year) + (1|solute) + (1|Stream_Name), data = df_cqslope, REML = TRUE)
	m2b_storage <- lmer(cqslope_sync_allpairs ~ storage_sd_c + (1|water_year) + (1|solute) + (1|Stream_Name), data = df_cqslope, REML = TRUE)
	m2c_cluster <- lmer(cqslope_sync_allpairs ~ storage_sd_c + Cluster + (1|water_year) + (1|solute) + (1|Stream_Name), data = df_cqslope, REML = TRUE)
	m2d_interaction <- lmer(cqslope_sync_allpairs ~ storage_sd_c * Cluster + (1|water_year) + (1|solute) + (1|Stream_Name), data = df_cqslope, REML = TRUE)
	anova_cqslope <- anova(m2a_null, m2b_storage, m2c_cluster, m2d_interaction)
	best_cqslope <- list(m2a_null, m2b_storage, m2c_cluster, m2d_interaction)[[which.min(anova_cqslope$AIC)]]
	fixed_cqslope <- broom.mixed::tidy(best_cqslope, effects = "fixed") %>% mutate(across(where(is.numeric), ~round(., 4)))
	random_cqslope <- broom.mixed::tidy(best_cqslope, effects = "ran_pars") %>% mutate(across(where(is.numeric), ~round(., 4)))
	r2_cqslope <- MuMIn::r.squaredGLMM(best_cqslope)
	readr::write_csv(fixed_cqslope, file.path(out_dir, "05_mixed_effects/model2_cqslope_sync_fixed_effects.csv"))
	readr::write_csv(random_cqslope, file.path(out_dir, "05_mixed_effects/model2_cqslope_sync_random_effects.csv"))
	readr::write_csv(as_tibble(anova_cqslope, rownames = "model"), file.path(out_dir, "05_mixed_effects/model2_cqslope_sync_comparison.csv"))
}

# Model 3: Storage × solute type
if ("conc_sync_allpairs" %in% names(model_data) && "solute_type" %in% names(model_data)) {
	df_soltype <- model_data %>% filter(!is.na(conc_sync_allpairs), !is.na(storage_sd_c), !is.na(solute_type)) %>% droplevels()
	m3_soltype <- lmer(conc_sync_allpairs ~ storage_sd_c * solute_type + (1|water_year) + (1|solute) + (1|Stream_Name), data = df_soltype, REML = TRUE)
	fixed_soltype <- broom.mixed::tidy(m3_soltype, effects = "fixed") %>% mutate(across(where(is.numeric), ~round(., 4)))
	r2_soltype <- MuMIn::r.squaredGLMM(m3_soltype)
	readr::write_csv(fixed_soltype, file.path(out_dir, "05_mixed_effects/model3_solute_type_fixed_effects.csv"))
}

# Summary table
model_summary <- tibble(
	Model = c("M1: Conc Sync (best)", "M2: CQ-Slope Sync (best)", "M3: Solute Type"),
	Response = c("conc_sync_allpairs", "cqslope_sync_allpairs", "conc_sync_allpairs")
)
readr::write_csv(model_summary, file.path(out_dir, "05_mixed_effects/model_summary_all.csv"))
