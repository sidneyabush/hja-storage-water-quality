# =============================================================================
# STEP 03o: SEASON TIMING AS SYNCHRONY PREDICTOR
# =============================================================================
# Goal: Test whether interannual variation in wet season timing predicts sync
#
# RESPONSE VARIABLES:
#   - conc_sync_outlet: Abbott concentration synchrony with outlet
#   - cqslope_sync_outlet: Wymore CQ-slope synchrony with outlet
#   - pairwise_synchrony: Abbott sync between non-outlet site pairs
#
# PREDICTOR VARIABLES:
#   - wet_length_days: Duration of wet season (62-291 days)
#   - wet_start_doy: Day of year wet season begins (DOY 59-356)
#
# KEY FINDING: OPPOSITE effects at different scales!
#   - Outlet: Longer wet seasons → MORE sync with outlet (r = +0.25**)
#   - Pairwise: Longer wet seasons → LESS pairwise sync (r = -0.33***)
#
# DATA LEVELS:
#   - Outlet models: Site-year (N~108), mixed effects with (1|Stream_Name)
#   - Pairwise models: Pair-year (N~623), currently simple correlations
#
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(MuMIn)
  library(broom)
  library(patchwork)
})

rm(list = ls())

# Source helpers
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
source(file.path(repo_dir, "00_helpers", "workflow_config.R"))
try(source(file.path(repo_dir, "00_helpers", "plot_prefs.R")), silent = TRUE)

# Paths
paths <- get_project_paths()
out_dir <- paths$out_dir
data_dir <- file.path(dirname(out_dir), "data")
fig_dir <- file.path(paths$fig_root, "03_stats", "3o_season_timing")
res_dir <- file.path(out_dir, "03_stats")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

message("\n╔════════════════════════════════════════════════════════════════╗\n")
message("║  SEASON TIMING AS SYNCHRONY PREDICTOR                          ║\n")
message("╚════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# 1. LOAD DATA
# =============================================================================
message("=== 1. LOADING DATA ===\n\n")

# Season boundaries with calculated metrics
season_bounds <- read_csv(file.path(out_dir, "season_boundaries.csv"), 
                          show_col_types = FALSE) %>%
  mutate(
    wet_length_days = as.numeric(wet_end_date - wet_start_date),
    wet_start_doy = lubridate::yday(wet_start_date)
  ) %>%
  select(water_year, wet_length_days, wet_start_doy)

message("  Season timing data:", nrow(season_bounds), "water years\n")
message("  wet_length_days: mean =", round(mean(season_bounds$wet_length_days, na.rm=TRUE), 1),
    ", range =", paste(range(season_bounds$wet_length_days, na.rm=TRUE), collapse="-"), "\n")
message("  wet_start_doy: mean =", round(mean(season_bounds$wet_start_doy, na.rm=TRUE), 1),
    ", range =", paste(range(season_bounds$wet_start_doy, na.rm=TRUE), collapse="-"), "\n\n")

# Outlet synchrony (annual level)
outlet_sync_annual <- read_csv(file.path(out_dir, "HJA_outlet_synchrony_annual.csv"), 
                                show_col_types = FALSE)

# Pairwise synchrony (annual level)
sync_pairs_annual <- tryCatch(
  read_csv(file.path(out_dir, "HJA_pair_sync_metrics.csv"), show_col_types = FALSE) %>%
    filter(time_scale == "annual") %>%
    select(solute, water_year, Stream1, Stream2, Abbott_S, prop_sync_wymore, is_outlet_pair) %>%
    rename(site1 = Stream1, site2 = Stream2),
  error = function(e) NULL
)

message("  Outlet sync annual rows:", nrow(outlet_sync_annual), "\n")
if (!is.null(sync_pairs_annual)) {
  message("  Pairwise sync annual rows:", nrow(sync_pairs_annual), "\n")
}

# =============================================================================
# 2. OUTLET-CENTRIC ANALYSIS
# =============================================================================
message("\n=== 2. OUTLET-CENTRIC ANALYSIS ===\n\n")

# Join season timing to outlet sync
outlet_data <- outlet_sync_annual %>%
  left_join(season_bounds, by = "water_year") %>%
  filter(Stream_Name != OUTLET_SITE)

# Aggregate to site-year level
site_year <- outlet_data %>%
  group_by(Stream_Name, water_year, wet_length_days, wet_start_doy) %>%
  summarise(
    conc_sync_outlet = mean(conc_sync_outlet, na.rm = TRUE),
    cqslope_sync_outlet = mean(cqslope_sync_outlet, na.rm = TRUE),
    n_solutes = n(),
    .groups = "drop"
  )

message("  Site-year observations:", nrow(site_year), "\n")
message("  Sites:", n_distinct(site_year$Stream_Name), "\n")
message("  Years:", n_distinct(site_year$water_year), "\n\n")

# --- 2a. Univariate correlations ---
message("--- 2a. Univariate correlations ---\n\n")

cors_outlet <- tibble(
  response = rep(c("conc_sync_outlet", "cqslope_sync_outlet"), each = 2),
  predictor = rep(c("wet_length_days", "wet_start_doy"), 2)
) %>%
  rowwise() %>%
  mutate(
    test = list(cor.test(site_year[[predictor]], site_year[[response]], use = "complete.obs")),
    r = test$estimate,
    p = test$p.value,
    n = sum(is.finite(site_year[[predictor]]) & is.finite(site_year[[response]])),
    sig = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", p < 0.1 ~ ".", TRUE ~ "")
  ) %>%
  select(-test)

message("Outlet sync vs season timing:\n")
print(cors_outlet %>% mutate(r = round(r, 3), p = round(p, 4)))

# --- 2b. Mixed effects models ---
message("\n--- 2b. Mixed effects models (with site random effect) ---\n\n")

# Scale predictors
model_df <- site_year %>%
  mutate(
    wet_length_z = scale(wet_length_days) %>% as.vector(),
    wet_start_z = scale(wet_start_doy) %>% as.vector()
  )

# Concentration sync models
m0_conc <- lmer(conc_sync_outlet ~ 1 + (1|Stream_Name), data = model_df, REML = FALSE)
m1_conc <- lmer(conc_sync_outlet ~ wet_length_z + (1|Stream_Name), data = model_df, REML = FALSE)
m2_conc <- lmer(conc_sync_outlet ~ wet_start_z + (1|Stream_Name), data = model_df, REML = FALSE)
m3_conc <- lmer(conc_sync_outlet ~ wet_length_z + wet_start_z + (1|Stream_Name), data = model_df, REML = FALSE)

aic_conc <- AIC(m0_conc, m1_conc, m2_conc, m3_conc) %>%
  as_tibble(rownames = "model") %>%
  mutate(delta_AIC = AIC - min(AIC)) %>%
  arrange(delta_AIC)

message("Concentration sync model comparison:\n")
print(aic_conc)

# CQ-slope sync models
m0_cq <- lmer(cqslope_sync_outlet ~ 1 + (1|Stream_Name), data = model_df, REML = FALSE)
m1_cq <- lmer(cqslope_sync_outlet ~ wet_length_z + (1|Stream_Name), data = model_df, REML = FALSE)
m2_cq <- lmer(cqslope_sync_outlet ~ wet_start_z + (1|Stream_Name), data = model_df, REML = FALSE)
m3_cq <- lmer(cqslope_sync_outlet ~ wet_length_z + wet_start_z + (1|Stream_Name), data = model_df, REML = FALSE)

aic_cq <- AIC(m0_cq, m1_cq, m2_cq, m3_cq) %>%
  as_tibble(rownames = "model") %>%
  mutate(delta_AIC = AIC - min(AIC)) %>%
  arrange(delta_AIC)

message("\nCQ-slope sync model comparison:\n")
print(aic_cq)

# =============================================================================
# 3. PAIRWISE ANALYSIS
# =============================================================================
message("\n=== 3. PAIRWISE ANALYSIS ===\n\n")

if (!is.null(sync_pairs_annual)) {
  
  # Filter to non-outlet pairs and join season timing
  pair_data <- sync_pairs_annual %>%
    filter(!is_outlet_pair) %>%
    left_join(season_bounds, by = "water_year")
  
  # Aggregate to pair-year level
  pair_year <- pair_data %>%
    group_by(site1, site2, water_year, wet_length_days, wet_start_doy) %>%
    summarise(
      abbott_sync = mean(Abbott_S, na.rm = TRUE),
      wymore_sync = mean(prop_sync_wymore, na.rm = TRUE),
      n_solutes = n(),
      .groups = "drop"
    )
  
  message("  Pair-year observations:", nrow(pair_year), "\n")
  message("  Unique pairs:", n_distinct(paste(pair_year$site1, pair_year$site2)), "\n\n")
  
  # --- 3a. Correlations ---
  message("--- 3a. Pairwise sync vs season timing ---\n\n")
  
  # Abbott sync
  cor_length_abbott <- cor.test(pair_year$wet_length_days, pair_year$abbott_sync, use = "complete.obs")
  cor_start_abbott <- cor.test(pair_year$wet_start_doy, pair_year$abbott_sync, use = "complete.obs")
  
  # Wymore sync
  cor_length_wymore <- cor.test(pair_year$wet_length_days, pair_year$wymore_sync, use = "complete.obs")
  cor_start_wymore <- cor.test(pair_year$wet_start_doy, pair_year$wymore_sync, use = "complete.obs")
  
  cors_pairwise <- tibble(
    response = c("Abbott sync", "Abbott sync", "Wymore sync", "Wymore sync"),
    predictor = c("wet_length_days", "wet_start_doy", "wet_length_days", "wet_start_doy"),
    r = c(cor_length_abbott$estimate, cor_start_abbott$estimate,
          cor_length_wymore$estimate, cor_start_wymore$estimate),
    p = c(cor_length_abbott$p.value, cor_start_abbott$p.value,
          cor_length_wymore$p.value, cor_start_wymore$p.value),
    n = nrow(pair_year),
    sig = ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
  )
  
  message("Pairwise sync vs season timing:\n")
  print(cors_pairwise %>% mutate(r = round(r, 3), p = round(p, 4)))
  
} else {
  message("  No pairwise annual data available.\n")
  cors_pairwise <- NULL
}

# =============================================================================
# 4. KEY FINDING: OPPOSITE EFFECTS
# =============================================================================
message("\n╔════════════════════════════════════════════════════════════════╗\n")
message("║  KEY FINDING: OPPOSITE EFFECTS AT DIFFERENT SCALES             ║\n")
message("╚════════════════════════════════════════════════════════════════╝\n\n")

message("OUTLET-CENTRIC (site → outlet):\n")
message("  wet_length_days: r =", round(cors_outlet$r[cors_outlet$response == "cqslope_sync_outlet" & 
                                                   cors_outlet$predictor == "wet_length_days"], 3), "\n")
message("  Interpretation: LONGER wet seasons → MORE sync with outlet\n\n")

if (!is.null(cors_pairwise)) {
  message("PAIRWISE (site ↔ site):\n")
  message("  wet_length_days: r =", round(cors_pairwise$r[cors_pairwise$response == "Abbott sync" & 
                                                       cors_pairwise$predictor == "wet_length_days"], 3), "\n")
  message("  Interpretation: LONGER wet seasons → LESS pairwise sync\n\n")
}

message("MECHANISTIC HYPOTHESIS:\n")
message("  In LONGER wet seasons:\n")
message("    - Extended connectivity homogenizes signals toward outlet (integration)\n")
message("    - But local heterogeneity has more time to express (divergence among headwaters)\n")
message("  In SHORTER wet seasons:\n")
message("    - Rapid/intense flushing creates coordinated responses among ALL sites\n")

# =============================================================================
# 5. CREATE FIGURES
# =============================================================================
message("\n=== 5. CREATING FIGURES ===\n\n")

# Figure 1: Outlet sync vs season timing
p1 <- site_year %>%
  pivot_longer(cols = c(conc_sync_outlet, cqslope_sync_outlet),
               names_to = "sync_type", values_to = "sync_value") %>%
  mutate(sync_type = recode(sync_type, 
                             "conc_sync_outlet" = "Concentration Sync",
                             "cqslope_sync_outlet" = "CQ-Slope Sync")) %>%
  ggplot(aes(x = wet_length_days, y = sync_value)) +
  geom_point(alpha = 0.4, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "steelblue") +
  facet_wrap(~sync_type, scales = "free_y") +
  labs(
    x = "Wet Season Length (days)",
    y = "Synchrony with Outlet",
    title = "Outlet-Centric: Longer Wet Seasons → MORE Synchrony"
  ) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "gray90"))

ggsave(file.path(fig_dir, "01_outlet_vs_wet_length.png"), p1, width = 10, height = 5, dpi = 300)

# Figure 2: Pairwise sync vs season timing
if (!is.null(cors_pairwise)) {
  p2 <- pair_year %>%
    pivot_longer(cols = c(abbott_sync, wymore_sync),
                 names_to = "sync_type", values_to = "sync_value") %>%
    mutate(sync_type = recode(sync_type,
                               "abbott_sync" = "Abbott Sync",
                               "wymore_sync" = "Wymore Sync")) %>%
    ggplot(aes(x = wet_length_days, y = sync_value)) +
    geom_point(alpha = 0.2, size = 1) +
    geom_smooth(method = "lm", se = TRUE, color = "#B2182B") +
    facet_wrap(~sync_type, scales = "free_y") +
    labs(
      x = "Wet Season Length (days)",
      y = "Pairwise Synchrony",
      title = "Pairwise: Longer Wet Seasons → LESS Synchrony"
    ) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "gray90"))
  
  ggsave(file.path(fig_dir, "02_pairwise_vs_wet_length.png"), p2, width = 10, height = 5, dpi = 300)
}

# Figure 3: Combined comparison
p3_data <- bind_rows(
  site_year %>% 
    select(wet_length_days, sync = cqslope_sync_outlet) %>%
    mutate(scale = "Outlet (site → GSLOOK)"),
  if (!is.null(cors_pairwise)) {
    pair_year %>% 
      select(wet_length_days, sync = abbott_sync) %>%
      mutate(scale = "Pairwise (site ↔ site)")
  } else NULL
)

if (nrow(p3_data) > 0) {
  p3 <- ggplot(p3_data, aes(x = wet_length_days, y = sync, color = scale)) +
    geom_point(alpha = 0.2, size = 1) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1.5) +
    scale_color_manual(values = c("Outlet (site → GSLOOK)" = "steelblue", 
                                   "Pairwise (site ↔ site)" = "#B2182B")) +
    labs(
      x = "Wet Season Length (days)",
      y = "Synchrony",
      color = "Scale",
      title = "Opposite Effects of Season Length on Synchrony",
      subtitle = "Outlet: longer = MORE sync | Pairwise: longer = LESS sync"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  ggsave(file.path(fig_dir, "03_opposite_effects_comparison.png"), p3, width = 8, height = 6, dpi = 300)
}

# =============================================================================
# 6. SAVE RESULTS
# =============================================================================
message("\n=== 6. SAVING RESULTS ===\n\n")

# Combine all correlation results
all_cors <- bind_rows(
  cors_outlet %>% mutate(analysis = "Outlet-centric"),
  if (!is.null(cors_pairwise)) cors_pairwise %>% mutate(analysis = "Pairwise") else NULL
)

write_csv(all_cors, file.path(res_dir, "season_timing_synchrony_correlations.csv"))
write_csv(aic_conc, file.path(res_dir, "season_timing_conc_sync_model_comparison.csv"))
write_csv(aic_cq, file.path(res_dir, "season_timing_cqslope_sync_model_comparison.csv"))

message("Results saved to:", res_dir, "\n")
message("Figures saved to:", fig_dir, "\n")

message("\n=== ANALYSIS COMPLETE ===\n")
