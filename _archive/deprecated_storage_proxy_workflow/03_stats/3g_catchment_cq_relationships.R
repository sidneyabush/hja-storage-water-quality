# =============================================================================
# 3g: Catchment Characteristics → CQ Behavior Relationships
# =============================================================================
# Comprehensive analysis of how catchment characteristics (geology, topography,
# land use, hydrology) relate to concentration-discharge behavior and synchrony
#
# Key findings integrated:
# - DSi CQ slope strongly related to Lava1 percentage (r=-0.89)
# - Ca, Mg CQ slopes related to storage (Q_dS_range_mm) and Ash percentage
# - PO4 strongly related to elevation, pyroclastic, RBI (flashiness)
# - Elevation strongly controls damping ratio (DR) and flashiness (RBI)
# =============================================================================

library(tidyverse)
library(corrplot)
library(ggrepel)
library(patchwork)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
base_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality"
out_dir <- file.path(base_dir, "outputs")
plot_dir <- file.path(base_dir, "exploratory_plots", "03_stats", "3g_catchment_cq")

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "03_stats"), showWarnings = FALSE, recursive = TRUE)

# Load plot preferences if available
if (file.exists("00_helpers/plot_prefs.R")) source("00_helpers/plot_prefs.R")

# Custom theme
theme_clean <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray95")
  )

# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------
message("Loading data...\n")

# Site-level CQ means with integrated catchment characteristics
site_cq <- read_csv(file.path(out_dir, "HJA_clean_site_means.csv"), show_col_types = FALSE)

# Full catchment characteristics (now generated in 01_data_prep/1h)
catchment_chars <- read_csv(file.path(out_dir, "HJA_master_site_means.csv"),
                            show_col_types = FALSE) %>%
  select(-any_of("...1"))

# DS Drawdown annual
ds_drawdown <- read_csv(file.path(base_dir, "data/DS_drawdown_annual.csv"), 
                        show_col_types = FALSE) %>%
  select(-any_of("...1")) %>%
  rename(site = SITECODE, year = waterYear)

# Storage metrics annual
storage_annual <- read_csv(file.path(base_dir, "data/HJA_StorageMetrics_Annual.csv"),
                           show_col_types = FALSE) %>%
  select(-any_of("...1"))

message("Data loaded successfully\n")
message("Sites in site_cq:", length(unique(site_cq$Stream_Name)), "\n")
message("Solutes:", unique(site_cq$solute), "\n")

# =============================================================================
# 1. Catchment Characteristics Summary
# =============================================================================
message("\n=== CATCHMENT CHARACTERISTICS SUMMARY ===\n")

# Key characteristics to analyze
geology_vars <- c("Lava1_per", "Lava2_per", "Ash_Per", "Pyro_per")
topo_vars <- c("Elevation_mean_m", "Slope_mean", "Area_km2")
landuse_vars <- c("Harvest", "Age", "Landslide_Total")
hydro_vars <- c("DR_Overall", "RBI_mean", "DS_sum_mean", "recession_curve_slope_mean", 
                "Q5norm_mean", "fdc_slope_mean")

# Summary table
char_summary <- catchment_chars %>%
  select(site, all_of(c(geology_vars, topo_vars, landuse_vars))) %>%
  pivot_longer(-site, names_to = "characteristic", values_to = "value") %>%
  group_by(characteristic) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    range = max - min
  )

print(char_summary)

# =============================================================================
# 2. CQ Slope ~ Catchment Characteristic Correlations
# =============================================================================
message("\n=== CQ SLOPE ~ CATCHMENT CORRELATIONS ===\n")

# Define characteristics to test
test_chars <- c("Elevation_mean_m", "Slope_mean", "Area_km2",
                "Lava1_per", "Lava2_per", "Ash_Per", "Pyro_per",
                "Harvest", "Landslide_Total", "DR_Overall", 
                "RBI", "Q_dS_range_mm", "MTT_final", "FYw_final")

# Compute correlations for each solute
solutes <- unique(site_cq$solute)

cq_char_cors <- map_dfr(solutes, function(sol) {
  df <- site_cq %>% 
    filter(solute == sol) %>% 
    filter(!is.na(cq_slope))
  
  if (nrow(df) < 5) return(NULL)
  
  map_dfr(test_chars, function(ch) {
    if (!ch %in% names(df)) return(NULL)
    x <- df[[ch]]
    y <- df$cq_slope
    valid <- !is.na(x) & !is.na(y)
    if (sum(valid) < 5) return(NULL)
    
    test <- cor.test(x[valid], y[valid])
    tibble(
      solute = sol,
      characteristic = ch,
      r = test$estimate,
      p = test$p.value,
      n = sum(valid),
      sig = ifelse(test$p.value < 0.05, "*", "")
    )
  })
})

# Summary table
message("\n--- Significant CQ slope ~ Catchment correlations (p < 0.05) ---\n")
sig_cors <- cq_char_cors %>%
  filter(p < 0.05) %>%
  arrange(p) %>%
  mutate(r = round(r, 3), p = round(p, 4))

print(sig_cors, n = 30)

# Create correlation heatmap matrix
cor_matrix <- cq_char_cors %>%
  select(solute, characteristic, r) %>%
  pivot_wider(names_from = characteristic, values_from = r) %>%
  column_to_rownames("solute") %>%
  as.matrix()

# Plot correlation heatmap
png(file.path(plot_dir, "cq_slope_catchment_correlation_heatmap.png"),
    width = 14, height = 8, units = "in", res = 300)

# Handle NA for corrplot
cor_matrix[is.na(cor_matrix)] <- 0

corrplot(cor_matrix, 
         method = "color",
         type = "full",
         tl.col = "black",
         tl.cex = 0.8,
         cl.cex = 0.8,
         addCoef.col = "black",
         number.cex = 0.6,
         col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         title = "CQ Slope ~ Catchment Characteristics Correlations",
         mar = c(0, 0, 2, 0))
dev.off()
message("Saved: cq_slope_catchment_correlation_heatmap.png\n")

# =============================================================================
# 3. Key Relationship Visualizations
# =============================================================================
message("\n=== GENERATING KEY RELATIONSHIP PLOTS ===\n")

# 3a. DSi vs Lava1 (strongest geology relationship)
p_dsi_lava <- site_cq %>%
  filter(solute == "DSi") %>%
  ggplot(aes(x = Lava1_per, y = cq_slope)) +
  geom_point(aes(color = Stream_Name), size = 4) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  geom_text_repel(aes(label = Stream_Name), size = 3) +
  scale_color_viridis_d() +
  labs(
    title = "DSi CQ Slope vs Lava-1 Percentage",
    subtitle = "r = -0.89, p = 0.001 | Basaltic bedrock controls silica release",
    x = "Lava-1 Percentage (%)",
    y = "DSi CQ Slope",
    color = "Site"
  ) +
  theme_clean +
  theme(legend.position = "none")

# 3b. Ca vs Storage (Q_dS_range_mm)
p_ca_storage <- site_cq %>%
  filter(solute == "Ca") %>%
  ggplot(aes(x = Q_dS_range_mm, y = cq_slope)) +
  geom_point(aes(color = Stream_Name), size = 4) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  geom_text_repel(aes(label = Stream_Name), size = 3) +
  scale_color_viridis_d() +
  labs(
    title = "Ca CQ Slope vs Storage Range",
    subtitle = "r = 0.88, p = 0.002 | Higher storage = more Ca enrichment",
    x = "Q-derived Storage Range (mm)",
    y = "Ca CQ Slope",
    color = "Site"
  ) +
  theme_clean +
  theme(legend.position = "none")

# 3c. PO4 vs Elevation
p_po4_elev <- site_cq %>%
  filter(solute == "PO4") %>%
  ggplot(aes(x = Elevation_mean_m, y = cq_slope)) +
  geom_point(aes(color = Stream_Name), size = 4) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  geom_text_repel(aes(label = Stream_Name), size = 3) +
  scale_color_viridis_d() +
  labs(
    title = "PO4 CQ Slope vs Mean Elevation",
    subtitle = "r = -0.85, p = 0.004 | Lower elevation = more PO4 enrichment",
    x = "Mean Elevation (m)",
    y = "PO4 CQ Slope",
    color = "Site"
  ) +
  theme_clean +
  theme(legend.position = "none")

# 3d. PO4 vs RBI (flashiness)
p_po4_rbi <- site_cq %>%
  filter(solute == "PO4") %>%
  ggplot(aes(x = RBI, y = cq_slope)) +
  geom_point(aes(color = Stream_Name), size = 4) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  geom_text_repel(aes(label = Stream_Name), size = 3) +
  scale_color_viridis_d() +
  labs(
    title = "PO4 CQ Slope vs Flashiness (RBI)",
    subtitle = "r = 0.87, p = 0.002 | Flashier = more PO4 enrichment",
    x = "Richards-Baker Index (RBI)",
    y = "PO4 CQ Slope",
    color = "Site"
  ) +
  theme_clean +
  theme(legend.position = "none")

# Combine key relationships
p_key_relationships <- (p_dsi_lava | p_ca_storage) / (p_po4_elev | p_po4_rbi) +
  plot_annotation(
    title = "Key CQ Slope ~ Catchment Characteristic Relationships",
    subtitle = "Strong controls: Geology (DSi), Storage (Ca, Mg), Topography (PO4)"
  )

ggsave(file.path(plot_dir, "cq_catchment_key_relationships.png"),
       p_key_relationships, width = 14, height = 12, dpi = 300)
message("Saved: cq_catchment_key_relationships.png\n")

# =============================================================================
# 4. Geology Controls on CQ Behavior
# =============================================================================
message("\n=== GEOLOGY CONTROLS ANALYSIS ===\n")

# Create geology summary for each site
geology_summary <- catchment_chars %>%
  select(site, Lava1_per, Lava2_per, Ash_Per, Pyro_per) %>%
  mutate(
    dominant_geology = case_when(
      Lava1_per >= 30 ~ "Lava-1 dominated",
      Lava2_per >= 30 ~ "Lava-2 dominated",
      Ash_Per >= 30 ~ "Ash dominated",
      Pyro_per >= 30 ~ "Pyroclastic dominated",
      TRUE ~ "Mixed"
    )
  )

print(geology_summary)

# Plot CQ slopes by geology type
p_geology_cq <- site_cq %>%
  left_join(geology_summary, by = c("Stream_Name" = "site")) %>%
  filter(!is.na(dominant_geology)) %>%
  ggplot(aes(x = dominant_geology, y = cq_slope, fill = dominant_geology)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~solute, scales = "free_y") +
  scale_fill_viridis_d() +
  labs(
    title = "CQ Slopes by Dominant Geology Type",
    x = "Dominant Geology",
    y = "CQ Slope"
  ) +
  theme_clean +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave(file.path(plot_dir, "cq_slope_by_geology.png"),
       p_geology_cq, width = 14, height = 10, dpi = 300)
message("Saved: cq_slope_by_geology.png\n")

# =============================================================================
# 5. DS Drawdown Analysis
# =============================================================================
message("\n=== DS DRAWDOWN (DRY SEASON STORAGE DEPLETION) ANALYSIS ===\n")

# Annual DS drawdown summary
ds_annual_summary <- ds_drawdown %>%
  group_by(site) %>%
  summarise(
    DS_mean = mean(DS_sum, na.rm = TRUE),
    DS_sd = sd(DS_sum, na.rm = TRUE),
    DS_cv = DS_sd / abs(DS_mean),
    n_years = n()
  ) %>%
  arrange(DS_mean)

print(ds_annual_summary)

# Join with catchment characteristics
ds_with_chars <- ds_annual_summary %>%
  left_join(catchment_chars, by = "site")

# Correlations of DS drawdown with catchment characteristics
message("\n--- DS Drawdown ~ Catchment Correlations ---\n")
ds_cors <- map_dfr(c("Elevation_mean_m", "Slope_mean", "Area_km2", 
                     "Lava1_per", "Lava2_per", "Ash_Per", "Pyro_per",
                     "Harvest", "Landslide_Total", "DR_Overall", "RBI_mean"), function(ch) {
  if (!ch %in% names(ds_with_chars)) return(NULL)
  test <- cor.test(ds_with_chars$DS_mean, ds_with_chars[[ch]])
  tibble(
    characteristic = ch,
    r = round(test$estimate, 3),
    p = round(test$p.value, 4)
  )
})
print(ds_cors)

# Plot DS drawdown vs elevation (shows storage capacity relationship)
p_ds_elev <- ds_with_chars %>%
  ggplot(aes(x = Elevation_mean_m, y = DS_mean)) +
  geom_point(aes(color = site), size = 4) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  geom_text_repel(aes(label = site), size = 3) +
  scale_color_viridis_d() +
  labs(
    title = "Dry Season Drawdown vs Mean Elevation",
    subtitle = "Higher elevation catchments show larger storage depletion",
    x = "Mean Elevation (m)",
    y = "DS Drawdown (mm, negative = depletion)"
  ) +
  theme_clean +
  theme(legend.position = "none")

ggsave(file.path(plot_dir, "ds_drawdown_vs_elevation.png"),
       p_ds_elev, width = 8, height = 6, dpi = 300)
message("Saved: ds_drawdown_vs_elevation.png\n")

# =============================================================================
# 6. Damping Ratio and Flashiness Controls
# =============================================================================
message("\n=== DAMPING RATIO & FLASHINESS CONTROLS ===\n")

# DR and RBI are key hydrologic metrics - what controls them?
dr_rbi_cors <- catchment_chars %>%
  select(site, DR_Overall, RBI_mean, 
         Elevation_mean_m, Slope_mean, Area_km2,
         Lava1_per, Lava2_per, Ash_Per, Pyro_per,
         Harvest, Landslide_Total) %>%
  drop_na()

# Correlations
message("\n--- What controls Damping Ratio (DR)? ---\n")
for (ch in c("Elevation_mean_m", "Slope_mean", "Lava1_per", "Lava2_per", 
             "Ash_Per", "Pyro_per", "Landslide_Total")) {
  if (ch %in% names(dr_rbi_cors)) {
    test <- cor.test(dr_rbi_cors$DR_Overall, dr_rbi_cors[[ch]])
    message(sprintf("%s: r = %.3f, p = %.4f\n", ch, test$estimate, test$p.value))
  }
}

message("\n--- What controls Flashiness (RBI)? ---\n")
for (ch in c("Elevation_mean_m", "Slope_mean", "Lava1_per", "Lava2_per", 
             "Ash_Per", "Pyro_per", "Landslide_Total")) {
  if (ch %in% names(dr_rbi_cors)) {
    test <- cor.test(dr_rbi_cors$RBI_mean, dr_rbi_cors[[ch]])
    message(sprintf("%s: r = %.3f, p = %.4f\n", ch, test$estimate, test$p.value))
  }
}

# Plot DR vs Elevation and Pyro
p_dr_controls <- dr_rbi_cors %>%
  pivot_longer(c(Elevation_mean_m, Pyro_per, Slope_mean), 
               names_to = "predictor", values_to = "value") %>%
  ggplot(aes(x = value, y = DR_Overall)) +
  geom_point(aes(color = site), size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  facet_wrap(~predictor, scales = "free_x") +
  scale_color_viridis_d() +
  labs(
    title = "Damping Ratio Controls",
    subtitle = "DR_Overall increases with pyroclastic content and slope",
    y = "Damping Ratio (DR)"
  ) +
  theme_clean

ggsave(file.path(plot_dir, "damping_ratio_controls.png"),
       p_dr_controls, width = 12, height = 5, dpi = 300)
message("Saved: damping_ratio_controls.png\n")

# =============================================================================
# 6B. SOLUTE-TYPE STRATIFIED ANALYSIS (Geogenic vs Biogenic vs Nutrient)
# =============================================================================
message("\n=== SOLUTE-TYPE STRATIFIED ANALYSIS ===\n")

# Define solute categories (matching plot_prefs.R)
GEOGENIC_SOLUTES <- c("Ca", "Mg", "Na", "K", "Cl", "SO4")
BIOGENIC_SOLUTES <- c("NO3", "PO4", "NH3", "DOC")
NUTRIENT_SOLUTES <- c("DSi")  # DSi behaves uniquely

# Add solute type to data
site_cq_typed <- site_cq %>%
  mutate(solute_type = case_when(
    solute %in% GEOGENIC_SOLUTES ~ "Geogenic",
    solute %in% BIOGENIC_SOLUTES ~ "Biogenic",
    solute %in% NUTRIENT_SOLUTES ~ "Nutrient (DSi)",
    TRUE ~ "Other"
  ))

# Compute correlations stratified by solute type
cq_char_cors_bytype <- site_cq_typed %>%
  group_by(solute_type) %>%
  group_modify(function(df, key) {
    map_dfr(test_chars, function(ch) {
      if (!ch %in% names(df)) return(NULL)
      x <- df[[ch]]
      y <- df$cq_slope
      valid <- !is.na(x) & !is.na(y)
      if (sum(valid) < 5) return(NULL)
      
      test <- cor.test(x[valid], y[valid])
      tibble(
        characteristic = ch,
        r = test$estimate,
        p = test$p.value,
        n = sum(valid),
        sig = ifelse(test$p.value < 0.05, "*", "")
      )
    })
  }) %>%
  ungroup()

message("\n--- Significant correlations by solute type (p < 0.05) ---\n")
sig_bytype <- cq_char_cors_bytype %>%
  filter(p < 0.05) %>%
  arrange(solute_type, p) %>%
  mutate(r = round(r, 3), p = round(p, 4))

print(sig_bytype, n = 40)

# Create heatmap for each solute type
for (stype in c("Geogenic", "Biogenic", "Nutrient (DSi)")) {
  type_cors <- cq_char_cors %>%
    filter(solute %in% switch(stype,
      "Geogenic" = GEOGENIC_SOLUTES,
      "Biogenic" = BIOGENIC_SOLUTES,
      "Nutrient (DSi)" = NUTRIENT_SOLUTES
    ))
  
  if (nrow(type_cors) < 2) next
  
  cor_mat <- type_cors %>%
    select(solute, characteristic, r) %>%
    pivot_wider(names_from = characteristic, values_from = r) %>%
    column_to_rownames("solute") %>%
    as.matrix()
  
  cor_mat[is.na(cor_mat)] <- 0
  
  fname <- paste0("cq_catchment_heatmap_", gsub(" \\(DSi\\)|[^a-zA-Z]", "", stype), ".png")
  png(file.path(plot_dir, fname), width = 12, height = 6, units = "in", res = 300)
  corrplot(cor_mat,
           method = "color",
           type = "full",
           tl.col = "black",
           tl.cex = 0.9,
           addCoef.col = "black",
           number.cex = 0.7,
           col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           title = paste(stype, "Solutes: CQ Slope ~ Catchment Correlations"),
           mar = c(0, 0, 2, 0))
  dev.off()
  message("Saved:", fname, "\n")
}

# Combined boxplot: CQ slopes by solute type
p_cq_bytype <- site_cq_typed %>%
  filter(solute_type != "Other") %>%
  ggplot(aes(x = solute_type, y = cq_slope, fill = solute_type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = solute), width = 0.2, alpha = 0.7, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = c("Geogenic" = "#1f78b4", "Biogenic" = "#33a02c", "Nutrient (DSi)" = "#ff7f00")) +
  labs(
    title = "CQ Slopes by Solute Type",
    subtitle = "Geogenic solutes tend toward dilution; Biogenic toward enrichment",
    x = "Solute Type",
    y = "Mean CQ Slope",
    fill = "Type"
  ) +
  theme_clean +
  theme(legend.position = "right")

ggsave(file.path(plot_dir, "cq_slopes_by_solute_type.png"), p_cq_bytype, 
       width = 10, height = 7, dpi = 300)
message("Saved: cq_slopes_by_solute_type.png\n")

# Catchment control comparison between geogenic and biogenic
p_type_comparison <- site_cq_typed %>%
  filter(solute_type %in% c("Geogenic", "Biogenic")) %>%
  select(Stream_Name, solute, solute_type, cq_slope, 
         Elevation_mean_m, Q_dS_range_mm, RBI, Lava1_per) %>%
  pivot_longer(cols = c(Elevation_mean_m, Q_dS_range_mm, RBI, Lava1_per),
               names_to = "predictor", values_to = "value") %>%
  ggplot(aes(x = value, y = cq_slope, color = solute_type)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  facet_wrap(~predictor, scales = "free_x") +
  scale_color_manual(values = c("Geogenic" = "#1f78b4", "Biogenic" = "#33a02c")) +
  labs(
    title = "Geogenic vs Biogenic: Different Catchment Controls",
    subtitle = "Geogenic controlled by geology; Biogenic by topography/hydrology",
    y = "CQ Slope",
    color = "Solute Type"
  ) +
  theme_clean

ggsave(file.path(plot_dir, "geogenic_vs_biogenic_controls.png"), p_type_comparison,
       width = 12, height = 8, dpi = 300)
message("Saved: geogenic_vs_biogenic_controls.png\n")

# =============================================================================
# 8. Summary: Key Catchment-CQ Relationships
# =============================================================================
message("\n" , strrep("=", 60), "\n")
message("SUMMARY: KEY CATCHMENT-CQ RELATIONSHIPS\n")
message(strrep("=", 60), "\n")

message("\n1. GEOLOGY CONTROLS:\n")
message("   - DSi CQ slope strongly controlled by Lava-1 % (r = -0.89, p = 0.001)\n")
message("   - K CQ slope controlled by Lava-1 % (r = -0.80, p = 0.009)\n")
message("   - Ca, Mg, Na CQ slopes correlate with Ash % (r = 0.75-0.85)\n")
message("   - PO4 strongly related to Pyroclastic % (r = 0.86, p = 0.003)\n")

message("\n2. TOPOGRAPHY CONTROLS:\n")
message("   - Elevation controls Damping Ratio (r = -0.86) and RBI (r = -0.86)\n")
message("   - Slope controls RBI (r = 0.72) - steeper = flashier\n")
message("   - PO4 CQ slope inversely related to elevation (r = -0.85)\n")
message("   - Ca CQ slope inversely related to slope (r = -0.81)\n")

message("\n3. STORAGE/HYDROLOGY CONTROLS:\n")
message("   - Ca, Mg CQ slopes strongly related to Q_dS storage (r = 0.82-0.88)\n")
message("   - PO4 related to RBI (r = 0.87) - flashier = more enrichment\n")
message("   - DOC related to RBI (r = 0.72)\n")

message("\n4. DS DRAWDOWN (Dry Season Storage Depletion):\n")
message("   - Ranges from -502 mm (GSWS02) to -425 mm (GSWS07)\n")
message("   - Moderately correlated with elevation (r = 0.47) and slope (r = -0.49)\n")
message("   - Higher elevation = more drawdown (larger storage reservoir)\n")

message("\n5. SOLUTE TYPE PATTERNS:\n")
message("   - GEOGENIC (Ca, Mg, Na, K, Cl, SO4): Primarily controlled by geology\n")
message("     → Lava-1, Ash, and storage (Q_dS) are key predictors\n")
message("     → Tend toward dilution behavior (negative CQ slopes)\n")
message("   - BIOGENIC (NO3, PO4, NH3, DOC): Primarily controlled by topography/hydrology\n")
message("     → Elevation, RBI, and pyroclastic content are key predictors\n")
message("     → Tend toward enrichment behavior (positive CQ slopes)\n")
message("   - NUTRIENT (DSi): Unique behavior - geology-controlled but nutrient-like\n")
message("     → Strongly controlled by Lava-1 (basaltic weathering)\n")

message("\n6. KEY MECHANISTIC INTERPRETATIONS:\n")
message("   - Basaltic (Lava-1) bedrock releases silica slowly → dilution behavior\n")
message("   - Ash layers store and release Ca, Mg, Na → enrichment at high flows\n")
message("   - Pyroclastic terrain = rapid flow paths → flashy, enrichment behavior\n")
message("   - High elevation = deep storage, buffered chemistry, lower synchrony\n")
message("   - Low elevation = shallow flow paths, rapid response, higher synchrony\n")

# =============================================================================
# 9. Save Summary Statistics
# =============================================================================

# Save correlation results
write_csv(cq_char_cors, file.path(out_dir, "CQ_catchment_correlations.csv"))
message("\nSaved: CQ_catchment_correlations.csv\n")

# Save significant correlations
write_csv(sig_cors, file.path(out_dir, "CQ_catchment_correlations_significant.csv"))
message("Saved: CQ_catchment_correlations_significant.csv\n")

# Save solute-type stratified correlations
write_csv(cq_char_cors_bytype, file.path(out_dir, "CQ_catchment_correlations_byType.csv"))
message("Saved: CQ_catchment_correlations_byType.csv\n")

write_csv(sig_bytype, file.path(out_dir, "CQ_catchment_correlations_byType_significant.csv"))
message("Saved: CQ_catchment_correlations_byType_significant.csv\n")

# Save DS drawdown summary
write_csv(ds_annual_summary, file.path(out_dir, "DS_drawdown_summary.csv"))
message("Saved: DS_drawdown_summary.csv\n")

message("\n=== Analysis complete! ===\n")
