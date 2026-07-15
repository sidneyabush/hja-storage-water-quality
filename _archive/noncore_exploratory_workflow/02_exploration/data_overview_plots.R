# =============================================================================
# Data Overview Plots
# =============================================================================
# Creates descriptive plots showing the range of solute concentrations and
# discharge values across sites - for giving readers an overview of the data

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

# Source helpers
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
source(file.path(repo_dir, "00_helpers", "plot_prefs.R"))

# Paths
base_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality"
output_dir <- file.path(base_dir, "outputs")
plot_dir <- file.path(base_dir, "exploratory_plots", "01_data_prep", "conc_sample_bias")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Load raw CQ master data
cq_master <- read_csv(file.path(output_dir, "HJA_CQ_master.csv"), show_col_types = FALSE)

# =============================================================================
# Prepare data
# =============================================================================

# Apply factor ordering
cq_data <- cq_master %>%
  mutate(
    Stream_Name = factor(Stream_Name, levels = site_order),
    variable = factor(variable, levels = solute_order)
  ) %>%
  filter(!is.na(variable), !is.na(value), !is.na(Q_cms))

# =============================================================================
# PLOT 1: Solute concentration ranges by site (boxplots)
# =============================================================================

p_solutes <- ggplot(cq_data, aes(x = Stream_Name, y = value, fill = Stream_Name)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
  facet_wrap(~variable, scales = "free_y", ncol = 3) +
  scale_fill_site() +
  scale_y_log10() +
  labs(
    x = "Site",
    y = "Concentration (mg/L, log scale)"
  ) +
  theme_hja() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.position = "none",
    strip.text = element_text(size = 14)
  )

# =============================================================================
# PLOT 2: Discharge ranges by site
# =============================================================================

# Get unique discharge values per site
discharge_data <- cq_data %>%
  select(Stream_Name, Date, Q_cms) %>%
  distinct()

p_discharge <- ggplot(discharge_data, aes(x = Stream_Name, y = Q_cms, fill = Stream_Name)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
  scale_fill_site() +
  scale_y_log10() +
  labs(
    x = "Site",
    y = "Discharge (cms, log scale)"
  ) +
  theme_hja() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.position = "none"
  )

# =============================================================================
# PLOT 3: Violin plots for better distribution visualization
# =============================================================================

p_solutes_violin <- ggplot(cq_data, aes(x = Stream_Name, y = value, fill = Stream_Name)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.size = 0.3, outlier.alpha = 0.3) +
  facet_wrap(~variable, scales = "free_y", ncol = 3) +
  scale_fill_site() +
  scale_y_log10() +
  labs(
    x = "Site",
    y = "Concentration (mg/L, log scale)"
  ) +
  theme_hja() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(size = 11)
  )

# =============================================================================
# PLOT 4: Summary statistics table
# =============================================================================

# Calculate summary statistics
summary_stats <- cq_data %>%
  group_by(Stream_Name, variable) %>%
  summarise(
    n = n(),
    median = median(value, na.rm = TRUE),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  )

# Save summary table
write_csv(summary_stats, file.path(plot_dir, "solute_summary_statistics.csv"))

discharge_stats <- discharge_data %>%
  group_by(Stream_Name) %>%
  summarise(
    n = n(),
    median = median(Q_cms, na.rm = TRUE),
    q25 = quantile(Q_cms, 0.25, na.rm = TRUE),
    q75 = quantile(Q_cms, 0.75, na.rm = TRUE),
    min = min(Q_cms, na.rm = TRUE),
    max = max(Q_cms, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(discharge_stats, file.path(plot_dir, "discharge_summary_statistics.csv"))

# =============================================================================
# PLOT 5: Solute ranges across all sites (compare solutes)
# =============================================================================

p_solute_comparison <- ggplot(cq_data, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
  scale_fill_solute() +
  scale_y_log10() +
  labs(
    x = "Solute",
    y = "Concentration (mg/L, log scale)"
  ) +
  theme_hja() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.position = "none"
  )

# =============================================================================
# COMBINED PLOT: Solutes + Discharge
# =============================================================================

# Create combined plot with solutes on top and discharge on bottom
p_combined_overview <- p_solutes / p_discharge +
  plot_layout(heights = c(3, 1))

# =============================================================================
# COMBINED PLOT: Solute Ranges by Site + Discharge by Site
# =============================================================================

# Create combined plot with solute ranges by site on top and discharge by site on bottom
p_solute_discharge_comparison <- p_solutes / p_discharge +
  plot_layout(heights = c(3, 1))

# =============================================================================
# PLOT 6: Facet by Site with Discharge as a Variable
# =============================================================================

# Combine solute data with discharge data (treat Q as another "variable")
combined_data <- bind_rows(
  # Solute concentrations
  cq_data %>%
    select(Stream_Name, Date, variable, value) %>%
    rename(Variable = variable, Value = value),
  # Discharge treated as a variable
  discharge_data %>%
    select(Stream_Name, Date, Q_cms) %>%
    mutate(Variable = "Q (cms)") %>%
    rename(Value = Q_cms)
) %>%
  mutate(
    Stream_Name = factor(Stream_Name, levels = site_order),
    # Create ordered variable factor with Q at the end
    Variable = factor(Variable, levels = c(solute_order, "Q (cms)"))
  )

# Create plot faceted by site
p_by_site_with_discharge <- ggplot(combined_data, aes(x = Variable, y = Value, fill = Variable)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
  facet_wrap(~Stream_Name, ncol = 3) +
  scale_y_log10() +
  scale_fill_manual(
    values = c(solute_colors, "Q (cms)" = "#D3D3D3"),  # Use solute_colors from plot_prefs + light gray for Q
    breaks = c(solute_order, "Q (cms)")
  ) +
  labs(
    x = "Variable",
    y = "Value (mg/L for solutes, cms for Q; log scale)",
    fill = "Variable"
  ) +
  theme_hja() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 10)
  )

# =============================================================================
# Save all plots
# =============================================================================

ggsave(
  file.path(plot_dir, "solute_ranges_by_site_boxplot.png"),
  plot = p_solutes,
  width = 12,
  height = 14,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(plot_dir, "discharge_ranges_by_site.png"),
  plot = p_discharge,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(plot_dir, "solute_ranges_by_site_violin.png"),
  plot = p_solutes_violin,
  width = 12,
  height = 14,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(plot_dir, "solute_comparison_across_sites.png"),
  plot = p_solute_comparison,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)

# Save combined plot
ggsave(
  file.path(plot_dir, "combined_solutes_discharge_overview.png"),
  plot = p_combined_overview,
  width = 12,
  height = 16,
  dpi = 300,
  bg = "white"
)

# Save solute ranges by site + discharge combined plot
ggsave(
  file.path(plot_dir, "solute_discharge_comparison.png"),
  plot = p_solute_discharge_comparison,
  width = 12,
  height = 16,
  dpi = 300,
  bg = "white"
)

# Save by-site plot with discharge
ggsave(
  file.path(plot_dir, "chemical_signatures_by_site.png"),
  plot = p_by_site_with_discharge,
  width = 14,
  height = 12,
  dpi = 300,
  bg = "white"
)

# Print summary
cat("\n=== Summary Statistics ===\n")
cat("\nSolute concentration ranges (median [Q25-Q75]):\n")
summary_stats %>%
  mutate(range_str = sprintf("%.2f [%.2f-%.2f]", median, q25, q75)) %>%
  select(Stream_Name, variable, range_str) %>%
  pivot_wider(names_from = variable, values_from = range_str) %>%
  print(n = Inf)

cat("\n\nDischarge ranges (median [Q25-Q75] cms):\n")
discharge_stats %>%
  mutate(range_str = sprintf("%.3f [%.3f-%.3f]", median, q25, q75)) %>%
  select(Stream_Name, range_str) %>%
  print(n = Inf)

message("\n✓ Data overview plots and summary statistics saved")
message("  - Chemical signatures by site: chemical_signatures_by_site.png")
message("  - Combined overview: combined_solutes_discharge_overview.png")
message("  - Solute/discharge comparison: solute_discharge_comparison.png")
message("  - Boxplots: solute_ranges_by_site_boxplot.png")
message("  - Violin plots: solute_ranges_by_site_violin.png")
message("  - Discharge: discharge_ranges_by_site.png")
message("  - Solute comparison: solute_comparison_across_sites.png")
message("  - Summary tables: *_summary_statistics.csv")
