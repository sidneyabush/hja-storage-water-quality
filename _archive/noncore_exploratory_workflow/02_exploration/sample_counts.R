# =============================================================================
# Sample Count Visualization
# =============================================================================
# Creates bar plots showing number of individual samples per solute per site

suppressPackageStartupMessages({
  library(tidyverse)
})

# Source helpers
repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
source(file.path(repo_dir, "00_helpers", "plot_prefs.R"))

# Paths
base_dir <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_Water_Quality"
output_dir <- file.path(base_dir, "outputs")
plot_dir <- file.path(base_dir, "exploratory_plots", "01_data_prep", "conc_sample_bias")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Load raw CQ master data (contains individual sample measurements)
cq_master <- read_csv(file.path(output_dir, "HJA_CQ_master.csv"),
                     show_col_types = FALSE)

# Count individual samples per site-solute combination
sample_counts <- cq_master %>%
  group_by(Stream_Name, variable) %>%
  summarise(
    n_samples = n(),  # number of individual sample measurements
    .groups = "drop"
  ) %>%
  rename(site = Stream_Name, solute = variable) %>%
  mutate(
    site = factor(site, levels = site_order),
    solute = factor(solute, levels = solute_order)
  ) %>%
  filter(!is.na(solute))

# Create bar plot: samples by site, faceted by solute
p_by_site <- ggplot(sample_counts, aes(x = site, y = n_samples, fill = site)) +
  geom_col() +
  geom_text(aes(label = n_samples), vjust = -0.3, size = 2.5) +
  facet_wrap(~solute, ncol = 4) +
  scale_fill_site() +
  labs(
    x = "Site",
    y = "Number of Samples"
  ) +
  theme_hja() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

# Create bar plot: samples by solute, faceted by site
p_by_solute <- ggplot(sample_counts, aes(x = solute, y = n_samples, fill = solute)) +
  geom_col() +
  geom_text(aes(label = n_samples), vjust = -0.3, size = 2.5) +
  facet_wrap(~site, ncol = 3) +
  scale_fill_solute() +
  labs(
    x = "Solute",
    y = "Number of Samples"
  ) +
  theme_hja() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

# Create single bar plot with all site-solute combinations
sample_counts_combined <- sample_counts %>%
  mutate(site_solute = paste(site, solute, sep = "_"))

p_combined <- ggplot(sample_counts_combined, aes(x = solute, y = n_samples, fill = site)) +
  geom_col(position = "dodge") +
  scale_fill_site() +
  labs(
    x = "Solute",
    y = "Number of Samples",
    fill = "Site"
  ) +
  theme_hja() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

# Save plots
ggsave(
  file.path(plot_dir, "sample_counts_by_site.png"),
  plot = p_by_site,
  width = 12,
  height = 10,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(plot_dir, "sample_counts_by_solute.png"),
  plot = p_by_solute,
  width = 12,
  height = 12,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(plot_dir, "sample_counts_combined.png"),
  plot = p_combined,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)

# Print summary table
cat("\nSample counts per site-solute combination:\n")
sample_counts %>%
  pivot_wider(names_from = solute, values_from = n_samples) %>%
  print(n = Inf)

message("\n✓ Sample count bar plots saved to: ", plot_dir)
