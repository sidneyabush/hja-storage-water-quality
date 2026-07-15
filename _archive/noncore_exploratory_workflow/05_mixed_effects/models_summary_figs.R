suppressPackageStartupMessages({
	library(tidyverse)
	library(ggrepel)
})

repo_dir <- "/Users/sidneybush/Documents/GitHub/hja-water-quality"
try(source(file.path(repo_dir, "00_helpers", "plot_prefs.R")), silent = TRUE)
try(source(file.path(repo_dir, "00_helpers", "plot_theme_set.R")), silent = TRUE)

base_dir    <- "/Users/sidneybush/Library/CloudStorage/Box-Box"
project_dir <- file.path(base_dir, "Sidney_Bush", "HJA_Water_Quality")
out_stats   <- file.path(project_dir, "outputs", "05_mixed_effects")
fig_root    <- file.path(project_dir, "exploratory_plots", "05_mixed_effects")
dir.create(fig_root, showWarnings = FALSE, recursive = TRUE)

summary_file <- file.path(out_stats, "model_summary_all.csv")
model1_file  <- file.path(out_stats, "model1_conc_sync_comparison.csv")
model2_file  <- file.path(out_stats, "model2_cqslope_sync_comparison.csv")

if (!file.exists(summary_file)) stop("model_summary_all.csv not found in outputs/05_mixed_effects")

summary_all <- readr::read_csv(summary_file, show_col_types = FALSE)
model1_cmp  <- if (file.exists(model1_file)) readr::read_csv(model1_file, show_col_types = FALSE) else NULL
model2_cmp  <- if (file.exists(model2_file)) readr::read_csv(model2_file, show_col_types = FALSE) else NULL

top_by_response <- summary_all %>%
	filter(!is.na(Response)) %>%
	group_by(Response) %>%
	arrange(AIC, .by_group = TRUE) %>%
	slice(1) %>%
	ungroup()

readr::write_csv(top_by_response, file.path(out_stats, "model_top_by_response.csv"))

figA <- top_by_response %>%
	mutate(Model = forcats::fct_inorder(Model)) %>%
	ggplot(aes(x = Model, y = -AIC, fill = Response)) +
	geom_col() +
	geom_text(aes(label = paste0("R2m=", sprintf("%.2f", R2_marginal))), vjust = -0.5, size = 3) +
	labs(title = "Top Mixed-Effects Models per Response",
			 subtitle = "Bars show -AIC (higher is better); labels show marginal R2",
			 x = "Model", y = "-AIC") +
	theme_hja() +
	theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggplot2::ggsave(file.path(fig_root, "figA_top_models_per_response.png"), figA, width = 10, height = 6, dpi = 300)

tbl_text <- top_by_response %>%
	mutate(line = sprintf("%s: %s (AIC=%.0f, R2m=%.2f)", Response, Model, AIC, R2_marginal)) %>%
	mutate(idx = dplyr::row_number())

figB <- ggplot(tbl_text, aes(y = rev(idx), x = 0, label = line)) +
	geom_text(hjust = 0, family = "mono") +
	labs(title = "Top Models Summary", x = NULL, y = NULL) +
	theme_void() +
	xlim(0, 1)

ggplot2::ggsave(file.path(fig_root, "figB_top_models_summary_text.png"), figB, width = 8, height = 6, dpi = 300)

message("Figures written to:", fig_root, "\n")
