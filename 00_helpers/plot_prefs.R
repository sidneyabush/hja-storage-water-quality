# =============================================================================
# plot_prefs.R - Central plotting preferences for HJA Stream Chemistry project
# =============================================================================

library(ggplot2)

# =============================================================================
# PLOT QUALITY SETTINGS
# =============================================================================
PLOT_DPI <- 300
PLOT_WIDTH <- 10
PLOT_HEIGHT <- 7
BASE_SIZE <- 13

# =============================================================================
# VARIABLE LABELS (for axis labels)
# =============================================================================
var_labels <- list(
  RCS_p = "Recession Constant (RCS_p)\nhigher = slower drainage",
  dynamic_storage_strength_z = "Dynamic storage strength (z)",
  mobile_mixing_no_bf_z = "Mobile mixing, tracer-only (z)",
  mobile_mixing_with_bf_z = "Mobile mixing + BF (z)",
  flow_path_partitioning_z = "Flow-path partitioning, BF (z)",
  unified_state_index = "Unified storage state index",
  geology_pc1 = "Geology/landslide PC1",
  geology_pc2 = "Geology/landslide PC2",
  storage_mean = "Mean Recession Constant (RCS_p)",
  D_storage = "Storage Divergence |delta RCS_p|",
  conc_sync = "Absolute subcatchment synchrony (Abbott)",
  conc_sync_outlet = "Outlet synchrony with GSLOOK (Abbott)",
  cqslope_sync = "Absolute CQ-slope synchrony (Abbott)",
  cqslope_sync_outlet = "Outlet CQ-slope synchrony (Abbott)",
  cq_slope = "CQ slope (log-log)",
  Slope_mean = "Basin slope"
)

get_label <- function(var) {
  if (var %in% names(var_labels)) return(var_labels[[var]])
  if (var %in% names(STORAGE_LABELS)) return(STORAGE_LABELS[[var]])
  if (var %in% names(SYNC_LABELS)) return(SYNC_LABELS[[var]])
  return(var)
}

# =============================================================================
# ORDERING
# =============================================================================
site_order_storage <- c("WS09", "WS10", "WS01", "Look",
                        "WS02", "WS03", "Mack", "WS06", "WS07", "WS08")

# CQ data only has 8 sites (no GSMACK) - use site_order_cq for CQ plots
site_order_cq <- c("GSWS09", "GSWS10", "GSWS01", "GSLOOK",
                   "GSWS02", "GSWS06", "GSWS07", "GSWS08")
# Storage/comprehensive data includes GSMACK - use site_order for full dataset plots
site_order <- c("GSWS09", "GSWS10", "GSWS01", "GSLOOK",
                "GSWS02", "GSWS03", "GSMACK", "GSWS06", "GSWS07", "GSWS08")
solute_order <- c("Ca", "Mg", "Na", "K", "DSi", "Cl", "SO4",
                  "DOC", "NH3", "NO3", "PO4")
cluster_levels <- c("1", "2", "3", "4")

# =============================================================================
# SOLUTE CATEGORIES: Geogenic vs Biogenic vs Nutrient
# =============================================================================
# DSi behaves differently from other geogenic solutes, classify as "Nutrient"
GEOGENIC_SOLUTES <- c("Ca", "Mg", "Na", "K", "Cl", "SO4")
BIOGENIC_SOLUTES <- c("NO3", "PO4", "NH3", "DOC")
NUTRIENT_SOLUTES <- c("DSi")  # DSi behaves uniquely - weathering-derived but nutrient-like

# Backward compatibility aliases (consolidated from deprecated solutes.R)
GEO_SOLUTES <- c("Ca", "Mg", "Na", "K", "DSi", "Cl", "SO4")  # includes DSi
BIO_SOLUTES <- c("DOC", "NH3", "NO3", "PO4")
ALL_SOLUTES <- c(GEO_SOLUTES, BIO_SOLUTES)

# Solute group tibble for legacy compatibility
SOLUTE_GROUPS <- tibble::tibble(
  solute = ALL_SOLUTES,
  group = dplyr::case_when(
    solute %in% GEO_SOLUTES ~ "Geogenic",
    solute %in% BIO_SOLUTES ~ "Biogenic",
    TRUE ~ "Unknown"
  )
)

# Legacy helper function (use categorize_solute() or add_solute_type() for new code)
get_solute_group <- function(x) {
  dplyr::recode(x,
    !!!setNames(SOLUTE_GROUPS$group, SOLUTE_GROUPS$solute),
    .default = "Unknown"
  )
}

# Helper to categorize solutes (3-way: geogenic, biogenic, nutrient)
categorize_solute <- function(solute) {
  case_when(
    solute %in% NUTRIENT_SOLUTES ~ "Nutrient",
    solute %in% GEOGENIC_SOLUTES ~ "Geogenic",
    solute %in% BIOGENIC_SOLUTES ~ "Biogenic",
    TRUE ~ "Other"
  )
}

# Helper for 2-way categorization (geogenic vs biogenic, DSi with geogenic)
categorize_solute_2way <- function(solute) {
  case_when(
    solute %in% c(GEOGENIC_SOLUTES, NUTRIENT_SOLUTES) ~ "Geogenic",
    solute %in% BIOGENIC_SOLUTES ~ "Biogenic",
    TRUE ~ "Other"
  )
}

# Add solute_type column to dataframe (3-way by default)
add_solute_type <- function(df, solute_col = "solute", three_way = TRUE) {
  if (!solute_col %in% names(df)) return(df)
  if (three_way) {
    df$solute_type <- categorize_solute(df[[solute_col]])
    df$solute_type <- factor(df$solute_type, levels = c("Geogenic", "Biogenic", "Nutrient"))
  } else {
    df$solute_type <- categorize_solute_2way(df[[solute_col]])
    df$solute_type <- factor(df$solute_type, levels = c("Geogenic", "Biogenic"))
  }
  df
}

# =============================================================================
# STORAGE METRIC LABELS (explicit naming)
# =============================================================================
STORAGE_LABELS <- list(
  dynamic_storage_strength_z = "Dynamic storage strength (z)",
  mobile_mixing_no_bf_z = "Mobile mixing, tracer-only (z)",
  mobile_mixing_with_bf_z = "Mobile mixing + baseflow fraction (z)",
  flow_path_partitioning_z = "Flow-path partitioning, baseflow fraction (z)",
  unified_state_index = "Unified storage state index",
  geology_pc1 = "Geology/landslide PC1",
  geology_pc2 = "Geology/landslide PC2",
  RBI = "Richards-Baker flashiness index",
  RCS = "Recession-curve storage metric",
  FDC = "Flow-duration storage metric",
  SD = "Storage-discharge storage metric",
  WB = "Water-balance deficit storage metric",
  BF = "Calcium-derived baseflow fraction",
  DR = "Isotope damping ratio",
  Fyw = "Young-water fraction",
  MTT = "Mean transit time",
  RBI_short = "RBI",
  RCS_short = "RCS",
  FDC_short = "FDC",
  SD_short = "SD",
  WB_short = "WB",
  BF_short = "BF",
  DR_short = "DR",
  Fyw_short = "Fyw",
  MTT_short = "MTT",
  dynamic_storage_strength_z_short = "Dynamic storage",
  mobile_mixing_no_bf_z_short = "Mobile mixing",
  mobile_mixing_with_bf_z_short = "Mobile mixing + BF",
  flow_path_partitioning_z_short = "Flow paths",
  unified_state_index_short = "Storage state"
)

# Synchrony metric labels
SYNC_LABELS <- list(
  conc_sync_allpairs = "Absolute subcatchment synchrony (Abbott)",
  conc_sync_outlet = "Outlet synchrony with GSLOOK (Abbott)",
  cqslope_sync_allpairs = "Absolute CQ-slope synchrony (Abbott)",
  cqslope_sync_outlet = "Outlet CQ-slope synchrony (Abbott)",
  wymore_crosssite_allpairs = "Cross-site CQ quadrant agreement (Wymore)",
  wymore_crosssite_outlet = "Outlet CQ quadrant agreement (Wymore)",
  wymore_cvcq_consistency = "CV(C)/CV(Q) consistency (Wymore)"
)

get_storage_label <- function(metric, short = FALSE) {
  key <- if(short) paste0(metric, "_short") else metric
  if (key %in% names(STORAGE_LABELS)) return(STORAGE_LABELS[[key]])
  if (metric %in% names(STORAGE_LABELS)) return(STORAGE_LABELS[[metric]])
  return(metric)
}

get_sync_label <- function(metric) {
  if (metric %in% names(SYNC_LABELS)) return(SYNC_LABELS[[metric]])
  return(metric)
}

# =============================================================================
# COLOR PALETTES
# =============================================================================
# Cluster colors: tan → green → teal → blue
# Cluster colors sampled from climatology reference figure
cluster_colors <- c(
  "1" = "#CFA980",
  "2" = "#98B89F",
  "3" = "#5E8AA1",
  "4" = "#526B8E"
)

# =============================================================================
# CLUSTER LABELS (based on DTW clustering of monthly concentration z-scores)
# =============================================================================
# Clusters are defined by SEASONAL CONCENTRATION PATTERNS, not CQ behavior!
# - DTW (Dynamic Time Warping) distance between 12-month normalized concentration patterns
# - DBA (DTW Barycenter Averaging) centroids
# - k=4 determined by cluster validity indices
#
# Actual seasonal patterns (wet = Dec-Feb high flow, dry = Jun-Aug low flow):
#   Cluster 1: Peak Sep, concentrations build during baseflow → "Baseflow Enriched"
#   Cluster 2: Flat pattern, stable year-round → "Chemostatic"  
#   Cluster 3: Peak Jun, elevated in spring/early summer → "Spring/Early Summer Enriched"
#   Cluster 4: Peak Jan, high conc during high flow → "Winter Flushing"
cluster_labels <- c(
  "1" = "1-Baseflow Enriched",
  "2" = "2-Chemostatic",
  "3" = "3-Spring/Early Summer Enriched", 
  "4" = "4-Winter Flushing"
)

# Named cluster colors with labels
cluster_colors_labeled <- setNames(
  cluster_colors,
  cluster_labels
)

# Function to add cluster labels to data
add_cluster_labels <- function(df, cluster_col = "Cluster") {
  if (!cluster_col %in% names(df)) return(df)
  df$cluster_label <- cluster_labels[as.character(df[[cluster_col]])]
  df$cluster_label <- factor(df$cluster_label, levels = cluster_labels)
  df
}
# Solute type colors
solute_type_colors <- c(
  "Geogenic" = "#F9D5A7",
  "Biogenic" = "#74B49B",
  "Nutrient" = "#A75D5D"
)
season_colors <- c("Wet" = "#4E8098", "Dry" = "#D7867E")

# Active paper figure colors. Site colors below match the storage paper; these
# additional colors keep ordination vectors and diverging summaries quieter.
ordination_solute_vector_color <- "#9CA3AF"
ordination_storage_vector_color <- "#374151"
ordination_storage_label_color <- "#374151"
diverging_low_color <- "#2F6B9A"
diverging_mid_color <- "#F7F7F2"
diverging_high_color <- "#C07F2C"

# CQ behavior colors:
# Mobilizing = darker (concentrations increase with flow)
# Diluting = lighter (concentrations diluted by flow)
# Chemostatic = neutral (stable concentrations)
cq_behavior_colors <- c(
  "mobilizing"  = "#2C5F7C",  # dark teal-blue
  "diluting"    = "#BBDDE6",  # light blue
  "chemostatic" = "#CCCCCC"   # medium gray
)

cq_behavior_order <- c("mobilizing", "diluting", "chemostatic")

# Site colors match the storage paper palette for continuity across papers.
site_colors_storage <- c(
  "WS09" = "#882255",
  "WS10" = "#AA4499",
  "WS01" = "#CC6677",
  "Look" = "#DDCC77",
  "WS02" = "#999933",
  "WS03" = "#117733",
  "Mack" = "#332288",
  "WS06" = "#44AA99",
  "WS07" = "#88CCEE",
  "WS08" = "#6699CC"
)

site_colors <- c(
  "GSWS09" = site_colors_storage[["WS09"]],
  "GSWS10" = site_colors_storage[["WS10"]],
  "GSWS01" = site_colors_storage[["WS01"]],
  "GSLOOK" = site_colors_storage[["Look"]],
  "GSWS02" = site_colors_storage[["WS02"]],
  "GSWS03" = site_colors_storage[["WS03"]],
  "GSMACK" = site_colors_storage[["Mack"]],
  "GSWS06" = site_colors_storage[["WS06"]],
  "GSWS07" = site_colors_storage[["WS07"]],
  "GSWS08" = site_colors_storage[["WS08"]]
)

# Reorder to match site_order
site_colors <- site_colors[site_order]

# For CQ plots (no GSMACK)
site_colors_cq <- site_colors[site_order_cq]

# Legacy gradient version (keep for backward compatibility)
site_gradient_fn <- grDevices::colorRampPalette(c("#08306B", "#C6DBEF"))
site_colors_gradient <- setNames(site_gradient_fn(length(site_order)), site_order)

solute_colors <- c(
  Ca  = "#F9D5A7",
  Mg  = "#F8A978",
  Na  = "#F4976C",
  K   = "#D7867E",
  Cl  = "#85586F",
  SO4 = "#E3B778",
  DSi = "#A75D5D",
  DOC = "#A7D0CD",
  NH3 = "#74B49B",
  NO3 = "#508CA4",
  PO4 = "#87A8A4"
)

# =============================================================================
# GGPLOT SCALES
# =============================================================================
scale_color_season <- function() scale_color_manual(values = season_colors, name = "Season")
scale_fill_season <- function() scale_fill_manual(values = season_colors, name = "Season")
scale_color_cluster <- function() scale_color_manual(values = cluster_colors, name = "Cluster")
scale_fill_cluster <- function() scale_fill_manual(values = cluster_colors, name = "Cluster")
scale_color_site <- function(name = "Site", ...) scale_color_manual(values = site_colors, name = name, ...)
scale_fill_site  <- function(name = "Site", ...) scale_fill_manual(values = site_colors, name = name, ...)
scale_color_site_cq <- function(name = "Site", ...) scale_color_manual(values = site_colors_cq, name = name, ...)
scale_fill_site_cq  <- function(name = "Site", ...) scale_fill_manual(values = site_colors_cq, name = name, ...)
scale_color_cq_behavior <- function(name = "CQ Behavior", ...) scale_color_manual(values = cq_behavior_colors, name = name, ...)
scale_fill_cq_behavior <- function(name = "CQ Behavior", ...) scale_fill_manual(values = cq_behavior_colors, name = name, ...)
scale_color_solute <- function(name = "Solute") scale_color_manual(values = solute_colors, name = name)
scale_fill_solute <- function(name = "Solute") scale_fill_manual(values = solute_colors, name = name)

# =============================================================================
# THEME
# =============================================================================
theme_hja <- function(base_size = BASE_SIZE) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 2),
      plot.subtitle = element_text(size = base_size),
      plot.caption = element_text(size = base_size - 3, hjust = 0),
      strip.text = element_text(size = base_size + 1),
      strip.background = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      axis.ticks.length = unit(0.15, "cm"),
      axis.text = element_text(color = "black")
    )
}

# Alternative clean theme (simpler, no grid)
theme_clean <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size + 1),
      legend.position = "bottom",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      axis.ticks.length = unit(0.15, "cm"),
      axis.text = element_text(color = "black")
    )
}

# =============================================================================
# MONTH LABELS
# =============================================================================
month_labels <- c("J","F","M","A","M","J","J","A","S","O","N","D")

# =============================================================================
# SOLUTE GROUP PALETTES (for backward compatibility)
# =============================================================================
# Separate geo and bio palettes (subsets of solute_colors above)
geo_palette <- c(
  Ca  = "#F9D5A7",
  Mg  = "#F8A978",
  Na  = "#F4976C",
  K   = "#D7867E",
  DSi = "#A75D5D",
  Cl  = "#85586F",
  SO4 = "#E3B778"
)

bio_palette <- c(
  DOC = "#A7D0CD",
  NH3 = "#74B49B",
  NO3 = "#508CA4",
  PO4 = "#87A8A4"
)

# =============================================================================
# SAVE FUNCTION
# =============================================================================
save_plot <- function(p, filename, dir, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(
    file.path(dir, filename),
    plot = p,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}

save_plot_pdf <- function(p, filename, dir, width = PLOT_WIDTH, height = PLOT_HEIGHT) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(
    file.path(dir, filename),
    plot = p,
    width = width,
    height = height,
    device = cairo_pdf,
    bg = "white"
  )
}

# =============================================================================
# FACTOR HELPER
# =============================================================================
apply_factor_orders <- function(df) {
  if ("Stream_Name" %in% names(df)) df$Stream_Name <- factor(df$Stream_Name, levels = site_order)
  if ("site" %in% names(df)) df$site <- factor(df$site, levels = site_order)
  if ("solute" %in% names(df)) df$solute <- factor(df$solute, levels = solute_order)
  if ("chemical" %in% names(df)) df$chemical <- factor(df$chemical, levels = solute_order)
  if ("cluster" %in% names(df)) df$cluster <- factor(df$cluster, levels = cluster_levels)
  if ("hydrologic_season" %in% names(df)) df$hydrologic_season <- factor(df$hydrologic_season, levels = c("Wet", "Dry"))
  df
}

scale_shape_season <- function() scale_shape_manual(values = c("Wet" = 19, "Dry" = 17), name = "Season")
legend_bottom <- function() theme(legend.position = "bottom")
legend_right  <- function() theme(legend.position = "right")
