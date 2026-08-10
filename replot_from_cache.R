# =============================================================================
# replot_from_cache.R
# =============================================================================
# Regenerate all key publication plots from cached RDS files.
# No heavy recomputation — only reads existing .rds caches and re-renders plots.
#
# Run this file whenever you want to adjust text, labels, or dimensions without
# re-running the full analysis.
#
# HOW TO ADJUST PLOT DIMENSIONS:
#   Each section below has a clearly marked "# ── DIMENSIONS ──" block.
#   Edit the width/height values there. Common reference sizes:
#     A4 landscape : 11.7 × 8.3  in
#     A4 portrait  : 8.3  × 11.7 in
#     Square panel : 9    × 9    in
#
#   Some plot functions have dimensions hard-coded in their source files.
#   Those locations are noted explicitly in each section below.
# =============================================================================

source("FULL_FUNCTION.R")              # provides regenerate_plots_from_cache()
source("plot_subsampling_comparison.R") # provides plot_subsampling_comparison()
source("plot_diagonal_comparison.R")   # provides plot_diagonal_comparison()
source("plot_marker_dotplots.R")       # provides plot_marker_dotplots()
# source("plot_normalized_sum_vs_pss.R") # provides plot_normalized_sum_vs_pss()
library(patchwork)                     # for combining plots into A4 layout


# =============================================================================
# Section 0: Load all caches
# =============================================================================
# All paths are relative to the working directory (your project root).
# If a file is missing, a warning is printed and the variable is set to NULL.

cat("\n===== Section 0: Loading caches =====\n")

.load_rds <- function(path, label) {
  if (file.exists(path)) {
    cat(sprintf("  [OK]      %s\n  <- %s\n", label, path))
    readRDS(path)
  } else {
    warning(sprintf("[MISSING] %s not found:\n  %s", label, path))
    NULL
  }
}

# -- Main flow results (from use_FULL_FUNCTION.R) --
Uchimura_full_flow <- .load_rds("six2gfp/12.3.26/Uchimura_full_flow10noise.rds",  "Uchimura_full_flow")
freedman_flow      <- .load_rds("six2gfp/12.3.26/freedman_flow10noise.rds",        "freedman_flow")
cell_atlas_flow    <- .load_rds("six2gfp/12.3.26/cell_atlas_flow10noise.rds",      "cell_atlas_flow")
Takasato_full_flow <- .load_rds("six2gfp/12.3.26/Takasato_full_flow10noise.rds",   "Takasato_full_flow")
Vanslambrouck_full_flow <- .load_rds("six2gfp/Vanslambrouck_full_flow10noise.rds", "Vanslambrouck_full_flow")
Vanslambrouck_d13_full_flow <- .load_rds("six2gfp/Vanslambrouck_d13_full_flow10noise.rds", "Vanslambrouck_d13_full_flow")

# -- Atlas subsampling: baseline --
cell_atlas_flow
# cell_atlas_flow_anchors <- .load_rds(
#   "six2gfp/subsampling/atlas_with_anchors/flow_result.rds",
#   "cell_atlas_flow_anchors (baseline)"
# )

# -- Atlas subsampling: only the 3 removal types we need --
# To add/remove types, edit this vector:
ATLAS_REMOVAL_TYPES <- c("LOH", "PROX_1", "DIST_CD")

atlas_removal_flows <- list()
for (.type in ATLAS_REMOVAL_TYPES) {
  atlas_removal_flows[[.type]] <- .load_rds(
    paste0("six2gfp/subsampling/atlas_no_", .type, "_anchors/flow_result.rds"),
    paste0("atlas_no_", .type)
  )
}

cat("\n===== All caches loaded =====\n\n")


# =============================================================================
# Section 1: Regenerate main_plots for each dataset
# =============================================================================
# Calls regenerate_plots_from_cache() which re-runs plot_and_save_stability_scatters()
# and copies the key plots into each dataset's main_plots/ subdirectory.
#
# ── DIMENSIONS ──
#   These plots are saved inside plot_and_save_stability_scatters() in FULL_FUNCTION.R.
#   To change dimensions, edit the ggsave() calls in that function:
#     - scatter "all_freq" plots: width = 12, height = 8  (~line 1138)
#     - scatter "stability_pss_linear_fit" plots: width = 8, height = 7  (~line 1172)
# =============================================================================

cat("===== Section 1: Regenerating main_plots =====\n")

# output_prefix_base must match the directory used when the flow was first run.
.main_flows <- list(
  list(cache  = "six2gfp/12.3.26/Uchimura_full_flow10noise.rds",
       output = "six2gfp/new/fixed_size/Uchimura_noised/"),
  list(cache  = "six2gfp/12.3.26/freedman_flow10noise.rds",
       output = "six2gfp/new/fixed_size/freedman_not_clean/"),
  list(cache  = "six2gfp/12.3.26/cell_atlas_flow10noise.rds",
       output = "six2gfp/new/fixed_size/kidney_cell_atlas_clean_no_fibroblast/"),
  list(cache  = "six2gfp/12.3.26/Takasato_full_flow10noise.rds",
       output = "six2gfp/new/fixed_size/Takasato_noised_tranpose/"),
  list(cache  = "six2gfp/Vanslambrouck_full_flow10noise.rds",
       output = "six2gfp/new/fixed_size/Vanslambrouck_noised/"),
  list(cache  = "six2gfp/Vanslambrouck_d13_full_flow10noise.rds",
       output = "six2gfp/new/fixed_size/Vanslambrouck_d13_noised/"),
  # -- Human-to-Human flows --
  list(cache  = "human_to_human/Uchimura_full_flow5noise.rds",
       output = "human_to_human/new/Uchimura_noised/")
)

# .main_flows <- list(
#   list(cache  = "six2gfp/12.3.26/cell_atlas_flow10noise.rds",
#        output = "six2gfp/new/kidney_cell_atlas_clean/")
# )

for (.entry in .main_flows) {
  cat(sprintf("\n  -> %s\n", .entry$cache))
  regenerate_plots_from_cache(
    cache_path         = .entry$cache,
    output_prefix_base = .entry$output
  )
}

cat("\n===== Section 1 complete =====\n\n")


# =============================================================================
# Section 2: Atlas single-removal comparison (LOH, PROX_1, DIST_CD vs Baseline)
# =============================================================================
# Calls plot_subsampling_comparison() which saves per-removal SVG files and a
# combined PDF to:
#   six2gfp/subsampling/atlas_single_removal_comparison*
#
# After that, the 3 bi-directional stability paired barplots (one per removal type) are
# arranged side-by-side into a single A4-landscape figure.
#
# ── DIMENSIONS ──
#   A) Individual paired SVG plots (Baseline / No X, stacked):
#      Edit width/height inside plot_subsampling_comparison.R, the ggsave()
#      call around the "paired_plot" variable (search "paired_barplot_dir").
#      Current: width = 10, height = 10 in.
#
#   B) Combined A4 row figure (all 3 side-by-side):
#      Edit A4_ROW_WIDTH / A4_ROW_HEIGHT directly below.
# =============================================================================

cat("===== Section 2: Atlas single-removal comparison =====\n")

# ── DIMENSIONS: combined A4 row ─────────────────────────────────────────────
A4_ROW_WIDTH  <- 11.7  # in — A4 landscape; use 8.3 for portrait
A4_ROW_HEIGHT <-  8.3  # in — A4 landscape; use 11.7 for portrait
# ────────────────────────────────────────────────────────────────────────────

# Build result list: Baseline first, then the 3 removals in order
.atlas_sub_results <- c(
  list("Baseline" = cell_atlas_flow),
  setNames(
    lapply(ATLAS_REMOVAL_TYPES, function(t) atlas_removal_flows[[t]]),
    paste0("No ", ATLAS_REMOVAL_TYPES)
  )
)
# Drop any NULLs (missing caches won't break the call)
.atlas_sub_results <- Filter(Negate(is.null), .atlas_sub_results)

atlas_subsampling_plots <- plot_subsampling_comparison(
  .atlas_sub_results,
  output_prefix = "six2gfp/new/atlas_single_removal_comparison"
)

# --- Combine the 3 bi-directional stability paired barplots in a row for A4 ---
.paired_bi <- list()
for (.type in ATLAS_REMOVAL_TYPES) {
  .key <- paste0("No ", .type, "_bi_stability")
  .p   <- atlas_subsampling_plots$paired_barplots[[.key]]
  if (!is.null(.p)) {
    .paired_bi[[.type]] <- .p
  } else {
    warning(paste0("Paired bi-directional stability plot not found for key: ", .key,
                   "\n  (check that bi-directional stability data exists in the cached flow)"))
  }
}

if (length(.paired_bi) > 0) {
  .combined_a4_row <- wrap_plots(.paired_bi, nrow = 1) +
    plot_annotation(
      title = "Bi-directional Stability by Cell Type: Effect of Atlas Cell-Type Removal",
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )

  .a4_out <- "six2gfp/new/atlas_single_removal_comparison_paired_bi_stability_A4row.svg"
  ggsave(
    .a4_out,
    plot   = .combined_a4_row,
    width  = A4_ROW_WIDTH,
    height = A4_ROW_HEIGHT,
    units  = "in"
  )
  cat(sprintf("  [SAVED] %s\n", .a4_out))
}

cat("\n===== Section 2 complete =====\n\n")


# =============================================================================
# Section 3: Diagonal comparison across all 4 datasets
# =============================================================================
# Calls plot_diagonal_comparison() and saves all metric bar plots to:
#   six2gfp/7.5.26/diagonal_comparison_*
#
# Key outputs (n=10 because all 4 datasets used noised_number = 10):
#   six2gfp/7.5.26/diagonal_comparison__n10_bi_stability_by_celltype_comparison.svg
#   six2gfp/7.5.26/diagonal_comparison__n10_stability_by_celltype_comparison.svg
#   six2gfp/7.5.26/diagonal_comparison__n10_f1_stability_by_celltype_comparison.svg
#   six2gfp/7.5.26/diagonal_comparison__n10_pss_by_celltype_comparison.svg
#
# ── DIMENSIONS ──
#   These are set by ggsave() calls inside plot_diagonal_comparison.R.
#   To change them, search for ggsave() in that file and edit width/height.
#   Current: width = 12, height = 8 in (for all bar plots in that function).
# =============================================================================

cat("===== Section 3: Diagonal comparison =====\n")

# Only include datasets whose cache loaded successfully
.diag_list <- list()
if (!is.null(cell_atlas_flow))     .diag_list[["Cell Atlas"]] <- cell_atlas_flow
if (!is.null(Takasato_full_flow))  .diag_list[["Takasato"]]   <- Takasato_full_flow
if (!is.null(Uchimura_full_flow))  .diag_list[["Uchimura"]]   <- Uchimura_full_flow
if (!is.null(freedman_flow))       .diag_list[["Harder"]]   <- freedman_flow # Freedman is Harder's paperW
if (!is.null(Vanslambrouck_d13_full_flow)) .diag_list[["Vanslambrouck d13"]] <- Vanslambrouck_d13_full_flow
if (!is.null(Vanslambrouck_full_flow)) .diag_list[["Vanslambrouck"]] <- Vanslambrouck_full_flow

if (length(.diag_list) < 2) {
  warning("Need at least 2 datasets for diagonal comparison — check caches above.")
} else {
  diagonal_plots <- plot_diagonal_comparison(
    .diag_list,
    output_prefix = "six2gfp/new/diagonal_comparison_"
  )
}

cat("\n===== Section 3 complete =====\n")

# =============================================================================
# Section 4: Marker Gene DotPlots (Optional)
# =============================================================================
# Calls plot_marker_dotplots() to generate dotplots showing expression of 
# top marker genes for each cell type in each dataset.
# 
# Uncomment the code below to run it. Note: FindAllMarkers can take a while.
# =============================================================================

cat("\n===== Section 4: Marker DotPlots =====\n")
if (length(.diag_list) > 0) {
  # cache_paths maps each dataset name (from .diag_list) to the _cache/ directory
  # where run_without_noise.rds was saved during the original analysis run.
  # These must match the output_prefix_base values used in use_FULL_FUNCTION.R.
  .dotplot_cache_paths <- list(
    "Cell Atlas"       = "six2gfp/12.3.26/kidney_cell_atlas_clean_no_fibroblast/_cache",
    "Takasato"         = "six2gfp/12.3.26/Takasato_noised_tranpose/_cache",
    "Uchimura"         = "six2gfp/12.3.26/Uchimura_noised/_cache",
    "Harder"           = "six2gfp/12.3.26/freedman_not_clean/_cache",
    "Vanslambrouck d13" = "six2gfp/Vanslambrouck_d13_noised/_cache",
    "Vanslambrouck"    = "six2gfp/Vanslambrouck_noised/_cache"
  )
  
  # Load original Six2GFP and add it to the list
  if (file.exists("six2gfp/six2gfp_seurat_full")) {
    cat("  -- Loading original Six2GFP dataset --\n")
    six_env <- new.env()
    load("six2gfp/six2gfp_seurat_full", envir = six_env) # /lab/six2gfp/SIX2GFP_as_ground.R

    # When only that needed, uncomment: 
    # .diag_list = c(list("Original Six2GFP" = list(test = six2gfp_seurat)))
    # .dotplot_cache_paths <- c(list("Original Six2GFP" = NULL))
    
    # Add to beginning of .diag_list
    .diag_list <- c(list("Original Six2GFP" = list(test = six_env$six2gfp_seurat)), .diag_list)
    .dotplot_cache_paths <- c(list("Original Six2GFP" = NULL), .dotplot_cache_paths)
  }
  
  cat("\n  -- Method 1: Dynamic Top Markers --\n")
  marker_dotplots_dynamic <- plot_marker_dotplots(
    .diag_list,
    output_prefix = "six2gfp/new/dotplots_dynamic/",
    n_markers = 3, # Top 3 markers per cell type
    cache_paths = .dotplot_cache_paths
  )
  
  cat("\n  -- Method 2: Hardcoded Marker List --\n")
  marker_dotplots_hardcoded <- plot_hardcoded_marker_dotplots(
    .diag_list,
    output_prefix = "six2gfp/new/dotplots_hardcoded/",
    cache_paths = .dotplot_cache_paths
  )
}
cat("\n===== Section 4 complete =====\n")

cat("\n===== replot_from_cache.R: ALL DONE =====\n")
