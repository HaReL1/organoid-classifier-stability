# =============================================================================
# Usage of Subsampling Pipeline
# =============================================================================
# This script shows how to use the subsampling_pipeline.R functions:
#   - Remove one or more cell types from a dataset
#   - Re-classify and compare stability/PSS changes
#   - Optionally export and analyze transfer anchors
#
# Prerequisites: 
#   - source("FULL_FUNCTION.R")   # for analyze_noise_impact_on_prediction
#   - source("subsampling_pipeline.R")
#   - Datasets and flow results loaded (from use_FULL_FUNCTION.R)
# =============================================================================

source("FULL_FUNCTION.R")
source("subsampling_pipeline.R")

# ============================================================
# Example 1: Remove a single type from Atlas — step by step
# ============================================================

# Step 1: Remove all LOH cells
no_loh <- create_subsampled_seurat(
  seurat_obj    = atlas_object,
  flow_result   = cell_atlas_flow,
  target_labels = "LOH",
  keep_fraction = 0.0,
  seed          = 42
)

# Step 2: Run full analysis on subsampled data
no_loh_flow <- analyze_noise_impact_on_prediction(
  no_loh$seurat,
  noised_number = 3,
  train_Six2GFP,
  test_Six2GFP,
  ref_cell_type_column = "type",
  dims = 1:30,
  train_title = "SIX2GFP",
  test_title_prefix = "Atlas_no_LOH",
  n_neighbors = 8,
  skip_neighbors = TRUE,
  output_prefix_base = "six2gfp/subsampling/atlas_no_LOH/",
  prediction_column_name = "predicted.type",
  colors_feature_plot_noise = c('grey', '#f03b20'),
  myColors_cell_types = NULL,
  return_all_suerats = FALSE,
  use_cache = FALSE
)

# Step 3: Compare — PSS, stability, and what gets predicted as LOH now
comparison <- compare_flows(
  original_flow  = cell_atlas_flow,
  subsampled_flow = no_loh_flow,
  removed_labels  = "LOH"
)

# Quick check: how many cells still get classified as LOH across noise iters?
for (col in 1:ncol(no_loh_flow$noised_prediction_matrix)) {
  cat(sprintf("  iter %d: %d cells predicted as LOH\n", col - 1,
              sum(no_loh_flow$noised_prediction_matrix[, col] == "LOH")))
}


# ============================================================
# Example 2: Remove multiple types from Atlas — all-in-one
# ============================================================

atlas_no_cmdiv_cm <- run_subsampling_analysis(
  seurat_obj       = atlas_object,
  flow_result      = cell_atlas_flow,
  target_labels    = c("CM_DIV", "CM"),
  keep_fraction    = 0.0,
  train_data       = train_Six2GFP,
  test_data        = test_Six2GFP,
  ref_cell_type_column = "type",
  dims             = 1:30,
  noised_number    = 3,
  output_prefix_base = "six2gfp/subsampling/atlas_no_CM_DIV_CM/",
  seed             = 42,
  train_title      = "SIX2GFP",
  test_title_prefix = "Atlas_no_CM_DIV_CM",
  n_neighbors      = 8,
  skip_neighbors   = TRUE,
  prediction_column_name = "predicted.type",
  return_all_suerats = FALSE,
  use_cache        = FALSE,
  return_anchors=TRUE
)


# ============================================================
# Example 3: Remove types from Uchimura
# ============================================================

# uchi_no_prox2 <- run_subsampling_analysis(
#   seurat_obj       = seurat_Uchimura_Humphreys_20,
#   flow_result      = Uchimura_full_flow,
#   target_labels    = "PROX_2",
#   keep_fraction    = 0.0,
#   train_data       = train_Six2GFP,
#   test_data        = test_Six2GFP,
#   ref_cell_type_column = "type",
#   dims             = 1:30,
#   noised_number    = 3,
#   output_prefix_base = "six2gfp/subsampling/uchimura_no_PROX2/",
#   seed             = 42,
#   train_title      = "SIX2GFP",
#   test_title_prefix = "Uchimura_no_PROX2",
#   n_neighbors      = 8,
#   skip_neighbors   = TRUE,
#   prediction_column_name = "predicted.type",
#   return_all_suerats = FALSE,
#   use_cache        = FALSE
# )


# ============================================================
# Example 4: Remove types from Freedman
# ============================================================

# freedman_no_um <- run_subsampling_analysis(
#   seurat_obj       = freedman_seurat_obj,
#   flow_result      = freedman_flow,
#   target_labels    = "UM",
#   keep_fraction    = 0.0,
#   train_data       = train_Six2GFP,
#   test_data        = test_Six2GFP,
#   ref_cell_type_column = "type",
#   dims             = 1:15,
#   noised_number    = 3,
#   output_prefix_base = "six2gfp/subsampling/freedman_no_UM/",
#   seed             = 42,
#   train_title      = "SIX2GFP",
#   test_title_prefix = "Freedman_no_UM",
#   n_neighbors      = 8,
#   skip_neighbors   = TRUE,
#   prediction_column_name = "predicted.type",
#   return_all_suerats = FALSE,
#   use_cache        = FALSE
# )


# ============================================================
# Example 5: Systematic removal — remove each type one-at-a-time
# ============================================================

# types_to_test <- c("UM", "CM", "PROX_1", "PROX_2", "LOH", "DIST_CD", "ENDO", "PODO")
# all_comparisons <- list()
#
# for (type_to_remove in types_to_test) {
#   cat(sprintf("\n\n######## Removing %s from Atlas ########\n", type_to_remove))
#
#   result <- run_subsampling_analysis(
#     seurat_obj       = atlas_object,
#     flow_result      = cell_atlas_flow,
#     target_labels    = type_to_remove,
#     keep_fraction    = 0.0,
#     train_data       = train_Six2GFP,
#     test_data        = test_Six2GFP,
#     ref_cell_type_column = "type",
#     dims             = 1:30,
#     noised_number    = 3,
#     output_prefix_base = paste0("six2gfp/subsampling/atlas_no_", type_to_remove, "/"),
#     seed             = 42,
#     train_title      = "SIX2GFP",
#     test_title_prefix = paste0("Atlas_no_", type_to_remove),
#     n_neighbors      = 8,
#     skip_neighbors   = TRUE,
#     prediction_column_name = "predicted.type",
#     return_all_suerats = FALSE,
#     use_cache        = FALSE
#   )
#
#   all_comparisons[[type_to_remove]] <- result$comparison
# }


# ============================================================
# Example 6: With anchor analysis (before vs after removal)
# ============================================================
# Run BOTH flows with return_anchors=TRUE. compare_flows() will
# automatically detect the anchors and run compare_anchors().
#
# Step 1: Original flow with anchors exported
cell_atlas_flow_anchors <- analyze_noise_impact_on_prediction(
  atlas_object, noised_number = 3,
  train_Six2GFP, test_Six2GFP,
  ref_cell_type_column = "type", dims = 1:30,
  train_title = "SIX2GFP", test_title_prefix = "Kidney Cell Atlas",
  n_neighbors = 8, skip_neighbors = TRUE,
  output_prefix_base = "six2gfp/subsampling/atlas_with_anchors/",
  prediction_column_name = "predicted.type",
  use_cache = FALSE,
  return_anchors = TRUE
)

# Step 2: Remove LOH and run with anchors
no_cm_cmdiv_prox1_info <- create_subsampled_seurat(atlas_object, cell_atlas_flow_anchors, c("CM_DIV", "CM", "PROX_1"), 0.0)
no_cm_cmdiv_prox1_flow_anchors <- analyze_noise_impact_on_prediction(
  no_cm_cmdiv_prox1_info$seurat, noised_number = 3,
  train_Six2GFP, test_Six2GFP,
  ref_cell_type_column = "type", dims = 1:30,
  train_title = "SIX2GFP", test_title_prefix = "Atlas_no_cm_cmdiv_prox1",
  n_neighbors = 8, skip_neighbors = TRUE,
  output_prefix_base = "six2gfp/subsampling/atlas_no_cm_cmdiv_prox1_anchors/",
  prediction_column_name = "predicted.type",
  use_cache = FALSE,
  return_anchors = TRUE
)

# Step 3: compare_flows handles anchor comparison automatically!
comparison_with_anchors <- compare_flows(
  original_flow   = cell_atlas_flow_anchors,
  subsampled_flow = no_cm_cmdiv_prox1_flow_anchors,
  removed_labels  = c("CM_DIV", "CM", "PROX_1"),
  train_data      = train_Six2GFP
)
# # Output includes:
# #   - PSS/Stability comparisons (as before)
# #   - Anchor Comparison:
# #     - Total anchor pairs before vs after
# #     - Anchors by reference cell type
# #     - Query cells whose dominant anchor type changed
# #     - Cells anchoring to the removed type's ref cells


# ============================================================
# Example 7: Full Atlas subsampling — multi-group removals
# ============================================================
# Prerequisites: cell_atlas_flow_anchors and no_cm_cmdiv_prox1_flow_anchors
# already computed in Example 6.

# --- Atlas: No LOH+CM ---
no_loh_cm_info <- create_subsampled_seurat(atlas_object, cell_atlas_flow_anchors, c("LOH", "CM"), 0.0)
no_loh_cm_flow <- analyze_noise_impact_on_prediction(
  no_loh_cm_info$seurat, noised_number = 3,
  train_Six2GFP, test_Six2GFP,
  ref_cell_type_column = "type", dims = 1:30,
  train_title = "SIX2GFP", test_title_prefix = "Atlas_no_LOH_CM",
  n_neighbors = 8, skip_neighbors = TRUE,
  output_prefix_base = "six2gfp/subsampling/atlas_no_LOH_CM_anchors/",
  prediction_column_name = "predicted.type",
  use_cache = FALSE,
  return_anchors = TRUE
)
comparison_no_loh_cm <- compare_flows(
  original_flow   = cell_atlas_flow_anchors,
  subsampled_flow = no_loh_cm_flow,
  removed_labels  = c("LOH", "CM"),
  train_data      = train_Six2GFP
)


# ============================================================
# Example 8: Atlas — remove EACH cell type alone
# ============================================================

atlas_single_type_flows <- list()
atlas_single_types <- c("UM", "CM", "CM_DIV", "PODO", "PROX_1", "PROX_2", "LOH", "DIST_CD", "ENDO", "MACROPHAG")

for (type_to_remove in atlas_single_types) {
  cat(sprintf("\n\n######## Atlas: Removing %s ########\n", type_to_remove))
  
  info <- create_subsampled_seurat(atlas_object, cell_atlas_flow_anchors, type_to_remove, 0.0)
  flow <- analyze_noise_impact_on_prediction(
    info$seurat, noised_number = 3,
    train_Six2GFP, test_Six2GFP,
    ref_cell_type_column = "type", dims = 1:30,
    train_title = "SIX2GFP", test_title_prefix = paste0("Atlas_no_", type_to_remove),
    n_neighbors = 8, skip_neighbors = TRUE,
    output_prefix_base = paste0("six2gfp/subsampling/atlas_no_", type_to_remove, "_anchors/"),
    prediction_column_name = "predicted.type",
    use_cache = FALSE,
    return_anchors = TRUE
  )
  comparison <- compare_flows(
    original_flow   = cell_atlas_flow_anchors,
    subsampled_flow = flow,
    removed_labels  = type_to_remove,
    train_data      = train_Six2GFP
  )
  
  atlas_single_type_flows[[type_to_remove]] <- list(
    info = info, flow = flow, comparison = comparison
  )
}


# ============================================================
# Example 9: Uchimura subsampling
# ============================================================

# Step 1: Uchimura baseline with anchors
Uchimura_flow_anchors <- analyze_noise_impact_on_prediction(
  seurat_Uchimura_Humphreys_20, noised_number = 3,
  train_Six2GFP, test_Six2GFP,
  ref_cell_type_column = "type", dims = 1:30,
  train_title = "SIX2GFP", test_title_prefix = "Uchimura",
  n_neighbors = 8, skip_neighbors = TRUE,
  output_prefix_base = "six2gfp/subsampling/uchimura_with_anchors/",
  prediction_column_name = "predicted.type",
  use_cache = FALSE,
  return_anchors = TRUE
)

# --- Uchimura: No DIST_CD+LOH ---
uchi_no_distcd_loh_info <- create_subsampled_seurat(seurat_Uchimura_Humphreys_20, Uchimura_flow_anchors, c("DIST_CD", "LOH"), 0.0)
uchi_no_distcd_loh_flow <- analyze_noise_impact_on_prediction(
  uchi_no_distcd_loh_info$seurat, noised_number = 3,
  train_Six2GFP, test_Six2GFP,
  ref_cell_type_column = "type", dims = 1:30,
  train_title = "SIX2GFP", test_title_prefix = "Uchimura_no_DISTCD_LOH",
  n_neighbors = 8, skip_neighbors = TRUE,
  output_prefix_base = "six2gfp/subsampling/uchimura_no_DISTCD_LOH_anchors/",
  prediction_column_name = "predicted.type",
  use_cache = FALSE,
  return_anchors = TRUE
)
comparison_uchi_no_distcd_loh <- compare_flows(
  original_flow   = Uchimura_flow_anchors,
  subsampled_flow = uchi_no_distcd_loh_flow,
  removed_labels  = c("DIST_CD", "LOH"),
  train_data      = train_Six2GFP
)

# --- Uchimura: No PROX_2 ---
uchi_no_prox2_info <- create_subsampled_seurat(seurat_Uchimura_Humphreys_20, Uchimura_flow_anchors, "PROX_2", 0.0)
uchi_no_prox2_flow <- analyze_noise_impact_on_prediction(
  uchi_no_prox2_info$seurat, noised_number = 3,
  train_Six2GFP, test_Six2GFP,
  ref_cell_type_column = "type", dims = 1:30,
  train_title = "SIX2GFP", test_title_prefix = "Uchimura_no_PROX2",
  n_neighbors = 8, skip_neighbors = TRUE,
  output_prefix_base = "six2gfp/subsampling/uchimura_no_PROX2_anchors/",
  prediction_column_name = "predicted.type",
  use_cache = FALSE,
  return_anchors = TRUE
)
comparison_uchi_no_prox2 <- compare_flows(
  original_flow   = Uchimura_flow_anchors,
  subsampled_flow = uchi_no_prox2_flow,
  removed_labels  = "PROX_2",
  train_data      = train_Six2GFP
)

# --- Uchimura: No CM+CM_DIV ---
uchi_no_cm_cmdiv_info <- create_subsampled_seurat(seurat_Uchimura_Humphreys_20, Uchimura_flow_anchors, c("CM", "CM_DIV"), 0.0)
uchi_no_cm_cmdiv_flow <- analyze_noise_impact_on_prediction(
  uchi_no_cm_cmdiv_info$seurat, noised_number = 3,
  train_Six2GFP, test_Six2GFP,
  ref_cell_type_column = "type", dims = 1:30,
  train_title = "SIX2GFP", test_title_prefix = "Uchimura_no_CM_CMDIV",
  n_neighbors = 8, skip_neighbors = TRUE,
  output_prefix_base = "six2gfp/subsampling/uchimura_no_CM_CMDIV_anchors/",
  prediction_column_name = "predicted.type",
  use_cache = FALSE,
  return_anchors = TRUE
)
comparison_uchi_no_cm_cmdiv <- compare_flows(
  original_flow   = Uchimura_flow_anchors,
  subsampled_flow = uchi_no_cm_cmdiv_flow,
  removed_labels  = c("CM", "CM_DIV"),
  train_data      = train_Six2GFP
)


# ============================================================
# Example 10: Plot subsampling comparisons
# ============================================================
source("plot_subsampling_comparison.R")

# --- Atlas: multi-group removals ---
atlas_multigroup_plots <- plot_subsampling_comparison(
  list(
    "Baseline"              = cell_atlas_flow_anchors,
    "No CM+CM_DIV+PROX_1"   = no_cm_cmdiv_prox1_flow_anchors,
    "No LOH+CM"             = no_loh_cm_flow
  ),
  output_prefix = "six2gfp/subsampling/atlas_multigroup_comparison"
)

# --- Atlas: each type removed alone ---
atlas_single_results <- list("Baseline" = cell_atlas_flow_anchors)
for (type_name in names(atlas_single_type_flows)) {
  atlas_single_results[[paste0("No ", type_name)]] <- atlas_single_type_flows[[type_name]]$flow
}
atlas_single_plots <- plot_subsampling_comparison(
  atlas_single_results,
  output_prefix = "six2gfp/subsampling/atlas_single_removal_comparison"
)

# --- Atlas: ALL removals combined ---
atlas_all_results <- c(
  list(
    "Baseline"              = cell_atlas_flow_anchors,
    "No CM+CM_DIV+PROX_1"   = no_cm_cmdiv_prox1_flow_anchors,
    "No LOH+CM"             = no_loh_cm_flow
  ),
  setNames(
    lapply(names(atlas_single_type_flows), function(t) atlas_single_type_flows[[t]]$flow),
    paste0("No ", names(atlas_single_type_flows))
  )
)
atlas_all_plots <- plot_subsampling_comparison(
  atlas_all_results,
  output_prefix = "six2gfp/subsampling/atlas_all_removal_comparison"
)

# --- Uchimura ---
uchimura_plots <- plot_subsampling_comparison(
  list(
    "Baseline"        = Uchimura_flow_anchors,
    "No DIST_CD+LOH"  = uchi_no_distcd_loh_flow,
    "No PROX_2"       = uchi_no_prox2_flow,
    "No CM+CM_DIV"    = uchi_no_cm_cmdiv_flow
  ),
  output_prefix = "six2gfp/subsampling/uchimura_removal_comparison"
)

# create braplot for each flavor 

# Atlas - CM+CM_DIV+PROX_1, LOH+CM, EACH celltype alone
# Uchimura - DICT_CD+LOH, PROX_2, CM+CM_DIV
# another try: calculate stability not only as what got out after noise - but also 
# what came in. 1-([out+in]/original) 
# dont remove old calculation, only add this