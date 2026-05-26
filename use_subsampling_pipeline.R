# =============================================================================
# Usage of Subsampling Pipeline
# =============================================================================
# This script runs the subsampling stability analysis for each dataset.
# For each dataset (Atlas, Uchimura, Freedman):
#   A. Compute baseline flow (with anchors) — or reuse existing
#   B. Remove EACH cell type individually (loop)
#   C. Remove specific multi-group combinations
#   D. Plot all comparisons at the end of each section
#
# Prerequisites:
#   - Datasets loaded from use_FULL_FUNCTION.R:
#     atlas_object, seurat_Uchimura_Humphreys_20, freedman_seurat_obj
#   - Baseline flows loaded: cell_atlas_flow, Uchimura_full_flow, freedman_flow
#   - train_Six2GFP, test_Six2GFP
# =============================================================================

source("FULL_FUNCTION.R")
source("subsampling_pipeline.R")
source("plot_subsampling_comparison.R")


# =============================================================================
# Helper: run or load cached flow result
# =============================================================================
# Saves the complete analyze_noise_impact_on_prediction() result as a single
# RDS file. On re-run, loads instantly instead of re-running everything.
# This is the most impactful optimization: if the script crashes mid-way,
# already-completed removals load in seconds, not hours.
run_or_load_flow <- function(cache_path, run_fn, label = "", force_replot = FALSE, expected_noised_number = NULL) {
  if (file.exists(cache_path)) {
    cached <- readRDS(cache_path)
    # Validate: if cached result is missing F1 (old code), recompute
    if (!is.null(cached$stability_pred$F1_Stability)) {
      # Validate: if noised_number changed, recompute
      cached_n <- ncol(cached$noised_prediction_matrix) - 1
      if (!is.null(expected_noised_number) && cached_n != expected_noised_number) {
        cat(sprintf("  [CACHE STALE] %s was run with n=%d but now requesting n=%s — recomputing\n",
                    label, cached_n, expected_noised_number))
        rm(cached); gc(verbose = FALSE)
      } else {
        cat(sprintf("  [CACHE HIT] Loading %s from %s\n", label, cache_path))
        if (force_replot) {
          regenerate_plots_from_cache(cache_path)
        }
        return(cached)
      }
    } else {
      cat(sprintf("  [CACHE STALE] %s missing F1_Stability — recomputing\n", label))
      rm(cached); gc(verbose = FALSE)
    }
  }
  t_start <- Sys.time()
  cat(sprintf("  [COMPUTING] %s ...\n", label))
  result <- run_fn()
  t_end <- Sys.time()
  cat(sprintf("  [DONE] %s — took %.1f min\n", label, as.numeric(difftime(t_end, t_start, units = "mins"))))
  
  # Ensure output directory exists
  cache_dir <- dirname(cache_path)
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
  
  saveRDS(result, cache_path)
  cat(sprintf("  [SAVED] %s\n", cache_path))
  return(result)
}


# #############################################################################
# Section A: Atlas Subsampling
# #############################################################################
cat("\n\n========== SECTION A: ATLAS SUBSAMPLING ==========\n\n")
section_a_start <- Sys.time()

# --- A1. Baseline flow with anchors ---
# Reuse existing flow ONLY if it has both anchors AND F1 (from updated FULL_FUNCTION.R) AND correct noised_number
has_atlas_anchors <- exists("cell_atlas_flow") && !is.null(cell_atlas_flow$initial_anchors)
has_atlas_f1 <- exists("cell_atlas_flow") && !is.null(cell_atlas_flow$stability_pred$F1_Stability)
has_atlas_correct_n <- exists("cell_atlas_flow") && (ncol(cell_atlas_flow$noised_prediction_matrix) - 1) == 10

if (has_atlas_anchors && has_atlas_f1 && has_atlas_correct_n) {
  cell_atlas_flow_anchors <- cell_atlas_flow
  cat("  [REUSE] Using cell_atlas_flow (has anchors + F1 + n=10)\n")
} else {
  if (!has_atlas_anchors) cat("  [INFO] cell_atlas_flow missing anchors — recomputing\n")
  if (!has_atlas_f1) cat("  [INFO] cell_atlas_flow missing F1_Stability — recomputing\n")
  if (!has_atlas_correct_n) cat("  [INFO] cell_atlas_flow has wrong noised_number — recomputing\n")
  cell_atlas_flow_anchors <- run_or_load_flow(
    "six2gfp/subsampling/atlas_with_anchors/flow_result.rds",
    function() {
      analyze_noise_impact_on_prediction(
        atlas_object, noised_number = 10,
        train_Six2GFP, test_Six2GFP,
        ref_cell_type_column = "type", dims = 1:30,
        train_title = "SIX2GFP", test_title_prefix = "Kidney Cell Atlas",
        n_neighbors = 8, skip_neighbors = TRUE,
        output_prefix_base = "six2gfp/subsampling/atlas_with_anchors/",
        prediction_column_name = "predicted.type",
        use_cache = TRUE,
        return_anchors = TRUE
      )
    },
    label = "Atlas baseline",
    expected_noised_number = 10
  )
}

# --- . Single-type removals (loop over ALL types) ---
atlas_single_types <- c("UM", "CM", "CM_DIV", "PODO", "PROX_1", "PROX_2", "LOH", "DIST_CD", "ENDO", "MACROPHAG")
atlas_single_type_flows <- list()

for (idx in seq_along(atlas_single_types)) {
  type_to_remove <- atlas_single_types[idx]
  cat(sprintf("\n######## Atlas: Removing %s (%d/%d) ########\n", type_to_remove, idx, length(atlas_single_types)))
  
  flow_cache <- paste0("six2gfp/subsampling/atlas_no_", type_to_remove, "_anchors/flow_result.rds")
  
  flow <- run_or_load_flow(flow_cache, function() {
    info <- create_subsampled_seurat(atlas_object, cell_atlas_flow_anchors, type_to_remove, 0.0)
    analyze_noise_impact_on_prediction(
      info$seurat, noised_number = 10,
      train_Six2GFP, test_Six2GFP,
      ref_cell_type_column = "type", dims = 1:30,
      train_title = "SIX2GFP", test_title_prefix = paste0("Atlas_no_", type_to_remove),
      n_neighbors = 8, skip_neighbors = TRUE,
      output_prefix_base = paste0("six2gfp/subsampling/atlas_no_", type_to_remove, "_anchors/"),
      prediction_column_name = "predicted.type",
      use_cache = TRUE,
      return_anchors = TRUE
    )
  }, label = paste0("Atlas no ", type_to_remove), expected_noised_number = 10)
  
  comparison <- compare_flows(
    original_flow   = cell_atlas_flow_anchors,
    subsampled_flow = flow,
    removed_labels  = type_to_remove,
    train_data      = train_Six2GFP
  )
  
  atlas_single_type_flows[[type_to_remove]] <- list(
    flow = flow, comparison = comparison
  )
  
  # Memory cleanup: remove the flow copy (it's saved in the list)
  rm(flow, comparison)
  gc(verbose = FALSE)
}

# --- A3. Multi-group combinations ---
cat("\n######## Atlas: Multi-group combinations ########\n")

# No CM+CM_DIV+PROX_1
no_cm_cmdiv_prox1_flow_anchors <- run_or_load_flow(
  "six2gfp/subsampling/atlas_no_cm_cmdiv_prox1_anchors/flow_result.rds",
  function() {
    info <- create_subsampled_seurat(atlas_object, cell_atlas_flow_anchors, c("CM_DIV", "CM", "PROX_1"), 0.0)
    analyze_noise_impact_on_prediction(
      info$seurat, noised_number = 3,
      train_Six2GFP, test_Six2GFP,
      ref_cell_type_column = "type", dims = 1:30,
      train_title = "SIX2GFP", test_title_prefix = "Atlas_no_cm_cmdiv_prox1",
      n_neighbors = 8, skip_neighbors = TRUE,
      output_prefix_base = "six2gfp/subsampling/atlas_no_cm_cmdiv_prox1_anchors/",
      prediction_column_name = "predicted.type",
      use_cache = TRUE,
      return_anchors = TRUE
    )
  },
  label = "Atlas no CM+CM_DIV+PROX_1",
  expected_noised_number = 3
)
gc(verbose = FALSE)

# No LOH+CM
no_loh_cm_flow <- run_or_load_flow(
  "six2gfp/subsampling/atlas_no_LOH_CM_anchors/flow_result.rds",
  function() {
    info <- create_subsampled_seurat(atlas_object, cell_atlas_flow_anchors, c("LOH", "CM"), 0.0)
    analyze_noise_impact_on_prediction(
      info$seurat, noised_number = 3,
      train_Six2GFP, test_Six2GFP,
      ref_cell_type_column = "type", dims = 1:30,
      train_title = "SIX2GFP", test_title_prefix = "Atlas_no_LOH_CM",
      n_neighbors = 8, skip_neighbors = TRUE,
      output_prefix_base = "six2gfp/subsampling/atlas_no_LOH_CM_anchors/",
      prediction_column_name = "predicted.type",
      use_cache = TRUE,
      return_anchors = TRUE
    )
  },
  label = "Atlas no LOH+CM",
  expected_noised_number = 3
)
gc(verbose = FALSE)

cat(sprintf("\n[ATLAS] Sections A1-A3 complete — total: %.1f min\n",
    as.numeric(difftime(Sys.time(), section_a_start, units = "mins"))))

# --- A4. Atlas plotting ---
cat("\n######## Atlas: Plotting ########\n")

# Single removals
atlas_single_results <- list("Baseline" = cell_atlas_flow_anchors)
for (type_name in names(atlas_single_type_flows)) {
  atlas_single_results[[paste0("No ", type_name)]] <- atlas_single_type_flows[[type_name]]$flow
}
atlas_single_plots <- plot_subsampling_comparison(
  atlas_single_results,
  output_prefix = "six2gfp/subsampling/atlas_single_removal_comparison"
)

# Multi-group removals
atlas_multigroup_plots <- plot_subsampling_comparison(
  list(
    "Baseline"              = cell_atlas_flow_anchors,
    "No CM+CM_DIV+PROX_1"   = no_cm_cmdiv_prox1_flow_anchors,
    "No LOH+CM"             = no_loh_cm_flow
  ),
  output_prefix = "six2gfp/subsampling/atlas_multigroup_comparison"
)

# ALL removals combined
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

cat(sprintf("\n[ATLAS] Section A complete (incl. plotting) — total: %.1f min\n",
    as.numeric(difftime(Sys.time(), section_a_start, units = "mins"))))


# #############################################################################
# Section B: Uchimura Subsampling
# #############################################################################
cat("\n\n========== SECTION B: UCHIMURA SUBSAMPLING ==========\n\n")
section_b_start <- Sys.time()

# --- B1. Baseline flow with anchors ---
has_uchi_anchors <- exists("Uchimura_full_flow") && !is.null(Uchimura_full_flow$initial_anchors)
has_uchi_f1 <- exists("Uchimura_full_flow") && !is.null(Uchimura_full_flow$stability_pred$F1_Stability)
has_uchi_correct_n <- exists("Uchimura_full_flow") && (ncol(Uchimura_full_flow$noised_prediction_matrix) - 1) == 3

if (has_uchi_anchors && has_uchi_f1 && has_uchi_correct_n) {
  Uchimura_flow_anchors <- Uchimura_full_flow
  cat("  [REUSE] Using Uchimura_full_flow (has anchors + F1 + n=3)\n")
} else {
  if (!has_uchi_anchors) cat("  [INFO] Uchimura_full_flow missing anchors — recomputing\n")
  if (!has_uchi_f1) cat("  [INFO] Uchimura_full_flow missing F1_Stability — recomputing\n")
  if (!has_uchi_correct_n) cat("  [INFO] Uchimura_full_flow has wrong noised_number — recomputing\n")
  Uchimura_flow_anchors <- run_or_load_flow(
    "six2gfp/subsampling/uchimura_with_anchors/flow_result.rds",
    function() {
      analyze_noise_impact_on_prediction(
        seurat_Uchimura_Humphreys_20, noised_number = 3,
        train_Six2GFP, test_Six2GFP,
        ref_cell_type_column = "type", dims = 1:30,
        train_title = "SIX2GFP", test_title_prefix = "Uchimura",
        n_neighbors = 8, skip_neighbors = TRUE,
        output_prefix_base = "six2gfp/subsampling/uchimura_with_anchors/",
        prediction_column_name = "predicted.type",
        use_cache = TRUE,
        return_anchors = TRUE
      )
    },
    label = "Uchimura baseline",
    expected_noised_number = 3
  )
}

# --- B2. Single-type removals (loop over ALL types) ---
uchimura_single_types <- c("UM", "CM", "CM_DIV", "PODO", "PROX_1", "PROX_2", "LOH", "DIST_CD", "ENDO")
uchimura_single_type_flows <- list()

for (idx in seq_along(uchimura_single_types)) {
  type_to_remove <- uchimura_single_types[idx]
  cat(sprintf("\n######## Uchimura: Removing %s (%d/%d) ########\n", type_to_remove, idx, length(uchimura_single_types)))
  
  flow_cache <- paste0("six2gfp/subsampling/uchimura_no_", type_to_remove, "_anchors/flow_result.rds")
  
  flow <- run_or_load_flow(flow_cache, function() {
    info <- create_subsampled_seurat(seurat_Uchimura_Humphreys_20, Uchimura_flow_anchors, type_to_remove, 0.0)
    analyze_noise_impact_on_prediction(
      info$seurat, noised_number = 3,
      train_Six2GFP, test_Six2GFP,
      ref_cell_type_column = "type", dims = 1:30,
      train_title = "SIX2GFP", test_title_prefix = paste0("uchimura_no_", type_to_remove),
      n_neighbors = 8, skip_neighbors = TRUE,
      output_prefix_base = paste0("six2gfp/subsampling/uchimura_no_", type_to_remove, "_anchors/"),
      prediction_column_name = "predicted.type",
      use_cache = TRUE,
      return_anchors = TRUE
    )
  }, label = paste0("Uchimura no ", type_to_remove), expected_noised_number = 3)
  
  comparison <- compare_flows(
    original_flow   = Uchimura_flow_anchors,
    subsampled_flow = flow,
    removed_labels  = type_to_remove,
    train_data      = train_Six2GFP
  )
  
  uchimura_single_type_flows[[type_to_remove]] <- list(
    flow = flow, comparison = comparison
  )
  
  rm(flow, comparison)
  gc(verbose = FALSE)
}

# --- B3. Multi-group combinations ---
cat("\n######## Uchimura: Multi-group combinations ########\n")

# No DIST_CD+LOH
uchi_no_distcd_loh_flow <- run_or_load_flow(
  "six2gfp/subsampling/uchimura_no_DISTCD_LOH_anchors/flow_result.rds",
  function() {
    info <- create_subsampled_seurat(seurat_Uchimura_Humphreys_20, Uchimura_flow_anchors, c("DIST_CD", "LOH"), 0.0)
    analyze_noise_impact_on_prediction(
      info$seurat, noised_number = 3,
      train_Six2GFP, test_Six2GFP,
      ref_cell_type_column = "type", dims = 1:30,
      train_title = "SIX2GFP", test_title_prefix = "Uchimura_no_DISTCD_LOH",
      n_neighbors = 8, skip_neighbors = TRUE,
      output_prefix_base = "six2gfp/subsampling/uchimura_no_DISTCD_LOH_anchors/",
      prediction_column_name = "predicted.type",
      use_cache = TRUE,
      return_anchors = TRUE
    )
  },
  label = "Uchimura no DIST_CD+LOH",
  expected_noised_number = 3
)
gc(verbose = FALSE)

# No CM+CM_DIV
uchi_no_cm_cmdiv_flow <- run_or_load_flow(
  "six2gfp/subsampling/uchimura_no_CM_CMDIV_anchors/flow_result.rds",
  function() {
    info <- create_subsampled_seurat(seurat_Uchimura_Humphreys_20, Uchimura_flow_anchors, c("CM", "CM_DIV"), 0.0)
    analyze_noise_impact_on_prediction(
      info$seurat, noised_number = 3,
      train_Six2GFP, test_Six2GFP,
      ref_cell_type_column = "type", dims = 1:30,
      train_title = "SIX2GFP", test_title_prefix = "Uchimura_no_CM_CMDIV",
      n_neighbors = 8, skip_neighbors = TRUE,
      output_prefix_base = "six2gfp/subsampling/uchimura_no_CM_CMDIV_anchors/",
      prediction_column_name = "predicted.type",
      use_cache = TRUE,
      return_anchors = TRUE
    )
  },
  label = "Uchimura no CM+CM_DIV",
  expected_noised_number = 3
)
gc(verbose = FALSE)

cat(sprintf("\n[UCHIMURA] Sections B1-B3 complete — total: %.1f min\n",
    as.numeric(difftime(Sys.time(), section_b_start, units = "mins"))))

# --- B4. Uchimura plotting ---
cat("\n######## Uchimura: Plotting ########\n")

# Single removals
uchimura_single_results <- list("Baseline" = Uchimura_flow_anchors)
for (type_name in names(uchimura_single_type_flows)) {
  uchimura_single_results[[paste0("No ", type_name)]] <- uchimura_single_type_flows[[type_name]]$flow
}
uchimura_single_plots <- plot_subsampling_comparison(
  uchimura_single_results,
  output_prefix = "six2gfp/subsampling/uchimura_single_removal_comparison"
)

# Multi-group removals
uchimura_multigroup_plots <- plot_subsampling_comparison(
  list(
    "Baseline"        = Uchimura_flow_anchors,
    "No DIST_CD+LOH"  = uchi_no_distcd_loh_flow,
    "No CM+CM_DIV"    = uchi_no_cm_cmdiv_flow
  ),
  output_prefix = "six2gfp/subsampling/uchimura_multigroup_comparison"
)

# ALL removals combined
uchimura_all_results <- c(
  list(
    "Baseline"        = Uchimura_flow_anchors,
    "No DIST_CD+LOH"  = uchi_no_distcd_loh_flow,
    "No CM+CM_DIV"    = uchi_no_cm_cmdiv_flow
  ),
  setNames(
    lapply(names(uchimura_single_type_flows), function(t) uchimura_single_type_flows[[t]]$flow),
    paste0("No ", names(uchimura_single_type_flows))
  )
)
uchimura_all_plots <- plot_subsampling_comparison(
  uchimura_all_results,
  output_prefix = "six2gfp/subsampling/uchimura_all_removal_comparison"
)

cat(sprintf("\n[UCHIMURA] Section B complete (incl. plotting) — total: %.1f min\n",
    as.numeric(difftime(Sys.time(), section_b_start, units = "mins"))))


# #############################################################################
# Section C: Freedman Subsampling (commented out)
# #############################################################################

# # --- C1. Baseline flow with anchors ---
# if (is.null(freedman_flow$initial_anchors)) {
#   freedman_flow_anchors <- run_or_load_flow(
#     "six2gfp/subsampling/freedman_with_anchors/flow_result.rds",
#     function() {
#       analyze_noise_impact_on_prediction(
#         freedman_seurat_obj, noised_number = 3,
#         train_Six2GFP, test_Six2GFP,
#         ref_cell_type_column = "type", dims = 1:15,
#         train_title = "SIX2GFP", test_title_prefix = "Freedman",
#         n_neighbors = 8, skip_neighbors = TRUE,
#         output_prefix_base = "six2gfp/subsampling/freedman_with_anchors/",
#         prediction_column_name = "predicted.type",
#         use_cache = TRUE,
#         return_anchors = TRUE
#       )
#     },
#     label = "Freedman baseline"
#   )
# } else {
#   freedman_flow_anchors <- freedman_flow
# }
#
# # --- C2. Single-type removals (loop over ALL types) ---
# freedman_single_types <- c("UM", "CM", "CM_DIV", "PODO", "PROX_1", "PROX_2", "LOH", "DIST_CD", "ENDO")
# freedman_single_type_flows <- list()
#
# for (idx in seq_along(freedman_single_types)) {
#   type_to_remove <- freedman_single_types[idx]
#   cat(sprintf("\n######## Freedman: Removing %s (%d/%d) ########\n", type_to_remove, idx, length(freedman_single_types)))
#   
#   flow_cache <- paste0("six2gfp/subsampling/freedman_no_", type_to_remove, "_anchors/flow_result.rds")
#   
#   flow <- run_or_load_flow(flow_cache, function() {
#     info <- create_subsampled_seurat(freedman_seurat_obj, freedman_flow_anchors, type_to_remove, 0.0)
#     analyze_noise_impact_on_prediction(
#       info$seurat, noised_number = 3,
#       train_Six2GFP, test_Six2GFP,
#       ref_cell_type_column = "type", dims = 1:15,
#       train_title = "SIX2GFP", test_title_prefix = paste0("Freedman_no_", type_to_remove),
#       n_neighbors = 8, skip_neighbors = TRUE,
#       output_prefix_base = paste0("six2gfp/subsampling/freedman_no_", type_to_remove, "_anchors/"),
#       prediction_column_name = "predicted.type",
#       use_cache = TRUE,
#       return_anchors = TRUE
#     )
#   }, label = paste0("Freedman no ", type_to_remove))
#   
#   comparison <- compare_flows(
#     original_flow   = freedman_flow_anchors,
#     subsampled_flow = flow,
#     removed_labels  = type_to_remove,
#     train_data      = train_Six2GFP
#   )
#   
#   freedman_single_type_flows[[type_to_remove]] <- list(
#     flow = flow, comparison = comparison
#   )
#   
#   rm(flow, comparison)
#   gc(verbose = FALSE)
# }
#
# # --- C3. Multi-group combinations (add as needed) ---
# # (No specific combinations for Freedman yet)
#
# # --- C4. Freedman plotting ---
# freedman_single_results <- list("Baseline" = freedman_flow_anchors)
# for (type_name in names(freedman_single_type_flows)) {
#   freedman_single_results[[paste0("No ", type_name)]] <- freedman_single_type_flows[[type_name]]$flow
# }
# freedman_single_plots <- plot_subsampling_comparison(
#   freedman_single_results,
#   output_prefix = "six2gfp/subsampling/freedman_single_removal_comparison"
# )
# freedman_all_plots <- plot_subsampling_comparison(
#   freedman_single_results,
#   output_prefix = "six2gfp/subsampling/freedman_all_removal_comparison"
# )


# =============================================================================
# Notes on Performance and Caching:
# =============================================================================
# 1. TOP-LEVEL CACHING: Each analyze_noise_impact_on_prediction() result is
#    saved as flow_result.rds in its output directory. On re-run, this loads
#    instantly instead of recomputing. If the script crashes, already-completed
#    removals are NOT recomputed.
#
# 2. INNER CACHING: use_cache = TRUE caches intermediate Seurat objects
#    (training data, test queries, no-noise run). These help when the flow
#    needs to actually run, but the top-level cache skips everything.
#
# 3. MEMORY CLEANUP: gc() is called after each removal to free memory.
#    With 26GB RAM + 15GB swap, this prevents OOM crashes.
#
# 4. TO FORCE RECOMPUTATION: Delete the flow_result.rds file for the
#    specific removal you want to recompute.
#
# Stability metrics computed:
#   - Stability:    diagonal of confusion matrix (recall / outflow only)
#   - Bi_Stability: 1 - (out + in) / original (penalizes both outflow & inflow)
#   - F1_Stability: 2TP / (2TP + FP + FN) (harmonic mean of precision & recall)