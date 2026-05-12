# =============================================================================
# Subsampling Pipeline for Stability Analysis
# =============================================================================
# This script tests whether Seurat's label transfer requires a critical mass
# of cells from a type to correctly classify them. Remove cell types (fully
# or partially) and re-classify to see how stability changes.
#
# Functions:
# 1. create_subsampled_seurat  — remove one or more cell types, keeping only
#                                 a specified fraction (e.g. 0 = full removal)
# 2. compare_flows             — compare original vs subsampled flow results
#                                 (PSS, stability, what gets predicted as removed types)
# 3. run_subsampling_analysis  — full workflow: subsample → classify → compare
#
# Usage: source("subsampling_pipeline.R") then call functions individually.
# See use_subsampling_pipeline.R for full workflow examples.
# =============================================================================

library(Seurat)
library(Matrix)


# =============================================================================
# 1. create_subsampled_seurat
# =============================================================================
#' Create a Seurat object where one or more target cell types are subsampled
#' (or fully removed). All other cells remain untouched.
#'
#' @param seurat_obj    Full Seurat object (e.g. atlas_object)
#' @param flow_result   Output of analyze_noise_impact_on_prediction()
#' @param target_labels Character vector of cell types to remove/subsample
#'                      (e.g. "PROX_2" or c("LOH", "CM"))
#' @param keep_fraction Numeric 0–1, fraction of target cells to keep.
#'                      0 = fully remove, 0.05 = keep 5%, etc.
#' @param seed          Integer, random seed for reproducibility
#'
#' @return List with:
#'   - seurat: Seurat object with subsampled/removed types
#'   - removal_log: data.frame with per-type removal statistics
#'   - kept_cells: named list of kept cell barcodes per type
#'   - removed_cells: named list of removed cell barcodes per type
#'   - target_labels: the target labels used
#'   - keep_fraction: the fraction used
create_subsampled_seurat <- function(seurat_obj,
                                     flow_result,
                                     target_labels,
                                     keep_fraction = 0.0,
                                     seed = 42) {
  set.seed(seed)
  
  # Ensure target_labels is a character vector
  target_labels <- as.character(target_labels)
  
  # Extract initial predictions from the flow result
  noised_pred <- flow_result$noised_prediction_matrix
  initial_predictions <- noised_pred[, 1]
  
  all_kept_cells <- list()
  all_removed_cells <- list()
  removal_log <- data.frame(
    type = character(),
    original_count = integer(),
    kept_count = integer(),
    removed_count = integer(),
    stringsAsFactors = FALSE
  )
  
  # Collect all cells to remove across all target types
  all_cells_to_remove <- character(0)
  
  for (label in target_labels) {
    # Find cells predicted as this label
    target_cells <- names(initial_predictions[initial_predictions == label])
    
    # Fallback: try rownames if names are empty
    if (length(target_cells) == 0) {
      matching_idx <- which(noised_pred[, 1] == label)
      target_cells <- rownames(noised_pred)[matching_idx]
    }
    
    if (length(target_cells) == 0) {
      cat(sprintf("Warning: No cells predicted as '%s' found, skipping.\n", label))
      next
    }
    
    # Validate cells exist in Seurat object
    target_cells <- intersect(target_cells, colnames(seurat_obj))
    if (length(target_cells) == 0) {
      cat(sprintf("Warning: No cells for '%s' found in seurat_obj, skipping.\n", label))
      next
    }
    
    original_count <- length(target_cells)
    n_keep <- round(original_count * keep_fraction)
    
    if (n_keep == 0) {
      kept <- character(0)
      removed <- target_cells
    } else {
      kept <- sample(target_cells, size = n_keep, replace = FALSE)
      removed <- setdiff(target_cells, kept)
    }
    
    all_kept_cells[[label]] <- kept
    all_removed_cells[[label]] <- removed
    all_cells_to_remove <- c(all_cells_to_remove, removed)
    
    removal_log <- rbind(removal_log, data.frame(
      type = label,
      original_count = original_count,
      kept_count = length(kept),
      removed_count = length(removed),
      stringsAsFactors = FALSE
    ))
    
    cat(sprintf("  '%s': %d total → keeping %d, removing %d\n",
                label, original_count, length(kept), length(removed)))
  }
  
  if (length(all_cells_to_remove) == 0) {
    stop("No cells to remove — check target_labels")
  }
  
  # Build subsampled Seurat
  cells_to_keep <- setdiff(colnames(seurat_obj), all_cells_to_remove)
  subsampled_seurat <- subset(seurat_obj, cells = cells_to_keep)
  
  cat(sprintf("\nSubsampled Seurat: %d → %d cells (removed %d total from %d type(s))\n",
              ncol(seurat_obj), ncol(subsampled_seurat),
              length(all_cells_to_remove), nrow(removal_log)))
  
  return(list(
    seurat = subsampled_seurat,
    removal_log = removal_log,
    kept_cells = all_kept_cells,
    removed_cells = all_removed_cells,
    target_labels = target_labels,
    keep_fraction = keep_fraction
  ))
}


# =============================================================================
# 2. compare_flows
# =============================================================================
#' Compare the original flow result with a subsampled flow result.
#' Reports per-type changes in PSS, stability, and what gets predicted as
#' the removed type(s).
#'
#' @param original_flow   Output of analyze_noise_impact_on_prediction() (full data)
#' @param subsampled_flow Output of analyze_noise_impact_on_prediction() (after removal)
#' @param removed_labels  Character vector of types that were removed
#' @param output_file     Optional path to save the comparison report as CSV
#' @param ref_cell_type_column Column name for cell types in the reference/training
#'                              data (needed for anchor comparison)
#' @param train_data       Training Seurat object (needed for anchor comparison
#'                          to look up cell types)
#'
#' @return List with:
#'   - pss_comparison: data.frame with PSS per type (original vs subsampled)
#'   - stability_comparison: data.frame with stability per type (original vs subsampled)
#'   - removed_type_counts: how many cells are predicted as each removed type
#'                          (in original vs subsampled), per noise iteration
#'   - prediction_shift: for all remaining cells, what changed
compare_flows <- function(original_flow,
                          subsampled_flow,
                          removed_labels,
                          ref_cell_type_column = "type",
                          output_file = NULL,
                          train_data = NULL) {
  
  cat("\n============================================================\n")
  cat("  Flow Comparison: Original vs After Removal\n")
  cat(sprintf("  Removed type(s): %s\n", paste(removed_labels, collapse = ", ")))
  cat("============================================================\n\n")
  
  # --- 1. PSS comparison (diagonal of PSS matrix) ---
  orig_pss <- original_flow$pss_matrix
  sub_pss <- subsampled_flow$pss_matrix
  
  # Use union so removed types appear with NA instead of disappearing
  all_types_pss <- union(rownames(orig_pss), rownames(sub_pss))
  
  pss_comparison <- data.frame(
    type = all_types_pss,
    pss_original = sapply(all_types_pss, function(t) {
      if (t %in% rownames(orig_pss) && t %in% colnames(orig_pss)) orig_pss[t, t] else NA
    }),
    pss_subsampled = sapply(all_types_pss, function(t) {
      if (t %in% rownames(sub_pss) && t %in% colnames(sub_pss)) sub_pss[t, t] else NA
    }),
    stringsAsFactors = FALSE
  )
  pss_comparison$pss_delta <- pss_comparison$pss_subsampled - pss_comparison$pss_original
  
  cat("--- PSS Comparison (per type, diagonal) ---\n")
  print(pss_comparison, row.names = FALSE)
  
  # --- 2. Stability comparison (diagonal of confusion matrix) ---
  orig_stab <- original_flow$stability_pred
  sub_stab <- subsampled_flow$stability_pred
  
  all_types_stab <- union(rownames(orig_stab), rownames(sub_stab))
  
  stability_comparison <- data.frame(
    type = all_types_stab,
    stability_original = sapply(all_types_stab, function(t) {
      if (t %in% rownames(orig_stab)) orig_stab[t, "Stability"] else NA
    }),
    stability_subsampled = sapply(all_types_stab, function(t) {
      if (t %in% rownames(sub_stab)) sub_stab[t, "Stability"] else NA
    }),
    stringsAsFactors = FALSE
  )
  stability_comparison$stability_delta <- stability_comparison$stability_subsampled - stability_comparison$stability_original
  
  cat("\n--- Stability Comparison (per type) ---\n")
  print(stability_comparison, row.names = FALSE)
  
  # --- 3. How many cells are predicted as the removed type(s)? ---
  orig_pred <- original_flow$noised_prediction_matrix
  sub_pred <- subsampled_flow$noised_prediction_matrix
  
  n_cols_orig <- ncol(orig_pred)
  n_cols_sub <- ncol(sub_pred)
  
  cat("\n--- Cells Predicted as Removed Type(s) ---\n")
  removed_type_counts <- list()
  for (label in removed_labels) {
    orig_counts <- sapply(1:n_cols_orig, function(col) sum(orig_pred[, col] == label))
    sub_counts <- sapply(1:n_cols_sub, function(col) sum(sub_pred[, col] == label))
    
    names(orig_counts) <- c("no_noise", paste0("noise_", 1:(n_cols_orig - 1)))
    names(sub_counts) <- c("no_noise", paste0("noise_", 1:(n_cols_sub - 1)))
    
    removed_type_counts[[label]] <- list(
      original = orig_counts,
      subsampled = sub_counts
    )
    
    cat(sprintf("\n  '%s' predictions:\n", label))
    cat(sprintf("    Original  : %s\n", paste(orig_counts, collapse = ", ")))
    cat(sprintf("    Subsampled: %s\n", paste(sub_counts, collapse = ", ")))
  }
  
  # --- 4. Prediction shift for remaining cells ---
  # Which cells are in both flows?
  common_cells <- intersect(rownames(orig_pred), rownames(sub_pred))
  
  orig_initial <- orig_pred[common_cells, 1]
  sub_initial <- sub_pred[common_cells, 1]
  
  changed_mask <- orig_initial != sub_initial
  n_changed <- sum(changed_mask)
  
  cat(sprintf("\n--- Prediction Shifts (no-noise, %d common cells) ---\n", length(common_cells)))
  cat(sprintf("  Unchanged: %d (%.1f%%)\n", 
              sum(!changed_mask), 100 * sum(!changed_mask) / length(common_cells)))
  cat(sprintf("  Changed:   %d (%.1f%%)\n",
              n_changed, 100 * n_changed / length(common_cells)))
  
  if (n_changed > 0) {
    shift_table <- table(
      from = orig_initial[changed_mask],
      to = sub_initial[changed_mask]
    )
    cat("\n  Transition table (from → to):\n")
    print(shift_table)
  }
  
  # Save to file if requested
  if (!is.null(output_file)) {
    write.csv(pss_comparison, paste0(output_file, "_pss.csv"), row.names = FALSE)
    write.csv(stability_comparison, paste0(output_file, "_stability.csv"), row.names = FALSE)
    cat(sprintf("\nSaved comparison to %s_pss.csv and %s_stability.csv\n",
                output_file, output_file))
  }
  
  # --- 5. Anchor comparison (if available) ---
  anchor_comparison <- NULL
  if (!is.null(original_flow$initial_anchors) && !is.null(subsampled_flow$initial_anchors)) {
    cat("\n--- Anchor Comparison ---\n")
    anchor_comparison <- compare_anchors(
      original_anchors = original_flow$initial_anchors,
      subsampled_anchors = subsampled_flow$initial_anchors,
      removed_labels = removed_labels,
      ref_cell_type_column = ref_cell_type_column,
      train_data = train_data,
      original_query_names = rownames(orig_pred),
      subsampled_query_names = rownames(sub_pred)
    )
  } else if (is.null(original_flow$initial_anchors) || is.null(subsampled_flow$initial_anchors)) {
    cat("\n--- Anchor Comparison: skipped (run with return_anchors=TRUE to enable) ---\n")
  }
  
  # --- 6. Ghost Prediction Analysis ---
  # When removed types reappear under noise (e.g., 0 CM → 1296 CM), trace
  # whether this is driven by pre-existing anchor connections.
  ghost_analysis <- NULL
  
  if (ncol(sub_pred) > 1) {
    cat("\n--- Ghost Prediction Analysis ---\n")
    cat("  (Cells predicted as a removed type under noise, but NOT in the no-noise run.)\n\n")
    
    sub_initial <- sub_pred[, 1]  # no-noise predictions
    ghost_analysis <- list()
    
    for (label in removed_labels) {
      # Find cells NOT predicted as this type without noise, but predicted as it WITH noise
      not_label_initial <- which(sub_initial != label)
      
      ghost_cells_per_iter <- list()
      for (col in 2:ncol(sub_pred)) {
        ghost_idx <- intersect(not_label_initial, which(sub_pred[, col] == label))
        ghost_cells_per_iter[[col - 1]] <- rownames(sub_pred)[ghost_idx]
      }
      
      # Union of ghost cells across all noise iterations
      all_ghost_cells <- unique(unlist(ghost_cells_per_iter))
      n_ghost <- length(all_ghost_cells)
      
      # How many were predicted as this type initially (already that type, not "ghost")
      n_initial_as_label <- sum(sub_initial == label)
      
      cat(sprintf("  '%s': %d cells initially predicted → %d ghost cells appear under noise\n",
                  label, n_initial_as_label, n_ghost))
      
      if (n_ghost == 0) {
        ghost_analysis[[label]] <- list(n_ghost = 0)
        next
      }
      
      # Level 1: What were ghost cells predicted as without noise?
      ghost_original_preds <- sub_initial[all_ghost_cells]
      from_table <- sort(table(ghost_original_preds), decreasing = TRUE)
      cat("    Where ghost cells came from (no-noise prediction → this type under noise):\n")
      for (i in seq_along(from_table)) {
        cat(sprintf("      %s: %d cells (%.1f%%)\n",
                    names(from_table)[i], from_table[i],
                    100 * from_table[i] / n_ghost))
      }
      
      label_ghost <- list(
        n_ghost = n_ghost,
        n_initial_as_label = n_initial_as_label,
        ghost_cells = all_ghost_cells,
        ghost_origin = from_table
      )
      
      # Level 2: Cross-reference with initial anchors (no-noise subsampled run)
      if (!is.null(subsampled_flow$initial_anchors) && !is.null(train_data)) {
        sub_anchor_profile <- anchor_comparison$sub_anchor_profile
        
        # Ghost cells that had at least one anchor to this ref type in the initial run
        ghost_with_anchor <- intersect(
          all_ghost_cells,
          unique(sub_anchor_profile$query_cell[which(sub_anchor_profile$ref_type == label)])
        )
        pct_with_anchor <- if (n_ghost > 0) 100 * length(ghost_with_anchor) / n_ghost else 0
        
        cat(sprintf("    Initial anchor link: %d/%d ghost cells (%.1f%%) had no-noise anchors to '%s' ref cells\n",
                    length(ghost_with_anchor), n_ghost, pct_with_anchor, label))
        
        # What were the dominant anchor types for ghost cells?
        ghost_in_anchors <- intersect(all_ghost_cells, unique(sub_anchor_profile$query_cell))
        if (length(ghost_in_anchors) > 0) {
          ghost_dominant <- tapply(
            sub_anchor_profile$ref_type[sub_anchor_profile$query_cell %in% ghost_in_anchors],
            sub_anchor_profile$query_cell[sub_anchor_profile$query_cell %in% ghost_in_anchors],
            function(types) names(sort(table(types), decreasing = TRUE))[1]
          )
          dominant_table <- sort(table(ghost_dominant), decreasing = TRUE)
          cat("    Dominant anchor types of ghost cells (no-noise run):\n")
          for (i in seq_along(dominant_table)) {
            cat(sprintf("      %s: %d cells\n", names(dominant_table)[i], dominant_table[i]))
          }
        }
        
        label_ghost$ghost_with_initial_anchor <- ghost_with_anchor
        label_ghost$pct_with_initial_anchor <- pct_with_anchor
      } else {
        cat("    (Run with return_anchors=TRUE and pass train_data for anchor cross-referencing)\n")
      }
      
      # Level 3: Cross-reference with noised anchors (per-iteration)
      if (!is.null(subsampled_flow$noised_anchors) && length(subsampled_flow$noised_anchors) > 0 &&
          !is.null(train_data)) {
        cat(sprintf("    Noised-run anchor analysis (checking if noise created new '%s' connections):\n", label))
        
        ghost_with_noised_anchor_counts <- integer(0)
        for (iter_name in names(subsampled_flow$noised_anchors)) {
          iter_anchors <- subsampled_flow$noised_anchors[[iter_name]]
          iter_profile <- get_anchor_profile(iter_anchors, ref_cell_type_column, train_data,
                                               original_query_names = rownames(sub_pred))
          
          iter_idx <- as.integer(sub("iter", "", iter_name))
          iter_ghost_cells <- ghost_cells_per_iter[[iter_idx]]
          
          if (length(iter_ghost_cells) > 0) {
            ghost_with_noised_anchor <- intersect(
              iter_ghost_cells,
              unique(iter_profile$query_cell[which(iter_profile$ref_type == label)])
            )
            pct <- 100 * length(ghost_with_noised_anchor) / length(iter_ghost_cells)
            cat(sprintf("      %s: %d/%d ghost cells (%.1f%%) have noised anchors to '%s' ref cells\n",
                        iter_name, length(ghost_with_noised_anchor),
                        length(iter_ghost_cells), pct, label))
            ghost_with_noised_anchor_counts[iter_name] <- length(ghost_with_noised_anchor)
          }
        }
        label_ghost$ghost_with_noised_anchor_counts <- ghost_with_noised_anchor_counts
      } else if (is.null(subsampled_flow$noised_anchors) || length(subsampled_flow$noised_anchors) == 0) {
        cat("    (Re-run with updated FULL_FUNCTION.R + return_anchors=TRUE for noised anchor analysis)\n")
      }
      
      ghost_analysis[[label]] <- label_ghost
    }
  }
  
  return(list(
    pss_comparison = pss_comparison,
    stability_comparison = stability_comparison,
    removed_type_counts = removed_type_counts,
    common_cells = common_cells,
    n_changed = n_changed,
    shift_table = if (n_changed > 0) shift_table else NULL,
    anchor_comparison = anchor_comparison,
    ghost_analysis = ghost_analysis
  ))
}


# =============================================================================
# 2b-helper. get_anchor_profile
# =============================================================================
#' Extract per-query-cell anchor composition from an AnchorSet.
#' Used by compare_anchors() and the ghost prediction analysis in compare_flows().
#'
#' @param anchor_set      A Seurat AnchorSet object
#' @param ref_cell_type_col Column name for cell types in the reference data
#' @param train_obj        Training Seurat object (for cell type lookup)
#'
#' @return data.frame with columns: cell1, cell2, score, ref_cell, query_cell, ref_type
get_anchor_profile <- function(anchor_set, ref_cell_type_col, train_obj = NULL,
                               original_query_names = NULL) {
  anchor_df <- as.data.frame(anchor_set@anchors)
  
  # Get reference and query cell names from the anchor set
  ref_cells <- anchor_set@reference.cells
  query_cells <- anchor_set@query.cells
  
  # Map anchor indices to cell names
  anchor_df$ref_cell <- ref_cells[anchor_df$cell1]
  anchor_df$query_cell <- query_cells[anchor_df$cell2]
  
  # --- Normalize query cell names back to original names ---
  # Seurat appends suffixes (_1 for ref, _2 for query) in FindTransferAnchors.
  # This maps them back so they match rownames(noised_prediction_matrix).
  if (!is.null(original_query_names)) {
    direct_match <- anchor_df$query_cell %in% original_query_names
    
    if (sum(direct_match) < nrow(anchor_df) * 0.5) {
      # Strategy 2: Strip Seurat-added suffixes
      stripped <- sub("_\\d+$", "", anchor_df$query_cell)
      if (sum(stripped %in% original_query_names) > sum(direct_match)) {
        anchor_df$query_cell <- stripped
      } else if (length(query_cells) == length(original_query_names)) {
        # Strategy 3: Positional mapping
        positional_map <- setNames(original_query_names, query_cells)
        anchor_df$query_cell <- as.character(positional_map[anchor_df$query_cell])
      }
    }
  }
  
  # --- Build a type lookup for reference cells ---
  if (!is.null(train_obj)) {
    train_types <- as.character(train_obj@meta.data[[ref_cell_type_col]])
    train_names <- colnames(train_obj)
    
    # Strategy 1: Direct name match
    type_lookup <- setNames(train_types, train_names)
    matched <- type_lookup[anchor_df$ref_cell]
    
    # Strategy 2: Strip Seurat-added suffixes (e.g., _1, _2) and retry
    if (sum(is.na(matched)) > length(matched) * 0.5) {
      stripped <- sub("_\\d+$", "", anchor_df$ref_cell)
      matched <- type_lookup[stripped]
    }
    
    # Strategy 3: Positional mapping (ref_cells should be in same order as train cells)
    if (sum(is.na(matched)) > length(matched) * 0.5 && length(ref_cells) == ncol(train_obj)) {
      positional_lookup <- setNames(train_types, ref_cells)
      matched <- positional_lookup[anchor_df$ref_cell]
    }
    
    anchor_df$ref_type <- as.character(matched)
  } else {
    # Fallback: try from the anchor set's object list
    ref_obj <- anchor_set@object.list[[1]]
    if (is.null(ref_obj) || is.null(ref_obj@meta.data[[ref_cell_type_col]])) {
      stop(paste0("Cannot look up cell types: train_data not provided and ",
                  "anchor_set@object.list[[1]] does not contain '",
                  ref_cell_type_col, "' metadata.\n",
                  "Pass train_data to compare_flows() or compare_anchors()."))
    }
    ref_types <- ref_obj@meta.data[[ref_cell_type_col]]
    names(ref_types) <- colnames(ref_obj)
    anchor_df$ref_type <- as.character(ref_types[anchor_df$ref_cell])
  }
  
  return(anchor_df)
}


# =============================================================================
# 2b. compare_anchors
# =============================================================================
#' Compare transfer anchors between original and subsampled flows.
#' Shows which reference cell types query cells anchor to, and how this changes.
#'
#' @param original_anchors   AnchorSet from original flow (flow$initial_anchors)
#' @param subsampled_anchors AnchorSet from subsampled flow
#' @param removed_labels     Character vector of removed cell types
#' @param ref_cell_type_column Column name for cell types in the reference data
#' @param train_data          Training Seurat object used to look up cell types
#'
#' @return List with anchor composition summaries
compare_anchors <- function(original_anchors,
                            subsampled_anchors,
                            removed_labels,
                            ref_cell_type_column = "type",
                            train_data = NULL,
                            original_query_names = NULL,
                            subsampled_query_names = NULL) {
  
  orig_profile <- get_anchor_profile(original_anchors, ref_cell_type_column, train_data,
                                     original_query_names = original_query_names)
  sub_profile <- get_anchor_profile(subsampled_anchors, ref_cell_type_column, train_data,
                                     original_query_names = subsampled_query_names)
  
  cat(sprintf("  Total anchor pairs — Original: %d, After removal: %d\n",
              nrow(orig_profile), nrow(sub_profile)))
  cat("  (Each pair links one query cell to one reference cell.\n")
  cat("   A single query cell can participate in multiple pairs.)\n")
  
  # --- 1. Anchor count by reference type ---
  orig_ref_counts <- table(orig_profile$ref_type)
  sub_ref_counts <- table(sub_profile$ref_type)
  
  all_ref_types <- union(names(orig_ref_counts), names(sub_ref_counts))
  ref_type_comparison <- data.frame(
    ref_type = all_ref_types,
    anchors_original = as.integer(orig_ref_counts[all_ref_types]),
    anchors_subsampled = as.integer(sub_ref_counts[all_ref_types]),
    stringsAsFactors = FALSE
  )
  ref_type_comparison[is.na(ref_type_comparison)] <- 0
  ref_type_comparison$delta <- ref_type_comparison$anchors_subsampled - ref_type_comparison$anchors_original
  ref_type_comparison <- ref_type_comparison[order(-abs(ref_type_comparison$delta)), ]
  
  cat("\n  Anchor PAIRS by reference cell type (how many pairs involve each ref type):\n")
  cat("  (Note: the reference is unchanged — only the query was subsampled.\n")
  cat("   Fewer query cells → fewer total pairs, especially to ref types\n")
  cat("   that the removed query cells were heavily anchoring to.)\n")
  print(ref_type_comparison, row.names = FALSE)
  
  # --- 2. Per-query-cell: dominant anchor type (most frequent ref type) ---
  get_dominant_type <- function(profile) {
    dominant <- tapply(profile$ref_type, profile$query_cell, function(types) {
      names(sort(table(types), decreasing = TRUE))[1]
    })
    return(dominant)
  }
  
  orig_dominant <- get_dominant_type(orig_profile)
  sub_dominant <- get_dominant_type(sub_profile)
  
  # Compare for common query cells
  common_query <- intersect(names(orig_dominant), names(sub_dominant))
  
  if (length(common_query) > 0) {
    dominant_changed <- orig_dominant[common_query] != sub_dominant[common_query]
    n_dominant_changed <- sum(dominant_changed, na.rm = TRUE)
    
    cat(sprintf("\n  Dominant anchor type changed: %d/%d common query cells (%.1f%%)\n",
                n_dominant_changed, length(common_query),
                100 * n_dominant_changed / length(common_query)))
    cat("  (For each query cell, the 'dominant' type = the most frequent ref type\n")
    cat("   among its anchors. This measures whether REMAINING cells changed\n")
    cat("   which ref type they are most strongly linked to.)\n")
    
    if (n_dominant_changed > 0 && n_dominant_changed <= 500) {
      shift <- table(
        from = orig_dominant[common_query][dominant_changed],
        to = sub_dominant[common_query][dominant_changed]
      )
      cat("\n  Dominant anchor type shifts (from → to):\n")
      print(shift)
    }
  }
  
  # --- 3. Cells that anchored to the removed type ---
  cat("\n  Unique QUERY CELLS with at least one anchor to each removed ref type:\n")
  cat("  (Unlike the pair count above, this counts each query cell only once,\n")
  cat("   regardless of how many anchor pairs it has to that ref type.)\n")
  for (label in removed_labels) {
    orig_anchored_to_removed <- unique(orig_profile$query_cell[which(orig_profile$ref_type == label)])
    sub_anchored_to_removed <- unique(sub_profile$query_cell[which(sub_profile$ref_type == label)])
    
    cat(sprintf("    ref '%s' — Original: %d query cells, After removal: %d query cells\n",
                label, length(orig_anchored_to_removed), length(sub_anchored_to_removed)))
  }
  
  return(list(
    ref_type_comparison = ref_type_comparison,
    orig_anchor_profile = orig_profile,
    sub_anchor_profile = sub_profile,
    common_query_cells = common_query,
    n_dominant_changed = if (exists("n_dominant_changed")) n_dominant_changed else NA
  ))
}


# =============================================================================
# 3. run_subsampling_analysis
# =============================================================================
#' Full subsampling workflow: remove type(s), re-classify, and compare.
#'
#' @param seurat_obj        Full Seurat object
#' @param flow_result       Original flow result
#' @param target_labels     Character vector of types to remove
#' @param keep_fraction     Numeric 0–1, fraction to keep (0 = full removal)
#' @param train_data        Labeled training Seurat
#' @param test_data         Labeled test Seurat
#' @param ref_cell_type_column Column name for cell types in training data
#' @param dims              PCA dimensions
#' @param noised_number     Number of noise iterations
#' @param output_prefix_base Output directory
#' @param seed              Random seed
#' @param ...               Additional args for analyze_noise_impact_on_prediction()
#'
#' @return List with subsample_info, flow, and comparison
run_subsampling_analysis <- function(seurat_obj,
                                     flow_result,
                                     target_labels,
                                     keep_fraction = 0.0,
                                     train_data,
                                     test_data,
                                     ref_cell_type_column = "type",
                                     dims = 1:30,
                                     noised_number = 1,
                                     output_prefix_base = "six2gfp/subsampling/",
                                     seed = 42,
                                     save_rds = TRUE,
                                     ...) {
  
  # If cache file exists, load it directly to save time (especially after reboot)
  rds_file <- paste0(output_prefix_base, "subsampling_result.rds")
  if (save_rds && file.exists(rds_file)) {
    cat(sprintf("\n========== Loading cached result from %s ==========\n", rds_file))
    return(readRDS(rds_file))
  }
  
  # Step 1: Create subsampled Seurat
  cat("\n========== Step 1: Subsampling ==========\n")
  subsample_info <- create_subsampled_seurat(
    seurat_obj = seurat_obj,
    flow_result = flow_result,
    target_labels = target_labels,
    keep_fraction = keep_fraction,
    seed = seed
  )
  
  # Step 2: Run classification on subsampled data
  cat("\n========== Step 2: Re-classifying ==========\n")
  subsampled_flow <- analyze_noise_impact_on_prediction(
    subsample_info$seurat,
    noised_number = noised_number,
    train_data,
    test_data,
    ref_cell_type_column = ref_cell_type_column,
    dims = dims,
    output_prefix_base = output_prefix_base,
    ...
  )
  
  # Step 3: Compare flows
  cat("\n========== Step 3: Comparing flows ==========\n")
  comparison <- compare_flows(
    original_flow = flow_result,
    subsampled_flow = subsampled_flow,
    removed_labels = target_labels,
    output_file = paste0(output_prefix_base, "comparison")
  )
  
  final_result <- list(
    subsample_info = subsample_info,
    flow = subsampled_flow,
    comparison = comparison
  )
  
  if (save_rds) {
    saveRDS(final_result, rds_file)
    cat(sprintf("\nSaved full subsampling result to %s\n", rds_file))
  }
  
  return(final_result)
}
