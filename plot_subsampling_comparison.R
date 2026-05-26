library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(patchwork)

#' Plot Subsampling Comparison Across Multiple Removal Groups
#' 
#' This function creates bar plots comparing the diagonal elements
#' (PSS self-similarity and stability) across multiple subsampling runs
#' where different cell type groups were removed from the SAME dataset.
#' 
#' @param run_results_list A named list of subsampling run results, where each
#'   element is the output from run_subsampling_analysis()$flow (i.e. the flow
#'   result after removal), OR directly from analyze_noise_impact_on_prediction().
#'   Names should describe what was removed (e.g., "No CM+CM_DIV", "No PROX_2").
#'   Include the original (no removal) run as a baseline — name it e.g. "Baseline".
#' @param output_prefix Path prefix for saving output plots
#' @param plot_title_pss Title for the PSS bar plot
#' @param plot_title_stability Title for the stability bar plot
#' 
#' @return A list containing both ggplot objects and combined data
#' 
#' @examples
#' # After running your subsampling analyses:
#' results_list <- list(
#'   "Baseline" = cell_atlas_flow_anchors,
#'   "No CM+CM_DIV+PROX_1" = no_cm_cmdiv_prox1_flow_anchors,
#'   "No PROX_2" = no_prox2_flow$flow,
#'   "No DIST_CD" = no_distcd_flow$flow
#' )
#' subsampling_plots <- plot_subsampling_comparison(
#'   results_list, 
#'   output_prefix = "six2gfp/subsampling/atlas_removal_comparison"
#' )
plot_subsampling_comparison <- function(
    run_results_list,
    output_prefix = "subsampling_comparison",
    plot_title_pss = "PSS Self-Similarity: Effect of Group Removal",
    plot_title_stability = "Stability: Effect of Group Removal"
) {
  
  # Validate input
  if (length(run_results_list) < 1) {
    stop("run_results_list must contain at least one run result")
  }

  # 1. Collect ALL cell types across runs (union, not intersection)
  # Using union ensures that cell types removed in some runs still appear
  # in the plot from runs where they existed.
  all_cell_types <- NULL
  
  for (run_name in names(run_results_list)) {
      run_result <- run_results_list[[run_name]]
      
      if (is.null(run_result$stability_pred)) {
          warning(paste("Run", run_name, "missing stability_pred, skipping"))
          next
      }
      
      current_types <- rownames(run_result$stability_pred)
      
      if (is.null(all_cell_types)) {
          all_cell_types <- current_types
      } else {
          all_cell_types <- union(all_cell_types, current_types)
      }
  }
  
  if (length(all_cell_types) == 0) {
      stop("No cell types found across the provided runs.")
  }
  print(paste("All cell types found:", paste(all_cell_types, collapse=", ")))

  
  # Extract diagonal data from all runs
  combined_data <- data.frame()
  
  for (run_name in names(run_results_list)) {
    run_result <- run_results_list[[run_name]]
    
    # Check if necessary data exists
    if (is.null(run_result$stability_pred)) {
      warning(paste("Run", run_name, "missing stability_pred, skipping"))
      next
    }
    
    # Extract cell types and values
    cell_types <- rownames(run_result$stability_pred)
    pss_values <- run_result$stability_pred$PSS
    stability_values <- run_result$stability_pred$Stability
    bi_stability_values <- if (!is.null(run_result$stability_pred$Bi_Stability)) {
      run_result$stability_pred$Bi_Stability
    } else {
      rep(NA, length(cell_types))  # backward compat with old flow results
    }
    f1_stability_values <- if (!is.null(run_result$stability_pred$F1_Stability)) {
      run_result$stability_pred$F1_Stability
    } else {
      rep(NA, length(cell_types))  # backward compat with old flow results
    }
    
    # Try to extract cell counts from noised_prediction_matrix
    cell_counts <- rep(NA, length(cell_types))
    
    if (!is.null(run_result$noised_prediction_matrix)) {
      prediction_vector <- run_result$noised_prediction_matrix[, 1]
      count_table <- table(prediction_vector)
      for (i in seq_along(cell_types)) {
        if (cell_types[i] %in% names(count_table)) {
          cell_counts[i] <- count_table[cell_types[i]]
        }
      }
      message(paste0("  ", run_name, ": Extracted counts from noised_prediction_matrix"))
    } else {
      message(paste0("  ", run_name, ": Could not find cell count data"))
    }
    
    # Extract per-run stability if available (for error bars)
    stability_sd_values <- rep(NA, length(cell_types))
    stability_mean_values <- stability_values  # default: use single-run values
    bi_stability_sd_values <- rep(NA, length(cell_types))
    bi_stability_mean_values <- bi_stability_values
    f1_stability_sd_values <- rep(NA, length(cell_types))
    f1_stability_mean_values <- f1_stability_values
    
    if (!is.null(run_result$stability_per_run) && ncol(run_result$stability_per_run) > 1) {
      for (i in seq_along(cell_types)) {
        if (cell_types[i] %in% rownames(run_result$stability_per_run)) {
          per_run_vals <- run_result$stability_per_run[cell_types[i], ]
          stability_mean_values[i] <- mean(per_run_vals, na.rm = TRUE)
          stability_sd_values[i] <- sd(per_run_vals, na.rm = TRUE)
        }
      }
      message(paste0("  ", run_name, ": Using mean/SD from ", ncol(run_result$stability_per_run), " noise runs"))
    }
    if (!is.null(run_result$bi_stability_per_run) && ncol(run_result$bi_stability_per_run) > 1) {
      for (i in seq_along(cell_types)) {
        if (cell_types[i] %in% rownames(run_result$bi_stability_per_run)) {
          per_run_vals <- run_result$bi_stability_per_run[cell_types[i], ]
          bi_stability_mean_values[i] <- mean(per_run_vals, na.rm = TRUE)
          bi_stability_sd_values[i] <- sd(per_run_vals, na.rm = TRUE)
        }
      }
    }
    if (!is.null(run_result$f1_stability_per_run) && ncol(run_result$f1_stability_per_run) > 1) {
      for (i in seq_along(cell_types)) {
        if (cell_types[i] %in% rownames(run_result$f1_stability_per_run)) {
          per_run_vals <- run_result$f1_stability_per_run[cell_types[i], ]
          f1_stability_mean_values[i] <- mean(per_run_vals, na.rm = TRUE)
          f1_stability_sd_values[i] <- sd(per_run_vals, na.rm = TRUE)
        }
      }
    }
    
    # Create data frame for this run
    run_data <- data.frame(
      cell_type = cell_types,
      pss = pss_values,
      stability = stability_mean_values,
      stability_sd = stability_sd_values,
      bi_stability = bi_stability_mean_values,
      bi_stability_sd = bi_stability_sd_values,
      f1_stability = f1_stability_mean_values,
      f1_stability_sd = f1_stability_sd_values,
      removed_group = run_name,
      n_cells = cell_counts,
      stringsAsFactors = FALSE
    )
    
    combined_data <- rbind(combined_data, run_data)
  }

  
  # Check if we have data
  if (nrow(combined_data) == 0) {
    stop("No valid data extracted from any run")
  }

  # Filter to keep only common cell types
  combined_data <- combined_data %>% 
      filter(cell_type %in% all_cell_types)

  # Enforce specific order for cell types
  desired_order <- c("UM","CM","CM_DIV","PODO","PROX_1","PROX_2","LOH","DIST_CD","ENDO","MACROPHAG")
  # Retrieve intersection of desired order and available (common) types to avoid NA levels
  final_levels <- intersect(desired_order, all_cell_types)
  
  if (length(final_levels) < length(all_cell_types)) {
      warning("Some common cell types are missing from the manual ordering list and will be excluded from the plot.")
      combined_data <- combined_data %>% filter(cell_type %in% final_levels)
  }

  combined_data$cell_type <- factor(combined_data$cell_type, levels = final_levels)
  
  # Preserve the input order for the legend (first entry = baseline, rest = removals)
  combined_data$removed_group <- factor(combined_data$removed_group, 
                                         levels = names(run_results_list))
  
  # Fill missing cell_type × removed_group combinations with 0
  # (removed cell types get 0-height bars instead of disappearing)
  combined_data <- combined_data %>%
      complete(cell_type, removed_group, 
               fill = list(pss = 0, stability = 0, stability_sd = 0,
                           bi_stability = 0, bi_stability_sd = 0,
                           f1_stability = 0, f1_stability_sd = 0, n_cells = 0))
  
  # Define colors for different removal groups
  num_groups <- length(unique(combined_data$removed_group))
  group_colors <- scales::hue_pal()(num_groups)
  names(group_colors) <- levels(combined_data$removed_group)
  
  # ===== Plot 1: PSS values comparison =====
  
  plot_pss_comparison <- ggplot(combined_data, aes(x = cell_type, y = pss, fill = removed_group)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.8) +
    geom_text(aes(label = ifelse(!is.na(n_cells), paste0("n=", n_cells), ""), y = 0), 
              position = position_dodge(width = 0.8), 
              vjust = 0.5, hjust = 0.5, size = 2.5, angle = 90) +
    scale_fill_manual(values = group_colors, name = "Removed Group") +
    labs(
      title = plot_title_pss,
      x = "Cell Type",
      y = "Prediction Specificity Score (PSS)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 13, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 13),
      plot.title = element_text(size = 15, face = "bold")
    )
  
  # Save plot 1
  ggsave(
    paste0(output_prefix, "_pss_by_celltype_comparison.svg"),
    plot = plot_pss_comparison,
    width = 12, height = 8, units = "in"
  )
  
  
  # ===== Plot 2: Stability values comparison =====
  
  # Check if we have SD data for error bars
  has_error_bars <- any(!is.na(combined_data$stability_sd))
  
  plot_stability_comparison <- ggplot(combined_data, aes(x = cell_type, y = stability, fill = removed_group)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.8) +
    geom_text(aes(label = ifelse(!is.na(n_cells), paste0("n=", n_cells), ""), y = 0), 
              position = position_dodge(width = 0.8), 
              vjust = 0.5, hjust = 0, size = 2.5, angle = 90) +
    { if (has_error_bars) 
        geom_errorbar(aes(ymin = stability - stability_sd, ymax = stability + stability_sd),
                      position = position_dodge(width = 0.8), width = 0.25, linewidth = 0.4)
    } +
    scale_fill_manual(values = group_colors, name = "Removed Group") +
    labs(
      title = plot_title_stability,
      x = "Cell Type",
      y = "Stability (Fraction of cells remaining same type)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 13, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 13),
      plot.title = element_text(size = 15, face = "bold")
    )
  
  # Save plot 2
  ggsave(
    paste0(output_prefix, "_stability_by_celltype_comparison.svg"),
    plot = plot_stability_comparison,
    width = 12, height = 8, units = "in"
  )
  
  # ===== Plot 3: Bidirectional Stability comparison =====
  # Bi_Stability = max(0, 1 - (out + in) / original)
  # Penalizes both cells leaving AND cells arriving from other types
  
  has_bi_stability <- any(!is.na(combined_data$bi_stability) & combined_data$bi_stability != 0)
  has_bi_error_bars <- any(!is.na(combined_data$bi_stability_sd))
  
  if (has_bi_stability) {
    plot_bi_stability_comparison <- ggplot(combined_data, aes(x = cell_type, y = bi_stability, fill = removed_group)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.8) +
      geom_text(aes(label = ifelse(!is.na(n_cells) & n_cells > 0, paste0("n=", n_cells), ""), y = 0), 
                position = position_dodge(width = 0.8), 
                vjust = 0.5, hjust = 0, size = 2.5, angle = 90) +
      { if (has_bi_error_bars) 
          geom_errorbar(aes(ymin = bi_stability - bi_stability_sd, ymax = bi_stability + bi_stability_sd),
                        position = position_dodge(width = 0.8), width = 0.25, linewidth = 0.4)
      } +
      scale_fill_manual(values = group_colors, name = "Removed Group") +
      labs(
        title = "Bidirectional Stability: Effect of Group Removal",
        subtitle = "Penalizes both outflow (cells leaving) and inflow (cells arriving from other types)",
        x = "Cell Type",
        y = "Bi-Stability: 1 - (out + in) / original"
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 13),
        plot.title = element_text(size = 15, face = "bold")
      )
    
    ggsave(
      paste0(output_prefix, "_bi_stability_by_celltype_comparison.svg"),
      plot = plot_bi_stability_comparison,
      width = 12, height = 8, units = "in"
    )
  } else {
    plot_bi_stability_comparison <- NULL
    message("Bi_Stability not available in flow results (re-run with updated FULL_FUNCTION.R)")
  }
  
  # ===== Plot 3b: F1 Stability comparison =====
  has_f1_stability <- any(!is.na(combined_data$f1_stability) & combined_data$f1_stability != 0)
  has_f1_error_bars <- any(!is.na(combined_data$f1_stability_sd))
  
  if (has_f1_stability) {
    plot_f1_stability_comparison <- ggplot(combined_data, aes(x = cell_type, y = f1_stability, fill = removed_group)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.8) +
      geom_text(aes(label = ifelse(!is.na(n_cells) & n_cells > 0, paste0("n=", n_cells), ""), y = 0), 
                position = position_dodge(width = 0.8), 
                vjust = 0.5, hjust = 0, size = 2.5, angle = 90) +
      { if (has_f1_error_bars) 
          geom_errorbar(aes(ymin = f1_stability - f1_stability_sd, ymax = f1_stability + f1_stability_sd),
                        position = position_dodge(width = 0.8), width = 0.25, linewidth = 0.4)
      } +
      scale_fill_manual(values = group_colors, name = "Removed Group") +
      labs(
        title = "F1 Stability: Effect of Group Removal",
        subtitle = "Harmonic mean of precision and recall: 2TP/(2TP+FP+FN)",
        x = "Cell Type",
        y = "F1 Stability"
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 13),
        plot.title = element_text(size = 15, face = "bold")
      )
    
    ggsave(
      paste0(output_prefix, "_f1_stability_by_celltype_comparison.svg"),
      plot = plot_f1_stability_comparison,
      width = 12, height = 8, units = "in"
    )
  } else {
    plot_f1_stability_comparison <- NULL
    message("F1_Stability not available in flow results (re-run with updated FULL_FUNCTION.R)")
  }
  
  # ===== Plot 4: Per-run stability barplots =====
  max_runs <- 0
  for (run_name in names(run_results_list)) {
    spr <- run_results_list[[run_name]]$stability_per_run
    if (!is.null(spr)) max_runs <- max(max_runs, ncol(spr))
  }
  
  per_run_plots <- list()
  if (max_runs > 1) {
    for (run_idx in 1:max_runs) {
      run_data_all <- data.frame()
      for (run_name in names(run_results_list)) {
        spr <- run_results_list[[run_name]]$stability_per_run
        if (!is.null(spr) && run_idx <= ncol(spr)) {
          cell_types_run <- rownames(spr)
          run_df <- data.frame(
            cell_type = cell_types_run,
            stability = spr[, run_idx],
            removed_group = run_name,
            stringsAsFactors = FALSE
          )
          run_data_all <- rbind(run_data_all, run_df)
        }
      }
      
      # Filter and order like the main plot
      run_data_all <- run_data_all %>% filter(cell_type %in% final_levels)
      run_data_all$cell_type <- factor(run_data_all$cell_type, levels = final_levels)
      run_data_all$removed_group <- factor(run_data_all$removed_group, 
                                            levels = names(run_results_list))
      
      p_run <- ggplot(run_data_all, aes(x = cell_type, y = stability, fill = removed_group)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.8) +
        scale_fill_manual(values = group_colors, name = "Removed Group") +
        labs(
          title = paste0("Stability by Cell Type - Run ", run_idx),
          x = "Cell Type",
          y = "Stability (Fraction of cells remaining same type)"
        ) +
        theme_minimal() +
        theme(
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 13, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 13),
          plot.title = element_text(size = 15, face = "bold")
        )
      
      ggsave(
        paste0(output_prefix, "_stability_run_", run_idx, ".svg"),
        plot = p_run,
        width = 12, height = 8, units = "in"
      )
      per_run_plots[[paste0("run_", run_idx)]] <- p_run
    }
  }
  
  
  # Print summary statistics
  cat("\n=== Summary Statistics ===\n")
  summary_stats <- combined_data %>%
    group_by(removed_group) %>%
    summarise(
      n_cell_types = n(),
      mean_pss = mean(pss, na.rm = TRUE),
      sd_pss = sd(pss, na.rm = TRUE),
      mean_stability = mean(stability, na.rm = TRUE),
      sd_stability = sd(stability, na.rm = TRUE),
      mean_bi_stability = mean(bi_stability, na.rm = TRUE),
      sd_bi_stability = sd(bi_stability, na.rm = TRUE),
      mean_f1_stability = mean(f1_stability, na.rm = TRUE),
      sd_f1_stability = sd(f1_stability, na.rm = TRUE),
      .groups = "drop"
    )
  print(summary_stats)
  
  # Compare stability ranking per cell type
  compare_stability_ranking <- function(data, run_results_list) {
    cat("\n=== Stability Ranking Per Cell Type ===\n")
    
    # Mean stability table
    ranking_df <- data %>%
      select(cell_type, removed_group, stability) %>%
      pivot_wider(names_from = removed_group, values_from = stability)
    
    group_names <- setdiff(colnames(ranking_df), "cell_type")
    
    ranking_df$max_group <- apply(ranking_df[, group_names], 1, function(row) {
      group_names[which.max(row)]
    })
    ranking_df$max_value <- apply(ranking_df[, group_names], 1, max, na.rm = TRUE)
    
    cat("\n--- Mean Stability ---\n")
    print(as.data.frame(ranking_df))
    
    # Summary: how many cell types each group "wins"
    cat("\n--- Stability wins per removal group ---\n")
    wins <- table(ranking_df$max_group)
    print(wins)
    
    # Per-run details for each group
    for (run_name in names(run_results_list)) {
      run_result <- run_results_list[[run_name]]
      if (!is.null(run_result$stability_per_run) && ncol(run_result$stability_per_run) > 1) {
        cat(paste0("\n--- ", run_name, ": Per-Run Stability ---\n"))
        per_run <- run_result$stability_per_run
        # Add mean and SD columns
        per_run_df <- as.data.frame(per_run)
        per_run_df$mean <- rowMeans(per_run, na.rm = TRUE)
        per_run_df$sd <- apply(per_run, 1, sd, na.rm = TRUE)
        print(round(per_run_df, 4))
      }
    }
    cat("\n")
    
    return(ranking_df)
  }
  
  stability_ranking <- compare_stability_ranking(combined_data, run_results_list)
  
  # ===== Plot 5: Delta plots — change from baseline =====
  # If the first entry is a baseline, compute deltas for PSS, Stability, and Bi_Stability
  baseline_name <- names(run_results_list)[1]
  if (length(run_results_list) > 1) {
    baseline_data <- combined_data %>% 
      filter(removed_group == baseline_name) %>%
      select(cell_type, pss_baseline = pss, stability_baseline = stability,
             bi_stability_baseline = bi_stability, f1_stability_baseline = f1_stability)
    
    delta_data <- combined_data %>%
      filter(removed_group != baseline_name) %>%
      left_join(baseline_data, by = "cell_type") %>%
      mutate(
        pss_delta = pss - pss_baseline,
        stability_delta = stability - stability_baseline,
        bi_stability_delta = bi_stability - bi_stability_baseline,
        f1_stability_delta = f1_stability - f1_stability_baseline
      )
    
    # PSS delta plot
    plot_pss_delta <- ggplot(delta_data, aes(x = cell_type, y = pss_delta, fill = removed_group)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
      scale_fill_manual(values = group_colors[names(group_colors) != baseline_name], 
                        name = "Removed Group") +
      labs(
        title = paste0("PSS Change vs Baseline (", baseline_name, ")"),
        x = "Cell Type",
        y = "ΔPSS (removal − baseline)"
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 13),
        plot.title = element_text(size = 15, face = "bold")
      )
    
    ggsave(
      paste0(output_prefix, "_pss_delta_vs_baseline.svg"),
      plot = plot_pss_delta,
      width = 12, height = 8, units = "in"
    )
    
    # Stability delta plot
    plot_stability_delta <- ggplot(delta_data, aes(x = cell_type, y = stability_delta, fill = removed_group)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
      scale_fill_manual(values = group_colors[names(group_colors) != baseline_name], 
                        name = "Removed Group") +
      labs(
        title = paste0("Stability Change vs Baseline (", baseline_name, ")"),
        x = "Cell Type",
        y = "ΔStability (removal − baseline)"
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 13),
        plot.title = element_text(size = 15, face = "bold")
      )
    
    ggsave(
      paste0(output_prefix, "_stability_delta_vs_baseline.svg"),
      plot = plot_stability_delta,
      width = 12, height = 8, units = "in"
    )
    
    # Bi-Stability delta plot
    if (has_bi_stability) {
      plot_bi_stability_delta <- ggplot(delta_data, aes(x = cell_type, y = bi_stability_delta, fill = removed_group)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.8) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
        scale_fill_manual(values = group_colors[names(group_colors) != baseline_name], 
                          name = "Removed Group") +
        labs(
          title = paste0("Bi-Stability Change vs Baseline (", baseline_name, ")"),
          x = "Cell Type",
          y = "ΔBi-Stability (removal − baseline)"
        ) +
        theme_minimal() +
        theme(
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 13, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 13),
          plot.title = element_text(size = 15, face = "bold")
        )
      
      ggsave(
        paste0(output_prefix, "_bi_stability_delta_vs_baseline.svg"),
        plot = plot_bi_stability_delta,
        width = 12, height = 8, units = "in"
      )
    } else {
      plot_bi_stability_delta <- NULL
    }
    
    # F1-Stability delta plot
    if (has_f1_stability) {
      plot_f1_stability_delta <- ggplot(delta_data, aes(x = cell_type, y = f1_stability_delta, fill = removed_group)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.8) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
        scale_fill_manual(values = group_colors[names(group_colors) != baseline_name], 
                          name = "Removed Group") +
        labs(
          title = paste0("F1-Stability Change vs Baseline (", baseline_name, ")"),
          x = "Cell Type",
          y = "ΔF1-Stability (removal − baseline)"
        ) +
        theme_minimal() +
        theme(
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 13, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 13),
          plot.title = element_text(size = 15, face = "bold")
        )
      
      ggsave(
        paste0(output_prefix, "_f1_stability_delta_vs_baseline.svg"),
        plot = plot_f1_stability_delta,
        width = 12, height = 8, units = "in"
      )
    } else {
      plot_f1_stability_delta <- NULL
    }
    
    # ===== Paired Barplots: Baseline (top) vs Removal (bottom) =====
    # For each individual removal, create a 2-panel plot per metric
    # showing baseline on top and removal on bottom, so the delta is visible
    
    paired_barplot_list <- list()
    paired_barplot_dir <- paste0(output_prefix, "_paired_barplots/")
    if (!dir.exists(paired_barplot_dir)) {
      dir.create(paired_barplot_dir, recursive = TRUE)
    }
    
    removal_names <- setdiff(names(run_results_list), baseline_name)
    
    # Define which metrics to plot
    metric_configs <- list(
      list(col = "stability",    sd_col = "stability_sd",    label = "Stability",    title_prefix = "Stability"),
      list(col = "bi_stability", sd_col = "bi_stability_sd", label = "Bi-Stability", title_prefix = "Bi-Stability"),
      list(col = "f1_stability", sd_col = "f1_stability_sd", label = "F1 Stability", title_prefix = "F1 Stability")
    )
    
    all_paired_plots <- list()  # for PDF
    
    for (removal_name in removal_names) {
      # Filter data for just baseline and this removal
      pair_data <- combined_data %>%
        filter(removed_group %in% c(baseline_name, removal_name))
      pair_data$removed_group <- factor(pair_data$removed_group, 
                                         levels = c(baseline_name, removal_name))
      
      # Clean name for filenames
      clean_name <- gsub("[^A-Za-z0-9_]", "_", removal_name)
      
      for (mc in metric_configs) {
        # Check if this metric has data
        metric_vals <- pair_data[[mc$col]]
        if (all(is.na(metric_vals) | metric_vals == 0)) next
        
        has_sd <- any(!is.na(pair_data[[mc$sd_col]]))
        
        # Top panel: baseline only
        baseline_panel_data <- pair_data %>% filter(removed_group == baseline_name)
        removal_panel_data <- pair_data %>% filter(removed_group == removal_name)
        
        common_theme <- theme_minimal() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
            axis.text.y = element_text(size = 12),
            axis.title = element_text(size = 13),
            plot.title = element_text(size = 14, face = "bold"),
            legend.position = "none"
          )
        
        p_top <- ggplot(baseline_panel_data, aes(x = cell_type, y = .data[[mc$col]])) +
          geom_bar(stat = "identity", fill = group_colors[baseline_name], alpha = 0.8) +
          { if (has_sd) 
              geom_errorbar(aes(ymin = .data[[mc$col]] - .data[[mc$sd_col]], 
                                ymax = .data[[mc$col]] + .data[[mc$sd_col]]),
                            width = 0.25, linewidth = 0.4)
          } +
          scale_y_continuous(limits = c(0, 1)) +
          labs(title = paste0(mc$title_prefix, ": ", baseline_name),
               x = NULL, y = mc$label) +
          common_theme
        
        p_bottom <- ggplot(removal_panel_data, aes(x = cell_type, y = .data[[mc$col]])) +
          geom_bar(stat = "identity", fill = group_colors[removal_name], alpha = 0.8) +
          { if (has_sd) 
              geom_errorbar(aes(ymin = .data[[mc$col]] - .data[[mc$sd_col]], 
                                ymax = .data[[mc$col]] + .data[[mc$sd_col]]),
                            width = 0.25, linewidth = 0.4)
          } +
          scale_y_continuous(limits = c(0, 1)) +
          labs(title = paste0(mc$title_prefix, ": ", removal_name),
               x = "Cell Type", y = mc$label) +
          common_theme
        
        # Stack vertically
        paired_plot <- p_top / p_bottom
        
        plot_filename <- paste0(paired_barplot_dir, clean_name, "_", mc$col, "_paired.svg")
        ggsave(plot_filename, plot = paired_plot, width = 10, height = 10, units = "in")
        
        key <- paste0(removal_name, "_", mc$col)
        paired_barplot_list[[key]] <- paired_plot
        all_paired_plots[[length(all_paired_plots) + 1]] <- paired_plot
      }
    }
    
    # Save combined PDF with all paired barplots
    if (length(all_paired_plots) > 0) {
      pdf_file <- paste0(output_prefix, "_paired_barplots_all.pdf")
      pdf(pdf_file, width = 4, height = 6)
      for (p in all_paired_plots) {
        print(p)
      }
      dev.off()
      message(paste0("Saved combined paired barplots PDF: ", pdf_file))
    }
    
  } else {
    plot_pss_delta <- NULL
    plot_stability_delta <- NULL
    plot_bi_stability_delta <- NULL
    plot_f1_stability_delta <- NULL
    paired_barplot_list <- list()
    delta_data <- NULL
  }
  
  # Return the plots and data
  return(list(
    plot_pss_comparison = plot_pss_comparison,
    plot_stability_comparison = plot_stability_comparison,
    plot_bi_stability_comparison = plot_bi_stability_comparison,
    plot_f1_stability_comparison = if (exists("plot_f1_stability_comparison")) plot_f1_stability_comparison else NULL,
    plot_pss_delta = plot_pss_delta,
    plot_stability_delta = plot_stability_delta,
    plot_bi_stability_delta = if (exists("plot_bi_stability_delta")) plot_bi_stability_delta else NULL,
    plot_f1_stability_delta = if (exists("plot_f1_stability_delta")) plot_f1_stability_delta else NULL,
    paired_barplots = if (exists("paired_barplot_list")) paired_barplot_list else list(),
    per_run_plots = per_run_plots,
    combined_data = combined_data,
    delta_data = delta_data,
    summary_stats = summary_stats,
    stability_ranking = stability_ranking
  ))
}

# ============================================================
# Example usage with Atlas subsampling results:
# ============================================================
# source("subsampling_pipeline.R")
# source("plot_subsampling_comparison.R")
#
# # Assuming you already ran the subsampling pipeline (example 6 style):
# atlas_plots <- plot_subsampling_comparison(
#   list(
#     "Baseline"           = cell_atlas_flow_anchors,
#     "No CM+CM_DIV+PROX_1" = no_cm_cmdiv_prox1_flow_anchors,
#     "No PROX_2"          = no_prox2_flow$flow,
#     "No DIST_CD"         = no_distcd_flow$flow,
#     "No ENDO"            = no_endo_flow$flow,
#     "No PODO"            = no_podo_flow$flow
#   ),
#   output_prefix = "six2gfp/subsampling/atlas_removal_comparison"
# )
