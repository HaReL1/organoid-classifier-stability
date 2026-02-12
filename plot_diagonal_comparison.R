library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)

#' Plot Diagonal Comparison Across Multiple Runs
#' 
#' This function creates scatter plots comparing the diagonal elements
#' (PSS self-similarity and stability) across multiple experimental runs.
#' 
#' @param run_results_list A named list of run results, where each element is the
#'   output from analyze_noise_impact_on_prediction(). Names should be the dataset
#'   names (e.g., "Uchimura", "Freedman", "Cell Atlas")
#' @param output_prefix Path prefix for saving output plots
#' @param plot_title_pss Title for the PSS diagonal plot
#' @param plot_title_stability Title for the stability diagonal plot
#' 
#' @return A list containing both ggplot objects
#' 
#' @examples
#' # After running your analysis:
#' results_list <- list(
#'   "Uchimura" = Uchimura_full_flow,
#'   "Freedman" = freedman_flow,
#'   "Cell Atlas" = cell_atlas_flow
#' )
#' diagonal_plots <- plot_diagonal_comparison(
#'   results_list, 
#'   output_prefix = "six2gfp/diagonal_comparison"
#' )
plot_diagonal_comparison <- function(
    run_results_list,
    output_prefix = "diagonal_comparison",
    plot_title_pss = "PSS Self-Similarity Across Datasets",
    plot_title_stability = "Stability Across Datasets"
) {
  
  # Validate input
  if (length(run_results_list) < 1) {
    stop("run_results_list must contain at least one run result")
  }

  # 1. Identify common cell types across ALL runs
  common_cell_types <- NULL
  
  for (run_name in names(run_results_list)) {
      run_result <- run_results_list[[run_name]]
      
      if (is.null(run_result$stability_pred)) {
          warning(paste("Run", run_name, "missing stability_pred, skipping for cell type intersection"))
          next
      }
      
      current_types <- rownames(run_result$stability_pred)
      
      if (is.null(common_cell_types)) {
          common_cell_types <- current_types
      } else {
          common_cell_types <- intersect(common_cell_types, current_types)
      }
  }
  
  if (length(common_cell_types) == 0) {
      stop("No common cell types found across the provided runs.")
  }
  print(paste("Common cell types found:", paste(common_cell_types, collapse=", ")))

  
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
    
    # DEBUG: Print available fields to understand why counts are missing
    print(paste("Analyzing", run_name, "- Available fields:"))
    print(names(run_result))
    
    # Try to extract cell counts from multiple potential sources
    cell_counts <- rep(NA, length(cell_types))
    
    # Method 1: Try noised_prediction_matrix
    if (!is.null(run_result$noised_prediction_matrix)) {
      prediction_vector <- run_result$noised_prediction_matrix[, 1]
      count_table <- table(prediction_vector)
      for (i in seq_along(cell_types)) {
        if (cell_types[i] %in% names(count_table)) {
          cell_counts[i] <- count_table[cell_types[i]]
        }
      }
      message(paste0("  ", run_name, ": SUCCEEDED - Extracted counts from noised_prediction_matrix"))
    } 
    # Method 2: Check if there's a test object with metadata
    else if (!is.null(run_result$test) && inherits(run_result$test, "Seurat")) {
        # Extract from Seurat object metadata
        predicted_types <- run_result$test@meta.data$predicted.type
        if (!is.null(predicted_types)) {
          count_table <- table(predicted_types)
          for (i in seq_along(cell_types)) {
            if (cell_types[i] %in% names(count_table)) {
              cell_counts[i] <- count_table[cell_types[i]]
            }
          }
          message(paste0("  ", run_name, ": SUCCEEDED - Extracted counts from Seurat test object"))
        }
    } else {
      message(paste0("  ", run_name, ": FAILED - Could not find cell count data in any known field"))
    }
    
    # Create data frame for this run
    run_data <- data.frame(
      cell_type = cell_types,
      pss = pss_values,
      stability = stability_values,
      dataset = run_name,
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
      filter(cell_type %in% common_cell_types)

  # Enforce specific order for cell types
  desired_order <- c("UM","CM","CM_DIV","PODO","PROX_1","PROX_2","LOH","DIST_CD","ENDO","MACROPHAG")
  # Retrieve intersection of desired order and available (common) types to avoid NA levels
  final_levels <- intersect(desired_order, common_cell_types)
  
  if (length(final_levels) < length(common_cell_types)) {
      warning("Some common cell types are missing from the manual ordering list and will be excluded from the plot.")
      combined_data <- combined_data %>% filter(cell_type %in% final_levels)
  }

  combined_data$cell_type <- factor(combined_data$cell_type, levels = final_levels)
  
  # Define colors and shapes for different datasets
  num_datasets <- length(unique(combined_data$dataset))
  dataset_colors <- scales::hue_pal()(num_datasets)
  names(dataset_colors) <- unique(combined_data$dataset)
  
  # Define different shape isn't strictly needed for bar plots but keeping logic if needed later
  display_shapes <- c(16, 17, 15, 18, 3, 4) 
  
  # ===== Plot 2a: PSS values comparison - NOT normalized =====
  
  plot_pss_comparison <- ggplot(combined_data, aes(x = cell_type, y = pss, fill = dataset)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.8) +
    # Removed max dashed line
    geom_text(aes(label = ifelse(!is.na(n_cells), paste0("n=", n_cells), "")), 
              position = position_dodge(width = 0.8), 
              vjust = -0.5, size = 2.5, angle = 0) +
    scale_fill_manual(values = dataset_colors, name = "Dataset") +
    labs(
      title = "PSS Self-Similarity by Cell Type Across Datasets",
      x = "Cell Type",
      y = "Prediction Specificity Score (PSS)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 13, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 13),
      plot.title = element_text(size = 15, face = "bold")
    )
  
  # Save plot 2a
  ggsave(
    paste0(output_prefix, "_pss_by_celltype_comparison.svg"),
    plot = plot_pss_comparison,
    width = 12, height = 8, units = "in"
  )
  
  
  # ===== Plot 3a: Stability values comparison - NOT normalized =====
  
  plot_stability_comparison <- ggplot(combined_data, aes(x = cell_type, y = stability, fill = dataset)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.8) +
    # Removed max dashed line
    geom_text(aes(label = ifelse(!is.na(n_cells), paste0("n=", n_cells), "")), 
              position = position_dodge(width = 0.8), 
              vjust = -0.5, size = 2.5, angle = 0) +
    scale_fill_manual(values = dataset_colors, name = "Dataset") +
    labs(
      title = "Stability by Cell Type Across Datasets",
      x = "Cell Type",
      y = "Stability (Fraction of cells remaining same type)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 13, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 13),
      plot.title = element_text(size = 15, face = "bold")
    )
  
  # Save plot 3a
  ggsave(
    paste0(output_prefix, "_stability_by_celltype_comparison.svg"),
    plot = plot_stability_comparison,
    width = 12, height = 8, units = "in"
  )
  
  
  # Print summary statistics
  cat("\n=== Summary Statistics ===\n")
  summary_stats <- combined_data %>%
    group_by(dataset) %>%
    summarise(
      n_cell_types = n(),
      mean_pss = mean(pss, na.rm = TRUE),
      sd_pss = sd(pss, na.rm = TRUE),
      mean_stability = mean(stability, na.rm = TRUE),
      sd_stability = sd(stability, na.rm = TRUE),
      .groups = "drop"
    )
  print(summary_stats)
  
  # Return the plots and data (Cleaned up list)
  return(list(
    plot_pss_comparison = plot_pss_comparison,
    plot_stability_comparison = plot_stability_comparison,
    combined_data = combined_data,
    summary_stats = summary_stats
  ))
}

a=plot_diagonal_comparison(list("Uchimura" = Uchimura_full_flow, 
                                "Freedman" = freedman_flow,
                                "Cell Atlas" = cell_atlas_flow,
                                "Takasato" = Takasato_full_flow),
                           output_prefix="six2gfp/26.1.26/refactored_diagonal_comparison_"
)
