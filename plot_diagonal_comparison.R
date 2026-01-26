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
    
    # Create data frame for this run
    run_data <- data.frame(
      cell_type = cell_types,
      pss = pss_values,
      stability = stability_values,
      dataset = run_name,
      stringsAsFactors = FALSE
    )
    
    combined_data <- rbind(combined_data, run_data)
  }
  
  # Check if we have data
  if (nrow(combined_data) == 0) {
    stop("No valid data extracted from any run")
  }
  
  # Define colors and shapes for different datasets
  num_datasets <- length(unique(combined_data$dataset))
  dataset_colors <- scales::hue_pal()(num_datasets)
  names(dataset_colors) <- unique(combined_data$dataset)
  
  # Define different shapes for datasets (max 6 different shapes)
  dataset_shapes <- c(16, 17, 15, 18, 3, 4)  # circle, triangle, square, diamond, plus, X
  dataset_shapes <- dataset_shapes[1:min(num_datasets, 6)]
  names(dataset_shapes) <- unique(combined_data$dataset)
  
  # ===== Plot 1: PSS vs Stability (same as the individual plots but with multiple datasets) =====
  plot_pss_stability <- ggplot(combined_data, aes(x = pss, y = stability, 
                                                    color = dataset, shape = dataset)) +
    geom_point(size = 4, alpha = 0.7) +
    geom_text_repel(aes(label = cell_type), size = 3, 
                    box.padding = 0.5, max.overlaps = 20) +
    scale_color_manual(values = dataset_colors, name = "Dataset") +
    scale_shape_manual(values = dataset_shapes, name = "Dataset") +
    labs(
      title = plot_title_pss,
      x = "Prediction Specificity Score (PSS)",
      y = "Stability (Fraction of cells remaining same type)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 13, face = "bold"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 13),
      plot.title = element_text(size = 15, face = "bold")
    ) +
    guides(
      color = guide_legend(override.aes = list(size = 4)),
      shape = guide_legend(override.aes = list(size = 4))
    )
  
  # Save plot 1
  ggsave(
    paste0(output_prefix, "_pss_vs_stability_comparison.svg"),
    plot = plot_pss_stability,
    width = 12, height = 8, units = "in"
  )
  
  # ===== Plot 2: Just PSS values comparison across datasets =====
  # This shows how PSS varies for each cell type across datasets
  plot_pss_comparison <- ggplot(combined_data, aes(x = cell_type, y = pss, 
                                                     color = dataset, shape = dataset, group = dataset)) +
    geom_point(size = 4, alpha = 0.7, position = position_dodge(width = 0.5)) +
    geom_line(aes(group = dataset), alpha = 0.3, position = position_dodge(width = 0.5)) +
    scale_color_manual(values = dataset_colors, name = "Dataset") +
    scale_shape_manual(values = dataset_shapes, name = "Dataset") +
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
  
  # Save plot 2
  ggsave(
    paste0(output_prefix, "_pss_by_celltype_comparison.svg"),
    plot = plot_pss_comparison,
    width = 12, height = 8, units = "in"
  )
  
  # ===== Plot 3: Just Stability values comparison across datasets =====
  plot_stability_comparison <- ggplot(combined_data, aes(x = cell_type, y = stability, 
                                                           color = dataset, shape = dataset, group = dataset)) +
    geom_point(size = 4, alpha = 0.7, position = position_dodge(width = 0.5)) +
    geom_line(aes(group = dataset), alpha = 0.3, position = position_dodge(width = 0.5)) +
    scale_color_manual(values = dataset_colors, name = "Dataset") +
    scale_shape_manual(values = dataset_shapes, name = "Dataset") +
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
  
  # Save plot 3
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
  
  # Return the plots and data
  return(list(
    plot_pss_stability = plot_pss_stability,
    plot_pss_comparison = plot_pss_comparison,
    plot_stability_comparison = plot_stability_comparison,
    combined_data = combined_data,
    summary_stats = summary_stats
  ))
}


#' Alternative version: Plot with facets instead of overlaying
#' 
#' This creates separate panels for each dataset which can make it easier
#' to see patterns within each dataset
plot_diagonal_comparison_faceted <- function(
    run_results_list,
    output_prefix = "diagonal_comparison_faceted",
    plot_title = "PSS vs Stability Comparison Across Datasets"
) {
  
  # Extract diagonal data (same as above)
  combined_data <- data.frame()
  
  for (run_name in names(run_results_list)) {
    run_result <- run_results_list[[run_name]]
    
    if (is.null(run_result$stability_pred)) {
      warning(paste("Run", run_name, "missing stability_pred, skipping"))
      next
    }
    
    cell_types <- rownames(run_result$stability_pred)
    pss_values <- run_result$stability_pred$PSS
    stability_values <- run_result$stability_pred$Stability
    
    run_data <- data.frame(
      cell_type = cell_types,
      pss = pss_values,
      stability = stability_values,
      dataset = run_name,
      stringsAsFactors = FALSE
    )
    
    combined_data <- rbind(combined_data, run_data)
  }
  
  if (nrow(combined_data) == 0) {
    stop("No valid data extracted from any run")
  }
  
  # Create faceted plot
  faceted_plot <- ggplot(combined_data, aes(x = pss, y = stability)) +
    geom_point(size = 3, alpha = 0.7, color = "red") +
    geom_text_repel(aes(label = cell_type), size = 2.5, 
                    box.padding = 0.3, max.overlaps = 15) +
    facet_wrap(~ dataset, ncol = 2) +
    labs(
      title = plot_title,
      x = "Prediction Specificity Score (PSS)",
      y = "Stability (Fraction of cells remaining same type)"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 14, face = "bold")
    )
  
  # Save faceted plot
  ggsave(
    paste0(output_prefix, "_faceted.svg"),
    plot = faceted_plot,
    width = 12, height = 8, units = "in"
  )
  
  return(list(
    faceted_plot = faceted_plot,
    combined_data = combined_data
  ))
}


a=plot_diagonal_comparison(list("Uchimura" = Uchimura_full_flow, 
                                "Freedman" = freedman_flow,
                                "Cell Atlas" = cell_atlas_flow),
                           output_prefix="six2gfp/26.1.26/diagonal_comparison"
)
