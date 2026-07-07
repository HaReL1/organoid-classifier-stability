library(Seurat)
library(ggplot2)
library(dplyr)

#' Create Marker DotPlots for multiple datasets
#'
#' @param run_results_list A named list of run results (e.g., from .diag_list)
#' @param output_prefix Path prefix for saving output plots
#' @param n_markers Number of top markers per cell type to display (if genes_to_plot is NULL)
#' @param cell_type_column The column in meta.data containing cell types
#' @param genes_to_plot Optional specific list of genes to plot for all datasets instead of finding markers dynamically
#' @param cache_paths Optional named list mapping run names to base paths to load 'run_without_noise.rds'
#'
#' @return A list of ggplot DotPlot objects
plot_marker_dotplots <- function(run_results_list, output_prefix = "six2gfp/new/dotplots/", 
                                 n_markers = 3, cell_type_column = "predicted.type",
                                 genes_to_plot = NULL, cache_paths = NULL) {
  
  if (!dir.exists(output_prefix)) {
    dir.create(output_prefix, recursive = TRUE, showWarnings = FALSE)
  }
  
  plot_list <- list()
  
  # Standard cell type order
  desired_order <- c("UM","CM","CM_DIV","PODO","PROX_1","PROX_2","LOH","DIST_CD","ENDO","MACROPHAG")
  
  for (run_name in names(run_results_list)) {
    print(paste("Generating DotPlot for dataset:", run_name))
    run_result <- run_results_list[[run_name]]
    
    # We use the 'test' Seurat object, which contains the predictions
    seurat_obj <- run_result$test
    if (is.null(seurat_obj) && !is.null(cache_paths[[run_name]])) {
      cache_file <- file.path(cache_paths[[run_name]], "run_without_noise.rds")
      if (file.exists(cache_file)) {
        print(paste("Loading test Seurat from sub-cache:", cache_file))
        run_without_noise <- readRDS(cache_file)
        seurat_obj <- run_without_noise$test
      }
    }
    
    if (is.null(seurat_obj)) {
      warning(paste("Run", run_name, "missing 'test' Seurat object and no valid cache found, skipping"))
      next
    }
    
    if (!cell_type_column %in% colnames(seurat_obj@meta.data)) {
      # Fallback to 'type' if 'predicted.type' is missing
      if (cell_type_column == "predicted.type" && "type" %in% colnames(seurat_obj@meta.data)) {
        seurat_obj[[cell_type_column]] <- seurat_obj[["type"]]
      } else {
        warning(paste("Column", cell_type_column, "not found in metadata for", run_name, "- skipping"))
        next
      }
    }
    
    # Enforce order on the groups
    available_levels <- intersect(desired_order, unique(seurat_obj[[cell_type_column]][[1]]))
    if (length(available_levels) < length(unique(seurat_obj[[cell_type_column]][[1]]))) {
      other_levels <- setdiff(unique(seurat_obj[[cell_type_column]][[1]]), desired_order)
      available_levels <- c(available_levels, other_levels)
    }
    seurat_obj[[cell_type_column]] <- factor(seurat_obj[[cell_type_column]][[1]], levels = available_levels)
    
    Idents(seurat_obj) <- cell_type_column
    
    # Switch to RNA assay for expression plotting if available
    if ("RNA" %in% Assays(seurat_obj)) {
      DefaultAssay(seurat_obj) <- "RNA"
    }
    
    plot_genes <- genes_to_plot
    
    if (is.null(plot_genes)) {
      print(paste("Finding top", n_markers, "markers for", run_name, "... (this may take a moment)"))
      # Find markers
      markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
      
      if (nrow(markers) == 0) {
        warning(paste("No markers found for dataset", run_name))
        next
      }
      
      # Ensure genes are ordered by cell type (cluster)
      markers$cluster <- factor(markers$cluster, levels = available_levels)
      top_markers <- markers %>%
        group_by(cluster) %>%
        slice_max(n = n_markers, order_by = avg_log2FC, with_ties = FALSE) %>%
        arrange(cluster)
        
      plot_genes <- unique(top_markers$gene)
    }
    
    if (length(plot_genes) == 0) {
      warning(paste("No genes to plot for dataset", run_name))
      next
    }
    
    # Check which genes are actually in the dataset
    plot_genes <- intersect(plot_genes, rownames(seurat_obj))
    if (length(plot_genes) == 0) {
      warning(paste("None of the specified genes were found in the dataset", run_name))
      next
    }
    
    # Generate DotPlot
    # By default, DotPlot puts features on X and idents on Y. We use coord_flip() 
    # to match the standard look (groups on X, genes on Y)
    p <- DotPlot(seurat_obj, features = plot_genes, group.by = cell_type_column) + 
      coord_flip() +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 10),
        axis.title = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        panel.grid.major = element_line(color = "grey90"),
        panel.border = element_rect(color = "black", fill = NA)
      ) +
      labs(title = paste("Marker Genes -", run_name)) +
      scale_color_gradient(low = "lightgrey", high = "red")
      
    # Dynamic sizing based on number of genes and cell types
    plot_width <- max(6, length(available_levels) * 0.4 + 2)
    plot_height <- max(5, length(plot_genes) * 0.15 + 2)
    
    out_file <- paste0(output_prefix, gsub(" ", "_", run_name), "_dotplot.svg")
    ggsave(out_file, plot = p, width = plot_width, height = plot_height, units = "in")
    
    plot_list[[run_name]] <- p
    print(paste("Saved DotPlot to", out_file))
  }
  
  return(plot_list)
}

#' Plot Hardcoded Marker DotPlots
#'
#' @param run_results_list A named list of run results (e.g., from .diag_list)
#' @param output_prefix Path prefix for saving output plots
#' @param cell_type_column The column in meta.data containing cell types
#' @param cache_paths Optional named list mapping run names to base paths to load 'run_without_noise.rds'
#'
#' @return A list of ggplot DotPlot objects
plot_hardcoded_marker_dotplots <- function(run_results_list, output_prefix = "six2gfp/new/dotplots_hardcoded/", 
                                           cell_type_column = "predicted.type",
                                           cache_paths = NULL) {
  
  hardcoded_genes <- c(
    "KDR", "SOX17", "COL1A2", "COL1A1", "PDGFRB", "DCN", "FOXD1", 
    "DLL1", "IRX5", "IRX3", "CENPF", "NOTCH2", "WNT4", "LYPD1", "MKI67","TOP2A",
    "MEG3", "DAPL1", "SIX2", "CITED1","EYA1","SALL2","CRYM", "LAMP5", "TRPM3", "HNF1B", 
    "LHX1", "JAG1", "CCND1", "PCNA", "TOP2A", "IRX1", "DKK1", 
    "UPK2", "TACSTD2", "AQP2", "RET", "GATA3", "CALB1", "SLC12A3", 
    "TMEM207", "PROX1", "SLC12A1", "UMOD", "CLDN10", "MAL", 
    "POU3F3", "AQP1", "SLC22A6", "SLC6A19", "SLC34A1", "LRP2", 
    "CUBN", "HNF4A", "CDH11","CDH6","CDH1", "SLC27A2", "CLEC18B", "CLDN1", 
    "CAV1", "OLFM3", "MAFB", "NPHS1", "NPHS2", "NTNG1", "DDN","KCNJ1","MUC1"
  )
  
  # Reverse the list so the first gene (KDR) appears at the top of the y-axis when coord_flip is used
  genes_to_plot <- rev(hardcoded_genes)
  
  plot_marker_dotplots(run_results_list, output_prefix, 
                       cell_type_column = cell_type_column,
                       genes_to_plot = genes_to_plot,
                       cache_paths = cache_paths)
}
