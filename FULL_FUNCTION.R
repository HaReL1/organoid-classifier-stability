# library(randomForest)
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)
library(RColorBrewer)
library(RANN)
library(gridExtra)
library(reshape2)
library(rlang)
library(ggrepel)
library(tidyr)
library(svglite)
library(Matrix)

analyze_noise_impact_on_prediction <- function(
    seurat_obj,
    noised_number = 3,
    train_Six2GFP,
    test_Six2GFP,
    ref_cell_type_column = "type",
    dims = 1:30,
    train_title = "SIX2GFP",
    test_title_prefix = "Uchimura",
    n_neighbors = 8,
    skip_neighbors = FALSE,
    output_prefix_base = "six2gfp/Uchimura_noised/",
    prediction_column_name = "predicted.type",
    colors_feature_plot_noise = c('grey', '#f03b20'),
    myColors_cell_types = NULL,
    return_all_suerats = FALSE,
    min_cell_count_for_type = 3,
    use_cache = TRUE,
    return_anchors = FALSE
) {
  # 0. Create output directory if it doesn't exist ####
  if (!dir.exists(output_prefix_base)) {
    dir.create(output_prefix_base, recursive = TRUE)
  }
  
  # 1. Add noise function ####
  add_poisson_noise <- function(seurat_obj) {
    counts <- GetAssayData(seurat_obj, layer = "counts")
    # A dgCMatrix has a slot 'x' that contains all the non-zero values.
    noisy_x_values <- rpois(n = length(counts@x), lambda = counts@x)
    # Reconstruct a sparse matrix directly, preserving dimensions and names
    noisy_counts_sparse <- Matrix::sparseMatrix(
      i = counts@i,         # Row indices of non-zero elements (from original)
      p = counts@p,         # Pointers to column starts (from original)
      x = noisy_x_values,   # The NEW noisy data
      dims = dim(counts),
      dimnames = dimnames(counts),
      index1 = FALSE        # Important: Seurat's slots are 0-indexed
    )
    
    # Create the new Seurat object
    noisy_seurat <- CreateSeuratObject(counts = noisy_counts_sparse, project = "noised")
    if ("orig.ident" %in% colnames(seurat_obj@meta.data)) {
      noisy_seurat@meta.data <- seurat_obj@meta.data[, "orig.ident", drop = FALSE]
    } else if ("cell_type" %in% colnames(seurat_obj@meta.data)) {
      noisy_seurat@meta.data <- seurat_obj@meta.data[, "cell_type", drop = FALSE]
    } else {
      print("Neither 'orig.ident' nor 'cell_type' column found in meta.data")
    }
    return(noisy_seurat)
  }
  
  # 2. Embedded process_seurat_data_iter function ####
  process_seurat_data_iter <- function(labeled_train_data,
                                       test_data,
                                       labeled_test_data = NULL,
                                       dims = 1:30,
                                       plot_on_this_UMAP = NULL,
                                       train_title = "SIX2GFP",
                                       ref_cell_type_column = "type",
                                       test_title = "Uchimura",
                                       n_neighbors = 8,
                                       skip_neighbors = FALSE,
                                       colors_for_feature = NULL,
                                       myColors = NULL,
                                       output_prefix = "output",
                                       min_cell_count_for_type = 3,
                                       use_cache = TRUE,
                                       return_anchors = FALSE)  {
    
    t1 <- Sys.time()
    
    # --- Caching Setup ---
    cache_dir <- paste0(output_prefix, "_cache/")
    if(use_cache & !dir.exists(cache_dir)) {
      dir.create(cache_dir, recursive = TRUE)
    }
    train_cache_file <- paste0(cache_dir, "train_labeled_seurat_pca_umap_",test_title,".rds")
    test_query_cache_file <- paste0(cache_dir, "test_query_initial_",test_title,".rds")
    labeled_test_query_cache_file <- paste0(cache_dir, "labeled_test_query_initial_",test_title,".rds")
    
    # Process training data. Input is Seurat object
    if (use_cache && file.exists(train_cache_file)) {
      print("Loading cached training data...")
      train_labeled_seurat <- readRDS(train_cache_file)
    } else {
      print("Processing training data...")
      train_labeled_seurat <- NormalizeData(labeled_train_data)
      train_labeled_seurat <- FindVariableFeatures(train_labeled_seurat, selection.method = "vst", nfeatures = 2000)
      train_labeled_seurat <- ScaleData(train_labeled_seurat)
      train_labeled_seurat <- RunPCA(train_labeled_seurat)
      train_labeled_seurat <- FindNeighbors(train_labeled_seurat, dims = dims)
      train_labeled_seurat <- FindClusters(train_labeled_seurat)
      train_labeled_seurat <- RunUMAP(train_labeled_seurat, dims = dims, return.model = TRUE)
      if (use_cache) {
        print("Saving cached training data...")
        saveRDS(train_labeled_seurat, train_cache_file)
      }
    }
    
    if(!is.null(labeled_test_data)){
      # Process labeled test data
      if (use_cache && file.exists(labeled_test_query_cache_file)) {
        print("Loading cached labeled test query...")
        labeled_test_query <- readRDS(labeled_test_query_cache_file)
      } else {
        print("Processing labeled test data...")
        labeled_test_data <- NormalizeData(labeled_test_data)
        
        labeled_test_anchors <- FindTransferAnchors(reference = train_labeled_seurat, query = labeled_test_data, dims = dims,
                                                    reference.reduction = "pca")
        labeled_test_query <- MapQuery(anchorset = labeled_test_anchors, reference = train_labeled_seurat, query = labeled_test_data,
                                       refdata = setNames(list(ref_cell_type_column),ref_cell_type_column),
                                       reference.reduction = "pca", reduction.model = "umap")
        if (use_cache) {
          print("Saving cached labeled test query...")
          saveRDS(labeled_test_query, labeled_test_query_cache_file)
        }
      }
    }
    else
    {
      labeled_test_query = NULL
    }

    test_anchors_result <- NULL
    if (use_cache && file.exists(test_query_cache_file)) {
      print(paste0("Loading cached test query (initial)... ", test_query_cache_file))
      test_query <- readRDS(test_query_cache_file)
    } else {
      print("Processing test data (initial)...")
      test_data <- NormalizeData(test_data)
      test_anchors <- FindTransferAnchors(reference = train_labeled_seurat, query = test_data, dims = dims,
                                          reference.reduction = "pca")
      if (return_anchors) {
        test_anchors_result <- test_anchors
      }
      test_query <- MapQuery(anchorset = test_anchors, reference = train_labeled_seurat, query = test_data,
                             refdata = setNames(list(ref_cell_type_column),ref_cell_type_column), reference.reduction = "pca",
                             reduction.model = "umap")
      if (use_cache) {
        print("Saving cached test query (initial)...")
        saveRDS(test_query, test_query_cache_file)
      }
    }
    
    #### decision boundaries ####
    visualize_decision_boundaries <- function(plot_ref, # x-y are from here
                                              query, # prediction are from here
                                              resolution = 100,
                                              ref_cell_type = ref_cell_type_column) {
    # Extract the coordinates from the dimensional reduction
    print("Running decision boundaries plot")
    umap_type <- NULL
    if ("ref.umap" %in% Reductions(plot_ref)) {
      umap_type <- "ref.umap"
    } else if ("umap" %in% Reductions(plot_ref)) {
      umap_type <- "umap"
    } else { stop("Neither 'umap' nor 'ref.umap' reductions found in the Seurat object") }
    coords <- Embeddings(plot_ref, reduction = umap_type)
    train_data <- data.frame(
      x = coords[, 1],
      y = coords[, 2],
      class = as.factor(query@meta.data[[paste0("predicted.",ref_cell_type)]])
    )
    # Create a grid for prediction #
    x_range <- range(coords[, 1])
    y_range <- range(coords[, 2])

    # Add some padding #
    x_padding <- 0.05 * (x_range[2] - x_range[1])
    y_padding <- 0.05 * (y_range[2] - y_range[1])
    x_seq <- seq(from = x_range[1] - x_padding, to = x_range[2] + x_padding, length.out = 100)
    y_seq <- seq(from = y_range[1] - y_padding, to = y_range[2] + y_padding, length.out = 100)
    
    grid <- expand.grid(x = x_seq, y = y_seq)
    
    # Fit a Random Forest classifier on the 2D embeddings
    rf_model <- randomForest(class ~ x + y, data = train_data, ntree = 100)
    
    # Predict on the grid
    grid$predicted <- predict(rf_model, newdata = grid)
    
    # Create the plot
    decision_boundaries_plot <- ggplot() +
      geom_raster(data = grid, aes(x = x, y = y, fill = predicted), alpha = 0.3) +
      geom_point(data = train_data, aes(x = x, y = y, color = class), size = 0.5) +
      scale_fill_discrete(name = "Predicted Class") +
      scale_color_discrete(name = "Original Label") +
      labs(x = paste0("umap", "_", dims[1]), 
           y = paste0("umap", "_", dims[2]),
           title = "Decision Boundaries") +
      theme_minimal()
    return(decision_boundaries_plot)
    }
    #####
    
    # Consistent Order in Heatmaps - Get cell type order from training data
    cell_type_order <- unique(train_labeled_seurat[[ref_cell_type_column]])[[ref_cell_type_column]]
    cell_type_order = factor(cell_type_order, levels = c("UM","CM","CM_DIV","PODO","PROX_1","PROX_2","LOH","DIST_CD","ENDO","MACROPHAG"))
    cell_type_order = sort(cell_type_order)
    if (is.null(myColors)) {
      myColors <- brewer.pal(length(cell_type_order), "Paired") # use length of cell_type_order
      names(myColors) <- cell_type_order # use names from cell_type_order
    }
    
    if (is.null(colors_for_feature)) {
      colors_for_feature <- c('grey', '#f03b20')
    }
    
    calculate_plot_limits <- function(train_data, labeled_test_data, test_data, padding = 0.05) {
      train_umap <- train_data@reductions$umap@cell.embeddings
      test_umap <- test_data@reductions$ref.umap@cell.embeddings
      
      if (!is.null(labeled_test_data)) {
        labeled_test_umap <- labeled_test_data@reductions$ref.umap@cell.embeddings
        all_umap <- rbind(train_umap, labeled_test_umap, test_umap)
      } else {
        all_umap <- rbind(train_umap, test_umap)
      }
      
      x_min <- min(all_umap[,1])
      x_max <- max(all_umap[,1])
      y_min <- min(all_umap[,2])
      y_max <- max(all_umap[,2])
      
      x_range <- x_max - x_min
      y_range <- y_max - y_min
      
      xlim <- c(x_min - padding * x_range, x_max + padding * x_range)
      ylim <- c(y_min - padding * y_range, y_max + padding * y_range)
      
      return(list(xlim = xlim, ylim = ylim))
    }
    
    # Calculate entropy function
    calculate_entropy <- function(query) {
      log_pred <- log2(query@assays[[paste0("prediction.score.",ref_cell_type_column)]]@data)
      log_pred[is.infinite(log_pred)] <- 0
      entropy_by_cell <- -colSums(query@assays[[paste0("prediction.score.",ref_cell_type_column)]]@data * log_pred)
      query <- AddMetaData(query, entropy_by_cell, col.name = "entropy")
      return(query)
    }
    
    # Calculate nearest neighbors
    calculate_neighbors <- function(train_data, query) {
      test_data <- query@reductions[["ref.pca"]]@cell.embeddings
      all_types = cell_type_order # Use ordered cell types from training data
      
      all_neighbors_dist <- data.frame(matrix(NA, nrow = nrow(test_data), ncol = length(all_types)))
      
      for (i in 1:nrow(test_data)) {
        for(t in all_types)
        {
          train_same_as_t <- (train_labeled_seurat[[ref_cell_type_column]] == t)
          training_relevant_data <- train_data[train_same_as_t, ]
          knn_closest <- nn2(training_relevant_data, as.data.frame(t(test_data[i, ])), k = n_neighbors)
          all_neighbors_dist[i,which(t==all_types)] <- mean(knn_closest$nn.dist[1, ])
          if (i %% 50 == 0) {
            print(paste0("Row:", i, " Time:", Sys.time()))
          }
        }
      }
      colnames(all_neighbors_dist)=all_types
      return(all_neighbors_dist)
    }
    
    # Create confusion matrix
    calculate_pairwise_prediction_overlap <- function(query, title, cell_type_order) {
      probs_tub <- query@assays[[paste0("prediction.score.",ref_cell_type_column)]]@data
      num_classes <- nrow(probs_tub)
      pairwise_confusion <- matrix(0, nrow = num_classes, ncol = num_classes)
      row.names(pairwise_confusion) <- gsub("-", "_", row.names(probs_tub))
      colnames(pairwise_confusion) <- row.names(pairwise_confusion)
      
      # 1. Consistent Order in Heatmaps - Reorder matrix rows and cols
      pairwise_confusion <- pairwise_confusion[cell_type_order, cell_type_order]
      
      for (i in 1:num_classes) {
        for (j in 1:num_classes) {
          if (i != j) {
            pairwise_confusion[i, j] <- mean(probs_tub[i, ] * probs_tub[j, ])
          }
        }
      }
      
      pairwise_confusion_df <- reshape2::melt(pairwise_confusion)
      p <- ggplot(data = pairwise_confusion_df, aes(Var1, Var2, fill = value)) +
        geom_tile() +
        geom_text(aes(label = round(value,5))) +
        scale_fill_gradient(low = "white", high = "red") +
        labs(x = "Class", y = "Class", title = paste("Pairwise Confusion Matrix -", title)) + # 2. Add Title
        theme_minimal()
      return(p)
    }
    
    # Calculate plot limits dynamically
    plot_limits <- calculate_plot_limits(train_labeled_seurat, labeled_test_query, test_query)
    xlim <- plot_limits$xlim
    ylim <- plot_limits$ylim
    
    # create plots
    train_plot <- DimPlot(train_labeled_seurat, reduction = "umap", group.by = ref_cell_type_column, label = TRUE, label.size = 3,
                          repel = TRUE) + ggtitle(paste0(train_title,"\n(train) annotations")) + # 2. Add Title
      scale_color_manual(values = myColors) +
      theme(legend.position = "left", legend.text = element_text(size=8)) + xlim(xlim) + ylim(ylim)
    
    test_plot <- DimPlot(test_query, reduction = "ref.umap", group.by = paste0("predicted.",ref_cell_type_column), label = TRUE,
                         label.size = 3, repel = TRUE) + ggtitle(paste(test_title, "\nTransferred Labels")) + # 2. Add Title
      scale_color_manual(values = myColors) + NoLegend() + xlim(xlim) + ylim(ylim)
    
    test_query <- calculate_entropy(test_query)
    quantile_80 = quantile(test_query$entropy, probs = 0.8)
    
    test_entropy_plot <- FeaturePlot(test_query, reduction = "ref.umap", features = "entropy",
                                     cols = colors_for_feature, keep.scale = "all") +
      ggtitle(paste(test_title, "\nEntropy")) + xlim(xlim) + ylim(ylim)
    
    train_data <- train_labeled_seurat@reductions[["pca"]]@cell.embeddings[, dims]
    test_query_neighbors=''
    if (!skip_neighbors){
      test_query_neighbors <- calculate_neighbors(train_data, test_query)
    }
    
    test_confusion <- calculate_pairwise_prediction_overlap(test_query, test_title, cell_type_order) # pass cell_type_order
    
    # decision_boundaries_mapped <- visualize_decision_boundaries(test_query,test_query)
    
    if (!is.null(labeled_test_query)) {
      labeled_test_query <- calculate_entropy(labeled_test_query)
      labeled_test_confusion <- calculate_pairwise_prediction_overlap(labeled_test_query, "Labeled Test", cell_type_order) # pass cell_type_order
      
      # plot labeled_test related plots
      labeled_test_plot <- DimPlot(labeled_test_query, reduction = "ref.umap", group.by = paste0("predicted.",ref_cell_type_column), label = TRUE,
                                   label.size = 3, repel = TRUE) + ggtitle(paste0(train_title," Labeled Test\nTransferred Labels")) + # 2. Add Title
        scale_color_manual(values = myColors) + NoLegend() + xlim(xlim) + ylim(ylim)
      
      labeled_test_entropy_plot <- FeaturePlot(labeled_test_query, reduction = "ref.umap", features = "entropy",
                                               cols = colors_for_feature, keep.scale = "all") +
        ggtitle(paste0(train_title," Labeled Test\nEntropy")) + xlim(xlim) + ylim(ylim)
      
      
      # Save combined train and labeled test plot
      combined_train_labeled_test_plot <- train_plot + labeled_test_plot
      ggsave(paste0(output_prefix, "_train_labeled_test_regular.svg"), plot = combined_train_labeled_test_plot,
             width = 16, height = 9, units = "in")
             # width = 1600, height = 900, units = "px", dpi = 100) when png
      
      combined_train_labeled_test_entropy_plot <- train_plot + labeled_test_plot + labeled_test_entropy_plot
      ggsave(paste0(output_prefix, "_train_labeled_test_entropy.svg"), plot = combined_train_labeled_test_entropy_plot,
             width = 16, height = 9, units = "in")
      
      
      ggsave(paste0(output_prefix, "_labeled_test_pairwise_confusion.svg"), labeled_test_confusion,
             width = 16, height = 9, units = "in")
    }
    
    # Save plots
    combined_train_test_plot <- train_plot + test_plot
    ggsave(paste0(output_prefix, "_train_test_regular.svg"), plot = combined_train_test_plot,
           width = 16, height = 9, units = "in")
    
    combined_train_test_entropy_plot <- train_plot + test_plot + test_entropy_plot
    ggsave(paste0(output_prefix, "_train_test_entropy.svg"), plot = combined_train_test_entropy_plot,
           width = 16, height = 9, units = "in")
    
    
    ggsave(paste0(output_prefix, train_title, "_test_pairwise_confusion.svg"), test_confusion,
           width = 16, height = 9, units = "in")
    
    # ggsave(paste0(output_prefix, "_decision_boundaries_mapped_plot.svg"), plot = decision_boundaries_mapped,
    #        width = 16, height = 9, units = "in")

    
    # if there is an original plot from first projection
    if (!is.null(plot_on_this_UMAP)) {
      
      current_cells=names(test_query$predicted.type)
      
      if ("predicted.type" %in% colnames(plot_on_this_UMAP@meta.data)) {
        if(!is.null(plot_on_this_UMAP@reductions$ref.umap))
        { reduction='ref.umap' } else { reduction='umap' }
        
        p1 = DimPlot(plot_on_this_UMAP, cells=current_cells, reduction = reduction,
                     group.by = paste0("predicted.",ref_cell_type_column), label = TRUE,
                     label.size = 3, repel = FALSE) + ggtitle("original")
        
        # change the prediction according to this run
        plot_on_this_UMAP[[prediction_column_name]][current_cells]=test_query[[paste0("predicted.",ref_cell_type_column)]]
        
        p2 = DimPlot(plot_on_this_UMAP, cells=current_cells, reduction = reduction,
                     group.by = paste0("predicted.",ref_cell_type_column), label = TRUE,
                     label.size = 3, repel = FALSE) + ggtitle("new iteration")
        p2_alone = DimPlot(plot_on_this_UMAP, cells=current_cells, reduction = reduction,
                group.by = paste0("predicted.",ref_cell_type_column), label = TRUE,
                label.size = 5, repel = FALSE) + theme(legend.position = "none")
        
        combined_original_view_plot <- p1+p2
        ggsave(paste0(output_prefix, "_original_view_vs_iter.svg"), plot = combined_original_view_plot,
               width = 16, height = 9, units = "in")
        
        ggsave(paste0(output_prefix, "_original_view_only_prediction_A.svg"), plot = p2_alone,
               width = 8, height = 8, units = "in")
        
        
      }
      else {
        p1 = DimPlot(plot_on_this_UMAP, reduction = "umap", label = TRUE,
                     label.size = 3, repel = FALSE)
        
        
        plot_on_this_UMAP=AddMetaData(plot_on_this_UMAP,test_query$predicted.type,col.name = "predicted.type")
        p2 = DimPlot(plot_on_this_UMAP, reduction = "umap", label = TRUE, group.by = "predicted.type",
                     label.size = 8, repel = FALSE)
        
        combined_original_view_plot <- p1+p2
        ggsave(paste0(output_prefix, "_original_view_vs_iter.svg"), plot = combined_original_view_plot,
               width = 16, height = 9, units = "in")
        
        p2_alone = DimPlot(plot_on_this_UMAP, reduction = "umap", label = TRUE, group.by = "predicted.type",
                           label.size = 5, repel = FALSE) + theme(legend.position = "none")
        
        ggsave(paste0(output_prefix, "_original_view_only_prediction_A.svg"), plot = p2_alone,
               width = 8, height = 8, units = "in")
        
      }
      
      # decision_boundaries_original <- visualize_decision_boundaries(plot_on_this_UMAP,test_query)
      # ggsave(paste0(output_prefix, "_decision_boundaries_original_view_plot.svg"), plot = decision_boundaries_original,
      #        width = 16, height = 9, units = "in")
      
    }
    
    t2 <- Sys.time()
    runtime_message = sprintf("Run took %.2f %s", t2-t1, units(difftime(t2, t1)))
    print(runtime_message)
    
    return(list(train = train_labeled_seurat, labeled_test = labeled_test_query,
                test = test_query, runtime = difftime(t2, t1),
                quantile_80 = quantile_80, test_query_neighbors = test_query_neighbors, runtime_message = runtime_message,
                cell_type_order = cell_type_order,
                test_anchors = test_anchors_result))
  }
  
  # 2.5 run process_seurat_data_iter for the first time and get prediction ####
  cache_dir <- paste0(output_prefix_base, "_cache/")
  if(use_cache & !dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }
  first_process_file <- paste0(cache_dir, "train_seurat_processed.rds")
  if (use_cache && file.exists(first_process_file)) {
    print("Loading cached processed data...")
    seurat_obj <- readRDS(first_process_file)
  } else {
    print("Processing training data...")
    seurat_obj <- NormalizeData(seurat_obj)
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
    seurat_obj <- ScaleData(seurat_obj)
    seurat_obj <- RunPCA(seurat_obj)
    seurat_obj <- FindNeighbors(seurat_obj, dims = dims)
    seurat_obj <- FindClusters(seurat_obj)
    seurat_obj <- RunUMAP(seurat_obj, dims = dims, reduction.name = "umap")
    if (use_cache) {
      print("Caching processed data...")
      saveRDS(seurat_obj, first_process_file)
    }
  }
  
  second_process_file <- paste0(cache_dir, "run_without_noise.rds")
  if (use_cache && file.exists(second_process_file)) {
    print("Loading cached 'without noise' data...")
    temp_seurat_obj <- readRDS(second_process_file)
  } else {
    print("Run once without noise...")
    temp_seurat_obj <- process_seurat_data_iter(
      labeled_train_data = train_Six2GFP,
      labeled_test_data = test_Six2GFP,
      test_data = seurat_obj,
      plot_on_this_UMAP = seurat_obj,
      ref_cell_type_column = ref_cell_type_column,
      dims = dims,
      train_title = train_title,
      test_title = test_title_prefix,
      n_neighbors = n_neighbors,
      skip_neighbors = skip_neighbors,
      output_prefix = paste0(output_prefix_base, "0_silent_"), # without noise...
      min_cell_count_for_type = min_cell_count_for_type,
      use_cache = use_cache,
      return_anchors = return_anchors
    )
    if (use_cache) {
      print("Caching without noise run...")
      saveRDS(temp_seurat_obj, second_process_file)
    }
  }
  # initial_test_query = temp_seurat_obj[["test"]] # save initial test_query object
  main_cell_type_order <- temp_seurat_obj[["cell_type_order"]] # Capture cell_type_order
  
  # 3. Run noise noised_number times and keep relevant data ####
  noised_prediction = matrix('0', nrow = ncol(seurat_obj), ncol = noised_number+1)
  seurat_noised_prediction_list <- list()
  noised_anchors_list <- list()
  rownames(noised_prediction) = colnames(seurat_obj)
  noised_prediction[,1] = temp_seurat_obj[["test"]]@meta.data[[prediction_column_name]]
  for( k in c(1:noised_number)) {
    print(paste0("running run #",k))
    seurat_obj_noised = add_poisson_noise(seurat_obj)
    results_noised <- process_seurat_data_iter(
      labeled_train_data = train_Six2GFP,
      labeled_test_data = test_Six2GFP,
      test_data = seurat_obj_noised,
      plot_on_this_UMAP=seurat_obj,
      dims = dims,
      train_title = train_title,
      ref_cell_type_column = ref_cell_type_column,
      test_title = paste0(test_title_prefix,"_noised",k),
      n_neighbors = n_neighbors,
      skip_neighbors = skip_neighbors,
      colors_for_feature = colors_feature_plot_noise,
      myColors = myColors_cell_types,
      output_prefix = paste0(output_prefix_base,"noised",k,"_"),
      min_cell_count_for_type = min_cell_count_for_type,
      use_cache = use_cache,
      return_anchors = return_anchors
    )
    # # restore initial test_query object to avoid re-run initial steps ??? I think gemini is wrong here
    # results_noised[["test"]]@reductions[["ref.umap"]] = initial_test_query@reductions[["ref.umap"]]
    # results_noised[["test"]]@reductions[["ref.pca"]] = initial_test_query@reductions[["ref.pca"]]
    # results_noised[["test"]]@neighbors[["refnn"]] = initial_test_query@neighbors[["refnn"]]
    # results_noised[["test"]]@misc[["refdr"]] = initial_test_query@misc[["refdr"]]
    
    noised_prediction[,k+1] = results_noised[["test"]]@meta.data[[prediction_column_name]]
    if (return_all_suerats){
      seurat_noised_prediction_list[[paste0("iter",k)]] = results_noised
    }
    if (return_anchors && !is.null(results_noised$test_anchors)) {
      noised_anchors_list[[paste0("iter",k)]] = results_noised$test_anchors
    }
  }
  
  # 4. Check data: count changes ####
  count_changes = apply(noised_prediction, 1, function(row) {
    sum(row[-1] != row[1], na.rm = TRUE)
  })
  change_counts_df = data.frame(cell_barcode = rownames(noised_prediction), changes = count_changes)
  
  count_number_of_changes = apply(noised_prediction, 1, function(row) {
    if(length(rle(row)$lengths)>2){
      return(TRUE)
    } else {return(FALSE)}
  })
  change_counts_df$complex_changes = count_number_of_changes
  
  # 5. Plot feature plot of those who changed ####
  seurat_obj_with_changes <- AddMetaData(seurat_obj, change_counts_df$changes, col.name="Changes")
  change_feature_plot <- FeaturePlot(seurat_obj_with_changes,
                                     reduction = "umap", # original view
                                     features = "Changes",
                                     cols = colors_feature_plot_noise,
                                     keep.scale = "all") +
    ggtitle("Number of Prediction Changes with Noise")
  ggsave(paste0(output_prefix_base, "change_feature_plot.svg"), plot = change_feature_plot,
         width = 16, height = 9, units = "in")
  
  # 5.5 and plot by prediction ####
  prediction_data_original = temp_seurat_obj[["test"]]@assays[["prediction.score.type"]]@data # Use Original prediction
  rownames(prediction_data_original)=gsub("-", "_", rownames(prediction_data_original))
  prediction_data = results_noised[["test"]]@assays[["prediction.score.type"]]@data # Use noised prediction
  rownames(prediction_data)=gsub("-", "_", rownames(prediction_data))
  # all_cell_types <- union(rownames(prediction_data), rownames(prediction_data_original))
  all_cell_types <- intersect(rownames(prediction_data), rownames(prediction_data_original))
  for (t in all_cell_types){
    # Prepare data for original prediction plot
    original_prediction_values <- if(t %in% rownames(prediction_data_original)) {
      prediction_data_original[t, , drop = FALSE]
    } else {
      print("NA_original")
      rep(NA, ncol(seurat_obj))
    }
    
    # Prepare data for noised prediction plot
    noised_prediction_values <- if(t %in% rownames(prediction_data)) {
      prediction_data[t, , drop = FALSE]
    } else {
      print("NA_noised")
      rep(NA, ncol(seurat_obj))
    }
    
    noised_prediction_plot = FeaturePlot(AddMetaData(seurat_obj, t(noised_prediction_values),
                                                     col.name=paste0(t,"_prediction")),
                                         reduction = "umap",
                                         features = paste0(t,"_prediction"),
                                         keep.scale = "all") + ggtitle(paste0("Noised ",t,"_prediction")) # 2. Add Title
    original_prediction_plot = FeaturePlot(AddMetaData(seurat_obj, t(original_prediction_values),
                                                       col.name=paste0(t,"_prediction")),
                                           reduction = "umap",
                                           features = paste0(t,"_prediction"),
                                           keep.scale = "all") + ggtitle(paste0("Original ",t,"_prediction")) # 2. Add Title
    
    combined_feature_plot <- noised_prediction_plot | original_prediction_plot
    ggsave(paste0(output_prefix_base,"noised_",t, "_prediction_original_view.svg"),plot = combined_feature_plot,
           width = 16, height = 9, units = "in")
  }
  
  # 6. df that show who changed to what (proportions) ####
  plot_data = as.data.frame(noised_prediction[,c(1, noised_number+1)]) %>%
    setNames(c("Initial_Prediction", "Noised_Prediction")) %>%
    group_by(Initial_Prediction) %>%
    count(Noised_Prediction) %>%
    mutate(proportion = n/sum(n)) %>%
    ungroup()
  
  label_threshold <- 0.01
  change_proportion_plot <- ggplot(plot_data, aes(x = factor(Initial_Prediction, levels = levels(main_cell_type_order)), 
                                                  y = proportion, fill = Noised_Prediction)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    geom_text(
      aes(label = ifelse(proportion > label_threshold, 
                         Noised_Prediction,""), # If the proportion is too small, the label is an empty string.
          size = proportion), 
      position = position_stack(vjust = 0.5),
      # size = 3.5, 
      color = 'white'
      # lineheight = .8 # Adjust line spacing if using newline "\n"
    ) +
    scale_size_continuous(range = c(2.5, 5)) +
    scale_y_continuous(labels = scales::percent_format(), expand = c(0,0)) +
    scale_x_discrete(expand = c(0,0)) + 
    labs(
      title = "Transitions with Noise",
      x = "Initial Prediction",
      y = "Percentage",
      fill = "Noised Prediction"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.ticks.y = element_line(color = "black", linewidth = 0.5), # Major ticks
      axis.text.x = element_text(size = 13, angle = 45, hjust = 1),  # X-axis labels
      axis.text.y = element_text(size = 13, angle = 0, hjust = 0.4), # Y-axis labels
      axis.title.x = element_text(size = 14),                        # X-axis title
      axis.title.y = element_text(size = 14),                        # Y-axis title
      plot.title = element_text(size = 16, face = "bold")           # Main title
    )
  
  ggsave(paste0(output_prefix_base, "change_proportion_plot.svg"), plot = change_proportion_plot,
         width = 9, height = 9, units = "in")
  
  
  # 7. Confusion matrix ####
  conf_matrix <- table(Original = noised_prediction[,1],
                       New = noised_prediction[,ncol(noised_prediction)])
  
  # Handling Small Groups - Filter ROWS by original count, keep ALL columns that appear in New
  types_counts_conf_matrix <- rowSums(conf_matrix)
  rows_to_keep <- names(types_counts_conf_matrix)[types_counts_conf_matrix >= min_cell_count_for_type]
  cols_to_keep <- colnames(conf_matrix)  # keep all types appearing in noised predictions
  # Only keep columns that also exist in rows (for consistent ordering), plus any new ones
  cols_to_keep <- union(rows_to_keep, colnames(conf_matrix)[colSums(conf_matrix[rows_to_keep, , drop=FALSE]) > 0])
  conf_matrix_filtered <- conf_matrix[rows_to_keep, cols_to_keep, drop=FALSE]
  
  conf_matrix_fraction <- prop.table(conf_matrix_filtered, margin = 1) # * 100
  conf_matrix_long <- as.data.frame(as.table(conf_matrix_fraction))
  names(conf_matrix_long) <- c("Original", "New", "Fraction")
  
  # Consistent Order in Heatmaps - Order the confusion matrix plot
  conf_matrix_long$Original <- factor(conf_matrix_long$Original, levels = main_cell_type_order)
  conf_matrix_long$New <- factor(conf_matrix_long$New, levels = main_cell_type_order)
  conf_matrix_long <- conf_matrix_long %>% drop_na() # remove rows with NA after filtering
  
  confusion_matrix_plot <- ggplot(conf_matrix_long, aes(x = New, y = Original, fill = Fraction)) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
              color = "black", fill = NA, linewidth = 1, inherit.aes = FALSE) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", Fraction)),
              size = 7, family = "Arial") +
    scale_fill_gradient(low = "white",
                        high = "blue",
                        name = "Fraction") +
    theme_minimal() +
    ggtitle("Transitions with noise") + # Confusion Matrix of Predictions
    theme(legend.position = "none",
          axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 14, angle = 0, hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  ggsave(paste0(output_prefix_base, "confusion_matrix.svg"), plot = confusion_matrix_plot,
         width = 9, height = 9, units = "in")
  
  # 8. mapping score https://www.nature.com/articles/s41467-021-25089-2 implemented by ChatGPT ####
  # compute_avg_expression <- function(data, labels) {
  #   data <- as.matrix(data)
  #   if(!is.numeric(data)) {
  #     data <- matrix(as.numeric(data), nrow=nrow(data), ncol=ncol(data))
  #   }
  #   unique_labels <- sort(unique(labels))
  #   avg_expression <- t(sapply(unique_labels, function(label) {
  #     colMeans(data[labels == label, , drop = FALSE])
  #   }))
  #   rownames(avg_expression) <- unique_labels
  #   return(avg_expression)
  # }
  # 
  # calculate_mapping_score <- function(query_avg, reference_avg) {
  #   cor_matrix <- cor(t(log2(1+query_avg)), t(log2(1+reference_avg)), method = "pearson")
  #   row_sums <- rowSums(cor_matrix)
  #   col_sums <- colSums(cor_matrix)
  #   total_sum <- sum(cor_matrix)
  #   mapping_score_matrix <- cor_matrix - (row_sums %*% t(col_sums)) / total_sum
  #   mapping_score <- mean(diag(cor_matrix))
  #   return(list(correlation_matrix = cor_matrix, mapping_score = mapping_score, mapping_score_matrix = mapping_score_matrix))
  # }
  # 
  # reference_data <- t(as.matrix(GetAssayData(seurat_obj, layer = "counts")))
  # labels_reference <- factor(noised_prediction[,1])
  # query_data <- t(as.matrix(GetAssayData(results_noised[["test"]], layer = "counts")))
  # labels_query <- factor(noised_prediction[,ncol(noised_prediction)])
  # 
  # query_avg <- compute_avg_expression(query_data, labels_query)
  # reference_avg <- compute_avg_expression(reference_data, labels_reference)
  # mapping_score_results <- calculate_mapping_score(query_avg, reference_avg)
 
  # 8. mapping score (Efficient Version) - disabled ####
  
  # # Add the prediction labels to the metadata to group by them
  # seurat_obj$labels_reference <- factor(noised_prediction[, 1])
  # # Note: The 'results_noised' object is from the *last* noise iteration
  # results_noised[["test"]]$labels_query <- factor(noised_prediction[, ncol(noised_prediction)])
  # 
  # print("Calculating average expression for mapping score (efficiently)...")
  # 
  # # Use Seurat's optimized function. It returns a list, we need the 'RNA' element.
  # reference_avg_mat <- AverageExpression(seurat_obj, 
  #                                        group.by = "labels_reference", 
  #                                        assays = "RNA")$RNA
  # 
  # query_avg_mat <- AverageExpression(results_noised[["test"]], 
  #                                    group.by = "labels_query", 
  #                                    assays = "RNA")$RNA
  # 
  # # Define the mapping score function
  # calculate_mapping_score <- function(query_avg, reference_avg) {
  #   # Align genes before correlation
  #   common_genes <- intersect(rownames(query_avg), rownames(reference_avg))
  #   query_avg <- query_avg[common_genes, ]
  #   reference_avg <- reference_avg[common_genes, ]
  #   
  #   query_avg_num <- as.matrix(query_avg[common_genes, ])
  #   reference_avg_num <- as.matrix(reference_avg[common_genes, ])
  #   storage.mode(query_avg_num) <- "numeric"
  #   storage.mode(reference_avg_num) <- "numeric"
  #   # The function expects (cell_types x genes), so we transpose
  #   cor_matrix <- cor(t(log2(1 + query_avg_num)), t(log2(1 + reference_avg_num)), method = "pearson")
  #   
  #   row_sums <- rowSums(cor_matrix)
  #   col_sums <- colSums(cor_matrix)
  #   total_sum <- sum(cor_matrix)
  #   mapping_score_matrix <- cor_matrix - (row_sums %*% t(col_sums)) / total_sum
  #   mapping_score <- mean(diag(cor_matrix))
  #   return(list(correlation_matrix = cor_matrix, mapping_score = mapping_score, mapping_score_matrix = mapping_score_matrix))
  # }
  # 
  # # Now call the function with the efficiently calculated average matrices
  # mapping_score_results <- calculate_mapping_score(query_avg_mat, reference_avg_mat)
  # 
  # mapping_score_matrix_long <- as.data.frame(as.table(mapping_score_results$mapping_score_matrix))
  # names(mapping_score_matrix_long) <- c("Original", "Noised", "MappingScore")
  # 
  # # Handling Small Groups - Filter out small groups from mapping score
  # types_to_keep_mapping <- intersect(types_to_keep_conf_matrix, main_cell_type_order) # use types from confusion matrix filtering
  # mapping_score_matrix_long_filtered <- mapping_score_matrix_long %>%
  #   filter(Original %in% types_to_keep_mapping, Noised %in% types_to_keep_mapping)
  # mapping_score_matrix_long_filtered$Original <- factor(mapping_score_matrix_long_filtered$Original, levels = main_cell_type_order)
  # mapping_score_matrix_long_filtered$Noised <- factor(mapping_score_matrix_long_filtered$Noised, levels = main_cell_type_order)
  # mapping_score_matrix_long_filtered <- mapping_score_matrix_long_filtered %>% drop_na() # remove rows with NA after filtering
  # 
  # mapping_score_heatmap <- ggplot(mapping_score_matrix_long_filtered, aes(x = Noised, y = Original, fill = MappingScore)) +
  #   geom_tile(color = "white") +
  #   geom_text(aes(label = sprintf("%.3f", MappingScore)),
  #             size = 3.5) +
  #   scale_fill_gradient(low = "white",
  #                       high = "blue",
  #                       name = "MappingScore") +
  #   theme_minimal() +
  #   ggtitle("Mapping Score Heatmap") + # 2. Add Title
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1),
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank())
  # 
  # ggsave(paste0(output_prefix_base, "mapping_score_heatmap.svg"), plot = mapping_score_heatmap,
  #        width = 12, height = 9, units = "in")
  
  # 9. JSD ####
  drop_nan_rows_cols <- function(input_matrix, dimension = c("both", "row", "col")) {
    dimension <- match.arg(dimension) # Ensure dimension is one of "both", "row", "col"
    
    if (dimension %in% c("both", "row")) {
      # Identify rows where ALL elements are NaN
      all_nan_rows <- apply(input_matrix, 1, function(row) all(is.nan(row)))
      rows_to_keep <- !all_nan_rows # Keep rows that are NOT all NaN
      input_matrix <- input_matrix[rows_to_keep, , drop = FALSE] # Subset rows
    }
    
    if (dimension %in% c("both", "col")) {
      # Identify columns where ALL elements are NaN
      all_nan_cols <- apply(input_matrix, 2, function(col) all(is.nan(col)))
      cols_to_keep <- !all_nan_cols # Keep columns that are NOT all NaN
      input_matrix <- input_matrix[, cols_to_keep, drop = FALSE] # Subset columns
    }
    
    return(input_matrix)
  }
  calculate_entropy_jsd <- function(query) {
    log_pred <- log2(query)
    log_pred[is.infinite(log_pred)] <- 0
    entropy_by_cell <- -colSums(query * log_pred)
    return(entropy_by_cell)
  }
  
  prediction_data = temp_seurat_obj[["test"]]@assays[["prediction.score.type"]]@data # Use original prediction
  rownames(prediction_data)=gsub("-", "_", rownames(prediction_data))
  # P^R, i
  norm_p_factor = rowSums(prediction_data)
  norm_prediction_table = sweep(prediction_data, 1, norm_p_factor, FUN = "/")
  # P^C, j
  # prediction_data = results_noised[["test"]]@assays[["prediction.score.type"]]@data # Use the last noised result
  # rownames(prediction_data)=gsub("-", "_", rownames(prediction_data))
  predicted_zeroes_mat <- (sapply(temp_seurat_obj[["test"]]@meta.data[["predicted.type"]], function(x) rownames(prediction_data) == x) * 1)
  dimnames(predicted_zeroes_mat) = dimnames(prediction_data)
  norm_one_hot_factor = rowSums(predicted_zeroes_mat)
  norm_one_hot = sweep(predicted_zeroes_mat, 1, norm_one_hot_factor, FUN = "/")
  
  JSD_mat=matrix(0, nrow = nrow(norm_prediction_table), ncol = nrow(norm_prediction_table))
  dimnames(JSD_mat) = list(row.names(norm_prediction_table),row.names(norm_prediction_table))
  for (i in row.names(norm_prediction_table)){
    for (j in row.names(norm_prediction_table)){
      JSD_i_j = calculate_entropy_jsd(as.matrix(norm_prediction_table[i,] + norm_one_hot[j,])/2) -
        (calculate_entropy_jsd(as.matrix(norm_prediction_table[i,]))+calculate_entropy_jsd(as.matrix(norm_one_hot[j,])))/2
      JSD_mat[i,j] = JSD_i_j
    }
  }
  RSS_mat = 1-sqrt(JSD_mat)
  
  # Remove NaN rows/cols before plotting
  JSD_mat_filtered <- drop_nan_rows_cols(JSD_mat, "both") # remove rows with NaN
  RSS_mat_filtered <- drop_nan_rows_cols(RSS_mat, "both") # also filter RSS_mat accordingly
  
  JSD_long <- as.data.frame(as.table(JSD_mat_filtered))
  names(JSD_long) <- c("P_R", "P_C", "JSD")
  
  # Consistent Order in Heatmaps - Order JSD heatmap
  JSD_long$P_R <- factor(JSD_long$P_R, levels = main_cell_type_order)
  JSD_long$P_C <- factor(JSD_long$P_C, levels = main_cell_type_order)
  JSD_long <- JSD_long %>% drop_na() # remove rows with NA after filtering
  
  jsd_heatmap <- ggplot(JSD_long, aes(x = P_C, y = P_R, fill = JSD)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", JSD)),
              size = 6) +
    scale_fill_gradient(low = "white",
                        high = "blue",
                        name = "JSD") +
    theme_minimal() +
    ggtitle("JSD Heatmap") + 
    theme(legend.position = "none",
          axis.text.x = element_text(size = 13, angle = 0, hjust = 1),
          axis.text.y = element_text(size = 13, angle = 90, hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  ggsave(paste0(output_prefix_base, "jsd_heatmap.svg"), plot = jsd_heatmap,
         width = 12, height = 9, units = "in")
  
  # filter out small groups
  types_to_run_on = intersect(colnames(RSS_mat_filtered), rownames(conf_matrix_fraction)) # use rownames from filtered conf_matrix
  types_to_run_on = intersect(types_to_run_on, main_cell_type_order) # ensure order is consisten
  
  RSS_long <- as.data.frame(as.table(RSS_mat_filtered[types_to_run_on, types_to_run_on]))
  names(RSS_long) <- c("P_R", "P_C", "RSS")
  
  # Consistent Order in Heatmaps - Order RSS heatmap
  RSS_long$P_R <- factor(RSS_long$P_R, levels = main_cell_type_order)
  RSS_long$P_C <- factor(RSS_long$P_C, levels = main_cell_type_order)
  RSS_long <- RSS_long %>% drop_na() # remove rows with NA after filtering
  RSS_long <- RSS_long %>%
    mutate(
      RSS_label = sprintf("%.2f", RSS), 
      RSS_label = ifelse(RSS_label == "-0.00", "0.00", RSS_label) # correct the "-0.00" issue
    )
  
  rss_heatmap <- ggplot(RSS_long, aes(x = P_C, y = P_R, fill = RSS)) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
              color = "black", fill = NA, linewidth = 1, inherit.aes = FALSE) +
    geom_tile(color = "white") +
    geom_text(aes(label = RSS_label),
              size = 7, family = "Arial") +
    scale_fill_gradient(low = "white",
                        high = "blue",
                        name = "RSS") +
    theme_minimal() +
    ggtitle("Jensen Shannon Similarity") +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 14, angle = 0, hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  ggsave(paste0(output_prefix_base, "pss_heatmap.svg"), plot = rss_heatmap,
         width = 9, height = 9, units = "in")
  
  # 10. PSS vs Stability ####
  
  # Handling Small Groups - Filter out small groups from stability plot
  # types_to_run_on = intersect(colnames(RSS_mat_filtered), rownames(conf_matrix_fraction)) # use rownames from filtered conf_matrix
  # types_to_run_on = intersect(types_to_run_on, main_cell_type_order) # ensure order is consistent
  # types_to_run_on_filtered <- types_to_run_on[types_to_run_on %in% types_to_keep_conf_matrix] # final filter based on confusion matrix filter
  # types_to_run_on_filtered <- factor(types_to_run_on_filtered, levels = main_cell_type_order) # keep order
  
  print(paste0("running on types:[",paste(types_to_run_on,sep =","),"]"))
  stability_pred = matrix(0, nrow = length(types_to_run_on),
                          ncol = 4,
                          dimnames = list(types_to_run_on,c("PSS","Stability","Bi_Stability","F1_Stability")))
  stability_pred = as.data.frame(stability_pred)
  for (i in seq_along(types_to_run_on)){
    type = types_to_run_on[i]
    # print(paste0("calculating stability for...", type))
    stability_pred[i,1] = RSS_mat_filtered[type, type] # from stage #9
    stability_pred[i,2] = conf_matrix_fraction[type, type] # from stage #5.5
    
    # Bidirectional stability: penalizes both outflow AND inflow
    # out = cells originally type X that changed to something else
    # in  = cells from other types that became type X after noise
    # Bi_Stability = max(0, 1 - (out + in) / original)
    stayed <- conf_matrix_filtered[type, type]           # diagonal (raw count)
    original_count <- sum(conf_matrix_filtered[type, ])  # row sum
    noised_count <- sum(conf_matrix_filtered[, type])    # column sum
    out_count <- original_count - stayed
    in_count <- noised_count - stayed
    stability_pred[i,3] = max(0, 1 - (out_count + in_count) / original_count)
    
    # F1 Stability: harmonic mean of precision and recall
    # TP = stayed (diagonal), FP = in_count (inflow), FN = out_count (outflow)
    # F1 = 2*TP / (2*TP + FP + FN)
    tp <- stayed
    fp <- in_count
    fn <- out_count
    stability_pred[i,4] = if ((2*tp + fp + fn) > 0) (2*tp) / (2*tp + fp + fn) else 0
  }
  
  # 10.1 Per-run stability (for variability / error bars) ####
  stability_per_run <- matrix(NA, nrow = length(types_to_run_on), ncol = noised_number,
                              dimnames = list(types_to_run_on, paste0("run_", 1:noised_number)))
  bi_stability_per_run <- matrix(NA, nrow = length(types_to_run_on), ncol = noised_number,
                                  dimnames = list(types_to_run_on, paste0("run_", 1:noised_number)))
  f1_stability_per_run <- matrix(NA, nrow = length(types_to_run_on), ncol = noised_number,
                                  dimnames = list(types_to_run_on, paste0("run_", 1:noised_number)))
  for (run_idx in 1:noised_number) {
    run_conf <- table(Original = noised_prediction[, 1],
                      New = noised_prediction[, run_idx + 1])
    # Filter to types we're working with
    common_types <- intersect(intersect(rownames(run_conf), colnames(run_conf)), types_to_run_on)
    if (length(common_types) > 0) {
      run_conf_filtered <- run_conf[common_types, common_types, drop = FALSE]
      run_conf_frac <- prop.table(run_conf_filtered, margin = 1)
      for (type in common_types) {
        stability_per_run[type, run_idx] <- run_conf_frac[type, type]
        # Bidirectional stability per run
        stayed_run <- run_conf_filtered[type, type]
        original_run <- sum(run_conf_filtered[type, ])
        noised_run <- sum(run_conf_filtered[, type])
        out_run <- original_run - stayed_run
        in_run <- noised_run - stayed_run
        bi_stability_per_run[type, run_idx] <- max(0, 1 - (out_run + in_run) / original_run)
        # F1 stability per run: 2*TP / (2*TP + FP + FN)
        tp_run <- stayed_run
        fp_run <- in_run
        fn_run <- out_run
        f1_stability_per_run[type, run_idx] <- if ((2*tp_run + fp_run + fn_run) > 0) (2*tp_run) / (2*tp_run + fp_run + fn_run) else 0
      }
    }
  }
  
  stability_plot = ggplot(stability_pred, aes(x=PSS, y=Stability)) + geom_point() + theme_minimal() +
    ggtitle("Stability vs PSS") + geom_text(label=rownames(stability_pred), vjust = 1.5)
  
  ggsave(paste0(output_prefix_base, "stability_vs_pss_scatter_plot.svg"), plot = stability_plot,
         width = 12, height = 9, units = "in")

  bi_stability_plot = ggplot(stability_pred, aes(x=PSS, y=Bi_Stability)) + geom_point() + theme_minimal() +
    ggtitle("Bidirectional Stability vs PSS") + geom_text(label=rownames(stability_pred), vjust = 1.5)
  
  ggsave(paste0(output_prefix_base, "bi_stability_vs_pss_scatter_plot.svg"), plot = bi_stability_plot,
         width = 12, height = 9, units = "in")

  f1_stability_plot = ggplot(stability_pred, aes(x=PSS, y=F1_Stability)) + geom_point() + theme_minimal() +
    ggtitle("F1 Stability vs PSS") + geom_text(label=rownames(stability_pred), vjust = 1.5)
  
  ggsave(paste0(output_prefix_base, "f1_stability_vs_pss_scatter_plot.svg"), plot = f1_stability_plot,
         width = 12, height = 9, units = "in")
  
  # plot also the full table
  # conf_matrix_fraction_ordered = conf_matrix_fraction[types_to_run_on,types_to_run_on, drop=FALSE]
  # RSS_mat_filtered_ordered = RSS_mat_filtered[types_to_run_on,types_to_run_on, drop=FALSE]
  # conf_matrix_long <- as.data.frame(as.table(conf_matrix_fraction_ordered))
  # RSS_long <- as.data.frame(as.table(RSS_mat_filtered_ordered)) 
  # names(RSS_long) <- c("P_R", "P_C", "PSS")
  # 
  # rss_and_stab = data.frame(conf_matrix_long,RSS_long[,3])
  # colnames(rss_and_stab) = c("Original","New","Freq","PSS")
  
  conf_matrix_fraction_ordered <- conf_matrix_fraction[types_to_run_on, types_to_run_on, drop=FALSE]
  RSS_mat_filtered_ordered <- RSS_mat_filtered[types_to_run_on, types_to_run_on, drop=FALSE]
  
  # Melt BOTH ordered matrices into long format
  conf_long <- as.data.frame(as.table(conf_matrix_fraction_ordered))
  names(conf_long) <- c("Original", "New", "Freq")
  
  pss_long <- as.data.frame(as.table(RSS_mat_filtered_ordered))
  names(pss_long) <- c("Original", "New", "PSS")
  
  # Safely JOIN the two data frames by their common keys
  rss_and_stab <- dplyr::left_join(conf_long, pss_long, by = c("Original", "New"))
  rss_and_stab$is_diagonal <- rss_and_stab$Original == rss_and_stab$New
  
  # Linear fit 
  stability_data <- rss_and_stab %>% 
    filter(is_diagonal == TRUE)
  fit <- lm(Freq ~ PSS, data = stability_data)
  fit_summary <- summary(fit)
  intercept <- coef(fit)[1]
  slope <- coef(fit)[2]
  r_squared <- fit_summary$r.squared
  p_value_pss=fit_summary$coefficients["PSS","Pr(>|t|)"]
  equation_string <- sprintf("y = %.2f x + %.2f", slope, intercept)
  r_squared_string <- sprintf("R² = %.3f ; Pval = %.3f", r_squared, p_value_pss)
  
  spearman_test_diagonal <- cor.test(x = stability_pred$PSS, 
                                     y = stability_pred$Stability, 
                                     method = "spearman",
                                     exact = FALSE) # Use exact=FALSE for robustness if there are ties
  spearman_corr_diagonal <- spearman_test_diagonal$estimate  # The rho coefficient
  spearman_pval_diagonal <- spearman_test_diagonal$p.value    # The p-value
  spearman_corr_all = cor.test(x=rss_and_stab[,3], y = rss_and_stab[,4], method = c("spearman"),exact = FALSE)
  # significant_threshold = 0.20
  # over_threshold <- rss_and_stab$Freq >= significant_threshold
  # spearman_corr_significant <- cor(x = rss_and_stab$Freq[over_threshold], 
  #                                  y = rss_and_stab$PSS[over_threshold], 
  #                                  method = "spearman")
  
  # Consistent Order in Heatmaps - Order full stability plot
  rss_and_stab$Original <- factor(rss_and_stab$Original, levels = main_cell_type_order)
  rss_and_stab$New <- factor(rss_and_stab$New, levels = main_cell_type_order)
  rss_and_stab <- rss_and_stab %>% drop_na()
  
  plot_all_stab =  ggplot(rss_and_stab, aes(x=PSS, y=Freq)) +
    # geom_rect(
    #   data = data.frame(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    #   aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    #   color = "black",
    #   fill = NA, linewidth = 1, inherit.aes = FALSE) +
    geom_point(aes(color = is_diagonal), size = 3.5, alpha = 0.8) + theme_minimal() + 
    scale_color_manual(
      name = "Transition Type", 
      values = c("TRUE" = "red", "FALSE" = "black"),
      labels = c("TRUE" = "Self-Transition (Stability)", "FALSE" = "Cross-Transition")
    ) +
    geom_text(label=paste0(rss_and_stab$Original,"-", rss_and_stab$New), vjust = 1.5, size = 12/.pt) + 
    # geom_hline(yintercept = significant_threshold, linetype = "dashed", color = "red", size = 0.7) +
    # annotate("text", x = -Inf, y = significant_threshold, label = paste0("Threshold = ",significant_threshold),
    #          hjust = -0.1, vjust = -0.5, size = 9, color = "red") +
    # annotate("text", x = Inf, y = -Inf,
    #          label = paste("Spearman above threshold:", round(spearman_corr_significant, 3)),
    #          hjust = 1.05, vjust = -2.5, size = 9) +
    annotate("text", x = Inf, y = -Inf,
             label = paste("Spearman Coefficient:", round(spearman_corr_diagonal, 3), "; pVal:", round(spearman_pval_diagonal, 3)),
             hjust = 1.05, vjust = -1.5, size = 9) +
    annotate("text", x = Inf, y = -Inf,
             label = paste("Spearman all:", round(spearman_corr_all$estimate, 3), "; pVal:", round(spearman_corr_all$p.value, 3)),
             hjust = 1.05, vjust = -0.5, size = 9) +
    xlab("Prediction specificity score") + 
    ylab("Fraction of transitions with noise")+
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 13, angle = 45),
          axis.text.y = element_text(size = 13, angle = 0, hjust = 0.5),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  ggsave(paste0(output_prefix_base, "all_freq_vs_pss_scatter_plot.svg"), plot = plot_all_stab,
         width = 12, height = 8, units = "in")
  
  plot_all_stab_with_fit <- plot_all_stab + 
    geom_abline(intercept = intercept, slope = slope, color = "blue", linetype = "dashed", size = 1) +
    annotate("text", x = Inf, y = Inf,
             label = paste(equation_string, r_squared_string, sep = "\n"),
             hjust = 1.05, vjust = 1.2, # Adjust to place it neatly inside the plot area
             size = 4, color = "blue")
  
  ggsave(paste0(output_prefix_base, "all_freq_vs_pss_scatter_plot_with_fit.svg"), 
         plot = plot_all_stab_with_fit,
         width = 12, height = 8, units = "in")
  
  stability_only_plot <- ggplot(stability_data, aes(x = PSS, y = Freq)) +
    geom_point(color = "red", size = 4, alpha = 0.7) +
    geom_abline(intercept = intercept, slope = slope, color = "blue", linetype = "dashed", size = 1) +
    annotate("text", 
             x = min(stability_data$PSS), # Position at the minimum x
             y = max(stability_data$Freq),  # Position at the maximum y
             label = paste(equation_string, r_squared_string, sep = "\n"), # '\n' creates a new line
             hjust = 0, vjust = 1, # Align text to top-left corner
             size = 5,
             parse = FALSE) +
    geom_text_repel(aes(label = Original), size = 3.5, box.padding = 0.5) +
    labs(
      title = "Stability vs. Prediction Specificity Score (PSS)",
      subtitle = "Analysis of self-transitions (Original = New Prediction)",
      x = "Prediction Specificity Score (PSS)",
      y = "Stability (% of cells remaining same type after noise)"
    ) +
    theme_minimal()
  ggsave(paste0(output_prefix_base, "stability_pss_linear_fit.svg"), plot = stability_only_plot,
         width = 8, height = 7, units = "in")
  
  
  # 11. set chosen plots in place ####
  source_dir = output_prefix_base
  target_dir = paste0(source_dir,"/main_plots")
  if (!dir.exists(target_dir)) {
    dir.create(target_dir, recursive = TRUE)
  }
  files_to_move = c("0_silent__original_view_only_prediction_A.svg",
                    "all_freq_vs_pss_scatter_plot_with_fit.svg",
                    "change_proportion_plot.svg",
                    "confusion_matrix.svg",
                    "pss_heatmap.svg",
                    "bi_stability_vs_pss_scatter_plot.svg",
                    "f1_stability_vs_pss_scatter_plot.svg")
  moved_files <- character(0)
  failed_files <- character(0)
  
  for (file in files_to_move) {
    source_path <- file.path(source_dir, file)
    target_path <- file.path(target_dir, file)
    
    # Check if source file exists
    if (!file.exists(source_path)) {
      cat("File not found, skipping:", source_path, "\n")
      failed_files <- c(failed_files, file)
      next
    }
    
    tryCatch({
      file.copy(source_path, target_path, overwrite=TRUE)
      moved_files <- c(moved_files, file)
      cat("Copied:", file, "\n")
    }, error = function(e) {
      cat("Failed to copy", file, ":", e$message, "\n")
      failed_files <- c(failed_files, file)
    })
  }
  
  #####
  # temp_seurat_obj saved as run_without_noise.rds
  # seurat_obj saved as train_seurat_processed.rds
  return(list(
    noised_prediction_matrix = noised_prediction,
    change_counts = change_counts_df,
    # mapping_score_results = mapping_score_results,
    jsd_matrix = JSD_mat,
    pss_matrix = RSS_mat,
    stability_pred = stability_pred,
    stability_per_run = stability_per_run,
    bi_stability_per_run = bi_stability_per_run,
    f1_stability_per_run = f1_stability_per_run,
    seurat_noised_prediction_list = seurat_noised_prediction_list,
    initial_anchors = temp_seurat_obj$test_anchors,
    noised_anchors = noised_anchors_list
  ))
}
