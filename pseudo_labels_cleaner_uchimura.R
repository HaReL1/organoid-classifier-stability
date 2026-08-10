# load uchimura
# clean it
# do flow

{
  Uchimura_Humphreys_20 = read.table("GSM3763147_Uchimura.dge.txt")
  seurat_Uchimura_Humphreys_20 <- CreateSeuratObject(counts = Uchimura_Humphreys_20, assay = "RNA")
  # remove(Uchimura_Humphreys_20)
} # load Full Uchimura (seurat_Uchimura_Humphreys_20)
# FeaturePlot(seurat_Uchimura_Humphreys_20, reduction = "umap", features = "MLANA",
#             keep.scale = "all") + ggtitle("MLANA")
# 
# FeaturePlot(seurat_Uchimura_Humphreys_20, reduction = "umap", features = "MAP2",
#             keep.scale = "all")
# 
# FeaturePlot(seurat_Uchimura_Humphreys_20, reduction = "umap", features = "MYOG",
#             keep.scale = "all")

# Done once and saved. commented out 
# MLANA_group = CellSelector(FeaturePlot(seurat_Uchimura_Humphreys_20, reduction = "umap", features = "MLANA",
#                          keep.scale = "all"))
# Neurons_group = CellSelector(FeaturePlot(seurat_Uchimura_Humphreys_20, reduction = "umap", features = "SOX2",
#                                          keep.scale = "all"))
# Muscle_group = CellSelector(FeaturePlot(seurat_Uchimura_Humphreys_20, reduction = "umap", features = "PITX2",
#                                         keep.scale = "all"))
# extra = CellSelector(FeaturePlot(seurat_Uchimura_Humphreys_20, reduction = "umap", features = "PITX2",
#                                  keep.scale = "all"))
# uchimura_off_target_cell_list = unique(c(MLANA_group,Neurons_group,Muscle_group,extra))
# save(uchimura_off_target_cell_list,file = "uchimura_off_target_cell_list")

{
  Uchimura_Humphreys_20 = read.table("GSM3763147_Uchimura.dge.txt")
  load(file = "uchimura_off_target_cell_list")
  seurat_Uchimura_Humphreys_20 <- CreateSeuratObject(counts = Uchimura_Humphreys_20[,!colnames(Uchimura_Humphreys_20) %in% uchimura_off_target_cell_list], 
                                                     assay = "RNA")
  remove(Uchimura_Humphreys_20, uchimura_off_target_cell_list)
} # load clean Uchimura (seurat_Uchimura_Humphreys_20)

{
  seurat_Uchimura_Humphreys_20 <- NormalizeData(seurat_Uchimura_Humphreys_20)
  seurat_Uchimura_Humphreys_20 <- FindVariableFeatures(seurat_Uchimura_Humphreys_20, selection.method = "vst", nfeatures = 2000)
  seurat_Uchimura_Humphreys_20 <- ScaleData(seurat_Uchimura_Humphreys_20)
  seurat_Uchimura_Humphreys_20 <- RunPCA(seurat_Uchimura_Humphreys_20)
  seurat_Uchimura_Humphreys_20 <- FindNeighbors(seurat_Uchimura_Humphreys_20, dims = 1:30)
  seurat_Uchimura_Humphreys_20 <- FindClusters(seurat_Uchimura_Humphreys_20)
  seurat_Uchimura_Humphreys_20 <- RunUMAP(seurat_Uchimura_Humphreys_20, dims = 1:30, reduction.name = "umap")
  DimPlot(seurat_Uchimura_Humphreys_20, reduction = "umap", label = TRUE,
          label.size = 3, repel = FALSE)
} # seurat flow

# FeaturePlot(seurat_Uchimura_Humphreys_20, reduction = "umap", features = "TMEM52B",
#             keep.scale = "all") # LOH
# FeaturePlot(seurat_Uchimura_Humphreys_20, reduction = "umap", features = "PODXL",
#             keep.scale = "all") # PODO
# FeaturePlot(seurat_Uchimura_Humphreys_20, reduction = "umap", features = "GATA3",
#             keep.scale = "all") # Collecting duct
# FeaturePlot(seurat_Uchimura_Humphreys_20, reduction = "umap", features = "GAP43",
#             keep.scale = "all") # neron (empty)
# FeaturePlot(seurat_Uchimura_Humphreys_20, reduction = "umap", features = "GATM",
#             keep.scale = "all") # PT
# FeaturePlot(seurat_Uchimura_Humphreys_20, reduction = "umap", features = "COL1A1",
#             keep.scale = "all") # Mes prog
# FeaturePlot(seurat_Uchimura_Humphreys_20, reduction = "umap", features = "MKI67",
#             keep.scale = "all") # Mes
FeaturePlot(seurat_Uchimura_Humphreys_20, reduction = "umap", features = "CDH11", keep.scale = "all")
FeaturePlot(seurat_Uchimura_Humphreys_20, reduction = "umap", features = "PAX2", keep.scale = "all") 


# now: mapquery from six2gfp
{
  load(file = "six2gfp/train_Six2GFP")
  load(file = "six2gfp/test_Six2GFP")
} # load six2gfp (test_Six2GFP, train_Six2GFP)

process_seurat_data_iter <- function(labeled_train_data, test_data, labeled_test_data = NULL, 
                                     dims = 1:30, 
                                     plot_on_this_UMAP = NULL,
                                     iteration = NULL,
                                     train_title = "SIX2GFP",
                                     ref_cell_type_column = "type",
                                     test_title = "Uchimura",
                                     n_neighbors = 8,
                                     skip_neighbors = FALSE,
                                     colors_for_feature = NULL, 
                                     myColors = NULL,
                                     output_prefix = "output") {
  
  library(Seurat)
  library(ggplot2)
  library(RANN)
  library(dplyr)
  library(gridExtra)
  library(reshape2)
  library(RColorBrewer)
  library(rlang)
  
  t1 <- Sys.time()
  
  # Process training data. Input is Seurat object
  train_labeled_seurat <- NormalizeData(labeled_train_data)
  train_labeled_seurat <- FindVariableFeatures(train_labeled_seurat, selection.method = "vst", nfeatures = 2000)
  train_labeled_seurat <- ScaleData(train_labeled_seurat)
  train_labeled_seurat <- RunPCA(train_labeled_seurat)
  train_labeled_seurat <- FindNeighbors(train_labeled_seurat, dims = dims)
  train_labeled_seurat <- FindClusters(train_labeled_seurat)
  train_labeled_seurat <- RunUMAP(train_labeled_seurat, dims = dims, return.model = TRUE)

  
  if(!is.null(labeled_test_data)){
    # Process labeled test data
    labeled_test_data <- NormalizeData(labeled_test_data)
    
    labeled_test_anchors <- FindTransferAnchors(reference = train_labeled_seurat, query = labeled_test_data, dims = dims,
                                                reference.reduction = "pca")
    labeled_test_query <- MapQuery(anchorset = labeled_test_anchors, reference = train_labeled_seurat, query = labeled_test_data,
                                 refdata = setNames(list(ref_cell_type_column),ref_cell_type_column), 
                                 reference.reduction = "pca", reduction.model = "umap")
  }
  else
  {
    labeled_test_query = NULL
  }
  
  # Process test data
  test_data <- NormalizeData(test_data)
  
  test_anchors <- FindTransferAnchors(reference = train_labeled_seurat, query = test_data, dims = dims,
                                      reference.reduction = "pca")
  test_query <- MapQuery(anchorset = test_anchors, reference = train_labeled_seurat, query = test_data,
                         refdata = setNames(list(ref_cell_type_column),ref_cell_type_column), reference.reduction = "pca", 
                         reduction.model = "umap")
  
  if (is.null(myColors)) {
    myColors <- brewer.pal(nrow(unique(train_labeled_seurat[[ref_cell_type_column]])), "Paired")
    names(myColors) <- unique(train_labeled_seurat[[ref_cell_type_column]])[[ref_cell_type_column]]
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
    all_types = unique(train_labeled_seurat[[ref_cell_type_column]])$type
    
    all_neighbors_dist <- data.frame(matrix(NA, nrow = nrow(test_data), ncol = length(all_types)))
    # all_neighbors_names <- data.frame(matrix(NA, nrow = nrow(test_data), ncol = n_neighbors))
    
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
      # winner_type <- query@meta.data[[paste0("predicted.",ref_cell_type_column)]][i]
      # train_same_as_winner <- (train_labeled_seurat[[ref_cell_type_column]] == winner_type)
      # training_relevant_data <- train_data[train_same_as_winner, ]
      # knn_closest <- nn2(training_relevant_data, as.data.frame(t(test_data[i, ])), k = n_neighbors)
      # all_neighbors_dist[i, ] <- knn_closest$nn.dist[1, ]
      # all_neighbors_names[i, ] <- rownames(training_relevant_data)[knn_closest$nn.idx[1, ]]
      # if (i %% 50 == 0) {
      #   print(paste0("Row:", i, " Time:", Sys.time()))
      # }
    }
    
    # query <- AddMetaData(query, rowMeans(all_neighbors_dist), col.name = "closest_neighbor")
    colnames(all_neighbors_dist)=all_types
    return(all_neighbors_dist)
  }
  
  # Create confusion matrix
  create_confusion_matrix <- function(query, title) {
    probs_tub <- query@assays[[paste0("prediction.score.",ref_cell_type_column)]]@data
    num_classes <- nrow(probs_tub)
    pairwise_confusion <- matrix(0, nrow = num_classes, ncol = num_classes)
    row.names(pairwise_confusion) <- gsub("-", "_", row.names(probs_tub))
    colnames(pairwise_confusion) <- row.names(pairwise_confusion)
    
    for (i in 1:num_classes) {
      for (j in 1:num_classes) {
        if (i != j) {
          pairwise_confusion[i, j] <- mean(probs_tub[i, ] * probs_tub[j, ])
        }
      }
    }
    
    pairwise_confusion_df <- melt(pairwise_confusion)
    p <- ggplot(data = pairwise_confusion_df, aes(Var1, Var2, fill = value)) +
      geom_tile() +
      geom_text(aes(label = round(value,5))) +
      scale_fill_gradient(low = "white", high = "red") +
      labs(x = "Class", y = "Class", title = paste("Pairwise Confusion Matrix -", title)) +
      theme_minimal()
    return(p)
  }
  
  # Calculate plot limits dynamically
  plot_limits <- calculate_plot_limits(train_labeled_seurat, labeled_test_query, test_query)
  xlim <- plot_limits$xlim
  ylim <- plot_limits$ylim
  
  # create plots
  train_plot <- DimPlot(train_labeled_seurat, reduction = "umap", group.by = ref_cell_type_column, label = TRUE, label.size = 3,
                        repel = TRUE) + ggtitle(paste0(train_title,"\n(train) annotations")) + 
    scale_color_manual(values = myColors) +
    theme(legend.position = "left", legend.text = element_text(size=8)) + xlim(xlim) + ylim(ylim)
  
  test_plot <- DimPlot(test_query, reduction = "ref.umap", group.by = paste0("predicted.",ref_cell_type_column), label = TRUE,
                       label.size = 3, repel = TRUE) + ggtitle(paste(test_title, "\nTransferred Labels")) + 
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
  # test_dist_plot <- FeaturePlot(test_query, reduction = "ref.umap", features = "closest_neighbor",
  #                               cols = colors_for_feature, keep.scale = "all") + 
  #   ggtitle(paste(test_title, "\nclosest neighbor")) + xlim(xlim) + ylim(ylim)
  
  test_confusion <- create_confusion_matrix(test_query, test_title)
  
  if (!is.null(labeled_test_query)) {
    # labeled_test_query <- calculate_neighbors(train_data, labeled_test_query)
    labeled_test_query <- calculate_entropy(labeled_test_query)
    labeled_test_confusion <- create_confusion_matrix(labeled_test_query, "Labeled Test")
    
    # plot labeled_test related plots
    labeled_test_plot <- DimPlot(labeled_test_query, reduction = "ref.umap", group.by = paste0("predicted.",ref_cell_type_column), label = TRUE,
                                 label.size = 3, repel = TRUE) + ggtitle(paste0(train_title," Labeled Test\nTransferred Labels")) + 
      scale_color_manual(values = myColors) + NoLegend() + xlim(xlim) + ylim(ylim)

    labeled_test_entropy_plot <- FeaturePlot(labeled_test_query, reduction = "ref.umap", features = "entropy",
                                             cols = colors_for_feature, keep.scale = "all") + 
      ggtitle(paste0(train_title," Labeled Test\nEntropy")) + xlim(xlim) + ylim(ylim)
    
    # labeled_test_dist_plot <- FeaturePlot(labeled_test_query, reduction = "ref.umap", features = "closest_neighbor",
    #                                       cols = colors_for_feature, keep.scale = "all") + 
    #   ggtitle("Labeled Test\nclosest neighbor") + xlim(xlim) + ylim(ylim)
    
    # Save combined train and labeled test plot
    train_plot + labeled_test_plot
    ggsave(paste0(output_prefix, "_train_labeled_test_regular.png"), 
           width = 1600, height = 900, units = "px", dpi = 100)
    
    train_plot + labeled_test_plot + labeled_test_entropy_plot
    ggsave(paste0(output_prefix, "_train_labeled_test_entropy.png"), 
           width = 1600, height = 900, units = "px", dpi = 100)
    
    # ggsave(paste0(output_prefix, "_train_labeled_test_dist.png"), 
    #        train_plot + labeled_test_plot + labeled_test_dist_plot, 
    #        width = 1600, height = 900, units = "px", dpi = 100)
    
    ggsave(paste0(output_prefix, "_labeled_test_pairwise_confusion.jpg"), labeled_test_confusion, 
           width = 1600, height = 900, units = "px", dpi = 100)
  }
  
  # Save plots
  train_plot + test_plot
  ggsave(paste0(output_prefix, "_train_test_regular.png"), 
         width = 1600, height = 900, units = "px", dpi = 100)

  train_plot + test_plot + test_entropy_plot
  ggsave(paste0(output_prefix, "_train_test_entropy.png"), 
         width = 1600, height = 900, units = "px", dpi = 100)
  
  # ggsave(paste0(output_prefix, "_train_test_distance.png"), 
  #        train_plot + test_plot + test_dist_plot, 
  #        width = 1600, height = 900, units = "px", dpi = 100)
  
  ggsave(paste0(output_prefix, train_title, "_test_pairwise_confusion.jpg"), test_confusion, 
         width = 1600, height = 900, units = "px", dpi = 100)
  
  # if there is an original plot from first projection
  if (!is.null(plot_on_this_UMAP)) {
    # plot the results_takasato_six2gfp_full with the new "type" from the new iteration
    # full_original_projection = results_takasato_six2gfp_full$test
    
    # TODO: support ref_cell_type_column instead of predicted.type
    
    current_cells=names(test_query$predicted.type)
    
    if ("predicted.type" %in% colnames(plot_on_this_UMAP@meta.data)) {
      if(!is.null(plot_on_this_UMAP@reductions$ref.umap))
      { reduction='ref.umap' } else { reduction='umap' }
      
      p1 = DimPlot(plot_on_this_UMAP, cells=current_cells, reduction = reduction, 
                   group.by = paste0("predicted.",ref_cell_type_column), label = TRUE,
                   label.size = 3, repel = FALSE) + ggtitle("original")
      
      # change the prediction according to this run
      plot_on_this_UMAP$predicted.type[current_cells]=test_query$predicted.type
      
      p2 = DimPlot(plot_on_this_UMAP, cells=current_cells, reduction = reduction, 
                   group.by = paste0("predicted.",ref_cell_type_column), label = TRUE, 
                   label.size = 3, repel = FALSE) + ggtitle("new iteration")
      
      p1+p2
      ggsave(paste0(output_prefix, "_original_view_vs_iter.png"), 
             width = 1600, height = 900, units = "px", dpi = 100)
      
    }
    else {
      p1 = DimPlot(plot_on_this_UMAP, reduction = "umap", label = TRUE,
              label.size = 3, repel = FALSE)
      
      ## if the object doesn't have predicted.type, we need to add
      plot_on_this_UMAP=AddMetaData(plot_on_this_UMAP,test_query$predicted.type,col.name = "predicted.type")
      p2 = DimPlot(plot_on_this_UMAP, reduction = "umap", label = TRUE, group.by = "predicted.type",
              label.size = 3, repel = FALSE)
      
      p1+p2
      ggsave(paste0(output_prefix, "_original_view_vs_iter.png"), 
             width = 1600, height = 900, units = "px", dpi = 100)

    }

  }
  
  t2 <- Sys.time()
  sprintf("Run took %.2f %s", t2-t1, units(difftime(t2, t1)))
  
  return(list(train = train_labeled_seurat, labeled_test = labeled_test_query, 
              test = test_query, runtime = difftime(t2, t1), 
              quantile_80 = quantile_80, test_query_neighbors = test_query_neighbors))
}

results_Uchimura_six2gfp_full <- process_seurat_data_iter(
  labeled_train_data = train_Six2GFP,
  labeled_test_data = test_Six2GFP,
  test_data = seurat_Uchimura_Humphreys_20,
  plot_on_this_UMAP=seurat_Uchimura_Humphreys_20,
  skip_neighbors=TRUE,
  ref_cell_type_column = "type",
  dims = 1:30,
  train_title = "SIX2GFP",
  test_title = "Uchimura",
  n_neighbors = 8,
  # colors_for_feature = c('grey', '#f03b20'),
  # myColors = myColors,
  output_prefix = "six2gfp/Uchimura_iter/full"
)

min_column_indices <- apply(results_Uchimura_six2gfp_full$test_query_neighbors, 1, which.min)

sum(unique(results_Uchimura_six2gfp_full$train[["type"]])$type[min_column_indices]==
      results_Uchimura_six2gfp_full$test$predicted.type)/length(results_Uchimura_six2gfp_full$test$predicted.type) # 90%

# results_Uchimura_six2gfp_full$test$predicted.type[results_Uchimura_six2gfp_full$test$predicted.type=="CM"]

# sum(min_column_indices[results_Uchimura_six2gfp_full$test$predicted.type=="CM"]==1)
# 
# colnames(results_Uchimura_six2gfp_full$test_query_neighbors)
# 
# colnames(results_Uchimura_six2gfp_full$test_query_neighbors)[min_column_indices[results_Uchimura_six2gfp_full$test$predicted.type=="CM"]]


FeaturePlot(AddMetaData(results_Uchimura_six2gfp_full$test,results_Uchimura_six2gfp_full$test_query_neighbors$CM ,col.name="CM_DISTANCE"),
            reduction = "ref.umap", features = "CM_DISTANCE",
                        keep.scale = "all")

FeaturePlot(AddMetaData(results_Uchimura_six2gfp_full$test,results_Uchimura_six2gfp_full$test_query_neighbors$UM ,col.name="UM_DISTANCE"),
            reduction = "ref.umap", features = "UM_DISTANCE",
            keep.scale = "all")
seurat_Uchimura_Humphreys_20= AddMetaData(seurat_Uchimura_Humphreys_20,results_Uchimura_six2gfp_full$test$predicted.type,col.name = "predicted.type")
DimPlot(seurat_Uchimura_Humphreys_20, reduction = "umap", label = TRUE, group.by = "predicted.type",
        label.size = 3, repel = FALSE)
seurat_Uchimura_Humphreys_20= AddMetaData(seurat_Uchimura_Humphreys_20,results_Uchimura_six2gfp_full$test$entropy,col.name = "entropy")
FeaturePlot(seurat_Uchimura_Humphreys_20, reduction = "umap", features = "entropy",
                        keep.scale = "all")



row.names(results_Uchimura_six2gfp_full$test@assays[["prediction.score.type"]]@data)
# results_Uchimura_six2gfp_full$test@assays[["prediction.score.type"]]@data["CM",]
FeaturePlot(AddMetaData(seurat_Uchimura_Humphreys_20,results_Uchimura_six2gfp_full$test@assays[["prediction.score.type"]]@data["CM",] ,col.name="CM_prediction"),
            reduction = "umap", features = "CM_prediction",
            keep.scale = "all")

FeaturePlot(AddMetaData(seurat_Uchimura_Humphreys_20,results_Uchimura_six2gfp_full$test@assays[["prediction.score.type"]]@data["UM",] ,col.name="UM_prediction"),
            reduction = "umap", features = "UM_prediction",
            keep.scale = "all")

# plot prediction on original view
for (t in row.names(results_Uchimura_six2gfp_full$test@assays[["prediction.score.type"]]@data)){
  FeaturePlot(AddMetaData(seurat_Uchimura_Humphreys_20,results_Uchimura_six2gfp_full$test@assays[["prediction.score.type"]]@data[t,] 
                          ,col.name=paste0(t,"_prediction")),
              reduction = "umap", features = paste0(t,"_prediction"),
              keep.scale = "all")
  ggsave(paste0("six2gfp/Uchimura_iter/full",t, "_prediction_original_view.png"), 
         width = 1600, height = 900, units = "px", dpi = 100)
}



DimPlot(seurat_Uchimura_Humphreys_20, reduction = "umap", label = FALSE, 
        cells = names(results_Uchimura_six2gfp_full$test$predicted.type[results_Uchimura_six2gfp_full$test$predicted.type=="CM"]),
        label.size = 2, repel = FALSE)
DimPlot(seurat_Uchimura_Humphreys_20, reduction = "umap", label = FALSE, 
        cells = names(results_Uchimura_six2gfp_full$test$predicted.type[results_Uchimura_six2gfp_full$test$predicted.type=="UM"]),
        label.size = 2, repel = FALSE)

ggsave(paste0(output_prefix, "_original_view_vs_iter.png"), 
       width = 1600, height = 900, units = "px", dpi = 100)


plot(results_Uchimura_six2gfp_full[["test"]]@assays[["prediction.score.type"]]@data["UM",],
     results_Uchimura_six2gfp_full[["test"]]@assays[["prediction.score.type"]]@data["CM",])

# FeaturePlot(results_Uchimura_six2gfp_full$test, reduction = "ref.umap", features = "TMEM52B",
#             keep.scale = "all") # LOH
# FeaturePlot(results_Uchimura_six2gfp_full$test, reduction = "ref.umap", features = "PODXL",
#             keep.scale = "all") # PODO
# FeaturePlot(results_Uchimura_six2gfp_full$test, reduction = "ref.umap", features = "GATA3",
#             keep.scale = "all") # Collecting duct -dist cd
# FeaturePlot(results_Uchimura_six2gfp_full$test, reduction = "ref.umap", features = "GAP43",
#             keep.scale = "all") # neron (empty)
# FeaturePlot(results_Uchimura_six2gfp_full$test, reduction = "ref.umap", features = "GATM",
#             keep.scale = "all") # PT
# FeaturePlot(results_Uchimura_six2gfp_full$test, reduction = "ref.umap", features = "LHX1",
#             keep.scale = "all") # Prox1
# FeaturePlot(results_Uchimura_six2gfp_full$test, reduction = "ref.umap", features = "COL1A1",
#             keep.scale = "all") # Mes prog
# FeaturePlot(results_Uchimura_six2gfp_full$test, reduction = "ref.umap", features = "MKI67",
#             keep.scale = "all") # Mes
# FeaturePlot(results_Uchimura_six2gfp_full$test, reduction = "ref.umap", features = "SIX2",
#             keep.scale = "all") # should divine between CM/UM, but doesn't
# FeaturePlot(results_Uchimura_six2gfp_full$test, reduction = "ref.umap", features = "CRYM",
#             keep.scale = "all") # 
# FeaturePlot(results_Uchimura_six2gfp_full$train, reduction = "umap", features = "COL1A1",
#             keep.scale = "all") # 
# 
# FeaturePlot(results_Uchimura_six2gfp_full$test, reduction = "ref.umap", features = "entropy",
#             keep.scale = "all") + ggtitle("original")
# 
# 
#   # Split to 3
# Uchimura_Humphreys_20 = read.table("GSM3763147_Uchimura.dge.txt")
# {set.seed(1)
# sample_vector <- sample(c(1, 2, 3), ncol(Uchimura_Humphreys_20), replace=TRUE, prob=c(0.33,0.33,0.34))
# third1  <- colnames(Uchimura_Humphreys_20)[sample_vector==1]
# third2  <- colnames(Uchimura_Humphreys_20)[sample_vector==2]
# third3  <- colnames(Uchimura_Humphreys_20)[sample_vector==3]
# remove(sample_vector)
# 
# Uchimura_third1 <- Uchimura_Humphreys_20[,third1]
# Uchimura_third2 <- Uchimura_Humphreys_20[,third2]
# Uchimura_third3 <- Uchimura_Humphreys_20[,third3]
# 
# seurat_Uchimura_Humphreys_20_third1 <- CreateSeuratObject(counts = Uchimura_third1, assay = "RNA")
# seurat_Uchimura_Humphreys_20_third2 <- CreateSeuratObject(counts = Uchimura_third2, assay = "RNA")
# seurat_Uchimura_Humphreys_20_third3 <- CreateSeuratObject(counts = Uchimura_third3, assay = "RNA")
# remove(third1,third2,third3,Uchimura_third1,Uchimura_third2,Uchimura_third3,Uchimura_Humphreys_20)
# }

