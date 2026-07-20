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
library(grid)

# ============================================================================
# Override ggsave globally to force the "real" figure panels (UMAP, heatmaps)
# to a fixed inches size. This ensures the panels are perfectly identical 
# in size across all SVGs for easy alignment in graphic programs, while the 
# overall file dimensions automatically adjust to accommodate varying text/legends.
# ============================================================================
ggsave <- function(filename, plot = last_plot(), width = NA, height = NA, units = c("in", "cm", "mm", "px"), ...) {
  # Do not apply the single-panel sizing to patchwork objects as it breaks their internal grid layout
  if (inherits(plot, "ggplot") && !inherits(plot, "patchwork")) {
    g <- ggplot2::ggplotGrob(plot)
    panels <- grep("panel", g$layout$name)
    if (length(panels) > 0) {
      panel_index_w <- unique(g$layout$l[panels])
      panel_index_h <- unique(g$layout$t[panels])
      
      # Force ZxY inches per panel
      g$widths[panel_index_w] <- rep(unit(8, "in"), length(panel_index_w))
      g$heights[panel_index_h] <- rep(unit(8, "in"), length(panel_index_h))
      
      # Calculate new bounding box to fit the fixed panels + all legends/text
      total_width <- grid::convertWidth(sum(g$widths), unitTo = "in", valueOnly = TRUE)
      total_height <- grid::convertHeight(sum(g$heights), unitTo = "in", valueOnly = TRUE)
      
      ggplot2::ggsave(filename, plot = g, width = total_width, height = total_height, units = "in", ...)
      return(invisible(NULL))
    }
  }
  
  # Fallback for other objects
  ggplot2::ggsave(filename, plot = plot, width = width, height = height, units = units, ...)
}

# ============================================================================
# Setup Roboto Font for all plots
# ============================================================================
font_dir <- file.path(getwd(), "Roboto_fonts")
if (requireNamespace("systemfonts", quietly = TRUE)) {
  systemfonts::register_font(
    name = "Roboto",
    plain = file.path(font_dir, "static", "Roboto-Regular.ttf"),
    bold = file.path(font_dir, "static", "Roboto-Bold.ttf"),
    italic = file.path(font_dir, "static", "Roboto-Italic.ttf"),
    bolditalic = file.path(font_dir, "static", "Roboto-BoldItalic.ttf")
  )
  systemfonts::register_font(
    name = "Roboto Medium",
    plain = file.path(font_dir, "static", "Roboto-Medium.ttf"),
    bold = file.path(font_dir, "static", "Roboto-Bold.ttf"),
    italic = file.path(font_dir, "static", "Roboto-MediumItalic.ttf"),
    bolditalic = file.path(font_dir, "static", "Roboto-BoldItalic.ttf")
  )
} else if (requireNamespace("sysfonts", quietly = TRUE)) {
  sysfonts::font_add(
    family = "Roboto",
    regular = file.path(font_dir, "static", "Roboto-Regular.ttf"),
    bold = file.path(font_dir, "static", "Roboto-Bold.ttf"),
    italic = file.path(font_dir, "static", "Roboto-Italic.ttf"),
    bolditalic = file.path(font_dir, "static", "Roboto-BoldItalic.ttf")
  )
  sysfonts::font_add(
    family = "Roboto Medium",
    regular = file.path(font_dir, "static", "Roboto-Medium.ttf"),
    bold = file.path(font_dir, "static", "Roboto-Bold.ttf"),
    italic = file.path(font_dir, "static", "Roboto-MediumItalic.ttf"),
    bolditalic = file.path(font_dir, "static", "Roboto-BoldItalic.ttf")
  )
  if (requireNamespace("showtext", quietly = TRUE)) showtext::showtext_auto()
}

# Apply Roboto universally to all ggplot and Seurat objects
ggplot2::theme_set(
  ggplot2::theme_get() + 
  ggplot2::theme(
    text = ggplot2::element_text(family = "Roboto Medium"),
    axis.text = ggplot2::element_text(family = "Roboto Medium"),
    legend.text = ggplot2::element_text(family = "Roboto Medium"),
    title = ggplot2::element_text(family = "Roboto", face = "bold", size = 20),
    plot.title = ggplot2::element_text(family = "Roboto", face = "bold", size = 20),
    axis.title = ggplot2::element_text(family = "Roboto", size = 20),
    legend.title = ggplot2::element_text(family = "Roboto", face = "bold", size = 20)
  )
)

# Override theme_minimal so its default base_family doesn't reset Roboto
default_theme_minimal <- ggplot2::theme_minimal
theme_minimal <- function(base_size = 12, base_family = "Roboto Medium", ...) {
  default_theme_minimal(base_size = base_size, base_family = base_family, ...) +
    ggplot2::theme(
      text = ggplot2::element_text(family = "Roboto Medium"),
      axis.text = ggplot2::element_text(family = "Roboto Medium"),
      legend.text = ggplot2::element_text(family = "Roboto Medium"),
      title = ggplot2::element_text(family = "Roboto", face = "bold", size = 20),
      plot.title = ggplot2::element_text(family = "Roboto", face = "bold", size = 20),
      axis.title = ggplot2::element_text(family = "Roboto", size = 20),
      legend.title = ggplot2::element_text(family = "Roboto", face = "bold", size = 20)
    )
}

# Update default fonts for geom_text and ggrepel
ggplot2::update_geom_defaults("text", list(family = "Roboto Medium"))
ggplot2::update_geom_defaults("label", list(family = "Roboto Medium"))
if (requireNamespace("ggrepel", quietly = TRUE)) {
  ggplot2::update_geom_defaults("text_repel", list(family = "Roboto Medium"))
  ggplot2::update_geom_defaults("label_repel", list(family = "Roboto Medium"))
}

# Fixed color palette: each cell type always gets the same color
CELL_TYPE_COLORS <- setNames(
  scales::hue_pal()(10),
  sort(c("CM", "CM_DIV", "DIST_CD", "ENDO", "LOH", "MACROPHAG", "PODO", "PROX_1", "PROX_2", "UM"))
)
# Resulting assignment (alphabetical → hue order):
#   CM="#F8766D"  CM_DIV="#E58700"  DIST_CD="#C99800"  ENDO="#A3A500"  LOH="#6BB100"
#   MACROPHAG="#00BA38"  PODO="#00C19F"  PROX_1="#00B8E7"  PROX_2="#619CFF"  UM="#FF61C3"

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
    return_anchors = FALSE,
    noise_model = "poisson", # options: "poisson", "dropout", "negative_binomial"
    dropout_rate = 0.1,      # for "dropout" model
    background_rate = 0.001, # for "dropout" model (fraction of matrix to add background)
    background_mean = 1,     # for "dropout" model (mean of background counts)
    nb_theta = 10            # for "negative_binomial" model (overdispersion parameter)
) {
  # 0. Create output directory if it doesn't exist ####
  if (!dir.exists(output_prefix_base)) {
    dir.create(output_prefix_base, recursive = TRUE)
  }
  
  # 1. Add noise function ####
  add_noise <- function(seurat_obj) {
    counts <- GetAssayData(seurat_obj, layer = "counts")
    
    if (noise_model == "poisson") {
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
    } else if (noise_model == "negative_binomial") {
      # Simulate biological overdispersion using Negative Binomial
      noisy_x_values <- rnbinom(n = length(counts@x), mu = counts@x, size = nb_theta)
      noisy_counts_sparse <- Matrix::sparseMatrix(
        i = counts@i,
        p = counts@p,
        x = noisy_x_values,
        dims = dim(counts),
        dimnames = dimnames(counts),
        index1 = FALSE
      )
    } else if (noise_model == "dropout") {
      # 1. Apply poisson noise to existing non-zeros
      noisy_x <- rpois(n = length(counts@x), lambda = counts@x)
      # 2. Apply dropout (set a fraction to 0)
      drop_mask <- runif(length(noisy_x)) < dropout_rate
      noisy_x[drop_mask] <- 0
      
      # 3. Add background contamination (drop-in)
      nrow_counts <- nrow(counts)
      ncol_counts <- ncol(counts)
      total_elements <- as.numeric(nrow_counts) * as.numeric(ncol_counts)
      # Approximate binomial with poisson for efficiency
      n_bg <- rpois(1, lambda = total_elements * background_rate)
      
      base_mat <- Matrix::sparseMatrix(
        i = counts@i,
        p = counts@p,
        x = noisy_x,
        dims = c(nrow_counts, ncol_counts),
        dimnames = dimnames(counts),
        index1 = FALSE
      )
      
      if (n_bg > 0) {
        # Generate random indices for background
        bg_i <- sample(1:nrow_counts, n_bg, replace = TRUE)
        bg_j <- sample(1:ncol_counts, n_bg, replace = TRUE)
        bg_x <- rpois(n_bg, lambda = background_mean)
        
        # Keep only > 0 values to maintain sparsity
        keep <- bg_x > 0
        if (sum(keep) > 0) {
          bg_mat <- Matrix::sparseMatrix(
            i = bg_i[keep],
            j = bg_j[keep],
            x = bg_x[keep],
            dims = c(nrow_counts, ncol_counts),
            dimnames = dimnames(counts)
          )
          noisy_counts_sparse <- base_mat + bg_mat
        } else {
          noisy_counts_sparse <- base_mat
        }
      } else {
        noisy_counts_sparse <- base_mat
      }
    } else {
      stop(paste("Unknown noise_model:", noise_model))
    }
    
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
    six2gfp_levels <- c("UM","CM","CM_DIV","PODO","PROX_1","PROX_2","LOH","DIST_CD","MACROPHAG","ENDO")
    if (any(cell_type_order %in% six2gfp_levels)) {
      # SIX2GFP reference: use hardcoded order
      cell_type_order = factor(cell_type_order, levels = six2gfp_levels)
      cell_type_order = cell_type_order[order(match(cell_type_order, levels(cell_type_order)))]
    } else {
      # Non-SIX2GFP reference: check for human cell atlas types
      atlas_levels <- c(
        # Nephron Lineage
        "Cap mesenchyme", "Proliferating cap mesenchyme",
        "Proximal renal vesicle", "Distal renal vesicle", "Proliferating distal renal vesicle",
        "Proximal S shaped body", "Medial S shaped body", "Distal S shaped body",
        "Podocyte", "Proximal tubule", "Loop of Henle",
        
        # Ureteric Bud / Collecting Duct
        "Proximal UB", "CNT/PC - proximal UB", "Pelvic epithelium - distal UB",
        
        # Stroma / Fibroblasts
        "Stroma progenitor", "Proliferating stroma progenitor",
        "Fibroblast 1", "Fibroblast 2",
        "Myofibroblast 1", "Myofibroblast 2", "Proliferating myofibroblast",
        
        # Endothelial
        "Endothelium",
        
        # Immune - Myeloid
        "Macrophage 1", "Macrophage 2", "Proliferating macrophage",
        "Monocyte", "Proliferating monocyte", "Mast cells",
        "cDC1", "cDC2", "pDC", "Neutrophil",
        
        # Immune - Lymphoid
        "B cell", "Proliferating B cell", "CD4 T cell", "CD8 T cell", "NK cell", "Innate like lymphocyte",
        
        # Other
        "Erythroid", "Megakaryocyte", "Neuron"
      )
      if (any(cell_type_order %in% atlas_levels)) {
        # Combine found atlas levels in the specified order with any other unknown cell types at the end
        found_levels <- atlas_levels[atlas_levels %in% cell_type_order]
        other_levels <- sort(setdiff(cell_type_order, atlas_levels))
        cell_type_order = factor(cell_type_order, levels = c(found_levels, other_levels))
        cell_type_order = cell_type_order[order(match(cell_type_order, levels(cell_type_order)))]
        cell_type_order = cell_type_order[!is.na(cell_type_order)]
      } else {
        # Unknown reference: fall back to sorted order
        cell_type_order = sort(unique(cell_type_order))
      }
    }
    if (is.null(myColors)) {
      if (any(cell_type_order %in% six2gfp_levels)) {
        # Use the fixed global palette so colors are consistent across plots
        # with different numbers of groups (subset of the full 10-type palette).
        myColors <- CELL_TYPE_COLORS
      } else {
        # Non-SIX2GFP reference: generate a palette from Paired or hue_pal
        n_types <- length(cell_type_order)
        if (n_types <= 12) {
          myColors <- setNames(brewer.pal(max(n_types, 3), "Paired")[1:n_types], cell_type_order)
        } else {
          myColors <- setNames(scales::hue_pal()(n_types), cell_type_order)
        }
      }
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
      pairwise_confusion <- pairwise_confusion[gsub("-", "_", cell_type_order), gsub("-", "_", cell_type_order)]
      
      for (i in 1:num_classes) {
        for (j in 1:num_classes) {
          if (i != j) {
            pairwise_confusion[i, j] <- mean(probs_tub[i, ] * probs_tub[j, ])
          }
        }
      }
      
      pairwise_confusion_df <- reshape2::melt(pairwise_confusion)
      pairwise_confusion_df <- pairwise_confusion_df %>%
        mutate(
          label_val = round(value, 5),
          label_str = ifelse(label_val == 0, "", as.character(label_val))
        )
      p <- ggplot(data = pairwise_confusion_df, aes(Var1, Var2, fill = value)) +
        geom_tile() +
        geom_text(aes(label = label_str)) +
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
      
      pred_col <- paste0("predicted.", ref_cell_type_column)
      current_cells=names(test_query@meta.data[[pred_col]])
      
      if (pred_col %in% colnames(plot_on_this_UMAP@meta.data)) {
        if(!is.null(plot_on_this_UMAP@reductions$ref.umap))
        { reduction='ref.umap' } else { reduction='umap' }
        
        p1 = DimPlot(plot_on_this_UMAP, cells=current_cells, reduction = reduction,
                     group.by = paste0("predicted.",ref_cell_type_column), label = TRUE,
                     label.size = 3, repel = FALSE) + ggtitle("original") +
          scale_color_manual(values = myColors)
        
        # change the prediction according to this run
        plot_on_this_UMAP[[prediction_column_name]][current_cells]=test_query[[paste0("predicted.",ref_cell_type_column)]]
        
        p2 = DimPlot(plot_on_this_UMAP, cells=current_cells, reduction = reduction,
                     group.by = paste0("predicted.",ref_cell_type_column), label = TRUE,
                     label.size = 3, repel = FALSE) + ggtitle("new iteration") +
          scale_color_manual(values = myColors)
        p2_alone = DimPlot(plot_on_this_UMAP, cells=current_cells, reduction = reduction,
                group.by = paste0("predicted.",ref_cell_type_column), label = TRUE,
                label.size = 5, repel = FALSE) + theme(legend.position = "none") +
          scale_color_manual(values = myColors)
        
        combined_original_view_plot <- p1+p2
        ggsave(paste0(output_prefix, "_original_view_vs_iter.svg"), plot = combined_original_view_plot,
               width = 16, height = 9, units = "in")
        
        ggsave(paste0(output_prefix, "_original_view_only_prediction_A.svg"), plot = p2_alone,
               width = 8, height = 8, units = "in")
               
        p2_alone_with_legend = DimPlot(plot_on_this_UMAP, cells=current_cells, reduction = reduction,
                group.by = paste0("predicted.",ref_cell_type_column), label = TRUE,
                label.size = 5, repel = FALSE) +
          scale_color_manual(values = myColors)
        ggsave(paste0(output_prefix, "_original_view_only_prediction_A_with_legend.svg"), plot = p2_alone_with_legend,
               width = 10, height = 8, units = "in")
        
        
      }
      else {
        p1 = DimPlot(plot_on_this_UMAP, reduction = "umap", label = TRUE,
                     label.size = 3, repel = FALSE) +
          scale_color_manual(values = myColors)
        
        
        plot_on_this_UMAP=AddMetaData(plot_on_this_UMAP,test_query@meta.data[[pred_col]],col.name = pred_col)
        p2 = DimPlot(plot_on_this_UMAP, reduction = "umap", label = TRUE, group.by = pred_col,
                     label.size = 8, repel = FALSE) +
          scale_color_manual(values = myColors)
        
        combined_original_view_plot <- p1+p2
        ggsave(paste0(output_prefix, "_original_view_vs_iter.svg"), plot = combined_original_view_plot,
               width = 16, height = 9, units = "in")
        
        p2_alone = DimPlot(plot_on_this_UMAP, reduction = "umap", label = TRUE, group.by = pred_col,
                           label.size = 5, repel = FALSE) + theme(legend.position = "none") +
          scale_color_manual(values = myColors)
        
        ggsave(paste0(output_prefix, "_original_view_only_prediction_A.svg"), plot = p2_alone,
               width = 8, height = 8, units = "in")
               
        p2_alone_with_legend = DimPlot(plot_on_this_UMAP, reduction = "umap", label = TRUE, group.by = pred_col,
                           label.size = 5, repel = FALSE) +
          scale_color_manual(values = myColors)
        ggsave(paste0(output_prefix, "_original_view_only_prediction_A_with_legend.svg"), plot = p2_alone_with_legend,
               width = 10, height = 8, units = "in")
        
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
  main_cell_type_order <- gsub("-", "_", main_cell_type_order) # Normalize names to match prediction score rownames
  
  # 3. Run noise noised_number times and keep relevant data ####
  noised_prediction = matrix('0', nrow = ncol(seurat_obj), ncol = noised_number+1)
  seurat_noised_prediction_list <- list()
  noised_anchors_list <- list()
  rownames(noised_prediction) = colnames(seurat_obj)
  noised_prediction[,1] = temp_seurat_obj[["test"]]@meta.data[[prediction_column_name]]
  for( k in c(1:noised_number)) {
    print(paste0("running run #",k))
    seurat_obj_noised = add_noise(seurat_obj)
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
  # Normalize all prediction values to match gsub'd names used in JSD/PSS matrices
  noised_prediction[,] <- gsub("-", "_", noised_prediction)
  
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
  plot_change_feature(seurat_obj, change_counts_df, colors_feature_plot_noise, output_prefix_base)
  
  # 5.5 and plot by prediction ####
  plot_noised_predictions_original_view(seurat_obj, temp_seurat_obj, results_noised, output_prefix_base)
  
  # 6. df that show who changed to what (proportions) ####
  plot_change_proportion(noised_prediction, main_cell_type_order, output_prefix_base)
  
  
  # 7. Confusion matrix ####
  conf_matrix <- table(Original = noised_prediction[,1],
                       New = noised_prediction[,ncol(noised_prediction)])
  
  # Handling Small Groups - Filter ROWS by original count, keep ALL columns that appear in New
  types_counts_conf_matrix <- rowSums(conf_matrix)
  rows_to_keep <- names(types_counts_conf_matrix)[types_counts_conf_matrix >= min_cell_count_for_type]
  
  # Identify removed types
  removed_types_present <- names(types_counts_conf_matrix)[types_counts_conf_matrix < min_cell_count_for_type]
  all_known_types <- levels(main_cell_type_order)
  if (is.null(all_known_types)) all_known_types <- as.character(main_cell_type_order)
  absent_types <- setdiff(all_known_types, names(types_counts_conf_matrix))
  removed_types <- unique(c(removed_types_present, absent_types))
  
  cols_to_keep <- colnames(conf_matrix)  # keep all types appearing in noised predictions
  # Only keep columns that also exist in rows (for consistent ordering), plus any new ones
  cols_to_keep <- union(rows_to_keep, colnames(conf_matrix)[colSums(conf_matrix[rows_to_keep, , drop=FALSE]) > 0])
  
  # Exclude removed types from columns entirely so they don't appear as empty columns
  cols_to_keep <- setdiff(cols_to_keep, removed_types)
  
  conf_matrix_filtered <- conf_matrix[rows_to_keep, cols_to_keep, drop=FALSE]
  
  conf_matrix_fraction <- prop.table(conf_matrix_filtered, margin = 1) # * 100
  
  plot_confusion_matrix_heatmap(conf_matrix_fraction, all_known_types, removed_types, min_cell_count_for_type, output_prefix_base)
  plot_confusion_matrix_raw_counts(conf_matrix_filtered, all_known_types, removed_types, min_cell_count_for_type, output_prefix_base)
  
  # 9. JSD ####
  calculate_entropy_jsd <- function(query) {
    log_pred <- log2(query)
    log_pred[is.infinite(log_pred)] <- 0
    entropy_by_cell <- -colSums(query * log_pred)
    return(entropy_by_cell)
  }
  
  prediction_data = temp_seurat_obj[["test"]]@assays[[paste0("prediction.score.",ref_cell_type_column)]]@data # Use original prediction
  rownames(prediction_data)=gsub("-", "_", rownames(prediction_data))
  # P^R, i
  norm_p_factor = rowSums(prediction_data)
  norm_prediction_table = sweep(prediction_data, 1, norm_p_factor, FUN = "/")
  # P^C, j
  # prediction_data = results_noised[["test"]]@assays[["prediction.score.type"]]@data # Use the last noised result
  # rownames(prediction_data)=gsub("-", "_", rownames(prediction_data))
  predicted_zeroes_mat <- (sapply(temp_seurat_obj[["test"]]@meta.data[[prediction_column_name]], function(x) rownames(prediction_data) == gsub("-", "_", x)) * 1)
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
  
  plot_jsd_heatmap(JSD_mat, main_cell_type_order, output_prefix_base)
  
  # filter out small groups
  types_to_run_on = intersect(colnames(RSS_mat_filtered), rownames(conf_matrix_fraction)) # use rownames from filtered conf_matrix
  types_to_run_on = intersect(types_to_run_on, main_cell_type_order) # ensure order is consisten
  
  plot_pss_heatmap(RSS_mat_filtered, types_to_run_on, main_cell_type_order, output_prefix_base)
  
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
    
    # Bi-directional stability: penalizes both outflow AND inflow
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
        # Bi-directional stability per run
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
  
  plot_stability_vs_pss_scatters(stability_pred, output_prefix_base)
  
  # plot also the full table
  # conf_matrix_fraction_ordered = conf_matrix_fraction[types_to_run_on,types_to_run_on, drop=FALSE]
  # RSS_mat_filtered_ordered = RSS_mat_filtered[types_to_run_on,types_to_run_on, drop=FALSE]
  # conf_matrix_long <- as.data.frame(as.table(conf_matrix_fraction_ordered))
  # RSS_long <- as.data.frame(as.table(RSS_mat_filtered_ordered)) 
  # names(RSS_long) <- c("P_R", "P_C", "PSS")
  # 
  # rss_and_stab = data.frame(conf_matrix_long,RSS_long[,3])
  # colnames(rss_and_stab) = c("Original","New","Freq","PSS")
  
  plot_and_save_stability_scatters(conf_matrix_fraction, RSS_mat_filtered, stability_pred, types_to_run_on, main_cell_type_order, output_prefix_base)
  
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
    noised_anchors = noised_anchors_list,
    main_cell_type_order = main_cell_type_order
  ))
}

plot_and_save_stability_scatters <- function(conf_matrix_fraction, RSS_mat_filtered, stability_pred, types_to_run_on, main_cell_type_order, output_prefix_base) {
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
  
  # Match Original column to stability_pred rownames
  match_idx <- match(as.character(rss_and_stab$Original), rownames(stability_pred))
  rss_and_stab$Bi_Stability <- ifelse(rss_and_stab$is_diagonal, stability_pred$Bi_Stability[match_idx], rss_and_stab$Freq)
  rss_and_stab$F1_Stability <- ifelse(rss_and_stab$is_diagonal, stability_pred$F1_Stability[match_idx], rss_and_stab$Freq)

  # Consistent Order in Heatmaps - Order full stability plot
  # If main_cell_type_order is factor, we extract levels, else we use it as levels directly if it's char array
  lvls <- if(is.factor(main_cell_type_order)) levels(main_cell_type_order) else main_cell_type_order
  
  rss_and_stab$Original <- factor(rss_and_stab$Original, levels = lvls)
  rss_and_stab$New <- factor(rss_and_stab$New, levels = lvls)
  rss_and_stab <- rss_and_stab %>% drop_na()

  # Define function to generate the 3 plots per metric
  generate_scatter_plots <- function(metric_col, output_suffix, title_prefix, df_all) {
    df_diag <- df_all %>% filter(is_diagonal == TRUE)
    
    fit <- tryCatch({
      lm(as.formula(paste(metric_col, "~ PSS")), data = df_diag)
    }, error = function(e) return(NULL))
    
    if (is.null(fit)) return()
    
    fit_summary <- summary(fit)
    intercept <- coef(fit)[1]
    slope <- coef(fit)[2]
    r_squared <- fit_summary$r.squared
    
    # Safely extract p-value
    p_value_pss <- NA
    if ("PSS" %in% rownames(fit_summary$coefficients)) {
      p_value_pss <- fit_summary$coefficients["PSS","Pr(>|t|)"]
    }
    if(is.na(p_value_pss) || is.null(p_value_pss)){
       p_value_pss <- 1
    }

    equation_string <- sprintf("y = %.2f x + %.2f", slope, intercept)
    r_squared_string <- sprintf("R² = %.3f ; Pval = %.3f", r_squared, p_value_pss)

    spearman_test_diagonal <- cor.test(x = df_diag$PSS, 
                                       y = df_diag[[metric_col]], 
                                       method = "spearman",
                                       exact = FALSE)
    spearman_corr_diagonal <- spearman_test_diagonal$estimate
    spearman_pval_diagonal <- spearman_test_diagonal$p.value
    
    spearman_corr_all <- cor.test(x = df_all$PSS, y = df_all[[metric_col]], method = "spearman", exact = FALSE)
    
    plot_all <- ggplot(df_all, aes(x = PSS, y = .data[[metric_col]])) +
      geom_point(aes(color = is_diagonal), size = 3.5, alpha = 0.8) + theme_minimal() + 
      scale_color_manual(
        name = "Transition Type", 
        values = c("TRUE" = "red", "FALSE" = "black"),
        labels = c("TRUE" = paste("Self-Transition (", title_prefix, ")", sep=""), "FALSE" = "Cross-Transition")
      ) +
      geom_text(label = paste0(df_all$Original,"-", df_all$New), vjust = 1.5, size = 12/.pt) + 
      annotate("text", x = Inf, y = -Inf,
               label = paste("Spearman Coefficient:", round(spearman_corr_diagonal, 3), "; pVal:", round(spearman_pval_diagonal, 3)),
               hjust = 1.05, vjust = -1.5, size = 9) +
      annotate("text", x = Inf, y = -Inf,
               label = paste("Spearman all:", round(spearman_corr_all$estimate, 3), "; pVal:", round(spearman_corr_all$p.value, 3)),
               hjust = 1.05, vjust = -0.5, size = 9) +
      xlab("Prediction specificity score") + 
      ylab(paste("Fraction of transitions /", title_prefix)) +
      scale_x_continuous(limits = c(0, 1)) +
      scale_y_continuous(limits = c(0, 1)) +
      theme(legend.position = "none",
            axis.text.x = element_text(size = 13, angle = 45),
            axis.text.y = element_text(size = 13, angle = 0, hjust = 0.5),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    
    ggsave(paste0(output_prefix_base, "all_freq_vs_pss_scatter_plot", output_suffix, ".svg"), plot = plot_all,
           width = 12, height = 8, units = "in")
           
    plot_all_with_fit <- plot_all + 
      geom_abline(intercept = intercept, slope = slope, color = "blue", linetype = "dashed", linewidth = 1) +
      annotate("text", x = Inf, y = Inf,
               label = paste(equation_string, r_squared_string, sep = "\n"),
               hjust = 1.05, vjust = 1.2,
               size = 4, color = "blue")
    
    ggsave(paste0(output_prefix_base, "all_freq_vs_pss_scatter_plot", output_suffix, "_with_fit.svg"), 
           plot = plot_all_with_fit,
           width = 12, height = 8, units = "in")
           
    stability_only_plot <- ggplot(df_diag, aes(x = PSS, y = .data[[metric_col]])) +
      geom_point(color = "red", size = 4, alpha = 0.7) +
      geom_abline(intercept = intercept, slope = slope, color = "blue", linetype = "dashed", linewidth = 1) +
      annotate("text", 
               x = 0.02,
               y = 0.98, 
               label = paste(equation_string, r_squared_string, sep = "\n"),
               hjust = 0, vjust = 1,
               size = 5,
               parse = FALSE) +
      geom_text_repel(aes(label = Original), size = 3.5, box.padding = 0.5) +
      labs(
        title = paste(title_prefix, "vs. Prediction Specificity Score (PSS)"),
        subtitle = "Analysis of self-transitions (Original = New Prediction)",
        x = "Prediction Specificity Score (PSS)",
        y = paste(title_prefix, "(% of cells remaining same type after noise)")
      ) +
      scale_x_continuous(limits = c(0, 1)) +
      scale_y_continuous(limits = c(0, 1)) +
      theme_minimal()
    ggsave(paste0(output_prefix_base, "stability_pss_linear_fit", output_suffix, ".svg"), plot = stability_only_plot,
           width = 8, height = 7, units = "in")
  }
  
  generate_scatter_plots("Freq", "", "Stability", rss_and_stab)
  generate_scatter_plots("Bi_Stability", "_bi", "Bi-directional Stability", rss_and_stab)
  generate_scatter_plots("F1_Stability", "_f1", "F1_Stability", rss_and_stab)
  
  
  # 11. set chosen plots in place ####
  source_dir = output_prefix_base
  target_dir = paste0(source_dir,"/main_plots")
  if (!dir.exists(target_dir)) {
    dir.create(target_dir, recursive = TRUE)
  }
  files_to_move = c("0_silent__original_view_only_prediction_A.svg",
                    "all_freq_vs_pss_scatter_plot_with_fit.svg",
                    "all_freq_vs_pss_scatter_plot_bi_with_fit.svg",
                    "all_freq_vs_pss_scatter_plot_f1_with_fit.svg",
                    "change_proportion_plot.svg",
                    "confusion_matrix_raw_counts.svg",
                    "confusion_matrix_raw_counts_sum.svg",
                    "confusion_matrix.svg",
                    "pss_heatmap.svg",
                    "bi_stability_vs_pss_scatter_plot.svg",
                    "f1_stability_vs_pss_scatter_plot.svg",
                    "stability_pss_linear_fit.svg",
                    "stability_pss_linear_fit_bi.svg",
                    "stability_pss_linear_fit_f1.svg")
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
}

# ============================================================================
# Shared helper functions for plotting (used by both main flow and regenerate)
# ============================================================================

plot_confusion_matrix_heatmap <- function(conf_matrix_fraction, all_known_types, removed_types, min_cell_count_for_type, output_prefix_base) {
  conf_matrix_long <- as.data.frame(as.table(conf_matrix_fraction))
  names(conf_matrix_long) <- c("Original", "New", "Fraction")
  
  kept_levels <- setdiff(all_known_types, removed_types)
  conf_matrix_long$Original <- factor(conf_matrix_long$Original, levels = kept_levels)
  conf_matrix_long$New <- factor(conf_matrix_long$New, levels = kept_levels)
  conf_matrix_long <- conf_matrix_long %>% drop_na()
  conf_matrix_long <- conf_matrix_long %>%
    mutate(
      Fraction_label = sprintf("%.2f", Fraction),
      Fraction_label = ifelse(Fraction_label == "-0.00", "0.00", Fraction_label),
      Fraction_label = ifelse(Fraction_label == "0.00", "0", Fraction_label)
    )
  
  caption_text <- if (length(removed_types) > 0) {
    paste0("Excluded due to 0 or insufficient initial cells (<", min_cell_count_for_type, "): ", paste(removed_types, collapse = ", "))
  } else {
    NULL
  }
  
  confusion_matrix_plot <- ggplot(conf_matrix_long, aes(x = New, y = Original, fill = Fraction)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
             color = "black", fill = NA, linewidth = 1) +
    geom_tile(color = "white") +
    geom_text(aes(label = Fraction_label),
              size = 7) +
    scale_fill_gradient(low = "white",
                        high = "blue",
                        name = "Fraction") +
    theme_minimal() +
    labs(title = "Transitions with noise", caption = caption_text) +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 14, angle = 0, hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.caption = element_text(hjust = 0, size = 10, color = "gray30"))
  
  ggsave(paste0(output_prefix_base, "confusion_matrix.svg"), plot = confusion_matrix_plot,
         width = 9, height = 9, units = "in")
}

plot_confusion_matrix_raw_counts <- function(conf_matrix_filtered, all_known_types, removed_types, min_cell_count_for_type, output_prefix_base) {
  # Convert raw counts matrix to long format
  conf_matrix_long <- as.data.frame(as.table(conf_matrix_filtered))
  names(conf_matrix_long) <- c("Original", "New", "Count")
  
  # Normalize counts
  row_sums <- rowSums(conf_matrix_filtered)
  col_sums <- colSums(conf_matrix_filtered)
  
  conf_matrix_long$NormalizedSum <- apply(conf_matrix_long, 1, function(row) {
    orig <- as.character(row["Original"])
    new_type <- as.character(row["New"])
    count <- as.numeric(row["Count"])
    
    r_sum <- row_sums[orig]
    
    if (orig == new_type) {
      # Diagonal cell: (out_count + in_count) / original_count
      stayed <- count
      original_count <- r_sum
      noised_count <- if (new_type %in% names(col_sums)) col_sums[new_type] else 0
      out_count <- original_count - stayed
      in_count <- noised_count - stayed
      
      if (is.na(original_count) || original_count == 0) {
        0
      } else {
        max(0, 1 - (out_count + in_count) / original_count)
      }
    } else {
      # Non-diagonal cell: [cell/(original col row sum)] + [cell/(new col row sum)]
      c_sum <- row_sums[new_type]
      
      term1 <- if (is.na(r_sum) || r_sum == 0) 0 else count / r_sum
      term2 <- if (is.na(c_sum) || c_sum == 0) 0 else count / c_sum
      
      term1 + term2
    }
  })
  
  kept_levels <- setdiff(all_known_types, removed_types)
  conf_matrix_long$Original <- factor(conf_matrix_long$Original, levels = kept_levels)
  conf_matrix_long$New <- factor(conf_matrix_long$New, levels = kept_levels)
  conf_matrix_long <- conf_matrix_long %>% drop_na()
  
  caption_text <- if (length(removed_types) > 0) {
    paste0("Excluded due to 0 or insufficient initial cells (<", min_cell_count_for_type, "): ", paste(removed_types, collapse = ", "))
  } else {
    NULL
  }
  
  conf_matrix_long <- conf_matrix_long %>%
    mutate(Count_label = ifelse(Count == 0, "0", as.character(Count)))
  
  # 1. Colored by raw counts
  confusion_matrix_raw_plot <- ggplot(conf_matrix_long, aes(x = New, y = Original, fill = Count)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
             color = "black", fill = NA, linewidth = 1) +
    geom_tile(color = "white") +
    geom_text(aes(label = Count_label),
              size = 7) +
    scale_fill_continuous(low = "white",
                          high = "blue",
                          name = "Count",
                          trans = "log1p") +
    theme_minimal() +
    labs(title = "Transitions with noise (colored by raw counts)", caption = caption_text) +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 14, angle = 0, hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.caption = element_text(hjust = 0, size = 10, color = "gray30"))
  
  ggsave(paste0(output_prefix_base, "confusion_matrix_raw_counts.svg"), plot = confusion_matrix_raw_plot,
         width = 9, height = 9, units = "in")

  # 2. Colored by normalized sum & diagonal stability-based instability
 
  if (!is.null(caption_text)) {
    caption_text <- paste0(caption_text, '\n')
  } 
  caption_text <- paste0(caption_text, 'Colored by: diagonal: Bi_Stability-based, \n\t\tnon-diagonal: cell/rowSum + cell/newColRowSum')

  confusion_matrix_sum_plot <- ggplot(conf_matrix_long, aes(x = New, y = Original, fill = NormalizedSum)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
             color = "black", fill = NA, linewidth = 1) +
    geom_tile(color = "white") +
    geom_text(aes(label = Count_label),
              size = 7) +
    scale_fill_continuous(low = "white",
                          high = "blue",
                          name = "NormalizedSum",
                          trans = "log1p") +
    theme_minimal() +
    labs(title = "Transitions with noise raw counts", caption = caption_text) +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 14, angle = 0, hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.caption = element_text(hjust = 0, size = 10, color = "gray30"))
  
  ggsave(paste0(output_prefix_base, "confusion_matrix_raw_counts_sum.svg"), plot = confusion_matrix_sum_plot,
         width = 9, height = 9, units = "in")
}

plot_change_proportion <- function(noised_prediction, main_cell_type_order, output_prefix_base) {
  lvls <- if(is.factor(main_cell_type_order)) levels(main_cell_type_order) else main_cell_type_order
  
  plot_data <- as.data.frame(noised_prediction[, c(1, ncol(noised_prediction))]) %>%
    setNames(c("Initial_Prediction", "Noised_Prediction")) %>%
    group_by(Initial_Prediction) %>%
    count(Noised_Prediction) %>%
    mutate(proportion = n/sum(n)) %>%
    ungroup()
  
  label_threshold <- 0.01
  change_proportion_plot <- ggplot(plot_data, aes(x = factor(Initial_Prediction, levels = lvls), 
                                                  y = proportion, fill = Noised_Prediction)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    geom_text(
      aes(label = ifelse(proportion > label_threshold, 
                         Noised_Prediction,""), 
          size = proportion), 
      position = position_stack(vjust = 0.5),
      color = 'white'
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
      axis.ticks.y = element_line(color = "black", linewidth = 0.5),
      axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 13, angle = 0, hjust = 0.4),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold")
    )
  
  ggsave(paste0(output_prefix_base, "change_proportion_plot.svg"), plot = change_proportion_plot,
         width = 9, height = 9, units = "in")
}

plot_stability_vs_pss_scatters <- function(stability_pred, output_prefix_base) {
  stability_plot = ggplot(stability_pred, aes(x=PSS, y=Stability)) + geom_point() + theme_minimal() +
    ggtitle("Stability vs PSS") + geom_text(label=rownames(stability_pred), vjust = 1.5)
  
  ggsave(paste0(output_prefix_base, "stability_vs_pss_scatter_plot.svg"), plot = stability_plot,
         width = 12, height = 9, units = "in")

  bi_stability_plot = ggplot(stability_pred, aes(x=PSS, y=Bi_Stability)) + geom_point() + theme_minimal() +
    ggtitle("Bi-directional Stability vs PSS") + geom_text(label=rownames(stability_pred), vjust = 1.5)
  
  ggsave(paste0(output_prefix_base, "bi_stability_vs_pss_scatter_plot.svg"), plot = bi_stability_plot,
         width = 12, height = 9, units = "in")

  f1_stability_plot = ggplot(stability_pred, aes(x=PSS, y=F1_Stability)) + geom_point() + theme_minimal() +
    ggtitle("F1 Stability vs PSS") + geom_text(label=rownames(stability_pred), vjust = 1.5)
  
  ggsave(paste0(output_prefix_base, "f1_stability_vs_pss_scatter_plot.svg"), plot = f1_stability_plot,
         width = 12, height = 9, units = "in")
}

plot_pss_heatmap <- function(RSS_mat_filtered, types_to_run_on, main_cell_type_order, output_prefix_base) {
  RSS_long <- as.data.frame(as.table(RSS_mat_filtered[types_to_run_on, types_to_run_on]))
  names(RSS_long) <- c("P_R", "P_C", "RSS")
  
  RSS_long$P_R <- factor(RSS_long$P_R, levels = main_cell_type_order)
  RSS_long$P_C <- factor(RSS_long$P_C, levels = main_cell_type_order)
  RSS_long <- RSS_long %>% drop_na()
  RSS_long <- RSS_long %>%
    mutate(
      RSS_label = sprintf("%.2f", RSS), 
      RSS_label = ifelse(RSS_label == "-0.00", "0.00", RSS_label),
      RSS_label = ifelse(RSS_label == "0.00", "0", RSS_label)
    )
  
  rss_heatmap <- ggplot(RSS_long, aes(x = P_C, y = P_R, fill = RSS)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
             color = "black", fill = NA, linewidth = 1) +
    geom_tile(color = "white") +
    geom_text(aes(label = RSS_label),
              size = 7) +
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
}

drop_nan_rows_cols <- function(input_matrix, dimension = c("both", "row", "col")) {
  dimension <- match.arg(dimension)
  
  if (dimension %in% c("both", "row")) {
    all_nan_rows <- apply(input_matrix, 1, function(row) all(is.nan(row)))
    rows_to_keep <- !all_nan_rows
    input_matrix <- input_matrix[rows_to_keep, , drop = FALSE]
  }
  
  if (dimension %in% c("both", "col")) {
    all_nan_cols <- apply(input_matrix, 2, function(col) all(is.nan(col)))
    cols_to_keep <- !all_nan_cols
    input_matrix <- input_matrix[, cols_to_keep, drop = FALSE]
  }
  
  return(input_matrix)
}

plot_change_feature <- function(seurat_obj, change_counts_df, colors_feature_plot_noise = c('grey', '#f03b20'), output_prefix_base) {
  seurat_obj_with_changes <- AddMetaData(seurat_obj, change_counts_df$changes, col.name="Changes")
  change_feature_plot <- FeaturePlot(seurat_obj_with_changes,
                                     reduction = "umap",
                                     features = "Changes",
                                     cols = colors_feature_plot_noise,
                                     keep.scale = "all") +
    ggtitle("Number of Prediction Changes with Noise")
  ggsave(paste0(output_prefix_base, "change_feature_plot.svg"), plot = change_feature_plot,
         width = 16, height = 9, units = "in")
}

plot_noised_predictions_original_view <- function(seurat_obj, temp_seurat_obj, results_noised, output_prefix_base) {
  assay_names <- names(temp_seurat_obj[["test"]]@assays)
  pred_assay <- assay_names[grepl("^prediction\\.score\\.", assay_names)][1]
  prediction_data_original = temp_seurat_obj[["test"]]@assays[[pred_assay]]@data
  rownames(prediction_data_original)=gsub("-", "_", rownames(prediction_data_original))
  prediction_data = results_noised[["test"]]@assays[[pred_assay]]@data
  rownames(prediction_data)=gsub("-", "_", rownames(prediction_data))
  
  all_cell_types <- intersect(rownames(prediction_data), rownames(prediction_data_original))
  for (t in all_cell_types){
    original_prediction_values <- if(t %in% rownames(prediction_data_original)) {
      prediction_data_original[t, , drop = FALSE]
    } else {
      print("NA_original")
      rep(NA, ncol(seurat_obj))
    }
    
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
                                         keep.scale = "all") + ggtitle(paste0("Noised ",t,"_prediction"))
    original_prediction_plot = FeaturePlot(AddMetaData(seurat_obj, t(original_prediction_values),
                                                       col.name=paste0(t,"_prediction")),
                                           reduction = "umap",
                                           features = paste0(t,"_prediction"),
                                           keep.scale = "all") + ggtitle(paste0("Original ",t,"_prediction"))
    
    combined_feature_plot <- noised_prediction_plot | original_prediction_plot
    ggsave(paste0(output_prefix_base,"noised_",t, "_prediction_original_view.svg"),plot = combined_feature_plot,
           width = 16, height = 9, units = "in")
  }
}

plot_jsd_heatmap <- function(JSD_mat, main_cell_type_order, output_prefix_base) {
  JSD_mat_filtered <- drop_nan_rows_cols(JSD_mat, "both")
  
  JSD_long <- as.data.frame(as.table(JSD_mat_filtered))
  names(JSD_long) <- c("P_R", "P_C", "JSD")
  
  lvls <- if(is.factor(main_cell_type_order)) levels(main_cell_type_order) else main_cell_type_order
  JSD_long$P_R <- factor(JSD_long$P_R, levels = lvls)
  JSD_long$P_C <- factor(JSD_long$P_C, levels = lvls)
  JSD_long <- JSD_long %>% drop_na()
  JSD_long <- JSD_long %>%
    mutate(
      JSD_label = sprintf("%.2f", JSD),
      JSD_label = ifelse(JSD_label == "-0.00", "0.00", JSD_label),
      JSD_label = ifelse(JSD_label == "0.00", "0", JSD_label)
    )
  
  jsd_heatmap <- ggplot(JSD_long, aes(x = P_C, y = P_R, fill = JSD)) +
    geom_tile(color = "white") +
    geom_text(aes(label = JSD_label),
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
}

plot_seurat_dependent_plots_from_cache <- function(cache_path, change_counts_df, seurat_noised_prediction_list, output_prefix_base) {
  clean_output_prefix <- sub("/$", "", output_prefix_base)
  original_dir <- file.path(dirname(cache_path), basename(clean_output_prefix))
  
  run_without_noise_path <- file.path(original_dir, "_cache", "run_without_noise.rds")
  train_seurat_processed_path <- file.path(original_dir, "_cache", "train_seurat_processed.rds")
  
  if (!file.exists(train_seurat_processed_path) || !file.exists(run_without_noise_path)) {
    cat("Could not find cache files. Skipping DimPlots, change features, and noised prediction original view plots.\n")
    return()
  }
  
  seurat_obj <- readRDS(train_seurat_processed_path)
  temp_seurat_obj <- readRDS(run_without_noise_path)
  test_query <- temp_seurat_obj[["test"]]
  
  pred_col <- names(test_query@meta.data)[grepl("^predicted\\.", names(test_query@meta.data))][1]
  
  current_cells <- names(test_query@meta.data[[pred_col]])
  if (is.null(current_cells)) current_cells <- colnames(test_query)
  
  # Reconstruct myColors
  cell_type_order <- if (!is.null(temp_seurat_obj[["cell_type_order"]])) {
    levels(temp_seurat_obj[["cell_type_order"]])
  } else {
    sort(unique(test_query[[pred_col]][[pred_col]]))
  }
  if (is.null(cell_type_order)) cell_type_order <- sort(unique(test_query[[pred_col]][[pred_col]]))
  
  six2gfp_levels <- c("UM","CM","CM_DIV","PODO","PROX_1","PROX_2","LOH","DIST_CD","MACROPHAG","ENDO")
  if (any(cell_type_order %in% six2gfp_levels)) {
    myColors <- CELL_TYPE_COLORS
  } else {
    n_types <- length(cell_type_order)
    if (n_types <= 12) {
      myColors <- setNames(brewer.pal(max(n_types, 3), "Paired")[1:n_types], cell_type_order)
    } else {
      myColors <- setNames(scales::hue_pal()(n_types), cell_type_order)
    }
  }

  # 1. DimPlots
  if (pred_col %in% colnames(seurat_obj@meta.data)) {
    if(!is.null(seurat_obj@reductions$ref.umap)) { 
      reduction <- 'ref.umap' 
    } else { 
      reduction <- 'umap' 
    }
    
    seurat_obj[[pred_col]][current_cells] <- test_query[[pred_col]]
    
    p2_alone = DimPlot(seurat_obj, cells=current_cells, reduction = reduction,
                       group.by = pred_col, label = TRUE,
                       label.size = 5, repel = FALSE) + theme(legend.position = "none") +
      scale_color_manual(values = myColors)
    ggsave(paste0(output_prefix_base, "0_silent__original_view_only_prediction_A.svg"), plot = p2_alone,
           width = 8, height = 8, units = "in")
           
    p2_alone_with_legend = DimPlot(seurat_obj, cells=current_cells, reduction = reduction,
                       group.by = pred_col, label = TRUE,
                       label.size = 5, repel = FALSE) +
      scale_color_manual(values = myColors)
    ggsave(paste0(output_prefix_base, "0_silent__original_view_only_prediction_A_with_legend.svg"), plot = p2_alone_with_legend,
           width = 10, height = 8, units = "in")
           
  } else {
    seurat_obj <- AddMetaData(seurat_obj, test_query@meta.data[[pred_col]], col.name = pred_col)
    p2_alone = DimPlot(seurat_obj, reduction = "umap", label = TRUE, group.by = pred_col,
                       label.size = 5, repel = FALSE) + theme(legend.position = "none") +
      scale_color_manual(values = myColors)
    ggsave(paste0(output_prefix_base, "0_silent__original_view_only_prediction_A.svg"), plot = p2_alone,
           width = 8, height = 8, units = "in")
           
    p2_alone_with_legend = DimPlot(seurat_obj, reduction = "umap", label = TRUE, group.by = pred_col,
                       label.size = 5, repel = FALSE) +
      scale_color_manual(values = myColors)
    ggsave(paste0(output_prefix_base, "0_silent__original_view_only_prediction_A_with_legend.svg"), plot = p2_alone_with_legend,
           width = 10, height = 8, units = "in")
  }
  
  # 2. Change Feature Plot
  if (!is.null(change_counts_df)) {
    plot_change_feature(seurat_obj, change_counts_df, c('grey', '#f03b20'), output_prefix_base)
  }
  
  # 3. Noised Predictions Original View Plot
  if (!is.null(seurat_noised_prediction_list) && length(seurat_noised_prediction_list) > 0) {
    results_noised <- seurat_noised_prediction_list[[length(seurat_noised_prediction_list)]]
    plot_noised_predictions_original_view(seurat_obj, temp_seurat_obj, results_noised, output_prefix_base)
  }
}

# ============================================================================

regenerate_plots_from_cache <- function(cache_path, min_cell_count_for_type = 3,
                                        output_prefix_base = NULL) {
  if (!file.exists(cache_path)) {
    cat("Cache file not found:", cache_path, "\n")
    return()
  }
  
  cat("Regenerating plots for:", cache_path, "\n")
  data <- readRDS(cache_path)
  noised_prediction <- data$noised_prediction_matrix
  stability_pred <- data$stability_pred
  RSS_mat_filtered <- data$pss_matrix
  JSD_mat <- data$jsd_matrix
  change_counts_df <- data$change_counts
  
  # Use the provided output dir, or fall back to the cache file's directory
  if (is.null(output_prefix_base)) {
    output_prefix_base <- paste0(dirname(cache_path), "/")
  }
  
  if (is.null(noised_prediction) || is.null(stability_pred) || is.null(RSS_mat_filtered)) {
      cat("Missing required matrices in cache. Skipping.\n")
      return()
  }
  
  # Detect cell type order: use cached value if present, otherwise infer from data
  if (!is.null(data$main_cell_type_order)) {
    main_cell_type_order <- data$main_cell_type_order
  } else {
    # Infer from the matrices — fall back to rownames of JSD/PSS matrix or unique prediction values
    six2gfp_levels <- c("UM", "CM", "CM_DIV", "PODO", "PROX_1", "PROX_2", "LOH", "DIST_CD", "MACROPHAG", "ENDO")
    candidate_types <- if (!is.null(JSD_mat)) rownames(JSD_mat) else rownames(RSS_mat_filtered)
    if (any(candidate_types %in% six2gfp_levels)) {
      main_cell_type_order <- factor(six2gfp_levels, levels = six2gfp_levels)
    } else {
      # Non-SIX2GFP: use atlas ordering for known types, append any others sorted
      atlas_levels <- gsub("-", "_", c(
        "Cap mesenchyme", "Proliferating cap mesenchyme",
        "Proximal renal vesicle", "Distal renal vesicle", "Proliferating distal renal vesicle",
        "Proximal S shaped body", "Medial S shaped body", "Distal S shaped body",
        "Podocyte", "Proximal tubule", "Loop of Henle",
        "Proximal UB", "CNT/PC - proximal UB", "Pelvic epithelium - distal UB",
        "Stroma progenitor", "Proliferating stroma progenitor",
        "Fibroblast 1", "Fibroblast 2",
        "Myofibroblast 1", "Myofibroblast 2", "Proliferating myofibroblast",
        "Endothelium",
        "Macrophage 1", "Macrophage 2", "Proliferating macrophage",
        "Monocyte", "Proliferating monocyte", "Mast cells",
        "cDC1", "cDC2", "pDC", "Neutrophil",
        "B cell", "Proliferating B cell", "CD4 T cell", "CD8 T cell", "NK cell", "Innate like lymphocyte",
        "Erythroid", "Megakaryocyte", "Neuron"
      ))
      found <- atlas_levels[atlas_levels %in% candidate_types]
      other <- sort(setdiff(candidate_types, atlas_levels))
      ordered_types <- c(found, other)
      main_cell_type_order <- factor(ordered_types, levels = ordered_types)
    }
  }
  
  # Normalize noised_prediction names to match gsub'd convention
  noised_prediction[,] <- gsub("-", "_", noised_prediction)
  
  conf_matrix <- table(Original = noised_prediction[,1],
                       New = noised_prediction[,ncol(noised_prediction)])
  types_counts_conf_matrix <- rowSums(conf_matrix)
  rows_to_keep <- names(types_counts_conf_matrix)[types_counts_conf_matrix >= min_cell_count_for_type]
  
  all_known_types <- levels(main_cell_type_order)
  removed_types_present <- names(types_counts_conf_matrix)[types_counts_conf_matrix < min_cell_count_for_type]
  absent_types <- setdiff(all_known_types, names(types_counts_conf_matrix))
  removed_types <- unique(c(removed_types_present, absent_types))
  
  cols_to_keep <- union(rows_to_keep, colnames(conf_matrix)[colSums(conf_matrix[rows_to_keep, , drop=FALSE]) > 0])
  cols_to_keep <- setdiff(cols_to_keep, removed_types)
  
  conf_matrix_filtered <- conf_matrix[rows_to_keep, cols_to_keep, drop=FALSE]
  conf_matrix_fraction <- prop.table(conf_matrix_filtered, margin = 1)
  
  types_to_run_on <- intersect(colnames(RSS_mat_filtered), rownames(conf_matrix_fraction))
  types_to_run_on <- intersect(types_to_run_on, levels(main_cell_type_order))
  
  if (length(types_to_run_on) < 2) {
      cat("Not enough valid cell types to plot. Skipping.\n")
      return()
  }
  
  # Reuse shared helper functions
  plot_change_proportion(noised_prediction, main_cell_type_order, output_prefix_base)
  plot_confusion_matrix_heatmap(conf_matrix_fraction, all_known_types, removed_types, min_cell_count_for_type, output_prefix_base)
  plot_confusion_matrix_raw_counts(conf_matrix_filtered, all_known_types, removed_types, min_cell_count_for_type, output_prefix_base)
  
  if (!is.null(JSD_mat)) {
    plot_jsd_heatmap(JSD_mat, main_cell_type_order, output_prefix_base)
  }
  
  plot_pss_heatmap(RSS_mat_filtered, types_to_run_on, main_cell_type_order, output_prefix_base)
  
  # DimPlots, change counts feature plot, and noised prediction feature plots from cache
  plot_seurat_dependent_plots_from_cache(cache_path, change_counts_df, data$seurat_noised_prediction_list, output_prefix_base)
  
  plot_stability_vs_pss_scatters(stability_pred, output_prefix_base)
  plot_and_save_stability_scatters(conf_matrix_fraction, RSS_mat_filtered, stability_pred, types_to_run_on, main_cell_type_order, output_prefix_base)
  # plot_normalized_sum_vs_pss(conf_matrix_filtered, RSS_mat_filtered, types_to_run_on, main_cell_type_order, output_prefix_base)
  cat("Success.\n")
}
