# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3204322	#	sample type: singly picked organoid (not whole well)
# https://pubmed.ncbi.nlm.nih.gov/30626756/ - paper, Harder et al
# sample type: singly picked organoid
# time point: Harvested day 20
library(data.table)
freedman_gene_matrix <- fread("freedman/GSM3204322_NLW39_GeneMatrix.txt.gz", header = TRUE, sep = "\t")
# freedman_gene_matrix = read.table("freedman/GSM3204322_NLW39_GeneMatrix.txt.gz")
# Convert to standard matrix format
# Assuming first column contains gene names/IDs
row_names <- freedman_gene_matrix[[1]]
freedman_gene_matrix <- as.matrix(freedman_gene_matrix[, -1])
rownames(freedman_gene_matrix) <- row_names
freedman_seurat_obj <- CreateSeuratObject(counts = freedman_gene_matrix)
remove(row_names,freedman_gene_matrix)

# train_labeled_seurat <- NormalizeData(labeled_train_data)
# train_labeled_seurat <- FindVariableFeatures(train_labeled_seurat, selection.method = "vst", nfeatures = 2000)
# train_labeled_seurat <- ScaleData(train_labeled_seurat)
# train_labeled_seurat <- RunPCA(train_labeled_seurat)
# train_labeled_seurat <- FindNeighbors(train_labeled_seurat, dims = dims)
# train_labeled_seurat <- FindClusters(train_labeled_seurat)
# train_labeled_seurat <- RunUMAP(train_labeled_seurat, dims = dims, return.model = TRUE)

# 1. Quality Control and Filtering
# Calculate QC metrics
freedman_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(freedman_seurat_obj, pattern = "^MT-")

# Visualize QC metrics
# VlnPlot(freedman_seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells based on QC metrics (thresholds as descrived in article - 1926 cells left)
freedman_seurat_obj <- subset(freedman_seurat_obj, 
                              subset = nFeature_RNA > 500 & 
                                nFeature_RNA < 4000 & 
                                percent.mt < 25)

freedman_seurat_obj <- NormalizeData(freedman_seurat_obj)
freedman_seurat_obj <- FindVariableFeatures(freedman_seurat_obj, 
                                            selection.method = "vst", 
                                            nfeatures = 2000)
# Plot variable features
# top10 <- head(VariableFeatures(freedman_seurat_obj), 10)
# plot1 <- VariableFeaturePlot(freedman_seurat_obj)
# LabelPoints(plot = plot1, points = top10, repel = TRUE)
# remove(top10,plot1)

freedman_seurat_obj <- ScaleData(freedman_seurat_obj, features = rownames(freedman_seurat_obj))
freedman_seurat_obj <- RunPCA(freedman_seurat_obj)

# Visualize PCA results
# ElbowPlot(freedman_seurat_obj)
# DimPlot(freedman_seurat_obj, reduction = "pca")

# 6. Determine dimensionality
# Based on ElbowPlot, choose appropriate number of dimensions (e.g., 15)
n_dims <- 15

freedman_seurat_obj <- FindNeighbors(freedman_seurat_obj, dims = 1:n_dims)
freedman_seurat_obj <- FindClusters(freedman_seurat_obj, resolution = 0.5)
freedman_seurat_obj <- RunUMAP(freedman_seurat_obj, dims = 1:n_dims)

# 9. Visualize clusters
DimPlot(freedman_seurat_obj, reduction = "umap")

# 10. Find cluster markers
all_markers <- FindAllMarkers(freedman_seurat_obj, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25)

# View top markers per cluster
top_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

# 11. Visualize markers
FeaturePlot(freedman_seurat_obj, 
            features = unique(top_markers$gene)[1:9], 
            ncol = 3)

FeaturePlot(freedman_seurat_obj, features=c("COL3A1","TOP2A", "AFP"), reduction = "umap")
FeaturePlot(freedman_seurat_obj, features=c("MLANA","SOX2", "MAP2", "PITX2"), reduction = "umap")

# 12. Save the processed object
save(freedman_seurat_obj,file = "freedman_seurat_processed")

load(file = "freedman_seurat_processed")

freedman_flow = analyze_noise_impact_on_prediction(
  freedman_seurat_obj,
  noised_number = 1,
  train_Six2GFP, 
  test_Six2GFP,
  ref_cell_type_column = "type",
  dims = 1:15,
  train_title = "SIX2GFP",
  test_title_prefix = "Freedman",
  n_neighbors = 8,
  skip_neighbors = TRUE,
  output_prefix_base = "six2gfp/newest_function/freedman_not_clean/",
  prediction_column_name = "predicted.type", # Column name for predictions in metadata
  colors_feature_plot_noise = c('grey', '#f03b20'),
  myColors_cell_types = NULL, # Colors for cell type plots, if NULL, Paired palette will be used
  return_all_suerats = FALSE
)
