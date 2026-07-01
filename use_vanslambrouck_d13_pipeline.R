library(Seurat)

# 1. Load FULL_FUNCTION and Training Data
source("FULL_FUNCTION.R")
load(file = "six2gfp/train_Six2GFP")
load(file = "six2gfp/test_Six2GFP")

# 2. Read the Vanslambrouck d13 data
# Note: We are using ReadMtx because the raw 10X files have a custom prefix.
data_dir <- "Vanslambrouck_d13"
expression_matrix <- ReadMtx(
  mtx = file.path(data_dir, "GSM5600482_ML203459_matrix.mtx.gz"),
  cells = file.path(data_dir, "GSM5600482_ML203459_barcodes.tsv.gz"),
  features = file.path(data_dir, "GSM5600482_ML203459_features.tsv.gz") 
)

seurat_Vanslambrouck_d13 <- CreateSeuratObject(counts = expression_matrix, assay = "RNA")
rm(expression_matrix) # clear up memory

# 3. Quality Control (QC)
seurat_Vanslambrouck_d13[["percent.mt"]] <- PercentageFeatureSet(seurat_Vanslambrouck_d13, pattern = "^MT-")
# nFeature_RNA > 500 & nFeature_RNA < 9000
# percent.mt < 15
seurat_Vanslambrouck_d13 <- subset(seurat_Vanslambrouck_d13, subset = nFeature_RNA > 500 & nFeature_RNA < 9000 & percent.mt < 15)

# 4. Run the Pipeline
Vanslambrouck_d13_full_flow = analyze_noise_impact_on_prediction(
  seurat_obj = seurat_Vanslambrouck_d13,
  noised_number = 10,
  train_Six2GFP = train_Six2GFP, 
  test_Six2GFP = test_Six2GFP,
  ref_cell_type_column = "type",
  dims = 1:30,
  train_title = "SIX2GFP",
  test_title_prefix = "Vanslambrouck_d13",
  n_neighbors = 8,
  skip_neighbors = TRUE,
  output_prefix_base = "six2gfp/Vanslambrouck_d13_noised/",
  prediction_column_name = "predicted.type",
  colors_feature_plot_noise = c('grey', '#f03b20'),
  myColors_cell_types = NULL,
  return_all_suerats = FALSE,
  use_cache = TRUE,
  return_anchors = TRUE,
  noise_model = "poisson"
)

# 5. Save the output
  saveRDS(Vanslambrouck_d13_full_flow, "six2gfp/Vanslambrouck_d13_full_flow10noise.rds")

# plot specific GENE

# 1. Load the "original view" Seurat object from cache
cache_file <- "six2gfp/Vanslambrouck_d13_noised/_cache/train_seurat_processed.rds"
Vanslambrouck_d13_init <- readRDS(cache_file)
# 2. Plot the specific gene (SIX2)
p <- 
  FeaturePlot(Vanslambrouck_d13_init, features = "SIX2", reduction = "umap") +
  ggtitle("SIX2 in Vanslambrouck d13 (Original View)")
# 3. Save the plot
ggsave("six2gfp/Vanslambrouck_d13_noised/SIX2_original_view.svg", plot = p, width = 10, height = 8, units = "in")

# working for the top group:
FeaturePlot(Vanslambrouck_d13_init, features = "CDH1", reduction = "umap")
FeaturePlot(Vanslambrouck_d13_init, features = "LHX1", reduction = "umap")
FeaturePlot(Vanslambrouck_d13_init, features = "EPCAM", reduction = "umap")

# FeaturePlot(Vanslambrouck_d13_init, features = "CDH11", reduction = "umap")
# FeaturePlot(Vanslambrouck_d13_init, features = "DLL1", reduction = "umap")
# FeaturePlot(Vanslambrouck_d13_init, features = "KCNJ1", reduction = "umap")
# FeaturePlot(Vanslambrouck_d13_init, features = "PAX2", reduction = "umap")
# FeaturePlot(Vanslambrouck_d13_init, features = "SCNN1B", reduction = "umap")
# FeaturePlot(Vanslambrouck_d13_init, features = "CDH2", reduction = "umap")
