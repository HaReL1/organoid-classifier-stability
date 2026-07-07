library(Seurat)

# 1. Load FULL_FUNCTION and Training Data
source("FULL_FUNCTION.R")
load(file = "six2gfp/train_Six2GFP")
load(file = "six2gfp/test_Six2GFP")

# 2. Read the Vanslambrouck data D13+14
# Note: We are using ReadMtx because the raw 10X files have a custom prefix ("GSM5600483_ML204060_").
# The Read10X function expects the files to be named exactly "matrix.mtx.gz", "features.tsv.gz", and "barcodes.tsv.gz".
# ReadMtx works exactly the same way but allows us to specify the exact file names.
data_dir <- "Vanslambrouck"
expression_matrix <- ReadMtx(
  mtx = file.path(data_dir, "GSM5600483_ML204060_matrix.mtx.gz"),
  cells = file.path(data_dir, "GSM5600483_ML204060_barcodes.tsv.gz"),
  features = file.path(data_dir, "GSM5600483_ML204060_features.tsv.gz") 
)

seurat_Vanslambrouck <- CreateSeuratObject(counts = expression_matrix, assay = "RNA")
rm(expression_matrix) # clear up memory

# 3. Quality Control (QC)
seurat_Vanslambrouck[["percent.mt"]] <- PercentageFeatureSet(seurat_Vanslambrouck, pattern = "^MT-")
# print(summary(seurat_Vanslambrouck$nFeature_RNA))
# print(summary(seurat_Vanslambrouck$percent.mt))
# nFeature_RNA > 500 & nFeature_RNA < 9000: This easily captures the interquartile range (2846 to 4894) and the mean (4048) while successfully cropping out extreme outliers at the very bottom (empty/dead droplets) and the very top (potential doublets >9000).
# percent.mt < 15: Since your 3rd quartile is ~6%, a cutoff of 15 safely preserves the vast majority of healthy cells while comfortably excluding the high-mitochondrial dead/dying cells (max was 95%).
seurat_Vanslambrouck <- subset(seurat_Vanslambrouck, subset = nFeature_RNA > 500 & nFeature_RNA < 9000 & percent.mt < 15)



# 4. Run the Pipeline
# This matches the parameters used in use_FULL_FUNCTION.R
Vanslambrouck_full_flow = analyze_noise_impact_on_prediction(
  seurat_obj = seurat_Vanslambrouck,
  noised_number = 10,
  train_Six2GFP = train_Six2GFP, 
  test_Six2GFP = test_Six2GFP,
  ref_cell_type_column = "type",
  dims = 1:30,
  train_title = "SIX2GFP",
  test_title_prefix = "Vanslambrouck",
  n_neighbors = 8,
  skip_neighbors = TRUE,
  output_prefix_base = "six2gfp/Vanslambrouck_noised/",
  prediction_column_name = "predicted.type",
  colors_feature_plot_noise = c('grey', '#f03b20'),
  myColors_cell_types = NULL,
  return_all_suerats = FALSE,
  use_cache = TRUE,
  return_anchors = TRUE,
  noise_model = "poisson"
)

# 5. Save the output
saveRDS(Vanslambrouck_full_flow, "six2gfp/Vanslambrouck_full_flow10noise.rds")


# plot specific GENE

# 1. Load the "original view" Seurat object from cache
cache_file <- "six2gfp/Vanslambrouck_noised/_cache/train_seurat_processed.rds"
Vanslambrouck_init <- readRDS(cache_file)
# 2. Plot the specific gene (SIX2)
p <- FeaturePlot(Vanslambrouck_init, features = "SIX2", reduction = "umap") +
  ggtitle("SIX2 in Vanslambrouck (Original View)")
# 3. Save the plot
ggsave("six2gfp/Vanslambrouck_noised/SIX2_original_view.svg", plot = p, width = 10, height = 8, units = "in")

FeaturePlot(Vanslambrouck_init, features = "CDH11", reduction = "umap")
