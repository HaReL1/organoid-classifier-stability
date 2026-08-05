# ==============================================================================
# Preprocessing All Datasets for Analysis Pipeline
# ==============================================================================
# This script consolidates the full, step-by-step preprocessing procedures 
# applied to each dataset evaluated in the manuscript.
# ==============================================================================

library(Seurat)
library(dplyr)
library(data.table)

# ==============================================================================
# Section 1: Six2GFP Mouse Reference Dataset
# ==============================================================================
# --- Full Preprocessing from Raw Data ---
load(file = "six2gfp/rawdata")
load(file = "six2gfp/metadata")

# 1. Filter out low quality and zero ACTB/GAPDH cells
Six2GFP_metadata <- Six2GFP_metadata[Six2GFP_metadata$Type != "LOWQUALITY", ]
Six2GFP_metadata <- Six2GFP_metadata[Six2GFP_metadata$Type != "ZERO_ACTB_GAPDH", ]
Six2GFP <- Six2GFP[, Six2GFP_metadata$CellName]

# 2. Partition 80:20 into training and test sets
set.seed(1)
sample_vector <- sample(c(TRUE, FALSE), ncol(Six2GFP), replace = TRUE, prob = c(0.8, 0.2))
train_cells <- colnames(Six2GFP)[sample_vector]
test_cells  <- colnames(Six2GFP)[!sample_vector]

# 3. Construct Seurat objects (min.cells = 10) and attach metadata
train_Six2GFP <- CreateSeuratObject(counts = Six2GFP[, train_cells], assay = "RNA", min.cells = 10)
train_Six2GFP <- AddMetaData(train_Six2GFP, "train_Six2GFP", col.name = "source")
train_Six2GFP <- AddMetaData(train_Six2GFP, Six2GFP_metadata[train_cells, 2], col.name = "type")

test_Six2GFP <- subset(CreateSeuratObject(counts = Six2GFP, assay = "RNA", min.cells = 10), cells = test_cells)
test_Six2GFP <- AddMetaData(test_Six2GFP, "test_Six2GFP", col.name = "source")

rm(Six2GFP, Six2GFP_metadata, sample_vector, train_cells, test_cells)

# --- Alternative Fast Load (pre-computed train/test objects) ---
# load(file = "six2gfp/train_Six2GFP")
# load(file = "six2gfp/test_Six2GFP")


# ==============================================================================
# Section 2: Uchimura et al. Organoid Dataset (GSM3763147)
# ==============================================================================
# 1. Load raw count matrix
Uchimura_Humphreys_20 <- read.table("GSM3763147_Uchimura.dge.txt")

# 2. Full workflow to identify off-target lineages (Run once to generate cell list):
# seurat_full <- CreateSeuratObject(counts = Uchimura_Humphreys_20, assay = "RNA")
# seurat_full <- NormalizeData(seurat_full)
# seurat_full <- FindVariableFeatures(seurat_full, selection.method = "vst", nfeatures = 2000)
# seurat_full <- ScaleData(seurat_full)
# seurat_full <- RunPCA(seurat_full)
# seurat_full <- FindNeighbors(seurat_full, dims = 1:30)
# seurat_full <- FindClusters(seurat_full)
# seurat_full <- RunUMAP(seurat_full, dims = 1:30, reduction.name = "umap")
#
# MLANA_group   <- CellSelector(FeaturePlot(seurat_full, features = "MLANA", keep.scale = "all")) # Melanocytes
# Neurons_group <- CellSelector(FeaturePlot(seurat_full, features = "SOX2", keep.scale = "all"))  # Neurons
# Muscle_group  <- CellSelector(FeaturePlot(seurat_full, features = "PITX2", keep.scale = "all")) # Muscle
# extra         <- CellSelector(FeaturePlot(seurat_full, features = "PITX2", keep.scale = "all"))
# uchimura_off_target_cell_list <- unique(c(MLANA_group, Neurons_group, Muscle_group, extra))
# save(uchimura_off_target_cell_list, file = "uchimura_off_target_cell_list")

# 3. Load pre-saved off-target cell barcode list
load(file = "uchimura_off_target_cell_list")

# 4. Construct cleaned Seurat object excluding off-target cells
seurat_Uchimura_Humphreys_20 <- CreateSeuratObject(
  counts = Uchimura_Humphreys_20[, !colnames(Uchimura_Humphreys_20) %in% uchimura_off_target_cell_list],
  assay = "RNA"
)

rm(Uchimura_Humphreys_20, uchimura_off_target_cell_list)


# ==============================================================================
# Section 3: Takasato et al. Organoid Dataset (GSM3763146)
# ==============================================================================
# 1. Load raw count matrix
Takasato_Humphreys_20 <- read.table("GSM3763146_Takasato.dge.txt")

# 2. Full workflow to identify off-target lineages (Run once to generate cell list):
# seurat_full <- CreateSeuratObject(counts = Takasato_Humphreys_20, assay = "RNA")
# seurat_full <- NormalizeData(seurat_full)
# seurat_full <- FindVariableFeatures(seurat_full, selection.method = "vst", nfeatures = 2000)
# seurat_full <- ScaleData(seurat_full)
# seurat_full <- RunPCA(seurat_full)
# seurat_full <- FindNeighbors(seurat_full, dims = 1:30)
# seurat_full <- FindClusters(seurat_full)
# seurat_full <- RunUMAP(seurat_full, dims = 1:30, reduction.name = "umap")
#
# MLANA_group    <- CellSelector(FeaturePlot(seurat_full, features = "MLANA", keep.scale = "all")) # Melanocytes
# Neurons_group  <- CellSelector(FeaturePlot(seurat_full, features = "MAP2", keep.scale = "all"))  # Neurons
# Neurons2_group <- CellSelector(FeaturePlot(seurat_full, features = "SOX2", keep.scale = "all"))  # Neurons 2
# extra          <- CellSelector(FeaturePlot(seurat_full, features = "COL1A1", keep.scale = "all"))# Stroma / Non-renal
# takasato_off_target_cell_list <- unique(c(MLANA_group, Neurons_group, Neurons2_group, extra))
# save(takasato_off_target_cell_list, file = "takasato_off_target_cell_list")

# 3. Load pre-saved off-target cell barcode list
load(file = "takasato_off_target_cell_list")

# 4. Construct cleaned Seurat object excluding off-target cells
seurat_Takasato_Humphreys_20 <- CreateSeuratObject(
  counts = Takasato_Humphreys_20[, !colnames(Takasato_Humphreys_20) %in% takasato_off_target_cell_list],
  assay = "RNA"
)

rm(Takasato_Humphreys_20, takasato_off_target_cell_list)


# ==============================================================================
# Section 4: Human Fetal Kidney Cell Atlas
# ==============================================================================
# Load CZ CELLxGENE RDS object
atlas_object <- readRDS("fcb6225a-b329-4ac9-8a3f-559ca9bac50e.rds")

# 1. Standardize Ensembl IDs to gene symbol names (hyphenated)
feature_map <- gsub("_", "-", atlas_object@assays[["RNA"]]@meta.features$feature_name)
names(feature_map) <- rownames(atlas_object@assays[["RNA"]]@meta.features)
current_features <- rownames(atlas_object)
new_features <- as.character(feature_map[current_features])

rownames(atlas_object@assays[["RNA"]]@counts) <- new_features
rownames(atlas_object@assays[["RNA"]]@data) <- new_features
if (!is.null(atlas_object@assays[["RNA"]]@scale.data) && nrow(atlas_object@assays[["RNA"]]@scale.data) > 0) {
  rownames(atlas_object@assays[["RNA"]]@scale.data) <- new_features
}
rownames(atlas_object@assays[["RNA"]]@meta.features) <- new_features
rm(feature_map, current_features, new_features)

# 2. Exclude non-renal, immune, and stromal lineages from metadata
cells_to_keep <- rownames(atlas_object@meta.data)[
  !atlas_object$cell_type %in% c(
    "B cell", "neuron", "lymphocyte", "natural killer cell", 
    "plasmacytoid dendritic cell", "erythroid lineage cell", 
    "CD4-positive, alpha-beta T cell", "megakaryocyte", "mast cell", 
    "CD8-positive, alpha-beta T cell", "neutrophil", 
    "conventional dendritic cell", "fibroblast"
  )
]
atlas_object <- subset(atlas_object, cells = cells_to_keep)
rm(cells_to_keep)


# ==============================================================================
# Section 5: Harder et al. Organoid Dataset (GSM3204322)
# ==============================================================================
# --- Full Preprocessing from Raw Matrix ---
freedman_gene_matrix <- fread("freedman/GSM3204322_NLW39_GeneMatrix.txt.gz", header = TRUE, sep = "\t")
row_names <- freedman_gene_matrix[[1]]
freedman_gene_matrix <- as.matrix(freedman_gene_matrix[, -1])
rownames(freedman_gene_matrix) <- row_names

freedman_seurat_obj <- CreateSeuratObject(counts = freedman_gene_matrix)
rm(row_names, freedman_gene_matrix)

# 1. Calculate mitochondrial gene percentage and apply published QC bounds
freedman_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(freedman_seurat_obj, pattern = "^MT-")
freedman_seurat_obj <- subset(
  freedman_seurat_obj, 
  subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 25
)

# 2. Full Seurat workflow (Normalization, HVG, Scaling, PCA, Clustering, UMAP)
freedman_seurat_obj <- NormalizeData(freedman_seurat_obj)
freedman_seurat_obj <- FindVariableFeatures(freedman_seurat_obj, selection.method = "vst", nfeatures = 2000)
freedman_seurat_obj <- ScaleData(freedman_seurat_obj, features = rownames(freedman_seurat_obj))
freedman_seurat_obj <- RunPCA(freedman_seurat_obj)
freedman_seurat_obj <- FindNeighbors(freedman_seurat_obj, dims = 1:15)
freedman_seurat_obj <- FindClusters(freedman_seurat_obj, resolution = 0.5)
freedman_seurat_obj <- RunUMAP(freedman_seurat_obj, dims = 1:15)

# Save processed object if desired:
# save(freedman_seurat_obj, file = "freedman_seurat_processed")

# --- Alternative Fast Load (pre-processed object) ---
# load(file = "freedman_seurat_processed")


# ==============================================================================
# Section 6: Vanslambrouck et al. Organoid Datasets (GSM5600483 & GSM5600482)
# ==============================================================================
# --- 6a. Vanslambrouck Day 13+14 (GSM5600483) ---
data_dir_v <- "Vanslambrouck"
expression_matrix_v <- ReadMtx(
  mtx = file.path(data_dir_v, "GSM5600483_ML204060_matrix.mtx.gz"),
  cells = file.path(data_dir_v, "GSM5600483_ML204060_barcodes.tsv.gz"),
  features = file.path(data_dir_v, "GSM5600483_ML204060_features.tsv.gz") 
)
seurat_Vanslambrouck <- CreateSeuratObject(counts = expression_matrix_v, assay = "RNA")
rm(expression_matrix_v, data_dir_v)

seurat_Vanslambrouck[["percent.mt"]] <- PercentageFeatureSet(seurat_Vanslambrouck, pattern = "^MT-")
seurat_Vanslambrouck <- subset(
  seurat_Vanslambrouck, 
  subset = nFeature_RNA > 500 & nFeature_RNA < 9000 & percent.mt < 15
)

# --- 6b. Vanslambrouck Day 13 (GSM5600482) ---
data_dir_v13 <- "Vanslambrouck_d13"
expression_matrix_v13 <- ReadMtx(
  mtx = file.path(data_dir_v13, "GSM5600482_ML203459_matrix.mtx.gz"),
  cells = file.path(data_dir_v13, "GSM5600482_ML203459_barcodes.tsv.gz"),
  features = file.path(data_dir_v13, "GSM5600482_ML203459_features.tsv.gz") 
)
seurat_Vanslambrouck_d13 <- CreateSeuratObject(counts = expression_matrix_v13, assay = "RNA")
rm(expression_matrix_v13, data_dir_v13)

seurat_Vanslambrouck_d13[["percent.mt"]] <- PercentageFeatureSet(seurat_Vanslambrouck_d13, pattern = "^MT-")
seurat_Vanslambrouck_d13 <- subset(
  seurat_Vanslambrouck_d13, 
  subset = nFeature_RNA > 500 & nFeature_RNA < 9000 & percent.mt < 15
)
