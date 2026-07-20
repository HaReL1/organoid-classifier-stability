# =============================================================================
# Human-to-Human Pipeline: Cell Atlas as Train/Test, Organoid as Reference
# =============================================================================
# This script reverses the usual flow:
#   CURRENT: SIX2GFP (mouse) -> train/test, Organoid -> query
#   THIS:    Cell Atlas (human fetal kidney) -> train/test, Organoid -> query
#
# This addresses the concern that the current pipeline uses mouse-to-human
# transfer and hasn't been validated with human-to-human.
# 
# mesenchymal stem cell — earliest progenitor (≈ UM)
# mesenchymal cell — committed mesenchyme (≈ CM)
# kidney cell — intermediate progenitor (≈ CM_DIV)
# podocyte — glomerular (≈ PODO)
# epithelial cell of proximal tubule — (≈ PROX)
# kidney loop of Henle epithelial cell — (≈ LOH)
# kidney epithelial cell — distal/collecting duct (≈ DIST_CD)
# myofibroblast cell — stromal
# macrophage — (≈ MACROPHAG)
# monocyte — immune
# endothelial cell — (≈ ENDO)
# =============================================================================

library(Seurat)
source("FULL_FUNCTION.R")

# =============================================================================
# 1. Load the Human Fetal Kidney Cell Atlas and prepare train/test
# =============================================================================
{
  # "Fetal kidney dataset: full" from https://cellxgene.cziscience.com/datasets
  # https://cellxgene.cziscience.com/e/d7dcfd8f-2ee7-4385-b9ac-e074c23ed190.cxg/
  atlas_object = readRDS("fcb6225a-b329-4ac9-8a3f-559ca9bac50e.rds")
  
  # Fix gene names (ENSG IDs -> gene symbols):
  {
    feature_map = gsub("_", "-",atlas_object@assays[["RNA"]]@meta.features$feature_name)
    names(feature_map) = rownames(atlas_object@assays[["RNA"]]@meta.features)
    current_features <- rownames(atlas_object)
    new_features <- as.character(feature_map[current_features])
    rownames(atlas_object@assays[["RNA"]]@counts) <- new_features
    rownames(atlas_object@assays[["RNA"]]@data) <- new_features
    if (!is.null(atlas_object@assays[["RNA"]]@scale.data) && nrow(atlas_object@assays[["RNA"]]@scale.data) > 0) {
      rownames(atlas_object@assays[["RNA"]]@scale.data) <- new_features
    }
    atlas_object@assays[["RNA"]]@meta.features <- atlas_object@assays[["RNA"]]@meta.features
    rownames(atlas_object@assays[["RNA"]]@meta.features) <- new_features
    rm(feature_map, current_features, new_features)
  }
  
  # Remove immune and irrelevant cell types (same as use_FULL_FUNCTION.R):
  cells_to_keep <- rownames(atlas_object@meta.data)[!atlas_object$cell_type %in% c("B cell", "neuron", "lymphocyte",
                                                                                   "natural killer cell","plasmacytoid dendritic cell",
                                                                                   "erythroid lineage cell", "CD4-positive, alpha-beta T cell",
                                                                                   "megakaryocyte", "mast cell", "CD8-positive, alpha-beta T cell",
                                                                                   "neutrophil", "conventional dendritic cell", "fibroblast")]
  atlas_object <- subset(atlas_object, cells = cells_to_keep)
  remove(cells_to_keep)
  
  cat("Atlas cell types after filtering:\n")
  print(table(atlas_object$author_cell_type))
  cat("\nTotal atlas cells:", ncol(atlas_object), "\n")
  
  # Change specific author_cell_type names
  atlas_object$author_cell_type <- as.character(atlas_object$author_cell_type)
  atlas_object$author_cell_type[atlas_object$author_cell_type == "CNT/PC - proximal UB"] <- "CNT PC proximal UB"
  atlas_object$author_cell_type[atlas_object$author_cell_type == "Pelvic epithelium - distal UB"] <- "Pelvic epithelium distal UB"
}

# =============================================================================
# 2. Split atlas 50/50 into train and test
# =============================================================================
{
  set.seed(42)
  all_cells <- colnames(atlas_object)
  sample_vector <- sample(c(TRUE, FALSE), length(all_cells), replace = TRUE, prob = c(0.8, 0.2))
  train_cells <- all_cells[sample_vector]
  test_cells  <- all_cells[!sample_vector]
  
  # Create train Seurat object with cell_type labels
  atlas_train <- subset(atlas_object, cells = train_cells)
  atlas_train <- AddMetaData(atlas_train, "atlas_train", col.name = "source")
  # The cell_type column already exists from the atlas metadata
  
  # Create test Seurat object
  atlas_test <- subset(atlas_object, cells = test_cells)
  atlas_test <- AddMetaData(atlas_test, "atlas_test", col.name = "source")
  
  cat("Train cells:", ncol(atlas_train), "\n")
  cat("Test cells:", ncol(atlas_test), "\n")
  rm(all_cells, sample_vector, train_cells, test_cells)
}

# =============================================================================
# 3. Load organoid datasets to use as query
# =============================================================================
{
  Uchimura_Humphreys_20 = read.table("GSM3763147_Uchimura.dge.txt")
  load(file = "uchimura_off_target_cell_list") # from pseudo_labels_cleaner_uchimura.R
  seurat_Uchimura_Humphreys_20 <- CreateSeuratObject(counts = Uchimura_Humphreys_20[,!colnames(Uchimura_Humphreys_20) %in% uchimura_off_target_cell_list],
                                                     assay = "RNA")
  remove(Uchimura_Humphreys_20, uchimura_off_target_cell_list)
} # load clean Uchimura

# load(file = "freedman_seurat_processed") # freedman_seurat_obj

# {
#   Takasato_Humphreys_20 = read.table("GSM3763146_Takasato.dge.txt")
#   load(file = "takasato_off_target_cell_list")
#   seurat_Takasato_Humphreys_20 <- CreateSeuratObject(counts = Takasato_Humphreys_20[,!colnames(Takasato_Humphreys_20) %in% takasato_off_target_cell_list],
#                                                      assay = "RNA")
#   remove(Takasato_Humphreys_20)
# } # load Takasato

# =============================================================================
# 4. Run the pipeline: Human atlas -> Uchimura (human-to-human)
# =============================================================================
cache_h2h_Uchimura <- "human_to_human/Uchimura_full_flow5noise.rds"
if (file.exists(cache_h2h_Uchimura)) {
  h2h_Uchimura_flow <- readRDS(cache_h2h_Uchimura)
} else {
  h2h_Uchimura_flow = analyze_noise_impact_on_prediction(
    seurat_Uchimura_Humphreys_20,
    noised_number = 5,
    train_Six2GFP = atlas_train,       # <-- atlas train instead of SIX2GFP
    test_Six2GFP  = atlas_test,        # <-- atlas test instead of SIX2GFP
    ref_cell_type_column = "author_cell_type",
    dims = 1:30,
    train_title = "Fetal Kidney Atlas",
    test_title_prefix = "Uchimura",
    n_neighbors = 8,
    skip_neighbors = TRUE,
    output_prefix_base = "human_to_human/Uchimura_noised/",
    prediction_column_name = "predicted.author_cell_type",
    colors_feature_plot_noise = c('grey', '#f03b20'),
    myColors_cell_types = NULL, # will use Paired palette (auto-generated for atlas types)
    return_all_suerats = FALSE,
    use_cache = TRUE,
    return_anchors = TRUE,
    noise_model = "poisson"
  )
  saveRDS(h2h_Uchimura_flow, cache_h2h_Uchimura)
}

# =============================================================================
# 5. Run the pipeline: Human atlas -> Freedman (human-to-human)
# =============================================================================
# cache_h2h_freedman <- "human_to_human/freedman_full_flow10noise.rds"
# if (file.exists(cache_h2h_freedman)) {
#   h2h_freedman_flow <- readRDS(cache_h2h_freedman)
# } else {
#   h2h_freedman_flow = analyze_noise_impact_on_prediction(
#     freedman_seurat_obj,
#     noised_number = 10,
#     train_Six2GFP = atlas_train,
#     test_Six2GFP  = atlas_test,
#     ref_cell_type_column = "author_cell_type",
#     dims = 1:15,
#     train_title = "Fetal Kidney Atlas",
#     test_title_prefix = "Freedman",
#     n_neighbors = 8,
#     skip_neighbors = TRUE,
#     output_prefix_base = "human_to_human/freedman_noised/",
#     prediction_column_name = "predicted.author_cell_type",
#     colors_feature_plot_noise = c('grey', '#f03b20'),
#     myColors_cell_types = NULL,
#     return_all_suerats = FALSE,
#     use_cache = TRUE,
#     return_anchors = TRUE,
#     noise_model = "poisson"
#   )
#   saveRDS(h2h_freedman_flow, cache_h2h_freedman)
# }

# =============================================================================
# 6. Run the pipeline: Human atlas -> Takasato (human-to-human)
# =============================================================================
# cache_h2h_Takasato <- "human_to_human/Takasato_full_flow10noise.rds"
# if (file.exists(cache_h2h_Takasato)) {
#   h2h_Takasato_flow <- readRDS(cache_h2h_Takasato)
# } else {
#   h2h_Takasato_flow = analyze_noise_impact_on_prediction(
#     seurat_Takasato_Humphreys_20,
#     noised_number = 10,
#     train_Six2GFP = atlas_train,
#     test_Six2GFP  = atlas_test,
#     ref_cell_type_column = "author_cell_type",
#     dims = 1:30,
#     train_title = "Fetal Kidney Atlas",
#     test_title_prefix = "Takasato",
#     n_neighbors = 8,
#     skip_neighbors = TRUE,
#     output_prefix_base = "human_to_human/Takasato_noised/",
#     prediction_column_name = "predicted.author_cell_type",
#     colors_feature_plot_noise = c('grey', '#f03b20'),
#     myColors_cell_types = NULL,
#     return_all_suerats = FALSE,
#     use_cache = TRUE,
#     return_anchors = TRUE,
#     noise_model = "poisson"
#   )
#   saveRDS(h2h_Takasato_flow, cache_h2h_Takasato)
# }

# =============================================================================
# 7. (Optional) Run atlas as its own query for self-validation
#    This tests atlas-to-atlas (human-to-human, same source)
# =============================================================================
# cache_h2h_atlas_self <- "human_to_human/atlas_self_flow10noise.rds"
# if (file.exists(cache_h2h_atlas_self)) {
#   h2h_atlas_self_flow <- readRDS(cache_h2h_atlas_self)
# } else {
#   h2h_atlas_self_flow = analyze_noise_impact_on_prediction(
#     atlas_object,  # full atlas as query
#     noised_number = 10,
#     train_Six2GFP = atlas_train,
#     test_Six2GFP  = atlas_test,
#     ref_cell_type_column = "author_cell_type",
#     dims = 1:30,
#     train_title = "Fetal Kidney Atlas",
#     test_title_prefix = "Atlas (self)",
#     n_neighbors = 8,
#     skip_neighbors = TRUE,
#     output_prefix_base = "human_to_human/atlas_self_noised/",
#     prediction_column_name = "predicted.author_cell_type",
#     colors_feature_plot_noise = c('grey', '#f03b20'),
#     myColors_cell_types = NULL,
#     return_all_suerats = FALSE,
#     use_cache = TRUE,
#     return_anchors = TRUE,
#     noise_model = "poisson"
#   )
#   saveRDS(h2h_atlas_self_flow, cache_h2h_atlas_self)
# }
