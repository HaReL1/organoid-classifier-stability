# using FULL_FUNCTION:
{
  load(file = "six2gfp/train_Six2GFP")
  load(file = "six2gfp/test_Six2GFP")
} # load six2gfp (test_Six2GFP, train_Six2GFP) - training data
{
  Uchimura_Humphreys_20 = read.table("GSM3763147_Uchimura.dge.txt")
  load(file = "uchimura_off_target_cell_list") # from pseudo_labels_cleaner_uchimura.R
  seurat_Uchimura_Humphreys_20 <- CreateSeuratObject(counts = Uchimura_Humphreys_20[,!colnames(Uchimura_Humphreys_20) %in% uchimura_off_target_cell_list],
                                                     assay = "RNA")
  remove(Uchimura_Humphreys_20, uchimura_off_target_cell_list)
} # load clean Uchimura (seurat_Uchimura_Humphreys_20) - Figure 2
load(file = "freedman_seurat_processed") # - Figure 3
{
  # install.packages("anndata") #ad <- anndata::read_h5ad('Fetal_full_v3.h5ad')
  # "Fetal kidney dataset: full" from https://cellxgene.cziscience.com/datasets
  # https://cellxgene.cziscience.com/e/d7dcfd8f-2ee7-4385-b9ac-e074c23ed190.cxg/
  atlas_object = readRDS("fcb6225a-b329-4ac9-8a3f-559ca9bac50e.rds")
  # fix gene names:
  {
    feature_map = gsub("_", "-",atlas_object@assays[["RNA"]]@meta.features$feature_name)
    names(feature_map) = rownames(atlas_object@assays[["RNA"]]@meta.features)
    current_features <- rownames(atlas_object)
    # Replace ENSG IDs with gene names
    new_features <- as.character(feature_map[current_features])
    # Rename features in assay data and metadata
    rownames(atlas_object@assays[["RNA"]]@counts) <- new_features
    rownames(atlas_object@assays[["RNA"]]@data) <- new_features
    if (!is.null(atlas_object@assays[["RNA"]]@scale.data) && nrow(atlas_object@assays[["RNA"]]@scale.data) > 0) {
      rownames(atlas_object@assays[["RNA"]]@scale.data) <- new_features
    }
    atlas_object@assays[["RNA"]]@meta.features <- atlas_object@assays[["RNA"]]@meta.features
    rownames(atlas_object@assays[["RNA"]]@meta.features) <- new_features
    
    rm(feature_map, current_features, new_features)
    # sum(rownames(atlas_object@assays[["RNA"]]@meta.features) %in% rownames(train_Six2GFP))
  }
  # remove irrelevent data:
  # table(atlas_object@meta.data[["cell_type"]][atlas_object@meta.data[["compartment"]] != "immune"])
  # atlas_object@meta.data[atlas_object@meta.data[["compartment"]] != "immune",]
  # atlas_object = subset(x = atlas_object, subset = compartment != "immune")
  cells_to_keep <- rownames(atlas_object@meta.data)[!atlas_object$cell_type %in% c("B cell", "neuron", "lymphocyte",
                                                                                   "natural killer cell","plasmacytoid dendritic cell",
                                                                                   "erythroid lineage cell", "CD4-positive, alpha-beta T cell",
                                                                                   "megakaryocyte", "mast cell", "CD8-positive, alpha-beta T cell",
                                                                                   "neutrophil", "conventional dendritic cell", "fibroblast")]
  atlas_object <- subset(atlas_object, cells = cells_to_keep)
  remove(cells_to_keep)
  # table(atlas_object@meta.data[["cell_type"]])
  
  ## run flow
  # atlas_object <- NormalizeData(atlas_object)
  # atlas_object <- FindVariableFeatures(atlas_object, selection.method = "vst", nfeatures = 2000)
  # atlas_object <- ScaleData(atlas_object)
  # atlas_object <- RunPCA(atlas_object)
  # atlas_object <- FindNeighbors(atlas_object, dims = 1:30)
  # atlas_object <- FindClusters(atlas_object)
  # atlas_object <- RunUMAP(atlas_object, dims = 1:30, reduction.name = "umap")
  # DimPlot(atlas_object, reduction = "umap", group.by = "author_cell_type", label = TRUE, label.size = 3,
  #         repel = TRUE)
  
} # kidney cell atlas - Figure 4
{
  Takasato_Humphreys_20 = read.table("GSM3763146_Takasato.dge.txt")
  load(file = "takasato_off_target_cell_list")
  seurat_Takasato_Humphreys_20 <- CreateSeuratObject(counts = Takasato_Humphreys_20[,!colnames(Takasato_Humphreys_20) %in% takasato_off_target_cell_list],
                                                     assay = "RNA")
  remove(Takasato_Humphreys_20)
} # load Takasato (seurat_Takasato_Humphreys_20)
{GSM2879360_Org1 <- read.table(gzfile("urine_based/GSM2879360_Org1_HMKJWBGXY_S1.gene.coutt.txt.gz"), header = TRUE, sep = "\t", row.names = 1)
  GSM2879361_Org2 <- read.table(gzfile("urine_based/GSM2879361_Org2_HMKJWBGXY_S2.gene.coutt.txt.gz"), header = TRUE, sep = "\t", row.names = 1)
  # rownames(GSM2879360_Org1)[grep("RNF5P1",rownames(GSM2879360_Org1))] # two with the sam name...
  # "RNF5P1" %in% rownames(fetal_train) # ...but not exisiting in fetal, so removing
  GSM2879360_Org1 = GSM2879360_Org1[-grep("RNF5P1",rownames(GSM2879360_Org1)), ]
  rownames(GSM2879360_Org1) = sapply(strsplit(rownames(GSM2879360_Org1),"__"), getElement, 1)
  rownames(GSM2879361_Org2) = sapply(strsplit(rownames(GSM2879361_Org2),"__"), getElement, 1)
  tubuloid = merge(GSM2879360_Org1,GSM2879361_Org2,by="row.names")
  row.names(tubuloid) = tubuloid$Row.names
  tubuloid = tubuloid[,-c(1)]
  remove(GSM2879360_Org1); remove(GSM2879361_Org2)
  tubuloid_seurat = CreateSeuratObject(counts = tubuloid)
  remove(tubuloid)
  } # Tubloid
{
  # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149134
  # GSM4648413
  expression_matrix <- Read10X(data.dir = "Nishinakamura")
  rownames(expression_matrix) = toupper(rownames(expression_matrix))
  nishinakamura_seurat = CreateSeuratObject(counts = expression_matrix)
  remove(expression_matrix)
} # Nishinakamura

Uchimura_full_flow = analyze_noise_impact_on_prediction(
  seurat_Uchimura_Humphreys_20,
  noised_number = 1,
  train_Six2GFP, # change to labeled_train_data
  test_Six2GFP, # change to labeled_test_data
  ref_cell_type_column = "type",
  dims = 1:30,
  train_title = "SIX2GFP",
  test_title_prefix = "Uchimura",
  n_neighbors = 8,
  skip_neighbors = TRUE,
  output_prefix_base = "six2gfp/26.1.26/Uchimura_noised/",
  prediction_column_name = "predicted.type", # Column name for predictions in metadata
  colors_feature_plot_noise = c('grey', '#f03b20'),
  myColors_cell_types = NULL, # Colors for cell type plots, if NULL, Paired palette will be used
  return_all_suerats = FALSE
)
# temp_seurat_obj <- readRDS("six2gfp/newest_function/Uchimura_noised/_cache/run_without_noise.rds") # debug
# main_cell_type_order <- temp_seurat_obj[["cell_type_order"]]
# noised_prediction = Uchimura_full_flow$noised_prediction_matrix

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
  output_prefix_base = "six2gfp/26.1.26/freedman_not_clean/",
  prediction_column_name = "predicted.type", # Column name for predictions in metadata
  colors_feature_plot_noise = c('grey', '#f03b20'),
  myColors_cell_types = NULL, # Colors for cell type plots, if NULL, Paired palette will be used
  return_all_suerats = FALSE
)

cell_atlas_flow = analyze_noise_impact_on_prediction(
  atlas_object,
  noised_number = 1,
  train_Six2GFP, 
  test_Six2GFP,
  ref_cell_type_column = "type",
  dims = 1:30,
  train_title = "SIX2GFP",
  test_title_prefix = "Kidney Cell Atlas",
  n_neighbors = 8,
  skip_neighbors = TRUE,
  output_prefix_base = "six2gfp/22.1.26/kidney_cell_atlas_clean_no_fibroblast/",
  prediction_column_name = "predicted.type", # Column name for predictions in metadata
  colors_feature_plot_noise = c('grey', '#f03b20'),
  myColors_cell_types = NULL, # Colors for cell type plots, if NULL, Paired palette will be used
  return_all_suerats = FALSE
)

###### Extra stuff ####

Takasato_full_flow = analyze_noise_impact_on_prediction(
  seurat_Takasato_Humphreys_20,
  noised_number = 1,
  train_Six2GFP, # change to labeled_train_data
  test_Six2GFP, # change to labeled_test_data
  ref_cell_type_column = "type",
  dims = 1:30,
  train_title = "SIX2GFP",
  test_title_prefix = "Takasato",
  n_neighbors = 8,
  skip_neighbors = TRUE,
  output_prefix_base = "six2gfp/22.1.26/Takasato_noised_tranpose/",
  prediction_column_name = "predicted.type", # Column name for predictions in metadata
  colors_feature_plot_noise = c('grey', '#f03b20'),
  myColors_cell_types = NULL, # Colors for cell type plots, if NULL, Paired palette will be used
  return_all_suerats = FALSE
)

Tubuloid_full_flow = analyze_noise_impact_on_prediction(
  tubuloid_seurat,
  noised_number = 2,
  train_Six2GFP, # change to labeled_train_data
  test_Six2GFP, # change to labeled_test_data
  ref_cell_type_column = "type",
  dims = 1:50, # with 30 there is an Error in FindWeights: Number of anchor cells is less than k.weight. Consider lowering k.weight to less than 47 or increase k.anchor
  train_title = "SIX2GFP",
  test_title_prefix = "Tubuloid",
  n_neighbors = 8,
  skip_neighbors = TRUE,
  output_prefix_base = "six2gfp/newest_function/Tubuloid_noised/",
  prediction_column_name = "predicted.type", # Column name for predictions in metadata
  colors_feature_plot_noise = c('grey', '#f03b20'),
  myColors_cell_types = NULL, # Colors for cell type plots, if NULL, Paired palette will be used
  return_all_suerats = FALSE
)

#six2gfp 50-50
{
load(file = "six2gfp/rawdata")
load(file = "six2gfp/metadata")
Six2GFP_metadata = Six2GFP_metadata[Six2GFP_metadata$Type != "LOWQUALITY",]
Six2GFP_metadata = Six2GFP_metadata[Six2GFP_metadata$Type != "ZERO_ACTB_GAPDH",]
Six2GFP = Six2GFP[,Six2GFP_metadata$CellName]
set.seed(1);
sample_vector <- sample(c(TRUE, FALSE), ncol(Six2GFP), replace=TRUE, prob=c(0.5,0.5))
train  <- colnames(Six2GFP)[sample_vector]
test   <- colnames(Six2GFP)[!sample_vector]

train_Six2GFP_50 <- CreateSeuratObject(counts = Six2GFP[,train], assay = "RNA", min.cells = 10)
train_Six2GFP_50 <- AddMetaData(train_Six2GFP_50, "train_Six2GFP", col.name = "source")
train_Six2GFP_50 <- AddMetaData(train_Six2GFP_50, Six2GFP_metadata[train,2], col.name = "type")
test_Six2GFP_50 <- subset(CreateSeuratObject(counts = Six2GFP, assay = "RNA", min.cells = 10), cells = test)
test_Six2GFP_50 <- AddMetaData(test_Six2GFP_50, "test_Six2GFP", col.name = "source")
remove(Six2GFP,Six2GFP_metadata,sample_vector)
}

test_full_flow = analyze_noise_impact_on_prediction(
  test_Six2GFP_50,
  noised_number = 2,
  train_Six2GFP_50, # change to labeled_train_data
  test_Six2GFP_50, # change to labeled_test_data
  ref_cell_type_column = "type",
  dims = 1:30,
  train_title = "SIX2GFP",
  test_title_prefix = "six2gfp_test",
  n_neighbors = 8,
  skip_neighbors = TRUE,
  output_prefix_base = "six2gfp/newest_function/six2gfp_test_noised_50/",
  prediction_column_name = "predicted.type", # Column name for predictions in metadata
  colors_feature_plot_noise = c('grey', '#f03b20'),
  myColors_cell_types = NULL, # Colors for cell type plots, if NULL, Paired palette will be used
  return_all_suerats = FALSE
)


nishinakamura_flow = analyze_noise_impact_on_prediction(
  nishinakamura_seurat,
  noised_number = 1,
  train_Six2GFP, 
  test_Six2GFP,
  ref_cell_type_column = "type",
  dims = 1:30,
  train_title = "SIX2GFP",
  test_title_prefix = "Nishinakamura",
  n_neighbors = 8,
  skip_neighbors = TRUE,
  output_prefix_base = "six2gfp/20.8.25/Nishinakamura_noised/",
  prediction_column_name = "predicted.type", # Column name for predictions in metadata
  colors_feature_plot_noise = c('grey', '#f03b20'),
  myColors_cell_types = NULL, # Colors for cell type plots, if NULL, Paired palette will be used
  return_all_suerats = FALSE
)

# fix titles so they will show alwayes the dataset name
# find another dataset of organoid (friedman or else, not Humphries) and run flow on it.
# freedman : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115986


