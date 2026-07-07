# another classified data as ground truth - SIX2GFP
# pre-process. can be loaded below 
install.packages("readxl")
library(readxl)
Six2GFP <- read_excel("six2gfp/SMART-Seq2_Six2mice_E18.5_300Six2pos_300Six2neg_cells_RAW+normalized.xlsx")
Six2GFP=data.frame(Six2GFP)
rownames(Six2GFP) = toupper(Six2GFP$...1) # mouse is low-case
Six2GFP = Six2GFP[,-c(1)]
Six2GFP_metadata = list()
for (sheet in c(1:12)) {
  tmp = read_excel("six2gfp/manually classified cell types.xlsx", sheet = sheet)
  Six2GFP_metadata[colnames(tmp)] = tmp
}; remove(sheet,tmp)
Six2GFP_metadata$CM = Six2GFP_metadata$CM_ALL[!Six2GFP_metadata$CM_ALL %in% Six2GFP_metadata$CM_DIV] # CM = CM_ALL-CM_DIV
Six2GFP_metadata = Six2GFP_metadata[-4] # remove CM_ALL
Six2GFP_metadata <- data.frame(CellName = unlist(Six2GFP_metadata, use.names = FALSE),
                                Type = rep(names(Six2GFP_metadata), sapply(Six2GFP_metadata, length)))
Six2GFP_metadata$CellName=gsub("_",".",Six2GFP_metadata$CellName)
rownames(Six2GFP_metadata) = Six2GFP_metadata$CellName
# save(Six2GFP, file = "six2gfp/rawdata")
# save(Six2GFP_metadata, file = "six2gfp/metadata")

load(file = "six2gfp/rawdata")
load(file = "six2gfp/metadata")
library(Seurat)

# 1st: run seurat on full data
Six2GFP_metadata = Six2GFP_metadata[Six2GFP_metadata$Type != "LOWQUALITY",]
Six2GFP_metadata = Six2GFP_metadata[Six2GFP_metadata$Type != "ZERO_ACTB_GAPDH",]
Six2GFP = Six2GFP[,Six2GFP_metadata$CellName]

six2gfp_seurat <- CreateSeuratObject(counts = Six2GFP, assay = "RNA", min.cells = 10)
# six2gfp_seurat <- subset(six2gfp_seurat, subset = Type != "LOWQUALITY")
# six2gfp_seurat <- AddMetaData(six2gfp_seurat, "six2gfp", col.name = "source")
six2gfp_seurat <- AddMetaData(six2gfp_seurat, Six2GFP_metadata$Type, col.name = "type")
six2gfp_seurat <- NormalizeData(six2gfp_seurat)
six2gfp_seurat <- FindVariableFeatures(six2gfp_seurat, selection.method = "vst", nfeatures = 2000)
six2gfp_seurat <- ScaleData(six2gfp_seurat)
six2gfp_seurat <- RunPCA(six2gfp_seurat)
ElbowPlot(six2gfp_seurat, ndims = 30)
six2gfp_seurat <- FindNeighbors(six2gfp_seurat, dims = 1:15)
six2gfp_seurat <- FindClusters(six2gfp_seurat)
six2gfp_seurat <- RunUMAP(six2gfp_seurat, dims = 1:15, return.model = TRUE)
save(six2gfp_seurat, file = "six2gfp/six2gfp_seurat_full")
DimPlot(six2gfp_seurat, reduction = "umap", group.by="type", label = TRUE, label.size = 3)
# train-test ####
# split to train-test
{set.seed(1);
sample_vector <- sample(c(TRUE, FALSE), ncol(Six2GFP), replace=TRUE, prob=c(0.8,0.2))
train  <- colnames(Six2GFP)[sample_vector]
test   <- colnames(Six2GFP)[!sample_vector]
remove(sample_vector)}

train_Six2GFP <- CreateSeuratObject(counts = Six2GFP[,train], assay = "RNA", min.cells = 10)
train_Six2GFP <- AddMetaData(train_Six2GFP, "train_Six2GFP", col.name = "source")
train_Six2GFP <- AddMetaData(train_Six2GFP, Six2GFP_metadata[train,2], col.name = "type")
test_Six2GFP <- subset(CreateSeuratObject(counts = Six2GFP, assay = "RNA", min.cells = 10), cells = test)
test_Six2GFP <- AddMetaData(test_Six2GFP, "test_Six2GFP", col.name = "source")
# test_Six2GFP <- AddMetaData(test_Six2GFP, Six2GFP_metadata[test,2], col.name = "type")
save(train_Six2GFP, file = "six2gfp/train_Six2GFP")
save(test_Six2GFP, file = "six2gfp/test_Six2GFP")


train_Six2GFP <- NormalizeData(train_Six2GFP)
train_Six2GFP <- FindVariableFeatures(train_Six2GFP, selection.method = "vst", nfeatures = 2000)
train_Six2GFP <- ScaleData(train_Six2GFP)
train_Six2GFP <- RunPCA(train_Six2GFP)
train_Six2GFP <- FindNeighbors(train_Six2GFP, dims = 1:30)
train_Six2GFP <- FindClusters(train_Six2GFP)
train_Six2GFP <- RunUMAP(train_Six2GFP, dims = 1:30, return.model = TRUE)

test_Six2GFP <- NormalizeData(test_Six2GFP)

test_Six2GFP.anchors <- FindTransferAnchors(reference = train_Six2GFP, query = test_Six2GFP, dims = 1:30,
                                          reference.reduction = "pca")
test_Six2GFP.query <- MapQuery(anchorset = test_Six2GFP.anchors, reference = train_Six2GFP, query = test_Six2GFP,
                             refdata = list(type = "type"), reference.reduction = "pca", reduction.model = "umap")
library(ggplot2)
library(RColorBrewer)
myColors <- brewer.pal(length(unique(Six2GFP_metadata[train,2])),"Paired"); names(myColors) <- unique(Six2GFP_metadata[train,2])
train_plot = DimPlot(train_Six2GFP, reduction = "umap", group.by = "type", label = TRUE, label.size = 3,
              repel = TRUE)  + ggtitle("Six2GFP (train) annotations") + scale_color_manual(values =  myColors)+
              theme(legend.position = "left", legend.text = element_text(size=8))+xlim(-10.5,7.1)+ylim(-8,10)

test_plot = DimPlot(test_Six2GFP.query, reduction = "ref.umap", group.by = "predicted.type", label = TRUE,
        label.size = 3, repel = TRUE) + ggtitle("Test - Transferred Labels")+scale_color_manual(values =  myColors)+NoLegend()
train_plot+test_plot
ggsave("six2gfp/train-test-regular.png",width = 1600, height = 900, units = "px",dpi=100)

# calculate entropy and distance mat
log_pred_test = log2(test_Six2GFP.query@assays[["prediction.score.type"]]@data)
log_pred_test[is.infinite(log_pred_test)] = 0
entropy_by_cell_test = -colSums(test_Six2GFP.query@assays[["prediction.score.type"]]@data*log_pred_test)
hist(entropy_by_cell_test)

colors_for_feature = c('grey','#f03b20')
test_Six2GFP.query <- AddMetaData(test_Six2GFP.query, entropy_by_cell_test, col.name = "entropy")
tets_entropy_plot = FeaturePlot(test_Six2GFP.query, reduction = "ref.umap", features = "entropy",
            cols = colors_for_feature, keep.scale="all") + ggtitle("Test - Entropy")
train_plot+test_plot+tets_entropy_plot
ggsave("six2gfp/train-test-entropy.png",width = 1600, height = 900, units = "px",dpi=100)

library(RANN)
train_data = train_Six2GFP@reductions[["pca"]]@cell.embeddings[,1:30]
test_data = test_Six2GFP.query@reductions[["ref.pca"]]@cell.embeddings # 30 PCs
all_neighbors_dist = data.frame(matrix(NA, nrow = nrow(test_data), ncol = 8))
all_neighbors_names = data.frame(matrix(NA, nrow = nrow(test_data), ncol = 8))

for (i in 1:nrow(test_data)){ # 2 minutes for 4100 rows!
  winner_type = test_Six2GFP.query@meta.data[["predicted.type"]][i]
  train_same_as_winner = (Six2GFP_metadata[train,2] == winner_type)
  training_relevant_data = train_data[train_same_as_winner,]
  # knn_closest = get.knn(rbind(test_data[i,], training_relevant_data), k=1)
  knn_closest = nn2(training_relevant_data, as.data.frame(t(test_data[i,])), k=8)
  all_neighbors_dist[i,] = knn_closest$nn.dist[1,]
  all_neighbors_names[i,] = rownames(training_relevant_data)[knn_closest$nn.idx[1,]]
  if(i%%50==0){print(paste0("Row:",i," Time:",Sys.time()))}
}
test_Six2GFP.query <- AddMetaData(test_Six2GFP.query, rowMeans(all_neighbors_dist), col.name = "closest_neighbor")
test_dist_plot = FeaturePlot(test_Six2GFP.query, reduction = "ref.umap", features = "closest_neighbor",
            cols = colors_for_feature, keep.scale="all") + ggtitle("Test - closest_neighbor")
train_plot+test_plot+test_dist_plot
ggsave("six2gfp/train-test-dist.png",width = 1600, height = 900, units = "px",dpi=100)

hist(rowMeans(all_neighbors_dist))
## Confusion Matrix ####
probs_tub=test_Six2GFP.query@assays[["prediction.score.type"]]@data
num_classes <- 10
pairwise_confusion <- matrix(0, nrow=num_classes, ncol=num_classes)
row.names(pairwise_confusion) = gsub("-","_",row.names(probs_tub))
colnames(pairwise_confusion) = row.names(pairwise_confusion)
for (i in 1:num_classes) {
  for (j in 1:num_classes) {
    if (i != j) {
      pairwise_confusion[i, j] <- mean(probs_tub[i, ] * probs_tub[j, ])
    }
  }
}
library(reshape2)
pairwise_confusion_df <- melt(pairwise_confusion)
ggplot(data = pairwise_confusion_df, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "Class", y = "Class", title = "Pairwise Confusion Matrix - Test") +
  theme_minimal()
ggsave("six2gfp/test_pairwise_confusion.jpg", width = 1200, height = 700, units = "px",dpi=100)


# UCHIMURA #######
#maybe needs better clean up?
Uchimura_Humphreys_20 = read.table("GSM3763147_Uchimura.dge.txt")
seurat_Uchimura_Humphreys_20 <- CreateSeuratObject(counts = Uchimura_Humphreys_20, assay = "RNA")
remove(Uchimura_Humphreys_20)

# seurat_Uchimura_Humphreys_20[["percent.mt"]] <- PercentageFeatureSet(seurat_Uchimura_Humphreys_20, pattern = "^MT-")
# VlnPlot(seurat_Uchimura_Humphreys_20, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter(seurat_Uchimura_Humphreys_20, feature1 = "nCount_RNA", feature2 = "percent.mt")
# FeatureScatter(seurat_Uchimura_Humphreys_20, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

seurat_Uchimura_Humphreys_20 <- AddMetaData(seurat_Uchimura_Humphreys_20, "Organoid", col.name = "type")
# norm
seurat_Uchimura_Humphreys_20 <- NormalizeData(seurat_Uchimura_Humphreys_20)

uchimura.anchorsSix2GFP <- FindTransferAnchors(reference = train_Six2GFP, query = seurat_Uchimura_Humphreys_20, dims = 1:30,
                                     reference.reduction = "pca")
uchimura.query <- MapQuery(anchorset = uchimura.anchorsSix2GFP, reference = train_Six2GFP, query = seurat_Uchimura_Humphreys_20,
                           refdata = list(type = "type"), reference.reduction = "pca", reduction.model = "umap")

uchimura_basic_plot = DimPlot(uchimura.query, reduction = "ref.umap", group.by = "predicted.type", label = TRUE,
                    label.size = 3, repel = TRUE) + ggtitle("Uchimura - Transferred Labels")+
                    scale_color_manual(values =  myColors)+NoLegend()+xlim(-10.5,7)+ylim(-8,10)

log_pred_uchimura = log2(uchimura.query@assays[["prediction.score.type"]]@data)
log_pred_uchimura[is.infinite(log_pred_uchimura)] = 0
entropy_by_cell_uchimura = -colSums(uchimura.query@assays[["prediction.score.type"]]@data*log_pred_uchimura)
hist(entropy_by_cell_uchimura)

organoid_data = uchimura.query@reductions[["ref.pca"]]@cell.embeddings # 30 PCs
all_neighbors_dist_uchimura = data.frame(matrix(NA, nrow = nrow(organoid_data), ncol = 8))
all_neighbors_names_uchimura = data.frame(matrix(NA, nrow = nrow(organoid_data), ncol = 8))

for (i in 1:nrow(organoid_data)){
  winner_type = uchimura.query@meta.data[["predicted.type"]][i]
  train_same_as_winner = (Six2GFP_metadata[train,2] == winner_type)
  training_relevant_data = train_data[train_same_as_winner,]
  knn_closest = nn2(training_relevant_data, as.data.frame(t(organoid_data[i,])), k=8)
  all_neighbors_dist_uchimura[i,] = knn_closest$nn.dist[1,]
  all_neighbors_names_uchimura[i,] = rownames(training_relevant_data)[knn_closest$nn.idx[1,]]
  if(i%%50==0){print(paste0("Row:",i," Time:",Sys.time()))}
}

uchimura.query <- AddMetaData(uchimura.query, entropy_by_cell_uchimura, col.name = "entropy")
uchimura.query <- AddMetaData(uchimura.query, rowMeans(all_neighbors_dist_uchimura), col.name = "closest_neighbor")

uchimura_neighbor_p <- FeaturePlot(uchimura.query, reduction = "ref.umap", features = "closest_neighbor", 
                                   cols = colors_for_feature) + ggtitle("Uchimura - closest neighbor")+xlim(-10.5,7)+ylim(-8,10)
uchimura_entropy_p <- FeaturePlot(uchimura.query, reduction = "ref.umap", features = "entropy",
                                  cols = colors_for_feature) + ggtitle("Uchimura - Entropy") +xlim(-10.5,7)+ylim(-8,10)
train_plot+uchimura_basic_plot+uchimura_entropy_p
ggsave("six2gfp/train-uchimura-entropy.png",width = 1600, height = 900, units = "px",dpi=100)
train_plot+uchimura_basic_plot+uchimura_neighbor_p
ggsave("six2gfp/train-uchimura-dist.png",width = 1600, height = 900, units = "px",dpi=100)

# VlnPlot(uchimura.query, features = c("closest_neighbor", "entropy"), ncol = 2)
violin_df = data.frame(predicted=c(as.factor(test_Six2GFP.query@meta.data[["predicted.type"]]), as.factor(uchimura.query@meta.data[["predicted.type"]])), 
           entropy=c(entropy_by_cell_test, entropy_by_cell_uchimura),
           distance=c(rowMeans(all_neighbors_dist),rowMeans(all_neighbors_dist_uchimura)),
           source=c(rep("Test",length(entropy_by_cell_test)),rep("Uchimura",length(entropy_by_cell_uchimura))))

levels(violin_df$predicted)
library(dplyr)
library(gridExtra)
violin_df %>%  group_by(predicted, source) %>% summarise(n=n(), avg = mean(entropy))-> violin_df.Summary
p = ggplot(violin_df,aes(x=source, y=entropy, fill=source))+ geom_violin(scale="area")+ xlab("ALL")+ ylim(0,2.5)
violin_plots_entropy = list(p)
p = ggplot(violin_df,aes(x=source, y=distance, fill=source))+ geom_violin(scale="area")+ xlab("ALL")
violin_plots_distance = list(p)
for (t in levels(violin_df$predicted)){
  numes = violin_df.Summary %>% filter(predicted == t) %>% pull (source, n)
  numes = paste(numes,names(numes),sep = ":")
  p = ggplot(violin_df %>% filter(predicted == t),
         aes(x=predicted, y=entropy, fill=source))+geom_violin(scale="width") +
          xlab(paste(numes[1],numes[2],sep = "|")) + ylim(0,2.5)
  violin_plots_entropy[[t]] = p
  p = ggplot(violin_df %>% filter(predicted == t),
             aes(x=predicted, y=distance, fill=source))+geom_violin(scale="width") + xlab(paste(numes[1],numes[2],sep = "|"))
  violin_plots_distance[[t]] = p
}
g=grid.arrange(grobs = violin_plots_entropy, ncol = 4)
ggsave("six2gfp/violin_plots_entropy_uchimura.png", g, width = 1600, height = 900, units = "px",dpi=100)
g=grid.arrange(grobs = violin_plots_distance, ncol = 4)
ggsave("six2gfp/violin_plots_distance_uchimura.png", g, width = 1600, height = 900, units = "px",dpi=100)

## Confusion Matrix ####
probs_tub=uchimura.query@assays[["prediction.score.type"]]@data
num_classes <- 10
pairwise_confusion <- matrix(0, nrow=num_classes, ncol=num_classes)
row.names(pairwise_confusion) = gsub("-","_",row.names(probs_tub))
colnames(pairwise_confusion) = row.names(pairwise_confusion)
for (i in 1:num_classes) {
  for (j in 1:num_classes) {
    if (i != j) {
      pairwise_confusion[i, j] <- mean(probs_tub[i, ] * probs_tub[j, ])
    }
  }
}
library(reshape2)
pairwise_confusion_df <- melt(pairwise_confusion)
ggplot(data = pairwise_confusion_df, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "Class", y = "Class", title = "Pairwise Confusion Matrix - Uchimura") +
  theme_minimal()
ggsave("six2gfp/uchimura_pairwise_confusion.jpg", width = 1200, height = 700, units = "px",dpi=100)


# TAKASATO #######
Takasato_Humphreys_20 = read.table("GSM3763146_Takasato.dge.txt")
seurat_Takasato_Humphreys_20 <- CreateSeuratObject(counts = Takasato_Humphreys_20, assay = "RNA")
remove(Takasato_Humphreys_20)

# seurat_Takasato_Humphreys_20[["percent.mt"]] <- PercentageFeatureSet(seurat_Takasato_Humphreys_20, pattern = "^MT-")
# VlnPlot(seurat_Takasato_Humphreys_20, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter(seurat_Takasato_Humphreys_20, feature1 = "nCount_RNA", feature2 = "percent.mt")
# FeatureScatter(seurat_Takasato_Humphreys_20, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

seurat_Takasato_Humphreys_20 <- AddMetaData(seurat_Takasato_Humphreys_20, "Organoid", col.name = "type")
# norm
seurat_Takasato_Humphreys_20 <- NormalizeData(seurat_Takasato_Humphreys_20)

Takasato.anchorsSix2GFP <- FindTransferAnchors(reference = train_Six2GFP, query = seurat_Takasato_Humphreys_20, dims = 1:30,
                                               reference.reduction = "pca")
Takasato.query <- MapQuery(anchorset = Takasato.anchorsSix2GFP, reference = train_Six2GFP, query = seurat_Takasato_Humphreys_20,
                           refdata = list(type = "type"), reference.reduction = "pca", reduction.model = "umap")

Takasato_basic_plot = DimPlot(Takasato.query, reduction = "ref.umap", group.by = "predicted.type", label = TRUE,
                              label.size = 3, repel = TRUE) + ggtitle("Takasato - Transferred Labels")+
  scale_color_manual(values =  myColors)+NoLegend()+xlim(-10.5,7)+ylim(-8,10)

log_pred_Takasato = log2(Takasato.query@assays[["prediction.score.type"]]@data)
log_pred_Takasato[is.infinite(log_pred_Takasato)] = 0
entropy_by_cell_Takasato = -colSums(Takasato.query@assays[["prediction.score.type"]]@data*log_pred_Takasato)
hist(entropy_by_cell_Takasato)

organoid_data = Takasato.query@reductions[["ref.pca"]]@cell.embeddings # 30 PCs
all_neighbors_dist_Takasato = data.frame(matrix(NA, nrow = nrow(organoid_data), ncol = 8))
all_neighbors_names_Takasato = data.frame(matrix(NA, nrow = nrow(organoid_data), ncol = 8))

for (i in 1:nrow(organoid_data)){
  winner_type = Takasato.query@meta.data[["predicted.type"]][i]
  train_same_as_winner = (Six2GFP_metadata[train,2] == winner_type)
  training_relevant_data = train_data[train_same_as_winner,]
  knn_closest = nn2(training_relevant_data, as.data.frame(t(organoid_data[i,])), k=8)
  all_neighbors_dist_Takasato[i,] = knn_closest$nn.dist[1,]
  all_neighbors_names_Takasato[i,] = rownames(training_relevant_data)[knn_closest$nn.idx[1,]]
  if(i%%50==0){print(paste0("Row:",i," Time:",Sys.time()))}
}

Takasato.query <- AddMetaData(Takasato.query, entropy_by_cell_Takasato, col.name = "entropy")
Takasato.query <- AddMetaData(Takasato.query, rowMeans(all_neighbors_dist_Takasato), col.name = "closest_neighbor")

Takasato_neighbor_p <- FeaturePlot(Takasato.query, reduction = "ref.umap", features = "closest_neighbor", 
                                   cols = colors_for_feature) + ggtitle("Takasato - closest neighbor")+xlim(-10.5,7)+ylim(-8,10)
Takasato_entropy_p <- FeaturePlot(Takasato.query, reduction = "ref.umap", features = "entropy",
                                  cols = colors_for_feature) + ggtitle("Takasato - Entropy") +xlim(-10.5,7)+ylim(-8,10)
train_plot+Takasato_basic_plot+Takasato_entropy_p
ggsave("six2gfp/train-Takasato-entropy.png",width = 1600, height = 900, units = "px",dpi=100)
train_plot+Takasato_basic_plot+Takasato_neighbor_p
ggsave("six2gfp/train-Takasato-dist.png",width = 1600, height = 900, units = "px",dpi=100)


violin_df_takastao = data.frame(predicted=c(as.factor(test_Six2GFP.query@meta.data[["predicted.type"]]), as.factor(Takasato.query@meta.data[["predicted.type"]])), 
                       entropy=c(entropy_by_cell_test, entropy_by_cell_Takasato),
                       distance=c(rowMeans(all_neighbors_dist),rowMeans(all_neighbors_dist_Takasato)),
                       source=factor(c(rep("Test",length(entropy_by_cell_test)),rep("Takasato",length(entropy_by_cell_Takasato))),levels=c("Test","Takasato")))
library(dplyr)
library(gridExtra)

violin_df_takastao %>%  group_by(predicted, source) %>% summarise(n=n(), avg = mean(entropy))-> violin_df_takastao.Summary
p = ggplot(violin_df_takastao,aes(x=source, y=entropy, fill=source))+ geom_violin(scale="area")+xlab("ALL")+ylim(-0.01,2.8)
violin_plots_entropy_takastao = list(p)
p = ggplot(violin_df_takastao,aes(x=source, y=distance, fill=source))+ geom_violin(scale="area")+ xlab("ALL")+ylim(-0.01,40)
violin_plots_distance_takastao = list(p)
for (t in levels(violin_df_takastao$predicted)){
  numes = violin_df_takastao.Summary %>% filter(predicted == t) %>% pull (source, n)
  numes = paste(numes,names(numes),sep = ":")
  p = ggplot(violin_df_takastao %>% filter(predicted == t),
             aes(x=predicted, y=entropy, fill=source))+geom_violin(scale="width") +
    xlab(paste(numes[1],numes[2],sep = "|")) + ylim(-0.01,2.8)
  violin_plots_entropy_takastao[[t]] = p
  p = ggplot(violin_df_takastao %>% filter(predicted == t),
             aes(x=predicted, y=distance, fill=source))+geom_violin(scale="width") +
    xlab(paste(numes[1],numes[2],sep = "|"))+ylim(-0.01,65)
  violin_plots_distance_takastao[[t]] = p
}
g = grid.arrange(grobs = violin_plots_entropy_takastao, ncol = 4)
ggsave("six2gfp/violin_plots_entropy_takastao.png", g, width = 1600, height = 900, units = "px",dpi=100)
g =grid.arrange(grobs = violin_plots_distance_takastao, ncol = 4)
ggsave("six2gfp/violin_plots_distance_takastao.png",g,width = 1600, height = 900, units = "px",dpi=100)

## Confusion Matrix ####
probs_tub=Takasato.query@assays[["prediction.score.type"]]@data
num_classes <- 10
pairwise_confusion <- matrix(0, nrow=num_classes, ncol=num_classes)
row.names(pairwise_confusion) = gsub("-","_",row.names(probs_tub))
colnames(pairwise_confusion) = row.names(pairwise_confusion)
for (i in 1:num_classes) {
  for (j in 1:num_classes) {
    if (i != j) {
      pairwise_confusion[i, j] <- mean(probs_tub[i, ] * probs_tub[j, ])
    }
  }
}
library(reshape2)
pairwise_confusion_df <- melt(pairwise_confusion)
ggplot(data = pairwise_confusion_df, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "Class", y = "Class", title = "Pairwise Confusion Matrix - Takasato") +
  theme_minimal()
ggsave("six2gfp/takasato_pairwise_confusion.jpg", width = 1200, height = 700, units = "px",dpi=100)

# WILMS #######
expression_matrix <- Read10X(data.dir = "filtered_feature_bc_matrix_wilms_tumor")
wilms_tumor_clean_seurat = CreateSeuratObject(counts = expression_matrix)
remove(expression_matrix)

wilms_tumor_clean_seurat <- NormalizeData(wilms_tumor_clean_seurat)
wilms_tumor_clean_seurat <- AddMetaData(wilms_tumor_clean_seurat, "Wilms", col.name = "type")

Wilms.anchorsSix2GFP <- FindTransferAnchors(reference = train_Six2GFP, query = wilms_tumor_clean_seurat, dims = 1:30,
                                               reference.reduction = "pca")
Wilms.query <- MapQuery(anchorset = Wilms.anchorsSix2GFP, reference = train_Six2GFP, query = wilms_tumor_clean_seurat,
                           refdata = list(type = "type"), reference.reduction = "pca", reduction.model = "umap")

Wilms_basic_plot = DimPlot(Wilms.query, reduction = "ref.umap", group.by = "predicted.type", label = TRUE,
                              label.size = 3, repel = TRUE) + ggtitle("Wilms - Transferred Labels")+
  scale_color_manual(values =  myColors)+NoLegend()+xlim(-10.5,7)+ylim(-8,10)

log_pred_Wilms = log2(Wilms.query@assays[["prediction.score.type"]]@data)
log_pred_Wilms[is.infinite(log_pred_Wilms)] = 0
entropy_by_cell_Wilms = -colSums(Wilms.query@assays[["prediction.score.type"]]@data*log_pred_Wilms)
hist(entropy_by_cell_Wilms)

organoid_data = Wilms.query@reductions[["ref.pca"]]@cell.embeddings # 30 PCs
all_neighbors_dist_Wilms = data.frame(matrix(NA, nrow = nrow(organoid_data), ncol = 8))
all_neighbors_names_Wilms = data.frame(matrix(NA, nrow = nrow(organoid_data), ncol = 8))

for (i in 1:nrow(organoid_data)){
  winner_type = Wilms.query@meta.data[["predicted.type"]][i]
  train_same_as_winner = (Six2GFP_metadata[train,2] == winner_type)
  training_relevant_data = train_data[train_same_as_winner,]
  knn_closest = nn2(training_relevant_data, as.data.frame(t(organoid_data[i,])), k=8)
  all_neighbors_dist_Wilms[i,] = knn_closest$nn.dist[1,]
  all_neighbors_names_Wilms[i,] = rownames(training_relevant_data)[knn_closest$nn.idx[1,]]
  if(i%%50==0){print(paste0("Row:",i," Time:",Sys.time()))}
}

Wilms.query <- AddMetaData(Wilms.query, entropy_by_cell_Wilms, col.name = "entropy")
Wilms.query <- AddMetaData(Wilms.query, rowMeans(all_neighbors_dist_Wilms), col.name = "closest_neighbor")

Wilms_neighbor_p <- FeaturePlot(Wilms.query, reduction = "ref.umap", features = "closest_neighbor", 
                                   cols = colors_for_feature) + ggtitle("Wilms - closest neighbor") +xlim(-10.5,7)+ylim(-8,10)
Wilms_entropy_p <- FeaturePlot(Wilms.query, reduction = "ref.umap", features = "entropy",
                                  cols = colors_for_feature) + ggtitle("Wilms - Entropy") +xlim(-10.5,7)+ylim(-8,10)
train_plot+Wilms_basic_plot+Wilms_entropy_p
ggsave("six2gfp/train-Wilms-entropy.png",width = 1600, height = 900, units = "px",dpi=100)
train_plot+Wilms_basic_plot+Wilms_neighbor_p
ggsave("six2gfp/train-Wilms-dist.png",width = 1600, height = 900, units = "px",dpi=100)


# violin plots: 
violin_df_wilms = data.frame(predicted=c(as.factor(test_Six2GFP.query@meta.data[["predicted.type"]]), as.factor(Wilms.query@meta.data[["predicted.type"]])), 
                       entropy=c(entropy_by_cell_test, entropy_by_cell_Wilms),
                       distance=c(rowMeans(all_neighbors_dist),rowMeans(all_neighbors_dist_Wilms)),
                       source=c(rep("Test",length(entropy_by_cell_test)),rep("Wilms",length(entropy_by_cell_Wilms))))

levels(violin_df_wilms$predicted)
library(dplyr)
library(gridExtra)
violin_df_wilms %>%  group_by(predicted, source) %>% summarise(n=n(), avg = mean(entropy))-> violin_df_wilms.Summary
p = ggplot(violin_df_wilms,aes(x=source, y=entropy, fill=source))+ geom_violin(scale="area")+ xlab("ALL")+ylim(-0.01,2.8)
violin_plots_entropy_wilms = list(p)
p = ggplot(violin_df_wilms,aes(x=source, y=distance, fill=source))+ geom_violin(scale="area")+ xlab("ALL")
violin_plots_distance_wilms = list(p)
for (t in levels(violin_df_wilms$predicted)){
  numes = violin_df_wilms.Summary %>% filter(predicted == t) %>% pull (source, n)
  numes = paste(numes,names(numes),sep = ":")
  p = ggplot(violin_df_wilms %>% filter(predicted == t),
             aes(x=predicted, y=entropy, fill=source))+geom_violin(scale="width") +
    xlab(paste(numes[1],numes[2],sep = "|"))+ylim(-0.01,2.8)
  violin_plots_entropy_wilms[[t]] = p
  p = ggplot(violin_df_wilms %>% filter(predicted == t),
             aes(x=predicted, y=distance, fill=source))+geom_violin(scale="width") + xlab(paste(numes[1],numes[2],sep = "|"))
  violin_plots_distance_wilms[[t]] = p
}
g =grid.arrange(grobs = violin_plots_entropy_wilms, ncol = 4)
ggsave("six2gfp/violin_plots_entropy_wilms.png",g,width = 1600, height = 900, units = "px",dpi=100)
g =grid.arrange(grobs = violin_plots_distance_wilms, ncol = 4)
ggsave("six2gfp/violin_plots_distance_wilms.png",g,width = 1600, height = 900, units = "px",dpi=100)

## Confusion Matrix ####
probs_tub=Wilms.query@assays[["prediction.score.type"]]@data
num_classes <- 10
pairwise_confusion <- matrix(0, nrow=num_classes, ncol=num_classes)
row.names(pairwise_confusion) = gsub("-","_",row.names(probs_tub))
colnames(pairwise_confusion) = row.names(pairwise_confusion)
for (i in 1:num_classes) {
  for (j in 1:num_classes) {
    if (i != j) {
      pairwise_confusion[i, j] <- mean(probs_tub[i, ] * probs_tub[j, ])
    }
  }
}
library(reshape2)
pairwise_confusion_df <- melt(pairwise_confusion)
ggplot(data = pairwise_confusion_df, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "Class", y = "Class", title = "Pairwise Confusion Matrix - Wilms") +
  theme_minimal()
ggsave("six2gfp/wilms_pairwise_confusion.jpg", width = 1200, height = 700, units = "px",dpi=100)

# Tubuloid ####
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
remove(tubuloid)}
sum(rownames(tubuloid_seurat) %in% rownames(six2gfp_seurat))

tubuloid_seurat <- NormalizeData(tubuloid_seurat)
tubuloid_seurat <- AddMetaData(tubuloid_seurat, "tubuloid", col.name = "type")
tubuloid_seurat <- AddMetaData(tubuloid_seurat, c(rep("Org1",96),rep("Org2",96)), col.name = "original_org")

tubuloid.anchorset <- FindTransferAnchors(reference = train_Six2GFP, query = tubuloid_seurat, dims = 1:8,
                                          reference.reduction = "pca", k.filter=NA, k.anchor=10, k.score=10)

tubuloid.query <- MapQuery(anchorset = tubuloid.anchorset, reference = train_Six2GFP, query = tubuloid_seurat,
                           refdata = list(type = "type"), reference.reduction = "pca",
                           reduction.model = "umap")

tubuloid_basic_plot = DimPlot(tubuloid.query, reduction = "ref.umap", group.by = "predicted.type", label = TRUE,
                           label.size = 3, repel = TRUE) + ggtitle("Tubuloid - Transferred Labels")+
  scale_color_manual(values =  myColors)+NoLegend()+xlim(-10.5,7)+ylim(-8,10)


log_pred_tubuloid = log2(tubuloid.query@assays[["prediction.score.type"]]@data)
log_pred_tubuloid[is.infinite(log_pred_tubuloid)] = 0
entropy_by_cell_tubuloid = -colSums(tubuloid.query@assays[["prediction.score.type"]]@data*log_pred_tubuloid)
hist(entropy_by_cell_tubuloid)
tubuloid.query <- AddMetaData(tubuloid.query, entropy_by_cell_tubuloid, col.name = "entropy")

tubuloid_entropy_p <- FeaturePlot(tubuloid.query, reduction = "ref.umap", features = "entropy", min.cutoff=0,
                                  cols = colors_for_feature, keep.scale="all") + ggtitle("tubuloid - Entropy") +xlim(-10.5,7)+ylim(-8,10)

# dist
train_data = fetal_train@reductions[["pca"]]@cell.embeddings[,1:8]
tubuloid_data = tubuloid.query@reductions[["ref.pca"]]@cell.embeddings # 8 PCs
all_neighbors_dist_tubuloid = data.frame(matrix(NA, nrow = nrow(tubuloid_data), ncol = 10))
all_neighbors_names_tubuloid = data.frame(matrix(NA, nrow = nrow(tubuloid_data), ncol = 10))
for (i in 1:nrow(tubuloid_data)){ 
  winner_type = tubuloid.query@meta.data[["predicted.type"]][i]
  train_same_as_winner = (Six2GFP_metadata[train,2] == winner_type)
  training_relevant_data = train_data[train_same_as_winner,]
  knn_closest = nn2(training_relevant_data, as.data.frame(t(tubuloid_data[i,])), k=10)
  all_neighbors_dist_tubuloid[i,] = knn_closest$nn.dist[1,]
  all_neighbors_names_tubuloid[i,] = rownames(training_relevant_data)[knn_closest$nn.idx[1,]]
  if(i%%50==0){print(paste0("Row:",i," Time:",Sys.time()))}
}
tubuloid.query <- AddMetaData(tubuloid.query, rowMeans(all_neighbors_dist_tubuloid), col.name = "closest_neighbor")
tubuloid_neighbor_p <- FeaturePlot(tubuloid.query, reduction = "ref.umap", features = "closest_neighbor", 
                                   cols = colors_for_feature) + 
  ggtitle("tubuloid - closest neighbor")+xlim(-10.5,7)+ylim(-8,10)

train_plot+tubuloid_basic_plot+tubuloid_entropy_p
ggsave("six2gfp/train-tubuloid-entropy.png",width = 1600, height = 900, units = "px",dpi=100)
train_plot+tubuloid_basic_plot+tubuloid_neighbor_p
ggsave("six2gfp/train-tubuloid-dist.png",width = 1600, height = 900, units = "px",dpi=100)

violin_df_tubuloid = data.frame(predicted=c(as.factor(test_Six2GFP.query@meta.data[["predicted.type"]]), as.factor(tubuloid.query@meta.data[["predicted.type"]])), 
                                entropy=c(entropy_by_cell_test, entropy_by_cell_tubuloid),
                                distance=c(rowMeans(all_neighbors_dist),rowMeans(all_neighbors_dist_tubuloid)),
                                source=factor(c(rep("Test",length(entropy_by_cell_test)),rep("Tubuloid",length(entropy_by_cell_tubuloid))),levels=c("Test","Tubuloid")))
violin_df_tubuloid %>%  group_by(predicted, source) %>% summarise(n=n(), avg = mean(entropy))-> violin_df_tubuloid.Summary
p = ggplot(violin_df_tubuloid,aes(x=source, y=entropy, fill=source))+ geom_violin(scale="area")+xlab("ALL")+ylim(-0.01,2)
violin_plots_entropy_tubuloid = list(p)
p = ggplot(violin_df_tubuloid,aes(x=source, y=distance, fill=source))+ geom_violin(scale="area")+ xlab("ALL")+ylim(-0.01,40)
violin_plots_distance_tubuloid = list(p)
for (t in levels(violin_df_tubuloid$predicted)){
  numes = violin_df_tubuloid.Summary %>% filter(predicted == t) %>% pull (source, n)
  numes = paste(numes,names(numes),sep = ":")
  p = ggplot(violin_df_tubuloid %>% filter(predicted == t),
             aes(x=predicted, y=entropy, fill=source))+geom_violin(scale="width") +
    xlab(paste(numes[1],numes[2],sep = "|")) + ylim(-0.01,2)
  violin_plots_entropy_tubuloid[[t]] = p
  p = ggplot(violin_df_tubuloid %>% filter(predicted == t),
             aes(x=predicted, y=distance, fill=source))+geom_violin(scale="width") +
    xlab(paste(numes[1],numes[2],sep = "|"))+ylim(-0.01,40)
  violin_plots_distance_tubuloid[[t]] = p
}
g = grid.arrange(grobs = violin_plots_entropy_tubuloid, ncol = 4)
ggsave("six2gfp/violin_plots_entropy_tubuloid.png", g, width = 1600, height = 900, units = "px",dpi=100)
g=grid.arrange(grobs = violin_plots_distance_tubuloid, ncol = 4)
ggsave("six2gfp/violin_plots_distance_tubuloid.png", g, width = 1600, height = 900, units = "px",dpi=100)

# PCA
tubuloid_PCA = DimPlot(tubuloid.query, reduction = "ref.pca", dims=c(1,2), group.by = "predicted.type", label = TRUE,
        label.size = 3, repel = TRUE) + ggtitle("Tubuloid - Transferred Labels")+
  scale_color_manual(values =  myColors)+NoLegend() + ylim(-75, 25) +xlim(-20,45)
DimPlot(tubuloid.query, reduction = "ref.pca", dims=c(1,2), group.by = "ori,inal_org", label = TRUE,
        label.size = 3, repel = TRUE) + ggtitle("Tubuloid - Transferred Labels")
trianing_pca_1_2 =DimPlot(train_Six2GFP, reduction = "pca", dims=c(1,2), group.by = "type", label = TRUE,
        label.size = 3, repel = TRUE) + ggtitle("train - Labels")+ ylim(-75, 25) +xlim(-20,45) + scale_color_manual(values =  myColors) + 
        theme(legend.position = "left", legend.text = element_text(size=8))

tubuloid_PCA_entropy= FeaturePlot(tubuloid.query, reduction = "ref.pca", dims = c(1,2), features = "entropy",
            cols = colors_for_feature) + ggtitle("tubuloid - Entropy") # + ylim(-75, 25) +xlim(-20,45)
tubuloid_PCA+tubuloid_PCA_entropy
ggsave("six2gfp/tubuloid_PCA.png", width = 1200, height = 700, units = "px",dpi=100)

trianing_pca_1_2+tubuloid_PCA+tubuloid_PCA_entropy
ggsave("six2gfp/tubuloid_PCA_entropy.png", width = 1200, height = 700, units = "px",dpi=100)

## Confusion Matrix ####
probs_tub=tubuloid.query@assays[["prediction.score.type"]]@data
# prediction_confusion <- matrix(0, nrow=10, ncol=10)
# row.names(prediction_confusion) = gsub("-","_",row.names(probs_tub))
# colnames(prediction_confusion) = gsub("-","_",row.names(probs_tub))
# # rows - predicted, columns - mean of scores
# for (prediction in row.names(probs_tub)) {
#   prediction = gsub("-","_",prediction)
#   prediction_confusion[prediction,] = rowMeans(probs_tub[,tubuloid.query@meta.data[["predicted.type"]]==prediction])
# }
# summary(t(probs_tub))

num_classes <- 10
pairwise_confusion <- matrix(0, nrow=num_classes, ncol=num_classes)
row.names(pairwise_confusion) = gsub("-","_",row.names(probs_tub))
colnames(pairwise_confusion) = row.names(pairwise_confusion)
for (i in 1:num_classes) {
  for (j in 1:num_classes) {
    if (i != j) {
      pairwise_confusion[i, j] <- mean(probs_tub[i, ] * probs_tub[j, ])
    }
  }
}
library(reshape2)
pairwise_confusion_df <- melt(pairwise_confusion)
ggplot(data = pairwise_confusion_df, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "Class", y = "Class", title = "Pairwise Confusion Matrix - Tubuloid") +
  theme_minimal()
ggsave("six2gfp/tubuloid_pairwise_confusion.jpg", width = 1200, height = 700, units = "px",dpi=100)
