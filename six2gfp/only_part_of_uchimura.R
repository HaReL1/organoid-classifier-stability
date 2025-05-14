six2gfp_seurat_part = subset(uchimura.query, subset=predicted.type %in% c("CM","UM","PODO","PROX_1","PROX_2"))

DimPlot(six2gfp_seurat_part, reduction = "ref.umap", group.by="predicted.type", label = TRUE, label.size = 3)+scale_color_manual(values =  myColors)+NoLegend()
DimPlot(uchimura.query, reduction = "ref.umap", group.by="predicted.type", label = TRUE, label.size = 3)+scale_color_manual(values =  myColors)+NoLegend()

selected_uchimura_pca = DimPlot(six2gfp_seurat_part, reduction = "ref.pca", dims = c(1,3), group.by="predicted.type", label = TRUE, label.size = 3)+scale_color_manual(values =  myColors)+NoLegend()
selected_uchimura_pca_entropy =FeaturePlot(six2gfp_seurat_part, reduction = "ref.pca", dims = c(1,3), features = "entropy",
            cols = colors_for_feature) + ggtitle("Uchimura - Entropy")
selected_uchimura_pca+selected_uchimura_pca_entropy
ggsave("six2gfp/Uchimura_CM_UM_Podo_PROXs_PCA.png", width = 1200, height = 700, units = "px",dpi=100)


# FeaturePlot(uchimura.query, reduction = "ref.pca", features = "entropy",
#             cols = colors_for_feature) + ggtitle("Uchimura - Entropy")
# 
# DimPlot(six2gfp_seurat_part, reduction = "ref.pca", 
#         dims = c(1, 2), group.by="predicted.type", 
#         label = TRUE, label.size = 3)+scale_color_manual(values =  myColors)+NoLegend()
# 
# DimPlot(uchimura.query, reduction = "ref.pca", group.by="predicted.type", label = TRUE, 
#         label.size = 3)+scale_color_manual(values =  myColors)+NoLegend()


