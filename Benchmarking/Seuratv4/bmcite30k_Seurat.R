library(Seurat)
library(SeuratData)
library(cowplot)
library(dplyr)
#InstallData("bmcite")
bm <- LoadData(ds = "bmcite")

#preprocessing RNA assay
DefaultAssay(bm) <- 'RNA'
bm <- NormalizeData(bm) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

ElbowPlot(bm, ndims=50)

#preprocessing ADT assay  
DefaultAssay(bm) <- 'ADT'
VariableFeatures(bm) <- rownames(bm[["ADT"]])
bm <- NormalizeData(bm, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')

ElbowPlot(bm, ndims=24, reduction = "apca")
  
#weighted nearest neighbour analysis
bm <- FindMultiModalNeighbors(
    bm, reduction.list = list("pca", "apca"), 
    dims.list = list(1:30, 1:18), modality.weight.name = list("RNA.weight", "ADT.weight")
)
  
#clustering
bm <- RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
bm <- FindClusters(bm, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)
#p1 <- DimPlot(bm, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
#p2 <- DimPlot(bm, reduction = 'wnn.umap', group.by = 'celltype.l2', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
#p1 + p2
  
DimPlot(bm, reduction = 'wnn.umap', group.by = 'donor', repel = TRUE, label.size = 2.5) #Check if batch correction is required
write.csv(bm@meta.data['seurat_clusters'], 'bmcite30k_labels_Seurat.csv')

  
