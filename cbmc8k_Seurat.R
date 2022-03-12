library(Seurat)
library(SeuratData)
library(cowplot)
library(dplyr)
#InstallData("cbmc")
cbmc <- LoadData(ds = "cbmc")

#preprocessing RNA assay
DefaultAssay(cbmc) <- 'RNA'
cbmc <- NormalizeData(cbmc) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

ElbowPlot(cbmc, ndims=50)

#preprocessing ADT assay
DefaultAssay(cbmc) <- 'ADT'
VariableFeatures(cbmc) <- rownames(cbmc[["ADT"]])
cbmc <- NormalizeData(cbmc, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')

ElbowPlot(cbmc, ndims=9, reduction = "apca")

#Weighted nearest neighbour analysis
cbmc <- FindMultiModalNeighbors(
  cbmc, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:9), modality.weight.name = list("RNA.weight", "ADT.weight")
)

#clustering
cbmc <- RunUMAP(cbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
cbmc <- FindClusters(cbmc, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)

write.csv(cbmc@meta.data['seurat_clusters'], 'cbmc8k_labels_Seurat.csv')

