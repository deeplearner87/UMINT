library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
# rename features.tsv to genes.tsv
pbmc.data <- Read10X(data.dir = "pbmc10K_MALT_filtered_feature_bc_matrix/")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data$`Gene Expression`)
pbmc

# create a new assay to store ADT information
adt_assay <- CreateAssayObject(counts = pbmc.data$`Antibody Capture`)

# add this assay to the previously created Seurat object
pbmc[["ADT"]] <- adt_assay

#preprocessing RNA assay
DefaultAssay(pbmc) <- 'RNA'
pbmc <- NormalizeData(pbmc) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

ElbowPlot(pbmc, ndims=50)

#preprocessing ADT assay
DefaultAssay(pbmc) <- 'ADT'
VariableFeatures(pbmc) <- rownames(pbmc[["ADT"]])
pbmc <- NormalizeData(pbmc, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')

ElbowPlot(pbmc, ndims=16, reduction = "apca")

#weighted nearest neighbour analysis
pbmc <- FindMultiModalNeighbors(
  pbmc, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:16), modality.weight.name = list("RNA.weight", "ADT.weight")
)

#clustering
pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc <- FindClusters(pbmc, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)


write.csv(pbmc@meta.data['seurat_clusters'], 'pbmc10k_MALT_labels_Seurat.csv')
