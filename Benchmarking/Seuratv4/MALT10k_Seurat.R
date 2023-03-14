library(dplyr)
library(Seurat)
library(patchwork)

# Load the malt dataset
# rename features.tsv to genes.tsv
malt.data <- Read10X(data.dir = "MALT10k_filtered_feature_bc_matrix/")

# Initialize the Seurat object with the raw (non-normalized data).
malt <- CreateSeuratObject(counts = malt.data$`Gene Expression`)
malt

# create a new assay to store ADT information
adt_assay <- CreateAssayObject(counts = malt.data$`Antibody Capture`)

# add this assay to the previously created Seurat object
malt[["ADT"]] <- adt_assay

#preprocessing RNA assay
DefaultAssay(malt) <- 'RNA'
malt <- NormalizeData(malt) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

ElbowPlot(malt, ndims=50)

#preprocessing ADT assay
DefaultAssay(malt) <- 'ADT'
VariableFeatures(malt) <- rownames(malt[["ADT"]])
malt <- NormalizeData(malt, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')

ElbowPlot(malt, ndims=16, reduction = "apca")

#weighted nearest neighbour analysis
malt <- FindMultiModalNeighbors(
  malt, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:16), modality.weight.name = list("RNA.weight", "ADT.weight")
)

#clustering
malt <- RunUMAP(malt, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
malt <- FindClusters(malt, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)


write.csv(malt@meta.data['seurat_clusters'], 'MALT10k_labels_Seurat.csv')
