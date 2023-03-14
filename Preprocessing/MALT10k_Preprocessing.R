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

#data preparation for UMINT, MOFA+
rna_data <- as.matrix(GetAssayData(malt, slot = "scale.data", assay="RNA"))
write.csv(rna_data, 'MALT10k_rna_scaled.csv')

#data preparation for TotalVI
malt_rna <- CreateSeuratObject(counts=GetAssay(object = malt, assay = "RNA"))
SaveH5Seurat(malt_rna, filename = "MALT10k_rna.h5Seurat")
Convert("MALT10k_rna.h5Seurat", dest = "h5ad")

#preprocessing ADT assay
DefaultAssay(malt) <- 'ADT'
VariableFeatures(malt) <- rownames(malt[["ADT"]])
malt <- NormalizeData(malt, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')

#data preparation for UMINT, MOFA+
adt_data <- as.matrix(GetAssayData(malt, slot = "scale.data", assay="ADT"))
write.csv(adt_data, 'MALT10k_adt_scaled.csv')

#data preparation for TotalVI
malt_adt <- CreateSeuratObject(counts=GetAssay(object = malt, assay = "ADT"))
SaveH5Seurat(malt_adt, filename = "MALT10k_adt.h5Seurat")
Convert("MALT10k_adt.h5Seurat", dest = "h5ad")
