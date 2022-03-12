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

#data preparation for UMINT, MOFA+
rna_data <- as.matrix(GetAssayData(pbmc, slot = "scale.data", assay="RNA"))
write.csv(rna_data, 'pbmc10k_MALT_rna_scaled.csv')

#data preparation for TotalVI
pbmc_rna <- CreateSeuratObject(counts=GetAssay(object = pbmc, assay = "RNA"))
SaveH5Seurat(pbmc_rna, filename = "pbmc10k_MALT_rna.h5Seurat")
Convert("pbmc10k_MALT_rna.h5Seurat", dest = "h5ad")

#preprocessing ADT assay
DefaultAssay(pbmc) <- 'ADT'
VariableFeatures(pbmc) <- rownames(pbmc[["ADT"]])
pbmc <- NormalizeData(pbmc, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')

#data preparation for UMINT, MOFA+
adt_data <- as.matrix(GetAssayData(pbmc, slot = "scale.data", assay="ADT"))
write.csv(adt_data, 'pbmc10k_MALT_adt_scaled.csv')

#data preparation for TotalVI
pbmc_adt <- CreateSeuratObject(counts=GetAssay(object = pbmc, assay = "ADT"))
SaveH5Seurat(pbmc_adt, filename = "pbmc10k_MALT_adt.h5Seurat")
Convert("pbmc10k_MALT_adt.h5Seurat", dest = "h5ad")
