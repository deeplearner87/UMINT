library(Seurat)
library(SeuratData)
library(cowplot)
library(dplyr)
library(SeuratDisk)
#InstallData("cbmc")

#Loading data
cbmc <- LoadData(ds = "cbmc")
#saving groundTruth
write.csv(cbmc@meta.data['rna_annotations'], 'cbmc8k_groundTruth.csv')

#preprocessing RNA assay
DefaultAssay(cbmc) <- 'RNA'
cbmc <- NormalizeData(cbmc) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

#data preparation for UMINT, MOFA+
rna_data <- as.matrix(GetAssayData(cbmc, slot = "scale.data", assay="RNA"))
write.csv(rna_data, 'cbmc8k_rna_scaled.csv')

#data preparation for TotalVI
cbmc_rna <- CreateSeuratObject(counts=GetAssay(object = cbmc, assay = "RNA"))
SaveH5Seurat(cbmc_rna, filename = "cbmc8k_rna.h5Seurat")
Convert("cbmc8k_rna.h5Seurat", dest = "h5ad")

#preprocessing ADT assay
DefaultAssay(cbmc) <- 'ADT'
VariableFeatures(cbmc) <- rownames(cbmc[["ADT"]])
cbmc <- NormalizeData(cbmc, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')

#data preparation for UMINT, MOFA+
adt_data <- as.matrix(GetAssayData(cbmc, slot = "scale.data", assay="ADT"))
write.csv(adt_data, 'cbmc8k_adt_scaled.csv')

#data preparation for TotalVI
cbmc_adt <- CreateSeuratObject(counts=GetAssay(object = cbmc, assay = "ADT"))
SaveH5Seurat(cbmc_adt, filename = "cbmc8k_adt.h5Seurat")
Convert("cbmc8k_adt.h5Seurat", dest = "h5ad")
