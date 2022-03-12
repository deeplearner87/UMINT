library(Seurat)
library(SeuratData)
library(cowplot)
library(dplyr)
library(SeuratDisk)
#InstallData("bmcite")

#Loading data
bm <- LoadData(ds = "bmcite")

#saving groundTruth
write.csv(bm@meta.data['celltype.l2'], 'bmcite30k_groundTruth.csv')

#preprocessing RNA assay
DefaultAssay(bm) <- 'RNA'
bm <- NormalizeData(bm) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

#data preparation for UMINT, MOFA+
rna_data <- as.matrix(GetAssayData(bm, slot = "scale.data", assay="RNA"))
write.csv(rna_data, 'bmcite30k_rna_scaled.csv')

#data preparation for TotalVI
bm_rna <- CreateSeuratObject(counts=GetAssay(object = bm, assay = "RNA"))
SaveH5Seurat(bm_rna, filename = "bmcite30k_rna.h5Seurat")
Convert("bmcite30k_rna.h5Seurat", dest = "h5ad")  

#preprocessing ADT assay  
DefaultAssay(bm) <- 'ADT'
VariableFeatures(bm) <- rownames(bm[["ADT"]])
bm <- NormalizeData(bm, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')

#data preparation for UMINT, MOFA+
adt_data <- as.matrix(GetAssayData(bm, slot = "scale.data", assay="ADT"))
write.csv(adt_data, 'bmcite30k_adt_scaled.csv')
  
#data preparation for TotalVI
bm_adt <- CreateSeuratObject(counts=GetAssay(object = bm, assay = "ADT"))
SaveH5Seurat(bm_adt, filename = "bmcite30k_adt.h5Seurat")
Convert("bmcite30k_adt.h5Seurat", dest = "h5ad")
