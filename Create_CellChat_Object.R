Aim: To construct a CellChat object from the processed Seurat single-cell RNA-seq dataset using SingleR-based cell type annotations.

library(CellChat)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)
library(patchwork)

#Load Seurat Object
seurat_data <- readRDS("....../Github project/GSE114725/seurat_data.rds")

# Checking required MetaData
colnames(seurat_data@meta.data)

#Setting Cell Identities
Idents(seurat_data) <- "SingleR_HCA"
levels(seurat_data)

#Creating CellChat Object
cellchat_seurat <- createCellChat(object = seurat_data, assay = "RNA", group.by = "SingleR_HCA")
cellchat_seurat
levels(cellchat_seurat@idents)
table(cellchat_seurat@idents)

#Adding Meta Data to CellChat object
cellchat_seurat@meta$tissue <- seurat_data$tissue

#Setting up ligand-receptor database
CellChatDB <- CellChatDB.human 
cellchat_seurat@DB <- CellChatDB

#Making subset expression data
cellchat_seurat <- subsetData(cellchat_seurat)
saveRDS(cellchat_seurat, "....../Github project/GSE114725/cellchat_seurat.rds")
