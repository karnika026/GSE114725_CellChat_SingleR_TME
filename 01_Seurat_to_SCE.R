# Aim: To convert preprocessed Seurat object to SingleCellExperiment
# This is required for downstream SingleR automated annotation and CellChat analysis

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")

library(Seurat)
library(SingleCellExperiment)

seurat_data <- readRDS("C:..../Github project/GSE114725/seurat_data.rds")
glimpse(seurat_data)
sce_seurat_data <- as.SingleCellExperiment(seurat_data)
sce_seurat_data
