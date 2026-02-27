# ─────────────────────────────────────────────────────────────
#         CellChat Cell-Cell Communication Analysis
# Input:  data/GSE114725/cellchat_seurat.rds
#         data/GSE114725/seurat_data.rds
# Output: data/GSE114725/cellchat_seurat_data.rds
#         plots/cellchat_circle.png
#         plots/cellchat_pathway_heatmap.png
#         plots/cellchat_pathway_circle.png
# ─────────────────────────────────────────────────────────────
library(CellChat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(Seurat)

base_dir <- "data/GSE114725"
plot_dir <- "plots"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

#  1. Load data 
cellchat_seurat <- readRDS(file.path(base_dir, "cellchat_seurat.rds"))


#  2. Run CellChat pipeline
cellchat_seurat@DB  <- CellChatDB.human
cellchat_seurat     <- subsetData(cellchat_seurat)
cellchat_seurat     <- identifyOverExpressedGenes(cellchat_seurat)
cellchat_seurat     <- identifyOverExpressedInteractions(cellchat_seurat)
cellchat_seurat     <- computeCommunProb(cellchat_seurat)
cellchat_seurat     <- filterCommunication(cellchat_seurat, min.cells = 10)
cellchat_seurat     <- computeCommunProbPathway(cellchat_seurat)
cellchat_seurat     <- aggregateNet(cellchat_seurat)

saveRDS(cellchat_seurat, file.path(base_dir, "cellchat_seurat_data.rds"))

#  3. Circle plot
groupSize <- as.numeric(table(cellchat_seurat@idents))

png(file.path(plot_dir, "cellchat_circle.png"), width = 2000, height = 2000, res = 300)
netVisual_circle(cellchat_seurat@net$count, vertex.weight = groupSize,
                 weight.scale = TRUE, label.edge = FALSE)
dev.off()

# Heatmap
png(file.path(plot_dir, "cellchat_pathway_heatmap.png"), width = 3000, height = 2000, res = 300)
netVisual_heatmap(cellchat_seurat, measure = "weight")
dev.off()

#Pathway Circle Plot
pathways.show <- cellchat_seurat@netP$pathways

png(file.path(plot_dir, "cellchat_pathway_circle.png"), width = 2000, height = 2000, res = 300)
netVisual_aggregate(cellchat_seurat, signaling = pathways.show[1], layout = "circle")
dev.off()
```


