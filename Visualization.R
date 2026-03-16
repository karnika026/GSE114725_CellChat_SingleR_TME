# 
# Generate CellChat communication network visualizations

# -------------------------------
# Load libraries
# -------------------------------
library(CellChat)
library(dplyr)
library(ggplot2)
library(patchwork)

# -------------------------------
# Load analyzed CellChat object
# -------------------------------
cellchat_seurat_data <- readRDS("cellchat/cellchat_seurat_data.rds")

# -------------------------------
# Inspect network
# -------------------------------
names(cellchat_seurat_data@net)
dim(cellchat_seurat_data@net$count)

# -------------------------------
# Define group size
# -------------------------------
groupSize <- as.numeric(table(cellchat_seurat_data@idents))

# -------------------------------
# 1. Global interaction network
# -------------------------------
png(
  "PLOTS/cellchat_circle_plot.png",
  width = 3000,
  height = 3000,
  res = 300
)

netVisual_circle(
  cellchat_seurat_data@net$count,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE
)

dev.off()

# -------------------------------
# 2. Interaction heatmap
# -------------------------------
png(
  "PLOTS/cellchat_heatmap_count.png",
  width = 3000,
  height = 2000,
  res = 300
)

p1 <- netVisual_heatmap(
  cellchatseurat_data,
  measure = "count"
)

print(p1)
dev.off()

# -------------------------------
# 3. Compute centrality scores
# -------------------------------
cellchat_seurat_data <- netAnalysis_computeCentrality(cellchat_seurat_data)

# -------------------------------
# 4. Outgoing signaling heatmap
# -------------------------------
png(
  "PLOTS/cellchat_outgoing_heatmap.png",
  width = 3000,
  height = 2000,
  res = 300
)

p2 <- netAnalysis_signalingRole_heatmap(
  cellchat_seurat_data,
  pattern = "outgoing"
)

print(p2)
dev.off()

# -------------------------------
# 5. Incoming signaling heatmap
# -------------------------------
png(
  "PLOTS/cellchat_incoming_heatmap.png",
  width = 3000,
  height = 2000,
  res = 300
)

p3 <- netAnalysis_signalingRole_heatmap(
  cellchat_seurat_data,
  pattern = "incoming"
)

print(p3)
dev.off()

# -------------------------------
# 6. Example signaling pathway
# -------------------------------
pathways.show <- head(cellchat_seurat_data@netP$pathways, 3)

png(
  "PLOTS/cellchat_pathway_example.png",
  width = 3000,
  height = 3000,
  res = 300
)

netVisual_aggregate(
  cellchat_seurat_data,
  signaling = pathways.show[1],
  layout = "circle"
)

dev.off()

# -------------------------------
# Save final object
# -------------------------------
saveRDS(
  cellchat,
  "cellchat/cellchat_analysis_final.rds"
)


sessionInfo()
