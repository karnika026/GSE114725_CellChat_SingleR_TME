
#Installing Packages
BiocManager::install("SingleR")
BiocManager::install("celldex")

#Installing libraries
library(SingleR)
library(celldex)
library(ggplot2)
library(dplyr)

theme_github <- theme_classic(base_size = 12) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )

#  Checking the SingleCellExperiment object
ncol(sce_seurat_data)          # number of cells
nrow(sce_seurat_data)          # number of genes
colnames(colData(sce_seurat_data))  # metadata
reducedDimNames(sce_seurat_data)   # PCA, UMAP

##  SingleR annotation
ref.set <- celldex::HumanPrimaryCellAtlasData()

##  Run SingleR
pred.cnts <- SingleR(test = sce_seurat_data, ref = ref.set, labels = ref.set$label.main,de.method = 'wilcox')
pred.cnts %>% head()

## Add labels to SCE and Seurat
colData(sce_seurat_data)$SingleR_label <- pred.cnts$pruned.labels
table(colData(sce_seurat_data)$SingleR_label)


##  Plot SingleR heatmap
singleR_heatmap <- plotScoreHeatmap(pred.cnts) + theme_github
png(
  "C:...../Github project/GSE114725/PLOTS/SingleR_heatmap.png",
  width = 12,
  height = 8
)
print(singleR_heatmap)
dev.off()


## Plot SingleR delta distribution
plot_delta <- plotDeltaDistribution(pred.cnts, ncol = 4, dots.on.top = FALSE)
png(
  "C:...../Github project/GSE114725/PLOTS/SingleR_delta_distribution.png",
  width = 3000,
  height = 2000,
  res = 300
)
print(plot_delta)
dev.off()


##  Plot Seurat UMAP colored by SingleR
Idents(seurat_data) <- "SingleR_HCA"
seurat_singleR_dimplot <- DimPlot(
  seurat_data,
  label = TRUE,
  repel = TRUE,
  label.size = 3
) + NoLegend() + theme_github

##  Save UMAP
png(
  "C:....../Github project/GSE114725/PLOTS/seurat_singleR_dimplot.png",
  width = 3000,
  height = 2000,
  res = 300
)
print(seurat_singleR_dimplot)
dev.off()
