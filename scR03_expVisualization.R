#----------------------------------------------------------------------------
# scR03_expVisualization.R
#----------------------------------------------------------------------------
#Last update: 2022-08-23
#Kohei Kumegawa, Japanese Foundation for Cancer Research
library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(viridisLite)
seu <- readRDS("rds/seu.rds")
markers <- intersect(as.character(read.table("ref/markerGenes.tsv")[,1]), rownames(seu))
Idents(seu) <- seu$Clusters

#---- visualize marker gene ----#
#UMAP overlay
p1 <- lapply(markers, function(x) FeaturePlot(seu, features = x, order = T) + 
               scale_colour_gradientn(colours = viridis(256, option = "A")) + ArchR::theme_ArchR())
pdf("output/Plots/GEx_markers_UMAPoverlay.pdf", height = 5, width = 4)
p1
dev.off()

#Heatmap for cluster
avg <- AverageExpression(seu, return.seurat = T, group.by = "Clusters")
mtx <- avg@assays$RNA@scale.data[markers,]

col_fun1 <- colorRamp2(c(-2,-1,0,1,2), viridis(5, option = "A"))
fh = function(x) hclust(dist(x), method="ward.D2")
ht1 <- Heatmap(mtx, name = "Relative expression", cluster_columns = fh, cluster_rows = fh, show_row_dend = T, show_row_names = T, 
               col = col_fun1, row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 5))
p2 <- draw(ht1)
pdf("output/Plots/GEx_markers_Heatmap.pdf", height = 10, width = 4)
p2
dev.off()

#dotplot - genes sorted by Heatmap above
p3 <- DotPlot(seu, features = markers[row_order(p2)], dot.scale = 7, dot.min = 0.05) + 
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke= 1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8), axis.text.y = element_text(size = 8)) + 
  scale_colour_gradientn(colours = viridis(256, option = "A"))
pdf("output/Plots/GEx_markers_DotPlot.pdf", height = 5, width = 20)
p3
dev.off()

#---- DEG ----#
DEG <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(DEG, "output/Tables/GEx_DEGClusters.csv", quote=F)
