#----------------------------------------------------------------------------
# scR_Clustering.R
#----------------------------------------------------------------------------
#Last update: 2022-08-26
#Kohei Kumegawa, Japanese Foundation for Cancer Research
library(Seurat)
library(dplyr)
library(ggplot2)
pre_seu <- readRDS("rds/pre_seu.rds")

#quality filter
seu <- subset(pre_seu, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 25)

#normalize data and run PCA using top 2000 variable features
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(., selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(., features = rownames(.)) %>% RunPCA(.)
p1 <- ElbowPlot(seu, ndims = 50)
pdf("output/Plots/GEx_ElbowPlot.pdf", width = 5, height = 5)
p1
dev.off()

#construct community graph and identify clusters
use_dimension <- 50
seu <- FindNeighbors(seu, dims = 1:use_dimension) %>% RunUMAP(., dims = 1:use_dimension)
for(i in seq(0.2, 1.2, 0.1)){seu <- FindClusters(seu, resolution = i)}
p2 <- lapply(grep("RNA_snn_res.", colnames(seu@meta.data), value = T),
             function(x) DimPlot(seu, reduction = "umap", label = TRUE, group.by = x) + ggtitle(x))
pdf("output/Plots/GEx_multiResolution_UMAP.pdf", width = 5, height = 5)
p2
dev.off()

#define clusters
use_cluster <- paste0("RNA_snn_res.", 0.6)
seu$Clusters <- factor(paste0("GEx_C", as.numeric(seu@meta.data[, use_cluster])), 
                       levels = paste0("GEx_C", c(1:max(as.numeric(seu@meta.data[, use_cluster])))))
cluster_cols <- ArchR::ArchRPalettes$calm %>% 
  c(., ArchR::ArchRPalettes$kelly[c(1:(length(levels(seu$Clusters)) - length(.)))]) %>%
  `names<-`(., levels(seu$Clusters))
sample_cols <- ArchR::ArchRPalettes$stallion[seq_along(unique(seu$orig.ident))] %>% 
  `names<-`(., unique(seu$orig.ident))

#plotting
p3 <- DimPlot(seu, reduction = "umap", label = T, label.size = 2, group.by = "Clusters", cols = cluster_cols) + 
  ggtitle("clusters") + ArchR::theme_ArchR(legendPosition = "right")
p4 <- DimPlot(seu, reduction = "umap", label = F, group.by = "orig.ident", cols = sample_cols) + 
  ggtitle("samples") + ArchR::theme_ArchR(legendPosition = "right")

pdf("output/Plots/GEx_UMAP.pdf", width = 6.5, height = 5)
p3
p4
dev.off()

#confusion matrix
cM <- table(seu$Clusters, seu$orig.ident)
write.csv(cM, "output/Tables/GEx_SampleClusterMatrix.csv", quote = F)

#save rds
saveRDS(seu, "rds/seu.rds")
