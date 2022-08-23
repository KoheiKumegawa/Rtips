#----------------------------------------------------------------------------
# RhapsodyMultiplexFilter.R
#----------------------------------------------------------------------------
#Last update: 2022-08-23
#Kohei Kumegawa, Japanese Foundation for Cancer Research
library(Seurat)
library(dplyr)
library(ggplot2)

#assign sample GEx count matrix
sampleName <- list.files(path = "data/", pattern = "_RSEC_MolsPerCell.csv")
names(sampleName) <- gsub(pattern = "_RSEC_MolsPerCell.csv", replacement = "", sampleName)

#tag info
tag_call <- read.csv(list.files("data/", pattern = "_Tag_Calls.csv", full.names = T), comment.char = "#") %>% `rownames<-`(., .[,1])
tag_read <- read.csv(list.files("data/", pattern = "_Tag_ReadsPerCell.csv", full.names = T), comment.char = "#") %>% 
  `rownames<-`(., .[,1]) %>% `colnames<-`(., gsub(".stAbO", "", colnames(.))) %>% 
  .[, which(colnames(.) %in% unique(tag_call$Sample_Tag))]
tag_read <- log2(tag_read+1)

#make seurat object
makeSeuratObject <- function(sampleName = NULL,
                             samplePath = NULL){
  d <- t(read.csv(samplePath,row.names=1,comment.char="#"))
  a <- CreateSeuratObject(counts = d, project = sampleName)
  #a <- RenameCells(a, new.names = paste0(sampleName, "#", colnames(a)))
  return(a)
}
seu <- makeSeuratObject(sampleName = "bc_org", samplePath = paste0("data/", sampleName))

seu$orig.ident <- tag_call[colnames(seu), 2] #sample assignment
seu <- AddMetaData(seu, metadata = tag_read) #sample tag numbers
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "MT\\.") #mtRNA

p1 <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
pdf("output/Plots/GEx_QC.pdf", width = 6, height = 6)
p1
dev.off()

#quality filter
seu <- subset(seu, subset = nFeature_RNA > 400 & nFeature_RNA < 9000 & percent.mt < 40)

#per sample analysis
SampleTag <- paste0("SampleTag0", c(1:7), "_hs")
seu_sample <- lapply(SampleTag, function(x){
  out <- seu[, which(seu$orig.ident == x)]
  out <- NormalizeData(out, normalization.method = "LogNormalize", scale.factor = 10000) %>%
    FindVariableFeatures(., selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(., features = rownames(.)) %>%
    RunPCA(.) %>% FindNeighbors(., dims = 1:20) %>% RunUMAP(., dims = 1:20) %>% FindClusters(., resolution = 1)
  return(out)
}) %>% `names<-`(., SampleTag)

#UMAP per sample
p2 <- lapply(SampleTag, function(x){
  a <- seu_sample[[x]]
  out <- DimPlot(a, reduction = "umap", label = TRUE, group.by = "seurat_clusters") + ArchR::theme_ArchR() + ggtitle(x)
  return(out)
})
pdf("output/Plots/GEx_QC_UMAPsample.pdf", width = 5, height = 6)
p2
dev.off()

#each sampleTag number in each cluster in each sample
p3 <- lapply(SampleTag, function(x){
  seu <- seu_sample[[x]]
  max_val <- max(seu@meta.data[, c(4:10)])
  out <- lapply(SampleTag, function(y){
    tmp <- VlnPlot(seu, features = y, group.by = "seurat_clusters") + ggtitle(x) + ylab(y) + ylim(0, max_val)
    return(tmp)})
  return(out)
})
pdf("output/Plots/GEx_QC_tag_violinPlot.pdf", width = 5, height = 5)
p3
dev.off()

#Density plot of each sampleTag in each cluster in each sample
p4 <- lapply(SampleTag, function(x){
  seu <- seu_sample[[x]]
  clusters <- levels(seu$seurat_clusters)
  out <- lapply(clusters, function(y){
    df <- seu@meta.data[which(seu$seurat_clusters == y), c(4:10)]
    df <- reshape2::melt(df)
    tmp <- ggplot(df, aes(x = value, fill = variable)) + geom_density(alpha = 0.5) + 
      ArchR::theme_ArchR() + ggtitle(paste0(x, " :cluster ", y)) +
      labs(x = "log2(TagCounts + 1)", y = "Density")
  })
  return(out)
})
pdf("output/Plots/GEx_QC_tag_density.pdf", width = 5, height = 3)
p4
dev.off()

#Density of each called tags in each cluster
p5 <- lapply(SampleTag, function(x){
  seu <- seu_sample[[x]]
  df <- seu@meta.data[, c(x, "seurat_clusters")] %>% `colnames<-`(., c("Tag", "Cluster"))
  out <- ggplot(df, aes(x = Tag, fill = Cluster)) + geom_density(alpha = 0.8) + ArchR::theme_ArchR() + ggtitle(x) +
    labs(x = "log2(AssignedTagCounts + 1)", y = "Density")
  return(out)
})
pdf("output/Plots/GEx_QC_tag_density_cluster.pdf", width = 5, height = 5)
p5
dev.off()

#save rds
saveRDS(seu, "rds/seu.rds")
saveRDS(seu_sample, "rds/seu_sample.rds")
