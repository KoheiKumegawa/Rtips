#----------------------------------------------------------------------------
# scR_Projection.R
#----------------------------------------------------------------------------
#Last update: 2022-09-05
#Kohei Kumegawa, Japanese Foundation for Cancer Research
library(Seurat)
library(ggplot2)
library(dplyr)
seu <- readRDS("rds/seu.rds")

#load and converted published data
load("data/RNAMagnetDataBundle/NicheData10x.rda")
p1 <- DimPlot(NicheData10x, reduction = "tsne", label = TRUE, label.size = 2) + 
      ggtitle("Baccin et al: tSNE plot") + ArchR::theme_ArchR()
pub_seu <- CreateSeuratObject(counts = NicheData10x@assays$RNA@counts)

#processing
pub_seu <- NormalizeData(pub_seu, normalization.method = "LogNormalize", scale.factor = 10000) %>%
           FindVariableFeatures(., selection.method = "vst", nfeatures = 2000) %>%
           ScaleData(., features = rownames(.)) %>%
           RunPCA(.) %>%
           FindNeighbors(., dims = 1:50) %>%
           FindClusters(., resolution = 1) %>%
           RunUMAP(., dims = 1:50, return.model = T)
pub_seu$cellType <-  NicheData10x@active.ident[colnames(pub_seu)]

#visualization
p2 <- DimPlot(pub_seu, reduction = "umap", label = TRUE, group.by = "cellType", label.size = 2) + 
      ggtitle("Baccin et al: Cell type") + ArchR::theme_ArchR()
p3 <- DimPlot(pub_seu, reduction = "umap", label = TRUE, group.by = "seurat_clusters", label.size = 2) + 
      ggtitle("Baccin et al: Seurat clusters") + ArchR::theme_ArchR()

#projection
anchors <- FindTransferAnchors(reference = pub_seu, query = seu, dims = 1:50, reference.reduction = "pca")
proj <- TransferData(anchorset = anchors, reference = pub_seu, query = seu, refdata = list(cellType = "cellType"))
proj <- IntegrateEmbeddings(anchorset = anchors, reference = pub_seu, query = proj, new.reduction.name = "refpca")
proj <- ProjectUMAP(query = proj, query.reduction = "refpca", reference = pub_seu, 
                    reference.reduction = "pca", reduction.model = "umap")

cluster_cols <- ArchR::ArchRPalettes$calm %>% 
  c(., ArchR::ArchRPalettes$kelly[c(1:(length(levels(seu$Clusters)) - length(.)))]) %>%
  `names<-`(., levels(seu$Clusters))

p4 <- DimPlot(proj, reduction = "umap", group.by = "predicted.cellType", label = TRUE, label.size = 2, repel = TRUE) + 
      ggtitle("Query transferred labels") + ArchR::theme_ArchR()
p5 <- DimPlot(proj, reduction = "umap", group.by = "Clusters", label = TRUE, label.size = 2, repel = TRUE, cols = cluster_cols) + 
      ggtitle("Query clusters") + ArchR::theme_ArchR()

pdf("output/Plots//GEx_projected_Baccin.pdf", width = 10, height = 7)
p1
p2+p3
p4+p5
dev.off()

cM1 <- table(proj$predicted.cellType, proj$Clusters)
cM2 <- table(proj$predicted.cellType, proj$orig.ident)
write.csv(cM1, "output/Tables/GEx_projected_Baccin_CelltypeClusterMatrix.csv", quote = F)
write.csv(cM2, "output/Tables/GEx_projected_Baccin_SampleClusterMatrix.csv", quote = F)
