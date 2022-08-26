#----------------------------------------------------------------------------
# scR_CellCycleRegression.R
#----------------------------------------------------------------------------
#Last update: 2022-08-26
#Kohei Kumegawa, Japanese Foundation for Cancer Research
library(Seurat)
library(org.Hs.eg.db)
seu <- readRDS("rds/seu.rds")
seu_ccreg <- seu

#cell cycle genes
cell_cycle_genes <- RCurl::getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Homo_sapiens.csv") %>% read.csv(text = .)
cell_cycle_genes <- data.frame(phase = cell_cycle_genes$phase,
                               gene  = AnnotationDbi::select(org.Hs.eg.db, 
                                                             keys = cell_cycle_genes$geneID, 
                                                             keytype = "ENSEMBL", 
                                                             columns = c("ENSEMBL", "SYMBOL"))[,2])
s_genes   <- cell_cycle_genes[which(cell_cycle_genes$phase == "S"), "gene"]
g2m_genes <- cell_cycle_genes[which(cell_cycle_genes$phase == "G2/M"), "gene"]

#cell cycle gene regression
seu_ccreg <- CellCycleScoring(seu_ccreg, s.features = s_genes, g2m.features = g2m_genes)
seu_ccreg <- ScaleData(seu_ccreg, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seu_ccreg)) %>%
             RunPCA(., features = VariableFeatures(.), nfeatures.print = 10)
p1 <- ElbowPlot(seu_ccreg, ndims = 50)
pdf("output/Plots/GEx_ElbowPlot_CCregress.pdf", width = 5, height = 5)
p1
dev.off()

#construct community graph and identify clusters
use_dimension <- 50
seu_ccreg <- FindNeighbors(seu_ccreg, dims = 1:use_dimension) %>% RunUMAP(., dims = 1:use_dimension)
for(i in seq(0.2, 1.2, 0.1)){seu_ccreg <- FindClusters(seu_ccreg, resolution = i)}
p2 <- lapply(grep("RNA_snn_res.", colnames(seu_ccreg@meta.data), value = T),
             function(x) DimPlot(seu_ccreg, reduction = "umap", label = TRUE, group.by = x) + ggtitle(x))
pdf("output/Plots/GEx_multiResolution_UMAP_CCregress.pdf", width = 5, height = 5)
p2
dev.off()

#define clusters
use_cluster <- paste0("RNA_snn_res.", 0.6)
seu_ccreg$Clusters <- factor(paste0("GEx_C", as.numeric(seu_ccreg@meta.data[, use_cluster])), 
                             levels = paste0("GEx_C", c(1:max(as.numeric(seu_ccreg@meta.data[, use_cluster])))))
cluster_cols <- ArchR::ArchRPalettes$calm[seq_along(levels(seu_ccreg$Clusters))] %>% `names<-`(., levels(seu_ccreg$Clusters))
sample_cols  <- ArchR::ArchRPalettes$stallion[seq_along(unique(seu_ccreg$Sample))] %>% `names<-`(., sampleName)

#plotting
p3 <- DimPlot(seu_ccreg, reduction = "umap", label = T, label.size = 2, group.by = "Clusters", cols = cluster_cols) + 
      ggtitle("clusters") + ArchR::theme_ArchR(legendPosition = "right")
p4 <- DimPlot(seu_ccreg, reduction = "umap", label = F, group.by = "Sample", cols = sample_cols) + 
      ggtitle("samples") + ArchR::theme_ArchR(legendPosition = "right")

pdf("output/Plots/GEx_UMAP_CCregress.pdf", width = 5.5, height = 5)
p3
p4
dev.off()

#confusion matrix
cM <- table(seu_ccreg$Sample, seu_ccreg$Clusters)
write.csv(cM, "output/Tables/GEx_SampleClusterMatrix_CCregress.csv", quote = F)

#save rds
saveRDS(seu_ccreg, "rds/seu_ccreg.rds")
