#----------------------------------------------------------------------------
# scR01r_makeSeuratObject_Rhapsody.R
#----------------------------------------------------------------------------
#Last update: 2022-08-23
#Kohei Kumegawa, Japanese Foundation for Cancer Research
library(Seurat)

#system("mkdir data code rds ref output output/Plots output/Tables")
#Move csv file into data directory

#assign sample GEx count matrix
sampleName <- list.files(path = "data/", pattern = "_RSEC_MolsPerCell.csv")
names(sampleName) <- gsub(pattern = "_RSEC_MolsPerCell.csv", replacement = "", sampleName)

#make seurat object
makeSeuratObject <- function(sampleName = NULL,
                             samplePath = NULL){
  d <- t(read.csv(samplePath,row.names=1,comment.char="#"))
  a <- CreateSeuratObject(counts = d, project = sampleName)
  
  a <- RenameCells(a, new.names = paste0(sampleName, "#", colnames(a)))
  return(a)
}
seu_ls <- lapply(names(sampleName), function(x) makeSeuratObject(sampleName = x, samplePath = paste0("data//", sampleName[x])))
pre_seu <- merge(x = seu_ls[[1]], y = seu_ls[[-1]])

#quality metrics
pre_seu[["percent.mt"]] <- PercentageFeatureSet(pre_seu, pattern = "mt\\.") # mm
#pre_seu[["percent.mt"]] <- PercentageFeatureSet(pre_seu, pattern = "Mt\\.") # hg

p1 <- VlnPlot(pre_seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
p2 <- FeatureScatter(pre_seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
p3 <- FeatureScatter(pre_seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf("output/Plots/GEx_QCplots.pdf")
p1
p2
p3
dev.off()

#save rds
saveRDS(pre_seu, "rds/pre_seu.rds")
