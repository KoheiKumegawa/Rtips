#----------------------------------------------------------------------------
# scA03_visualize.R
#----------------------------------------------------------------------------
#Last update: 2022-08-22
#Kohei Kumegawa, Japanese Foundation for Cancer Research
library(ArchR)
arc <- readRDS("rds/arc.rds")

#parameter settings
dotsize <- 1
markers <- c("Epcam", "Vim")

#UMAP colored by cluster or sample
p1 <- plotEmbedding(arc, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", plotAs = "points", size = dotsize) 
p2 <- plotEmbedding(arc, colorBy = "cellColData", name = "Sample", embedding = "UMAP", plotAs = "points", size = dotsize)
plotPDF(p1, p2, name = "CA_UMAP.pdf", ArchRProj = arc, addDOC = FALSE, width = 5, height = 5)

#Sample x Cluster matrix
cM <- table(arc$Sample, arc$Clusters)[, paste0("C", seq_along(unique(arc$Clusters)))]
write.csv(cM, "output/Tables/CA_SampleClusterMatrix.csv")

#marker gene activity UMAP overlay
arc <- addImputeWeights(arc)
p3 <- plotEmbedding(arc, 
                    colorBy = "GeneScoreMatrix", 
                    name = markers, 
                    embedding = "UMAP", 
                    imputeWeights = getImputeWeights(arc), 
                    plotAs = "points", size = dotsize)
plotPDF(p3, name = "CA_UMAPoverlay_markerGS.pdf", ArchRProj = arc, addDOC = FALSE, width = 5, height = 5)

#motif score variance
varHomer <- getVarDeviations(arc, name = "homerMatrix", plot = F)
varCisbp <- getVarDeviations(arc, name = "cisbpMatrix", plot = F)

varHomer$motif <- "homer"
varCisbp$motif <- "cisbp"

varMotif <- rbind(varHomer, varCisbp)
write.csv(varMotif, "output/Tables/CA_MotifVariance.csv")

#motif score UMAP overlay, top50 variable motifs
p4 <- plotEmbedding(arc, 
                    colorBy = "homerMatrix", 
                    name = paste0("z:", varHomer$name[c(1:50)]), 
                    embedding = "UMAP", 
                    imputeWeights = getImputeWeights(arc), 
                    plotAs = "points", size = dotsize)
p5 <- plotEmbedding(arc, 
                    colorBy = "cisbpMatrix", 
                    name = paste0("z:", varCisbp$name[c(1:50)]), 
                    embedding = "UMAP", 
                    imputeWeights = getImputeWeights(arc), 
                    plotAs = "points", size = dotsize)
plotPDF(p4, name = "CA_UMAPoverlay_homerMotif_var50.pdf", ArchRProj = arc, addDOC = FALSE, width = 5, height = 5)
plotPDF(p5, name = "CA_UMAPoverlay_cisbpMotif_var50.pdf", ArchRProj = arc, addDOC = FALSE, width = 5, height = 5)
