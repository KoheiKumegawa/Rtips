#----------------------------------------------------------------------------
# scA_makeArchRProj.R
#----------------------------------------------------------------------------
#Last update: 2022-08-30
#Kohei Kumegawa, Japanese Foundation for Cancer Research
#system("mkdir data code rds ref output output/Plots output/Tables")
#Move tsv file into data directory
library(ArchR)
addArchRThreads(threads = 8)
addArchRGenome("mm10")

#assign sample tsv files
sampleName <- list.files(path = "data/", pattern = "tsv.gz")
names(sampleName) <- gsub(pattern = ".fragments.tsv.gz", replacement = "", sampleName)
outFile <- as.character(sampleName)

#make ArrowFiles
ArrowFiles <- character(length(sampleName))
ArrowFiles <- createArrowFiles(inputFiles = paste0("data//", outFile), 
                               sampleNames = names(sampleName),
                               minTSS = 4, 
                               minFrags = 1000, 
                               addTileMat = TRUE, addGeneScoreMat = TRUE, 
                               force = TRUE)

#infer doublets
doubScores <- addDoubletScores(ArrowFiles, k = 10, knnMethod = "UMAP", LSIMethod = 1)

#make pre-filtered ArchRProject
pre_arc <- ArchRProject(ArrowFiles, 
                        outputDirectory = "output", 
                        copyArrows = F)
pre_arc <- filterDoublets(pre_arc)

#quality metrics
df <- data.frame(nFrags = log10(pre_arc$nFrags), 
                 TSSe = pre_arc$TSSEnrichment,
                 Sample = pre_arc$Sample)
p1 <- lapply(names(sampleName), function(x){
  dfx <- df[which(df$Sample == x),]
  out <- ggplot(data = dfx, aes(x = nFrags, y = TSSe)) + geom_hex(bins = 100) + 
    theme_ArchR() + scale_fill_viridis_c(trans = "log") + lims(x = c(0, NA), y = c(0, NA)) + 
    labs(x = "Log10 Unique Fragments", y = "TSS Enrichment") +
    geom_vline(xintercept = c(log10(1000), log10(2000), log10(3000)), lty = "dashed") +
    ggtitle(x)
})
plotPDF(p1, name = "CA_scatter_nFrag_TSSe.pdf", ArchRProj = pre_arc, addDOC = FALSE, width = 5, height = 5)

mtx <- lapply(c(1000, 1500, 2000, 2500, 3000), 
              function(x) table(df$Sample[which(df$nFrags > log10(x) & df$TSSe > 8)])) %>% 
  do.call(rbind, .) %>% `rownames<-`(., paste0("nFrag_", c(1000, 1500, 2000, 2500, 3000)))
write.csv(mtx, file = "output/Tables/CA_QC_cellNumberSample_nFrag_TSSe8.csv", quote = F)

#save rds
saveRDS(pre_arc, "rds/pre_arc.rds")
