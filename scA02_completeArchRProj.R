#----------------------------------------------------------------------------
# scA02_completeArchRProj.R
#----------------------------------------------------------------------------
#Last update: 2022-08-22
#Kohei Kumegawa, Japanese Foundation for Cancer Research
library(ArchR)
library(dplyr)
addArchRThreads(threads = 8)

#cutoff setting
cutoff_TSSe   <- 8
cutoff_nFrags <- 1500 

#filter cells
pre_arc <- readRDS("rds/pre_arc.rds")
arc <- pre_arc[which(pre_arc$TSSEnrichment > cutoff_TSSe & pre_arc$nFrags > cutoff_nFrags)]

#clustering
arc <- addIterativeLSI(arc, useMatrix = "TileMatrix", name = "IterativeLSI") %>%
       addClusters(., reducedDims = "IterativeLSI", maxClusters = 30, force = T) %>%
       addUMAP(., reducedDims = "IterativeLSI", force = T)

#call peaks
pathToMacs2 <- findMacs2()
arc <- addGroupCoverages(arc, groupBy = "Clusters") %>%
       addReproduciblePeakSet(., groupBy = "Clusters", 
                              pathToMacs2 = pathToMacs2, 
                              method = "q", cutOff = 0.05) %>% 
       addPeakMatrix(.)

#calculate motif deviation
arc <- addMotifAnnotations(arc, motifSet = "homer", name = "homer")
arc <- addMotifAnnotations(arc, motifSet = "cisbp", name = "cisbp")
arc <- addBgdPeaks(arc) %>% addDeviationsMatrix(., peakAnnotation = "homer") %>% addDeviationsMatrix(., peakAnnotation = "cisbp")

#save rds
saveRDS(arc, "rds/arc.rds")
