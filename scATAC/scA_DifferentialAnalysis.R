#----------------------------------------------------------------------------
# scA_DifferentialAnalysis.R
#----------------------------------------------------------------------------
#Last update: 2022-08-22
#Kohei Kumegawa, Japanese Foundation for Cancer Research
library(ArchR)
arc <- readRDS("rds/arc.rds")

#differential test for all features
diff_se <- list(
  GeneScore = getMarkerFeatures(
    ArchRProj = arc, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"),
  Peak = getMarkerFeatures(
      ArchRProj = arc, 
      useMatrix = "PeakMatrix", 
      groupBy = "Clusters",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      testMethod = "wilcoxon"),
  HomerMotif = getMarkerFeatures(
    ArchRProj = arc, 
    useMatrix = "homerMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon", useSeqnames = "z"),
  CisbpMotif = getMarkerFeatures(
    ArchRProj = arc, 
    useMatrix = "cisbpMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon", useSeqnames = "z")
)
saveRDS(diff_se, "rds/diff_se.rds")

#mean score for each clusters
summarize_means <- lapply(c(c("GeneScore", "HomerMotif", "CisbpMotif")), function(x){
  out <- assays(diff_se[[x]])$Mean %>% `rownames<-`(., mcols(diff_se[[x]])$name)
  write.csv(out, paste0("output/Tables/CA_ClusterMeans_", x, ".csv"))
  return(NULL)
})

#summarize test for genescore, homermotif, cisbpmotif
summarize_test <- lapply(c("GeneScore", "HomerMotif", "CisbpMotif"), function(x){
  res <- lapply(names(diff_se[[x]]@assays@data), function(y){
    res <- round(diff_se[[x]]@assays@data[[y]], digits = 3) %>% `colnames<-`(., colnames(diff_se[[x]]))
    res$feature <- diff_se[[x]]@elementMetadata$name
    res <- reshape2::melt(res) %>% `colnames<-`(., c("feature", "cluster", y))
  })
  out <- lapply(res, function(z) z[,3]) %>% do.call(cbind, .) %>% `colnames<-`(., names(diff_se[[x]]@assays@data))
  out <- cbind(res[[1]][,c(1:2)], out)
  return(out)
}) %>% `names<-`(., c("GeneScore", "HomerMotif", "CisbpMotif"))
#summarize test for peak
res <- lapply(names(diff_se[["Peak"]]@assays@data), function(x){
  res <- round(diff_se[["Peak"]]@assays@data[[x]], digits = 3) %>% `colnames<-`(., colnames(diff_se[["Peak"]]))
  res$feature <- paste0(diff_se[["Peak"]]@elementMetadata$seqnames, ":", 
                        diff_se[["Peak"]]@elementMetadata$start, "-", diff_se[["Peak"]]@elementMetadata$end)
  res <- reshape2::melt(res) %>% `colnames<-`(., c("feature", "cluster", x))
  return(res)
})
df <- lapply(res, function(z) z[,3]) %>% do.call(cbind, .) %>% `colnames<-`(., names(diff_se[["Peak"]]@assays@data))
df <- cbind(res[[1]][,c(1:2)], df)
summarize_test[["Peak"]] <- df

#output csv
lapply(names(summarize_test), function(x) write.csv(summarize_test[[x]], paste0("output/Tables/CA_DiffTest_", x, ".csv")))
