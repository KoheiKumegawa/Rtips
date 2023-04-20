#-------------------------------------------------------
# edgeR_PairwiseFunction.R
#-------------------------------------------------------
#this code originally introduced by Jeff Granja
#modified by Kohei Kumegawa
library(edgeR)
library(dplyr)
'%ni%' = Negate('%in%')

#edgeR pairwise test function
edgeR_pairwise <- function(se, compareCol = "Group", topGroup, bottomGroup, comparisonName = "Differential Test"){
  
  stopifnot("counts" %in% assayNames(se))
  assays(se) <- assays(se)[c(which(assayNames(se) %in% "counts"), which(assayNames(se) %ni% "counts"))]
  
  #Identify Top and Bottom Idx
  stopifnot(!any(topGroup %in% bottomGroup)) #topGroup and bottomGroup can not intersect
  topGroupIdx <- which(paste0(colData(se)[,compareCol]) %in% topGroup)
  bottomGroupIdx <- which(paste0(colData(se)[,compareCol]) %in% bottomGroup)
  stopifnot(length(topGroupIdx) > 1) #topGroup needs 2 samples
  stopifnot(length(bottomGroupIdx) > 1) #bottomGroup needs 2 samples
  message(sprintf("Comparing\nTop Samples: %s\nvs\nBottom Samples: %s", paste(colnames(se)[topGroupIdx],collapse=" "), paste(colnames(se)[bottomGroupIdx],collapse=" ")))
  
  #Keep only Columns important for comparison
  se <- se[,c(topGroupIdx, bottomGroupIdx)]
  
  #Design
  design <- as.formula(paste0("~",compareCol))
  
  #Run Differential Testing
  diffObj <- edgeR_pairwise_helper(se = se, design = design)
  
  #Pairwise Test
  diffTest <- diffObj$diffTestResults
  method <- diffObj$diffTestMethod
  rm(diffObj)
  
  suppressPackageStartupMessages(require(edgeR))
  message(sprintf("Running Pairwise test using edgeR on comparison %s...", comparisonName))
  
  #Check input
  stopifnot(all(topGroup %in% colnames(diffTest$design)))
  stopifnot(all(bottomGroup %in% colnames(diffTest$design)))
  
  #Make Contrasts
  contrast <- rep(0, ncol(diffTest$design))
  contrast[which(colnames(diffTest$design) %in% topGroup)] <- 1/length(which(colnames(diffTest$design) %in% topGroup))
  contrast[which(colnames(diffTest$design) %in% bottomGroup)] <- -1/length(which(colnames(diffTest$design) %in% bottomGroup))
  
  #Test
  e <- glmQLFTest(diffTest, contrast = contrast) %>% topTags(n = nrow(diffTest$counts), sort.by = "none") %>% data.frame()
  
  message("Organizing Output...")
  out <- data.frame(row.names = row.names(e),
                    log2Mean = e$logCPM,
                    log2FC = e$logFC,
                    pval = e$PValue,
                    FDR = e$FDR)
  
  message("Returning Pairwise Summarized Experiment...")
  #if a padj value is NA set to 1 ie non significant and if 0 convert to the least non-zero because -log10(fdr) plots become unhappy
  out[is.na(out$FDR),"FDR"] <- 1
  out[which(out$FDR==0),"FDR"] <- min(out[which(out$FDR!=0),"FDR"])
  seDiff <- SummarizedExperiment(assays = SimpleList(differential = as.matrix(out)), rowRanges = rowRanges(se))
  seDiff@metadata$compName <- comparisonName
  seDiff@metadata$compGroups <- list(topGroup = topGroup, bottomGroup = bottomGroup)
  
  return(seDiff)
  
}

edgeR_pairwise_helper <- function(se, design = ~Group){
  
  #Make Model Matrix
  design0 <- as.formula(paste0(paste0(as.character(design),collapse="")," + 0"))
  modelMatrix <- model.matrix(design0, data = colData(se))
  
  if(length(all.vars(design)) == 1){
    colnames(modelMatrix) <- gsub(all.vars(design),"",colnames(modelMatrix))
  }
  
  mainMethod <- "edgeR"
  message("Running EdgeR TMM Pipeline...")
  suppressPackageStartupMessages(require(edgeR))
  
  # GLM QLF
  out <- DGEList(assay(se)) %>% calcNormFactors(., method = "TMM") %>% 
    estimateDisp(., design = modelMatrix, robust = TRUE) %>%
    glmQLFit(., design = modelMatrix)
  
  return(list(design = design, diffTestMethod = mainMethod, diffTestResults = out))
}
