#------------------------------------------------------------------------------
# MakeSE_RNAseq.R
#------------------------------------------------------------------------------
#Kohei Kumegawa, Japanese Foundation for Cancer Research
library(SummarizedExperiment)
library(data.table)
library(dplyr)
library(ggplot2)
countToTpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

#make summrized experiment
counts <- fread("data/counts_hg38.txt")
len <- counts$Length
gene <- counts$Geneid
counts <- counts[, -c(1:6)] %>% 
  `colnames<-`(., gsub("_umi_hg38.dedup.bam", "", colnames(counts[, -c(1:6)]))) %>%
  as.matrix %>% `rownames<-`(., gene) 
tpms <- countToTpm(counts, len)
se <- SummarizedExperiment(assays = list(counts = counts, log2tpm = log2(tpms+1)),
                           rowData = DataFrame(symbol = gene),
                           colData = DataFrame(sample = colnames(counts)))
saveRDS(se, "rds/se.rds")
