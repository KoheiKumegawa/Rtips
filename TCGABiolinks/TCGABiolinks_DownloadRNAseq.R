#----------------------------------------------------------------------------
# TCGABiolinks_DownloadRNAseq.R
#----------------------------------------------------------------------------
#Last update: 2022-08-26
#Kohei Kumegawa, Japanese Foundation for Cancer Research
library(TCGAbiolinks)
setwd("data/")

# download TCGA-BRCA RNA-seq Counts data
query <- GDCquery(project="TCGA-BRCA",
                  data.category="Transcriptome Profiling",
                  data.type="Gene Expression Quantification",
                  workflow.type="STAR - Counts")
GDCdownload(query)
data <- GDCprepare(query)

# RangedSummarizedExperiment
saveRDS(data, "../rds/TCGA_BRCA_RNAExp.rds")
