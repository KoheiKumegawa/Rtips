#------------------------------------------------------------------------------
# scA_CreatCustomArchRGenome_mm10.R
# Creating a custom ArchRGenome: mm10
#------------------------------------------------------------------------------
#Last update: 2022-08-22
#Kohei Kumegawa, Japanese Foundation for Cancer Research

library(ArchR)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(BSgenome.Mmusculus.UCSC.mm10)
"%ni%" <- Negate("%in%")

#import reference genes to add
add_ref_list <- "ref/filename.csv"
add_ref <- read.csv(add_ref_list)
add_gene <- GRanges(seqnames = add_ref$chr, 
                     IRanges(start = add_ref$start, end = add_ref$end), 
                     strand = add_ref$strand, 
                     symbol = paste0("ADD_", add_ref$X.4))
add_tss <- resize(add_gene, width = 1, fix = "start")

#mm10
geneAnno <- createGeneAnnotation(TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene, 
                                 OrgDb = org.Mm.eg.db)
gene <- c(geneAnno$genes, add_gene)
exon <- c(geneAnno$exons, add_gene[add_gene %ni% geneAnno$genes])
tss <- unique(c(geneAnno$TSS, add_tss))

# custom annotation
customGenomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Mmusculus.UCSC.mm10)
customGeneAnnotation <- createGeneAnnotation(TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene, 
                                             OrgDb = org.Mm.eg.db,
                                             genes = gene, exons = exon, TSS = tss)

saveRDS(customGenomeAnnotation, "rds/customGenomeAnnotation.rds")
saveRDS(customGeneAnnotation, "rds/customGeneAnnotation.rds")
