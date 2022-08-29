#------------------------------------------------------------------------------
# ConstructATACseqCountMatrix.R
#------------------------------------------------------------------------------
#last update: 2022-08-29
#Kohei Kumegawa, Japanese Foundation for Cancer Research
#method and code from Corces et al Science 2018 
library(Rsamtools)
library(data.table)
library(dplyr)
library(rtracklayer)
library(SummarizedExperiment)
'%ni%' <- Negate("%in%")

#--------- functions ---------#
FragmentGRToInsert <- function(
  gr = NULL,
  name = NULL,
  outputDirBED = NULL
){
  inserts <- c(
    GRanges(seqnames = seqnames(gr), ranges = IRanges(start(gr), start(gr))),
    GRanges(seqnames = seqnames(gr), ranges = IRanges(end(gr), end(gr)))
  )
  #save bed file (2bp)
  d <- data.frame(seqnames = seqnames(inserts), start = start(inserts)-1, end = end(inserts))
  write.table(d, paste0(outputDirBED, "/", name, "_inserts.bed"), row.names = F, col.names = F, quote = F, sep = "\t")
  return(NULL)
}

RunMacs2 <- function(
  inputBedPATH = NULL,
  inputBedName = NULL,
  genome = NULL,
  outputDir = NULL,
  shift = NULL,
  extsize = NULL,
  method = c("p", "q"),
  cutoff = 0.05
){
  if(genome %ni% c("hg19", "hg38", "mm10")){
    stop("Please set genome as hg19, hg38 and mm10!")
  }
  if(genome %in% c("hg19", "hg38")){
    gen <- "hs"
  }
  if(genome == "mm10"){
    gen <- "mm"
  }
  
  commandPeaks <- sprintf(
    "macs2 callpeak -g %s --name %s --treatment %s --outdir %s --format BED --nomodel --call-summits --nolambda --keep-dup all",
    gen, inputBedName, inputBedPATH, outputDir)
  
  if (!is.null(shift) & !is.null(extsize)) {
    commandPeaks <- sprintf("%s --shift %s --extsize %s", commandPeaks, shift, extsize)
  }
  if (tolower(method) == "p") {
    commandPeaks <- sprintf("%s -p %s", commandPeaks, cutoff)
  } else {
    commandPeaks <- sprintf("%s -q %s", commandPeaks, cutoff)
  }
  message("Running Macs2...")
  message(commandPeaks)
  system(commandPeaks, intern = TRUE)
  
  return(NULL)
}

MakeSamplePeakSet <- function(gr, by = "score"){
  #nonOverlappingGRanges
  stopifnot(by %in% colnames(mcols(gr)))
  
  #function for picking up most significant peaks
  clusterGRanges <- function(gr, by = "score"){
    gr <- sort(sortSeqlevels(gr))
    r <- GenomicRanges::reduce(gr, min.gapwidth=0L, ignore.strand=TRUE)
    o <- findOverlaps(gr,r)
    mcols(gr)$cluster <- subjectHits(o)
    gr <- gr[order(mcols(gr)[,by], decreasing = TRUE),]
    gr <- gr[!duplicated(mcols(gr)$cluster),]
    gr <- sort(sortSeqlevels(gr))
    mcols(gr)$cluster <- NULL
    return(gr)
  }
  
  #iteration of filtering overlapping peaks
  i <-  0
  gr_converge <- gr
  while(length(gr_converge) > 0){
    i <-  i + 1
    gr_selected <- clusterGRanges(gr = gr_converge, by = by)
    gr_converge <- subsetByOverlaps(gr_converge, gr_selected, invert=TRUE) #blacklist selected gr
    if(i == 1){ #if i=1 then set gr_all to clustered
      gr_all <- gr_selected
    }else{
      gr_all <- c(gr_all, gr_selected)
    }
  }
  gr_all <- sort(sortSeqlevels(gr_all))
  return(gr_all)
}

MakeATACSummarizedExperiment <- function(
  fragmentGRangesList = NULL,
  unionPeaks = NULL,
  blacklist = NULL,
  sampleName = NULL,
  prior.count = 5,
  by = "sample"
){
  fragments <- unlist(as(fragmentGRangesList, "GRangesList"))
  inserts <- c(
    GRanges(seqnames = seqnames(fragments), ranges = IRanges(start(fragments), start(fragments)), sample = mcols(fragments)[,by]),
    GRanges(seqnames = seqnames(fragments), ranges = IRanges(end(fragments), end(fragments)), sample = mcols(fragments)[,by])
  )
  overlapDF <- DataFrame(findOverlaps(unionPeaks, inserts, ignore.strand = TRUE, maxgap=-1L, minoverlap=0L, type = "any"))
  overlapDF$name <- mcols(inserts)[overlapDF[, 2], by]
  overlapTDF <- transform(overlapDF, id = match(name, unique(name)))
  #Summarize
  sparseM <- Matrix::sparseMatrix(
    i = overlapTDF[, 1], 
    j = overlapTDF[, 4],
    x = rep(1, nrow(overlapTDF)), 
    dims = c(length(unionPeaks), length(unique(overlapDF$name))))
  colnames(sparseM) <- unique(overlapDF$name)
  sparseM <- sparseM[, sampleName]
  rownames(sparseM) <- unionPeaks$name
  sparseM.cpm <- edgeR::cpm(sparseM, log = TRUE, prior.count = prior.count)
  sparseM.norm <- preprocessCore::normalize.quantiles(sparseM.cpm)
  colnames(sparseM.norm) <- colnames(sparseM.cpm)
  rownames(sparseM.norm) <- unionPeaks$name
  #SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = sparseM, normcounts = sparseM.norm),
                             rowRanges = unionPeaks, 
                             colData = DataFrame(sample = sampleName))
  return(se)
}

#--------- analysis ---------#
#read fragment bed as GRanges
FragFiles   <- list.files("output/fragments_bed/") %>% `names<-`(., gsub("-fragments.bed", "", .))
fragments <- parallel::mclapply(FragFiles, function(x){
  fr <- fread(paste0("output/fragments_bed/", x)) %>% data.frame(.)
  out <- GRanges(seqnames = fr[,1], IRanges(start = fr[,2]+1, end = fr[,3]))
  return(out)
}, mc.cores = 10)

#fragments to inserts
parallel::mclapply(names(fragments), function(x){
  FragmentGRToInsert(gr = fragments[[x]], 
                     name = x, 
                     outputDirBED = "output/inserts_bed/")
}, mc.cores = 10)

#Peak call for inserts
insertFiles <- list.files("output/inserts_bed/") %>% `names<-`(., gsub("_inserts.bed", "", .))
parallel::mclapply(names(insertFiles), function(x){
  RunMacs2(inputBedPATH = paste0("output/inserts_bed/", insertFiles[x]),
           inputBedName = x,
           genome = "hg38",
           outputDir = "output/sample_peaks_row/",
           shift = -75,
           extsize = 150,
           method = "p",
           cutoff = 0.01)
}, mc.cores = 10)

#Make sample peak set
BSgenome   <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
chromSizes <- GRanges(names(seqlengths(BSgenome)), IRanges(1, seqlengths(BSgenome))) %>% 
              GenomeInfoDb::keepStandardChromosomes(., pruning.mode = "coarse")
blacklist  <- rtracklayer::import.bed("ref/ENCFF356LFX.bed")
gr_ls <- GenomicRanges::GRangesList(lapply(list.files("output/sample_peaks_row/", pattern = "summits.bed", full.names = T), function(x) import.bed(x)))
names(gr_ls) <- gsub("_summits.bed", "", list.files("output/sample_peaks_row/", pattern = "summits.bed"))
  
gr_ls_proc <- parallel::mclapply(gr_ls, function(x){
  gr <- resize(x, width = 501, fix = "center") %>%
    subsetByOverlaps(., chromSizes, type = "within") %>%
    subsetByOverlaps(., blacklist, invert=TRUE) %>%
    MakeSamplePeakSet(., by = "score")
  mcols(gr)$scorePerMillion <- mcols(gr)$score / (sum(mcols(gr)$score) / 1000000)
  return(gr)
}, mc.cores = 10) %>% GenomicRanges::GRangesList(.)

parallel::mclapply(names(gr_ls_proc), function(x){
  gr <- gr_ls_proc[[x]]
  d <- data.frame(seqnames = seqnames(gr), start = start(gr)-1, end = end(gr), name = mcols(gr)$name, score = mcols(gr)$score, scorePerMillion = mcols(gr)$scorePerMillion)
  write.table(d, paste0("output/sample_peaks_processed/", x, "_samplePeaks_processed.bed"), row.names = F, col.names = F, quote = F, sep = "\t")
  return(NULL)
}, mc.cores = 10)

#cancer type specific peakset
gr_cumulative <- MakeSamplePeakSet(unlist(gr_ls_proc), by = "scorePerMillion")
mcols(gr_cumulative)$sampleOverlap <- countOverlaps(gr_cumulative, gr_ls_proc)
reproduciblePeaks <- gr_cumulative[which(mcols(gr_cumulative)$scorePerMillion >= 5 &
                                mcols(gr_cumulative)$sampleOverlap   >= 2 &
                                seqnames(gr_cumulative) %ni% "chrY")]
names(reproduciblePeaks) <- paste0("repPeaks_", c(1:length(reproduciblePeaks)))
mcols(reproduciblePeaks)$origName <- mcols(reproduciblePeaks)$name
mcols(reproduciblePeaks)$name <- names(reproduciblePeaks)

#construct ATAC count matrix
fragments_named <- lapply(names(fragments), function(x){
  fr <- fragments[[x]]
  mcols(fr)$sample <- x
  return(fr)
})

se <- MakeATACSummarizedExperiment(fragmentGRangesList = fragments_named,
                                   unionPeaks = reproduciblePeaks,
                                   blacklist = blacklist,
                                   sampleName = sampleName,
                                   by = "sample")
saveRDS(se, "rds/se.rds")
