#------------------------------------------------------------------------------
# ATACPreprocess.R
#------------------------------------------------------------------------------
#last update: 2022-08-29
#this code is modified by Kohei Kumegawa
#original code by ArchR package(Granja et al)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(dplyr)
library(Rcpp)
library(Rsamtools)
library(data.table)

#--------- functions ---------#
bamToFragmentGR <- function(
  bamPATH = NULL,
  bamNAME = NULL,
  offsetPlus = 4,
  offsetMinus = -5,
  bamFlag = NULL
){
  if(is.null(bamPATH)){
    stop("Please set PATH to bam files")
  }
  if(is.null(bamNAME)){
    stop("No input bamNAME; please recheck your input")
  }
  if(is.null(bamFlag)){
    stop("Please set bamFlag using Rsamtools's scanBamFlag!")
  }
  
  #1. Read In Bam File
  sF <- scanBam(bamPATH, param = ScanBamParam(flag = bamFlag, what = c("rname","pos", "isize")))[[1]]
  
  #2. Make Fragment Table
  dt <- data.table(seqnames = sF$rname, start = sF$pos + offsetPlus, end = sF$pos + abs(sF$isize) - 1 + offsetMinus)
  
  #3. Make fragment Granges and remove unwanted chromosomes
  gr <- GRanges(seqnames = dt$seqnames, IRanges(start = dt$start, end = dt$end))
  idy = which(seqnames(gr) %in% seqlevels(gr)[grep("random|chrM|chrUn|chrEBV", seqlevels(gr))])
  gr <- gr[-idy]
  gr <- dropSeqlevels(gr, seqlevels(gr)[grep("random|chrM|chrUn|chrEBV", seqlevels(gr))])
  mcols(gr) <- DataFrame(sample = bamNAME)
  
  #4. output Granges List
  return(gr)
}

FragmentGRToBED <- function(
  gr = NULL,
  name = NULL,
  outputDir = NULL
){
  d <- data.frame(seqnames = seqnames(gr), start = start(gr)-1, end = end(gr))
  write.table(d, paste0(outputDir, "/", name, "-fragments.bed"), row.names = F, col.names = F, quote = F, sep = "\t")
  return(NULL)
}

sourceCpp(code='
          #include <Rcpp.h>
          using namespace Rcpp;
          using namespace std;
          // [[Rcpp::export]]
          IntegerMatrix tabulate2dCpp(IntegerVector x1, int xmin, int xmax, IntegerVector y1, int ymin, int ymax){
          if(x1.size() != y1.size()){
          stop("width must equal size!");
          }
          IntegerVector x = clone(x1);
          IntegerVector y = clone(y1);
          int n = x.size();
          IntegerVector rx = seq(xmin,xmax);
          IntegerVector ry = seq(ymin,ymax);
          IntegerMatrix mat( ry.size() , rx.size() );
          int xi,yi;
          for(int i = 0; i < n; i++){
          xi = (x[i] - xmin);
          yi = (y[i] - ymin);
          if(yi >= 0 && yi < ry.size()){
          if(xi >= 0 && xi < rx.size()){
          mat( yi , xi ) = mat( yi , xi ) + 1; 
          }
          }
          }
          return mat;
          }')

#function
TSSenrich <- function(
  fragments = NULL,
  TSSgranges = NULL,
  flank = 2000,
  norm = 100,
  smooth = 51,
  range = 50,
  by = NULL
){
  message(sprintf("Convert %s fragments to inserts", length(fragments)))
  inserts <- c(
    GRanges(seqnames = seqnames(fragments), ranges = IRanges(start(fragments), start(fragments)), sample = mcols(fragments)[,by]),
    GRanges(seqnames = seqnames(fragments), ranges = IRanges(end(fragments), end(fragments)), sample = mcols(fragments)[,by]))
  
  message(sprintf("center the features"))
  center <- unique(resize(TSSgranges, width = 1, fix = "center", ignore.strand = FALSE))
  
  message(sprintf("get overlaps between the feature and insertions only up to %s bp", flank))
  overlap <- DataFrame(findOverlaps(query = center, subject = inserts, maxgap = flank, ignore.strand = TRUE))
  overlap$strand <- strand(center)[overlap[,1]]
  overlap$name <- mcols(inserts)[overlap[,2],by]
  overlap <- transform(overlap, id=match(name, unique(name)))
  ids <- length(unique(overlap$name))
  
  message(sprintf("calculate distance from TSS"))
  overlap$dist <- NA
  minus <- which(overlap$strand == "-")
  other <- which(overlap$strand != "-")
  overlap$dist[minus] <- start(center[overlap[minus,1]]) - start(inserts[overlap[minus,2]])
  overlap$dist[other] <- start(inserts[overlap[other,2]]) - start(center[overlap[other,1]])
  
  message(sprintf("make insertion matrix"))
  profile_mat <- tabulate2dCpp(x1 = overlap$id, y1 = overlap$dist, xmin = 1, xmax = ids, ymin = -flank, ymax = flank)
  colnames(profile_mat) <- unique(overlap$name)
  profile <- rowSums(profile_mat)
  
  message(sprintf("calculate TSS enrichments"))
  profile_mat_norm <- apply(profile_mat, 2, function(x) x/mean(x[c(1:norm,(flank*2-norm+1):(flank*2+1))]))
  profile_norm <- profile/mean(profile[c(1:norm,(flank*2-norm+1):(flank*2+1))])
  profile_mat_norm_smooth <- apply(profile_mat_norm, 2, function(x) zoo::rollmean(x, smooth, fill = 1))
  profile_norm_smooth <- zoo::rollmean(profile_norm, smooth, fill = 1)
  max_finite <- function(x){
    suppressWarnings(max(x[is.finite(x)], na.rm=TRUE))
  }
  
  #TSS enrichment score                                 
  e <- max_finite(profile_norm_smooth[(flank-range):(flank+range)])
  
  #Insersion profiles                                  
  rownames(profile_mat_norm_smooth) <- c(-2000:2000)
  plotDF <- reshape2::melt(profile_mat_norm_smooth)
  colnames(plotDF) <- c("distance", "sample", "normInsert")
  
  #fragment width                                   
  w <- as.numeric(width(fragments))
  
  out <- list(enrichScore = e, plotDF = plotDF, fragmentWidth = w)
  message(sprintf("ATAC QC metrics successfully calculated!"))
  return(out)
}


#--------- analysis ---------#
#import fragments
bamFiles  <- list.files("data/", pattern = ".rmdup.bam") %>% `names<-`(., gsub(".rmdup.bam", "", .))
fragments <- parallel::mclapply(names(bamFiles), 
                                function(x){
                                  out <- bamToFragmentGR(bamPATH = paste0("../data/", bamFiles[x]), 
                                                         bamNAME = x, 
                                                         bamFlag = scanBamFlag(isMinusStrand = FALSE, isProperPair  = TRUE))
                                  return(out)
                                  }, mc.cores = 10) %>% `names<-`(., names(bamFiles))

#quality metrics
tss <- TxDb.Hsapiens.UCSC.hg38.knownGene %>% transcripts(.) %>% resize(., width = 1, fix = "start") %>% unique
QualityMetricsATAC <- parallel::mclapply(names(fragments), FUN = function(x){
  message(sprintf("Data processing : %s", x))
  a <- fragments[[x]]
  b <- TSSenrich(fragments = a, TSSgranges = tss, by = "sample")
  return(b)
}, mc.cores = 10)
names(QualityMetricsATAC) <- names(fragments)

#check TSS enrichment score
lapply(QualityMetricsATAC, function(x) x$enrichScore) %>% unlist(.) %>% sort(.)

#output
saveRDS(QualityMetricsATAC, "rds/QualityMetricsATAC.rds")
parallel::mclapply(names(fragments), function(x) FragmentGRToBED(gr = fragments[[x]], name = x, outputDir = "output/fragments_bed/"), mc.cores = 10)
