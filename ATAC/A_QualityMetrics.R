#------------------------------------------------------------------------------
# A_QualityMetrics.R
#------------------------------------------------------------------------------
#last update: 2022-08-26
#this code is modified by Kohei Kumegawa
#original code by ArchR package(Granja et al)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(dplyr)
library(Rcpp)

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
  
  e <- max_finite(profile_norm_smooth[(flank-range):(flank+range)])
  
  rownames(profile_mat_norm_smooth) <- c(-2000:2000)
  plotDF <- reshape2::melt(profile_mat_norm_smooth)
  colnames(plotDF) <- c("distance", "sample", "normInsert")
  
  w <- as.numeric(width(fragments))

  out <- list(enrichScore = e, plotDF = plotDF, fragmentWidth = w)
  message(sprintf("ATAC QC metrics successfully calculated!"))
  return(out)
}

#Analysis
fragment_ls <- readRDS("rds/fragment_ls.rds")
tss <- TxDb.Hsapiens.UCSC.hg38.knownGene %>% transcripts(.) %>% resize(., width = 1, fix = "start") %>% unique

out <- parallel::mclapply(names(fragment_ls), FUN = function(x){
  message(sprintf("Data processing : %s", x))
  a <- fragment_ls[[x]]
  b <- TSSenrich(fragments = a, TSSgranges = tss, by = "sample")
  return(b)
}, mc.cores = 8)

saveRDS(out, "rds/ATAC_QualityMetrics.rds")
