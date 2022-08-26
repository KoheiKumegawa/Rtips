

FragmentGRtoBED <- function(
  gr = NULL,
  name = NULL,
  outputDir = NULL
){
  d <- data.frame(seqnames = seqnames(gr), start = start(gr)-1, end = end(gr))
  write.table(d, paste0(outputDir, "/", name, "-fragments.bed"), row.names = F, col.names = F, quote = F, sep = "\t")
  
  return(NULL)
}
