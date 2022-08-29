#------------------------------------------------------------------------------
# FragmentGRToBED.R
#------------------------------------------------------------------------------
#last update: 2022-08-26
#Kohei Kumegawa, Japanese Foundation for Cancer Research
library(data.table)

FragmentGRToBED <- function(
  gr = NULL,
  name = NULL,
  outputDir = NULL
){
  d <- data.frame(seqnames = seqnames(gr), start = start(gr)-1, end = end(gr))
  write.table(d, paste0(outputDir, "/", name, "-fragments.bed"), row.names = F, col.names = F, quote = F, sep = "\t")
  return(NULL)
}
