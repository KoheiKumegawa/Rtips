#----------------------------------------------------------------------------
# 00_BamToTsv.R
#----------------------------------------------------------------------------
library(Rsamtools)
library(data.table)
library(parallel)
`%ni%` = Negate(`%in%`)

#----- function -----#
MakeArchRImput <- function(inputBam = NULL, inputName = NULL, inputPATH = NULL, outputPATH = NULL){
  bamFlag <- scanBamFlag(isMinusStrand = FALSE, isProperPair  = TRUE)
  outFile <- paste0(outputPATH, "/", inputName, ".fragments.tsv.gz")
  
  #1. Read In Bam File
  scanFragments <- scanBam(paste0(inputPATH,"/", inputBam),
                           param = ScanBamParam(
                             flag = bamFlag,
                             what = c("rname","pos", "isize"),
                             tag = c("DB") #### CAREFULLY CHECKED
                           ))[[1]]
  
  #2. Make RGTag
  RGTag <- gsub("alignments.possorted.tagged_", "", scanFragments$tag[[1]])
  
  #3. Make Fragment Table
  dt <- data.table(
    seqnames = scanFragments$rname,
    start = scanFragments$pos,
    end = scanFragments$pos + abs(scanFragments$isize) - 1,
    RG = RGTag
  )
  
  #5. Sort de-duplicate
  dt <- unique(dt[order(dt$seqnames, dt$start, dt$end), ])
  dt <- dt[dt$end > dt$start,] #this can sometimes be a weird aligner error!
  
  #6. Write To Fragments
  data.table::fwrite(dt, stringr::str_split(outFile, pattern="\\.gz", simplify=TRUE)[,1], sep = "\t", col.names = FALSE)
  
  #7. Bgzip
  Rsamtools::bgzip(stringr::str_split(outFile, pattern="\\.gz", simplify=TRUE)[,1], outFile)
  
  return(NULL)
}

#----- analysis -----#
bamFiles <- list.files("data/bam/", pattern = ".bam")
names(bamFiles) <- gsub(".bam", "", bamFiles)
mclapply(names(bamFiles), function(i){
  inputBam <- bamFiles[i]
  inputName <- i
  MakeArchRImput(inputBam = inputBam, inputName = inputName, inputPATH = "data/bam/", outputPATH = "data/")
})
