#------------------------------------------------------------------------------
# RunMacs2.R
#------------------------------------------------------------------------------
#last update: 2022-08-26
#Kohei Kumegawa, Japanese Foundation for Cancer Research
library(data.table)

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
