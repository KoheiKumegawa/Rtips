#------------------------------------------------------------------------------
# homerFreq.R
#------------------------------------------------------------------------------
library(data.table)
library(dplyr)
library(ggplot2)

listFiles <- list.files("output/homer_motif/")
motifDF <- lapply(listFiles, function(x){
  out <- fread(paste0("output/homer_motif/CA/", x, "/knownResults.txt")) %>% data.frame()
  out[,1] <- paste0(stringr::str_split(out[,1], pattern = "/", simplify = T)[,1], "#", 
                    stringr::str_split(stringr::str_split(out[,1], pattern = "\\/", simplify = T)[,2], pattern = "-", simplify = T)[,1])
  return(out)
}) %>% `names<-`(., listFiles)

sig_motifs <- motifDF$CA_totalPeaks_48h$Motif.Name[motifDF$CA_totalPeaks_48h$q.value..Benjamini. < 0.01]
TgTFreqDF <- lapply(motifDF, function(x){
  df <- x
  df <- df[df[,1] %in% sig_motifs,]
  df <- df[order(df[,1]), ]
  df <- df[, c(1,7)] %>% `colnames<-`(., c("Motif", "Freq"))
  df$Freq <- gsub("%", "", df$Freq) %>% as.numeric()
  return(df)
})

FreqDF <- lapply(TgTFreqDF[c(1:2)], function(x){
  df <- x
  out <- data.frame(Motif = df$Motif,
                    Freq  = log2(df$Freq / TgTFreqDF$CA_totalPeaks_48h$Freq))
  out <- out[order(out$Freq, decreasing = T), ]
  return(out)
})
