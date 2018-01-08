## to use this script, just change the logpath to the directory where all of the clipLengthCounts.txt files are for an alignment. 
## also change the name of the file to be output (last line of the script, while running plotting function)
## the functions will not be screwed up by other files present as it looks for " clipLengthCounts" in the name of the file

library(lattice)
library(ggplot2)
library(gplots)
logpath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/alignment_11-28-17/FisherPools/logs"
fig_dir = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/figures/readLengthDistributions/"
dir.create(fig_dir)

importReadLengthsDists = function(logpath = logpath){
  DF = matrix(ncol = 3, nrow = 0)
  colnames(DF)<- c("length", "count", "sample")
  sample = 0
  sampleNames = c()
  for (readDist in Sys.glob(file.path(logpath, "*.txt"))) {
    file.name = strsplit(strsplit(readDist, "/")[[1]][length(strsplit(readDist, "/")[[1]])], '\\.')[[1]][1]
    if (grepl("clipLengthCounts", file.name)){
      sample = sample+1
      cat(file.name,"\n")
      sampleNames = c(sampleNames, file.name)
      lengthDF = read.table(file = paste(logpath,"/", file.name, ".txt", sep = ""), sep = "")
      lengthDF = cbind(lengthDF, sample)
      DF = rbind(DF, lengthDF)
  }
  }
  return(list(DF, sampleNames))
}

test = importReadLengthsDists(logpath = logpath)

batchPlotReadLengths = function(lengthDat  = test, filename = "test.pdf"){
  lenDat = lengthDat[[1]]
  sampleNames = lengthDat[[2]]
  for (i in 1:length(sampleNames)){
    lenDat$sample[lenDat$sample == i] <- sampleNames[i]
  }
  ggplot(lenDat, aes(x = V1, y = V2)) + 
         geom_bar(stat="identity") +  # stat="identity" allows you to use y column to set bar heights
         facet_wrap(~ sample) +
         xlab("Read Length") + ylab("Count")
  
  ggsave(file = paste(fig_dir, filename, sep = ""),plot = last_plot(), width=15, height=15)
 }


batchPlotReadLengths(lengthDat  = test, filename = "alignment_11-28-17_FisherPools_ReadLengthDist.pdf")



