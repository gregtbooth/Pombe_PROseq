## Charles' Wave calling scripts were run on each biological replicate data for the cdk9as time course:
## it was also run using different size windows (larger windows were suggested to improve performance)
fig_dir = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/figures/cdk9asTC_waves/10-05-17/BothWaves/TryTSmooth_6Kbgenes/"
dir.create(fig_dir)
infopath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/pombeinfo/"
filteredGL = read.table(paste(infopath, "SP_CompleteFilteredGenes.bed", sep = ""))
BaseCountDat = read.table("/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/countData/AllGenesObsTSS/ePRO_Cdk9as_0min_combinedPIcountData.txt")
GL = filteredGL[,c(1,2,3,6,4,5)]
idx = (GL[,3] - GL[,2] >= 6000)
BasCountDat_6kbGenes = merge(filteredGL[idx,], BaseCountDat, by.x = "V4", by.y = 0)[,c(1,2,3,4,5,6,15)]
#####
# 10-06-17 Note: new runs give 4 more genes passing filter (n = 32), but these four extra genes do not have great wave calls): results in slightly different rate estimate.
## try running script using only original set of 28 
Original_clearWaves = read.table("/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/figures/cdk9asTC_waves/06-01-17/BothWaves/TryTSmooth_6Kbgenes/5_TSmooth/WaveLocations_50bpWin/30secClearWaves.bed")
OriginalGL = merge(GL[idx,], Original_clearWaves, by.x = "V4", by.y = "V4")[, c(2,3,4,5,1,6)]
colnames(OriginalGL) = c("V1", "V2","V3","V6","V4","V5")
######
# Below is a list of values used for offsetting the TSS and upStream window used for initializing the HMM when I ran the algorithm 
# to search for the region between advancing and clearing waves.
## These will be used in the function "CleanWaves" to back calculate the positions of wave fronts relative to actual TSS. 
TSSoffSet_List = list("30sec" = 1000, "1min" = 1000, "2.5min" = 1000)#, "5min" = 2500)
upStream_List = list("30sec" = 1000, "1min" = 1000, "2.5min" = 1000)#, "5min" = 1000)
# load Combined replicate data, rep1 data, rep2 data for comparison of wave calls (when using 50 bp windows)
LoadWaveRepData = function(r1File, r2File, combinedFile){
  result = list()
  load(file = r1File)
  load(file = r2File)
  load(file = combinedFile)
  # creaet new variables only for simplified wave call output = last (225th) list element (nomenclature: rep_time_window)
  # 30 second
  r1_30sec = data.frame(pw_30sec_r1[length(pw_30sec_r1)])
  r2_30sec = data.frame(pw_30sec_r2[length(pw_30sec_r2)])
  combined_30sec= data.frame(pw_30sec[length(pw_30sec)])
  # 1 min
  r1_1min = data.frame(pw_1min_r1[length(pw_1min_r1)])
  r2_1min = data.frame(pw_1min_r2[length(pw_1min_r2)])
  combined_1min = data.frame(pw_1min[length(pw_1min)])
  # 2.5 min
  r1_2.5min = data.frame(pw_2_5min_r1[length(pw_2_5min_r1)])
  r2_2.5min = data.frame(pw_2_5min_r2[length(pw_2_5min_r2)])
  combined_2.5min = data.frame(pw_2_5min[length(pw_2_5min)])
  # 5 min
  #r1_5min = data.frame(pw_5min_r1[length(pw_5min_r1)])
  #r2_5min = data.frame(pw_5min_r2[length(pw_5min_r2)])
  #combined_5min = data.frame(pw_5min[length(pw_5min)])
  # place all dataframes into results list
  result[["r1_30sec"]] = r1_30sec
  result[["r2_30sec"]] = r2_30sec
  result[["combined_30sec"]] = combined_30sec
  result[["r1_1min"]] = r1_1min
  result[["r2_1min"]] = r2_1min
  result[["combined_1min"]] = combined_1min
  result[["r1_2.5min"]] = r1_2.5min
  result[["r2_2.5min"]] = r2_2.5min
  result[["combined_2.5min"]] = combined_2.5min
  #result[["r1_5min"]] = r1_5min
  #result[["r2_5min"]] = r2_5min
  #result[["combined_5min"]] = combined_5min
  return(result)
}
# run loading script for each window length used (loop over all Tsmooth values)
## NOTE: now the data are accessible as nested lists (i.e. List_Waves25bp[["NA"]][[1]])
List_Waves25bp = list()
for (Tsmooth in c("NA", 2, 5, 10, 20)){
  List_Waves25bp[[Tsmooth]] = LoadWaveRepData(r1File = paste(fig_dir, "Rep1/", 25, "bpWindowing/", Tsmooth, "_TSmooth/","cdk9_TCr1_BothWavesCalldat.RData", sep=""), 
                                            r2File = paste(fig_dir, "Rep2/", 25, "bpWindowing/", Tsmooth, "_TSmooth/", "cdk9_TCr2_BothWavesCalldat.RData", sep=""), 
                                            combinedFile = paste(fig_dir, "CombinedReps/", 25, "bpWindowing/", Tsmooth, "_TSmooth/", "cdk9_TC_BothWavesCalldat.RData", sep=""))
}
List_Waves50bp = list()
for (Tsmooth in c("NA", 2, 5, 10, 20)){
  List_Waves50bp[[Tsmooth]] = LoadWaveRepData(r1File = paste(fig_dir, "Rep1/", 50, "bpWindowing/", Tsmooth, "_TSmooth/", "cdk9_TCr1_BothWavesCalldat.RData", sep=""), 
                            r2File = paste(fig_dir, "Rep2/", 50, "bpWindowing/", Tsmooth, "_TSmooth/", "cdk9_TCr2_BothWavesCalldat.RData", sep=""), 
                            combinedFile = paste(fig_dir, "CombinedReps/", 50, "bpWindowing/", Tsmooth, "_TSmooth/", "cdk9_TC_BothWavesCalldat.RData", sep=""))
}
List_Waves100bp = list()
for (Tsmooth in c("NA", 2, 5, 10, 20)){
  List_Waves100bp[[Tsmooth]] = LoadWaveRepData(r1File = paste(fig_dir, "Rep1/", 100, "bpWindowing/", Tsmooth, "_TSmooth/", "cdk9_TCr1_BothWavesCalldat.RData", sep=""), 
                             r2File = paste(fig_dir, "Rep2/", 100, "bpWindowing/", Tsmooth, "_TSmooth/", "cdk9_TCr2_BothWavesCalldat.RData", sep=""), 
                             combinedFile = paste(fig_dir, "CombinedReps/", 100, "bpWindowing/", Tsmooth, "_TSmooth/", "cdk9_TC_BothWavesCalldat.RData", sep=""))
}
List_Waves250bp = list()
for (Tsmooth in c("NA", 2, 5, 10, 20)){
  List_Waves250bp[[Tsmooth]] = LoadWaveRepData(r1File = paste(fig_dir, "Rep1/", 250, "bpWindowing/", Tsmooth, "_TSmooth/", "cdk9_TCr1_BothWavesCalldat.RData", sep=""), 
                             r2File = paste(fig_dir, "Rep2/", 250, "bpWindowing/", Tsmooth, "_TSmooth/", "cdk9_TCr2_BothWavesCalldat.RData", sep=""), 
                             combinedFile = paste(fig_dir, "CombinedReps/", 250, "bpWindowing/", Tsmooth, "_TSmooth/", "cdk9_TC_BothWavesCalldat.RData", sep=""))
}
List_Waves500bp = list()
for (Tsmooth in c("NA", 2, 5, 10, 20)){
  List_Waves500bp[[Tsmooth]] = LoadWaveRepData(r1File = paste(fig_dir, "Rep1/", 500, "bpWindowing/", Tsmooth, "_TSmooth/", "cdk9_TCr1_BothWavesCalldat.RData", sep=""), 
                             r2File = paste(fig_dir, "Rep2/", 500, "bpWindowing/", Tsmooth, "_TSmooth/", "cdk9_TCr2_BothWavesCalldat.RData", sep=""), 
                             combinedFile = paste(fig_dir, "CombinedReps/", 500, "bpWindowing/", Tsmooth, "_TSmooth/", "cdk9_TC_BothWavesCalldat.RData", sep=""))
}
# Clearing Wave calls will be filtered by the following criteria: 
## All time points must have a successful wave call (KL divergence > 0 ) 
## Clearing Wave must not move toward TSS at any subsequent time point 
filterWaves = function(WaveDatList, GeneList, TSSoffsetList = TSSoffSet_List, upStreamList = upStream_List, Filt = "conservative"){
  result = list()
  allDat = cbind(data.frame(WaveDatList["combined_30sec"]), data.frame(WaveDatList["combined_1min"]), data.frame(WaveDatList["combined_2.5min"]))
  # make a row index (qualFilt) that determines all rows that satisfy the wave filtering criteria. 
  qualFilt =  (allDat$combined_30sec.EndWave <= allDat$combined_1min.EndWave & allDat$combined_1min.EndWave <= allDat$combined_2.5min.EndWave & 
               allDat$combined_30sec.KLdiv > 1 & allDat$combined_1min.KLdiv > 1 & allDat$combined_2.5min.KLdiv > 1)
  Clean_combined_30sec_Dat = allDat[qualFilt, c("combined_30sec.ID", "combined_30sec.StartWave", "combined_30sec.EndWave", "combined_30sec.Rate", "combined_30sec.KLdiv")]
  Clean_combined_1min_Dat = allDat[qualFilt, c("combined_1min.ID", "combined_1min.StartWave", "combined_1min.EndWave", "combined_1min.Rate", "combined_1min.KLdiv")]
  Clean_combined_2.5min_Dat = allDat[qualFilt, c("combined_2.5min.ID", "combined_2.5min.StartWave", "combined_2.5min.EndWave", "combined_2.5min.Rate", "combined_2.5min.KLdiv")]
  #Clean_combined_5min_Dat = allDat[qualFilt5min, c("combined_5min.ID", "combined_5min.StartWave", "combined_5min.EndWave", "combined_5min.Rate", "combined_5min.KLdiv")]
  ## rather than just having the lengths of "spacer" between advancing and clearing waves, we also want wave positions in the genome. 
  ## function below returns information about the genomic positions of both advancing and clearing waves. 
  returnWaveLocs <- function(Dat, genes, timePoint = "30sec"){
    mergeDat = merge(Dat, genes, by.x = 1, by.y = 5) # note that we have reordered the genelist columns
    plusStrand = mergeDat$V6 == "+"
    minusStrand = mergeDat$V6 == "-"
    TSSoff = TSSoffsetList[[timePoint]]
    upstreamSearch = upStreamList[[timePoint]]
    mergeDat$AdvWaveEnd = 0 ## placeholder so I can compute on plus and minus strand genes separately
    mergeDat$ClearWaveStart = 0 ## placeholder
    mergeDat$AdvWaveDist = 0 ## placeholder so I can compute on plus and minus strand genes separately
    mergeDat$ClearWaveDist = 0 ## placeholder
    mergeDat[plusStrand,]$AdvWaveEnd = (mergeDat[plusStrand,]$V2 + TSSoff + mergeDat[plusStrand, 2] - upstreamSearch)
    mergeDat[minusStrand,]$AdvWaveEnd = (mergeDat[minusStrand,]$V3 - TSSoff - mergeDat[minusStrand, 2] + upstreamSearch)
    mergeDat[plusStrand,]$ClearWaveStart = (mergeDat[plusStrand,]$V2 + TSSoff + mergeDat[plusStrand, 3] - upstreamSearch)
    mergeDat[minusStrand,]$ClearWaveStart = (mergeDat[minusStrand,]$V3 - TSSoff - mergeDat[minusStrand, 3] + upstreamSearch)
    # Get Wave boundary distances relative to TSS
    mergeDat[plusStrand,]$AdvWaveDist = mergeDat[plusStrand,]$AdvWaveEnd - mergeDat[plusStrand,]$V2
    mergeDat[minusStrand,]$AdvWaveDist = mergeDat[minusStrand,]$V3 - mergeDat[minusStrand,]$AdvWaveEnd
    mergeDat[plusStrand,]$ClearWaveDist = mergeDat[plusStrand,]$ClearWaveStart - mergeDat[plusStrand,]$V2
    mergeDat[minusStrand,]$ClearWaveDist = mergeDat[minusStrand,]$V3 - mergeDat[minusStrand,]$ClearWaveStart
    res = mergeDat[,c(1,2,3,4,5,6,7,8,9,11,12,13,14)]
    return(res)
  }
  Clean_combined_30sec_Dat_exp = returnWaveLocs(Dat = Clean_combined_30sec_Dat, genes = GeneList, timePoint = "30sec")
  Clean_combined_1min_Dat_exp = returnWaveLocs(Dat = Clean_combined_1min_Dat, genes = GeneList, timePoint = "1min")
  Clean_combined_2.5min_Dat_exp = returnWaveLocs(Dat = Clean_combined_2.5min_Dat, genes = GeneList, timePoint = "2.5min")
  #Clean_combined_5min_Dat_exp = returnWaveLocs(Dat = Clean_combined_5min_Dat, genes = GeneList, timePoint = "5min")
  result[["30sec"]] = Clean_combined_30sec_Dat_exp[complete.cases(Clean_combined_30sec_Dat_exp),]
  result[["1min"]] = Clean_combined_1min_Dat_exp[complete.cases(Clean_combined_1min_Dat_exp),]
  result[["2.5min"]] = Clean_combined_2.5min_Dat_exp[complete.cases(Clean_combined_2.5min_Dat_exp),]
  #result[["5min"]] = Clean_combined_5min_Dat_exp[complete.cases(Clean_combined_5min_Dat_exp),]
  return(result)
}

###################################################################################################################
## function for plotting histograms of wave front calls for each time point
## input is the CleanWaves output
require(ggplot2)
library(plyr) # for calculating the group means in function below
histWaveFronts <- function(FiltWaveDat, filename = "HistWaveDists.pdf", Bins = 30, waveType = "advancing"){
  print("waveType must be set to: advancing, clearing, or between")
  if (waveType == "advancing"){pos = 12}
  else if (waveType == "clearing"){pos = 13}
  else if (waveType == "between"){pos = 4}
  dists_30sec = cbind(FiltWaveDat[["30sec"]][,pos], 0.5)
  dists_1min = cbind(FiltWaveDat[["1min"]][,pos], 1)
  dists_2.5min = cbind(FiltWaveDat[["2.5min"]][,pos], 2.5)
  #dists_5min = cbind(FiltWaveDat[["5min"]][,pos], 5)
  result = data.frame(rbind(dists_30sec, dists_1min, dists_2.5min))
  colnames(result) <- c("WaveDist", "TimePoint")
  mu <- ddply(result, "TimePoint", summarise, grp.mean=mean(WaveDist)) # get mean wave distance for each time point
  p<-ggplot(result, aes(x=WaveDist, fill=TimePoint, color=TimePoint, alpha = .5)) + 
    geom_histogram(bins = Bins) +
    facet_grid(TimePoint ~ .) +
    theme(legend.position="none") +
    geom_vline(data = mu, aes(xintercept=grp.mean, color=TimePoint),
               linetype="dashed") + 
    xlab("Wave Distance (bp)") + 
    ylab("Number of Genes")
  ggsave(filename = paste(fig_dir, filename, sep = ""), plot = p, width = 5, height = 8)
  return(result)
}
# test above plotting 
#x = histWaveFronts(Clean50bpWaves, filename = "test.pdf")

### Determining the average rate of Pol II in advancing wave: 
## function for regressing distance traveled and time elapsed.
getPolIIrate <- function(FiltWaveDat, filename = "scatterWaveDistsVsTime.pdf", waveType = "advancing"){
  print("waveType should be either 'advancing' or 'clearing'")
  if (waveType == "advancing"){pos = 12}
  if (waveType == "clearing"){pos = 13}
  dists_30sec = cbind(FiltWaveDat[["30sec"]][,pos], 0.5)
  dists_1min = cbind(FiltWaveDat[["1min"]][,pos], 1)
  dists_2.5min = cbind(FiltWaveDat[["2.5min"]][,pos], 2.5)
  #dists_5min = cbind(FiltWaveDat[["5min"]][,4], 5)
  result = data.frame(rbind(dists_30sec, dists_1min, dists_2.5min))
  colnames(result) <- c("WaveDist", "TimePoint")
  lm_eqn = function(result){ # function for grabbing the equation of the line from regression of WaveDist and TimePoint
    m = lm(WaveDist ~ TimePoint, result);
    eq <- substitute(italic(y) == b %.% italic(x) + a,#*","~~italic(r)^2~"="~r2, 
                     list(a = format(coef(m)[1], digits = 2), 
                          b = format(coef(m)[2], digits = 2))) 
    #r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));                 
  }
  eq <- data.frame(V1 = lm_eqn(result))
  p <- ggplot(result, aes(x=TimePoint, y=WaveDist)) +
    geom_point(shape=1, alpha = 0.2) + 
    geom_smooth(method=lm) + # adds linear regression line with 95% confidence interval
    xlab("Time (mins)") + 
    ylab("Wave Distance from TSS (bp)")
  p1 = p + geom_text(data=eq,aes(x = 1, y = 4000,label=V1), parse = TRUE, inherit.aes=FALSE)
  ggsave(filename = paste(fig_dir, filename, sep = ""), plot = p1, width = 3, height = 3)
  return(result)
}
# test above plotting 
#y = getPolIIrate(Clean50bpWaves, filename = "testScatter.pdf")

## Function for returning a geneList with the Coordinates of the Wave calls (in .bed format) for visualization in IGV
WaveCoordinates <- function(genes = filteredGL, FiltWaveData, filePrefix = "WaveLocations/"){ ## writes a separate file of bed-formatted wave coordinates for each time point
  dir.create(paste(fig_dir, filePrefix, sep = ""))
  dat_30sec = FiltWaveData[[1]]
  dat_1min = FiltWaveData[[2]]
  dat_2.5min =  FiltWaveData[[3]]
  #dat_5min  = FiltWaveData[[4]]
  getBothWaveCoords <- function(Dat){ ## returns a list of dataframes.  one gives TSS to advancing Wave end, the other give Clear wave start to CPS, both formatted for geneLists
    plusStrand = Dat[,"V6"] == "+"
    minusStrand = Dat[,"V6"] == "-"
    AdvCoords = Dat
    ClrCoords = Dat
    AdvCoords$Advstart = 0 # set placeholder columns (col 14)
    AdvCoords$AdvEnd = 0 # (col 15)
    ClrCoords$Clrstart = 0 # (col 14)
    ClrCoords$ClrEnd = 0 # (col 15)
    AdvCoords[plusStrand,]$Advstart = AdvCoords[plusStrand,]$V2
    AdvCoords[minusStrand,]$Advstart = AdvCoords[minusStrand,]$AdvWaveEnd ## switch order of TSS and Wave end on minus strand
    AdvCoords[plusStrand,]$AdvEnd = AdvCoords[plusStrand,]$AdvWaveEnd
    AdvCoords[minusStrand,]$AdvEnd = AdvCoords[minusStrand,]$V3 ## switch order of TSS and Wave end on minus strand
    ClrCoords[plusStrand,]$Clrstart = ClrCoords[plusStrand,]$ClearWaveStart
    ClrCoords[minusStrand,]$Clrstart = ClrCoords[minusStrand,]$V2 ## switch order of CPS and Wave end on minus strand
    ClrCoords[plusStrand,]$ClrEnd = ClrCoords[plusStrand,]$V3
    ClrCoords[minusStrand,]$ClrEnd = ClrCoords[minusStrand,]$ClearWaveStart ## switch order of CPS and Wave end on minus strand
    AdvTable = AdvCoords[,c(6,14,15,1,12,9)]
    ClrTable = ClrCoords[,c(6,14,15,1,13,9)]
    return(list(AdvTable, ClrTable))
  }
  waves_30sec = getBothWaveCoords(Dat = dat_30sec)
  waves_1min = getBothWaveCoords(Dat = dat_1min)
  waves_2.5min = getBothWaveCoords(Dat = dat_2.5min)
  #waves_5min = getBothWaveCoords(Dat = dat_5min)
  write.table(waves_30sec[[1]], file = paste(fig_dir, filePrefix, "30secAdvWaves.bed", sep = ""), quote = F, col.names = T, row.names = F, sep = "\t")
  write.table(waves_30sec[[2]], file = paste(fig_dir, filePrefix, "30secClearWaves.bed", sep = ""), quote = F, col.names = T, row.names = F, sep = "\t")
  write.table(waves_1min[[1]], file = paste(fig_dir, filePrefix, "1minAdvWaves.bed", sep = ""), quote = F, col.names = T, row.names = F, sep = "\t")
  write.table(waves_1min[[2]], file = paste(fig_dir, filePrefix, "1minClearWaves.bed", sep = ""), quote = F, col.names = T, row.names = F, sep = "\t")
  write.table(waves_2.5min[[1]], file = paste(fig_dir, filePrefix, "2.5minAdvWaves.bed", sep = ""), quote = F, col.names = T, row.names = F, sep = "\t")
  write.table(waves_2.5min[[2]], file = paste(fig_dir, filePrefix, "2.5minClearWaves.bed", sep = ""), quote = F, col.names = T, row.names = F, sep = "\t")
  #write.table(waves_5min[[1]], file = paste(fig_dir, filePrefix, "5minAdvWaves.bed", sep = ""), quote = F, col.names = T, row.names = F, sep = "\t")
  #write.table(waves_5min[[2]], file = paste(fig_dir, filePrefix, "5minClearWaves.bed", sep = ""), quote = F, col.names = T, row.names = F, sep = "\t")
}

###################################################################################################################
## Run the processing and plotting functions for waves called using different window sizes: 
## Loop over all Tsmooth values: 
for (Tsmooth in c("NA", 2, 5, 10, 20)){
  # 25 bp window calls
  dir.create(paste(fig_dir, Tsmooth, "_TSmooth/", sep = ""))
  dir.create(paste(fig_dir, Tsmooth, "_TSmooth/", "figs_25bp_window_calls/", sep = ""))
  Clean25bpWaves = filterWaves(List_Waves25bp[[Tsmooth]], GeneList = GL, TSSoffsetList = TSSoffSet_List, upStreamList = upStream_List, Filt = "none") ## Note: Relaxed reduces the filtering of waves called based on individual replicates
  hist_25bpAdvWaves = histWaveFronts(Clean25bpWaves, filename = paste(Tsmooth, "_TSmooth/","figs_25bp_window_calls/HistDistsAdvancingWaves_25bpWin.pdf", sep = ""), Bins = 20, waveType = "advancing")
  hist_25bpClearingWaves = histWaveFronts(Clean25bpWaves, filename = paste(Tsmooth, "_TSmooth/","figs_25bp_window_calls/HistDistsClearingWaves_25bpWin.pdf", sep = ""), Bins = 20, waveType = "clearing")
  hist_25bpBetweenWaves = histWaveFronts(Clean25bpWaves, filename = paste(Tsmooth, "_TSmooth/","figs_25bp_window_calls/HistDistsBetweenWaves_25bpWin.pdf", sep = ""), Bins = 20, waveType = "between")
  scatter_25bpAdvWaves = getPolIIrate(Clean25bpWaves, filename = paste(Tsmooth, "_TSmooth/","figs_25bp_window_calls/scatterAdvWaveDistsVsTime_25bpWin.pdf", sep = ""), waveType = "advancing")
  scatter_25bpClrWaves = getPolIIrate(Clean25bpWaves, filename = paste(Tsmooth, "_TSmooth/","figs_25bp_window_calls/scatterClearWaveDistsVsTime_25bpWin.pdf", sep = ""), waveType = "clearing")
  WaveCoordinates(genes = filteredGL, Clean25bpWaves, filePrefix = paste(Tsmooth, "_TSmooth/", "WaveLocations_25bpWin/", sep = "")) # writes the bed files with "high confidence" wave coords for all time points 
  # 50 bp window calls
  dir.create(paste(fig_dir, Tsmooth, "_TSmooth/", "figs_50bp_window_calls/", sep = ""))
  Clean50bpWaves = filterWaves(List_Waves50bp[[Tsmooth]], GeneList = GL, TSSoffsetList = TSSoffSet_List, upStreamList = upStream_List, Filt = "none")
  hist_50bpAdvWaves = histWaveFronts(Clean50bpWaves, filename = paste(Tsmooth, "_TSmooth/","figs_50bp_window_calls/HistDistsAdvancingWaves_50bpWin.pdf", sep = ""), Bins = 20, waveType = "advancing")
  hist_50bpClearingWaves = histWaveFronts(Clean50bpWaves, filename = paste(Tsmooth, "_TSmooth/","figs_50bp_window_calls/HistDistsClearingWaves_50bpWin.pdf", sep = ""), Bins = 20, waveType = "clearing")
  hist_50bpBetweenWaves = histWaveFronts(Clean50bpWaves, filename = paste(Tsmooth, "_TSmooth/","figs_50bp_window_calls/HistDistsBetweenWaves_50bpWin.pdf", sep = ""), Bins = 20, waveType = "between")
  scatter_50bpAdvWaves = getPolIIrate(Clean50bpWaves, filename = paste(Tsmooth, "_TSmooth/","figs_50bp_window_calls/scatterAdvWaveDistsVsTime_50bpWin.pdf", sep = ""), waveType = "advancing")
  scatter_50bpClrWaves = getPolIIrate(Clean50bpWaves, filename = paste(Tsmooth, "_TSmooth/","figs_50bp_window_calls/scatterClearWaveDistsVsTime_50bpWin.pdf", sep = ""), waveType = "clearing")
  WaveCoordinates(genes = filteredGL, Clean50bpWaves, filePrefix = paste(Tsmooth, "_TSmooth/", "WaveLocations_50bpWin/", sep = "")) # writes the bed files with "high confidence" wave coords for all time points 
  # 100 bp window calls
  dir.create(paste(fig_dir, Tsmooth, "_TSmooth/", "figs_100bp_window_calls/", sep = ""))
  Clean100bpWaves = filterWaves(List_Waves100bp[[Tsmooth]], GeneList = GL, TSSoffsetList = TSSoffSet_List, upStreamList = upStream_List, Filt = "none")
  hist_100bpAdvancingWaves = histWaveFronts(Clean100bpWaves, filename = paste(Tsmooth, "_TSmooth/","figs_100bp_window_calls/HistDistsAdvancingWaves_100bpWin.pdf", sep = ""), Bins = 20, waveType = "advancing")
  hist_100bpClearingWaves = histWaveFronts(Clean100bpWaves, filename = paste(Tsmooth, "_TSmooth/","figs_100bp_window_calls/HistDistsClearingWaves_100bpWin.pdf", sep = ""), Bins = 20, waveType = "clearing")
  hist_100bpBetweenWaves = histWaveFronts(Clean100bpWaves, filename = paste(Tsmooth, "_TSmooth/","figs_100bp_window_calls/HistDistsBetweenWaves_100bpWin.pdf", sep = ""), Bins = 20, waveType = "between")
  scatter_100bpAdvWaves = getPolIIrate(Clean100bpWaves, filename = paste(Tsmooth, "_TSmooth/","figs_100bp_window_calls/scatterAdvWaveDistsVsTime_100bpWin.pdf", sep = ""), waveType = "advancing")
  scatter_100bpClrWaves = getPolIIrate(Clean100bpWaves, filename = paste(Tsmooth, "_TSmooth/","figs_100bp_window_calls/scatterClearWaveDistsVsTime_100bpWin.pdf", sep = ""), waveType = "clearing")
  WaveCoordinates(genes = filteredGL, Clean100bpWaves, filePrefix = paste(Tsmooth, "_TSmooth/", "WaveLocations_100bpWin/", sep = ""))
  # 250 bp window calls
  dir.create(paste(fig_dir, Tsmooth, "_TSmooth/", "figs_250bp_window_calls/", sep = ""))
  Clean250bpWaves = filterWaves(List_Waves250bp[[Tsmooth]], GeneList = GL, TSSoffsetList = TSSoffSet_List, upStreamList = upStream_List, Filt = "none")
  hist_250bpAdvancingWaves = histWaveFronts(Clean250bpWaves, filename = paste(Tsmooth, "_TSmooth/","figs_250bp_window_calls/HistDistsAdvancingWaves_250bpWin.pdf", sep = ""), Bins = 20, waveType = "advancing")
  hist_250bpClearingWaves = histWaveFronts(Clean250bpWaves, filename = paste(Tsmooth, "_TSmooth/","figs_250bp_window_calls/HistDistsClearingWaves_250bpWin.pdf", sep = ""), Bins = 20, waveType = "clearing")
  hist_250bpBetweenWaves = histWaveFronts(Clean250bpWaves, filename = paste(Tsmooth, "_TSmooth/","figs_250bp_window_calls/HistDistsBetweenWaves_250bpWin.pdf", sep = ""), Bins = 20, waveType = "between")
  scatter_250bpAdvWaves = getPolIIrate(Clean250bpWaves, filename = paste(Tsmooth, "_TSmooth/","figs_250bp_window_calls/scatterAdvWaveDistsVsTime_250bpWin.pdf", sep = ""), waveType = "advancing")
  scatter_250bpClrWaves = getPolIIrate(Clean250bpWaves, filename = paste(Tsmooth, "_TSmooth/","figs_250bp_window_calls/scatterClearWaveDistsVsTime_250bpWin.pdf", sep = ""), waveType = "clearing")
  WaveCoordinates(genes = filteredGL, Clean250bpWaves, filePrefix = paste(Tsmooth, "_TSmooth/", "WaveLocations_250bpWin/", sep = "")) 
  # 500 bp window calls
  dir.create(paste(fig_dir, Tsmooth, "_TSmooth/", "figs_500bp_window_calls/", sep = ""))
  Clean500bpWaves = filterWaves(List_Waves500bp[[Tsmooth]], GeneList = GL, TSSoffsetList = TSSoffSet_List, upStreamList = upStream_List, Filt = "none")
  hist_500bpAdvancingWaves = histWaveFronts(Clean500bpWaves, filename = paste(Tsmooth, "_TSmooth/","figs_500bp_window_calls/HistDistsAdvancingWaves_500bpWin.pdf", sep = ""), Bins = 10, waveType = "advancing")
  hist_500bpClearingWaves = histWaveFronts(Clean500bpWaves, filename = paste(Tsmooth, "_TSmooth/","figs_500bp_window_calls/HistDistsClearingWaves_500bpWin.pdf", sep = ""), Bins = 10, waveType = "clearing")
  hist_500bpBetweenWaves = histWaveFronts(Clean500bpWaves, filename = paste(Tsmooth, "_TSmooth/","figs_500bp_window_calls/HistDistsBetweenWaves_500bpWin.pdf", sep = ""), Bins = 10, waveType = "between")
  scatter_500bpAdvWaves = getPolIIrate(Clean500bpWaves, filename = paste(Tsmooth, "_TSmooth/","figs_500bp_window_calls/scatterAdvWaveDistsVsTime_500bpWin.pdf", sep = ""), waveType = "advancing")
  scatter_500bpClrWaves = getPolIIrate(Clean500bpWaves, filename = paste(Tsmooth, "_TSmooth/","figs_500bp_window_calls/scatterClearWaveDistsVsTime_500bpWin.pdf", sep = ""), waveType = "clearing")
  WaveCoordinates(genes = filteredGL, Clean500bpWaves, filePrefix = paste(Tsmooth, "_TSmooth/", "WaveLocations_500bpWin/", sep = ""))
}
