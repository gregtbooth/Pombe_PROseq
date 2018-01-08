## Charles' Wave calling scripts were run on each biological replicate data for the cdk9as time course:
## it was also run using different size windows (larger windows were suggested to improve performance)
fig_dir = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/figures/cdk9asTC_waves/10-05-17/fixedApproxDist_2kb/"
wave_dir = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/figures/cdk9asTC_waves/10-05-17/fixedApproxDist_2kb/"
dir.create(fig_dir)
infopath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/pombeinfo/"
filteredGL = read.table(paste(infopath, "SP_CompleteFilteredGenes.bed", sep = ""))
GL = filteredGL[,c(1,2,3,6,4,5)]
idx = (GL[,3] - GL[,2] >= 4000)
# try using only the genes used in initial paper (n = )

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
  r1_2.5min = data.frame(pw_2.5min_r1[length(pw_2.5min_r1)])
  r2_2.5min = data.frame(pw_2.5min_r2[length(pw_2.5min_r2)])
  combined_2.5min = data.frame(pw_2.5min[length(pw_2.5min)])
  # 1 min
  r1_5min = data.frame(pw_5min_r1[length(pw_5min_r1)])
  r2_5min = data.frame(pw_5min_r2[length(pw_5min_r2)])
  combined_5min = data.frame(pw_5min[length(pw_5min)])
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
  result[["r1_5min"]] = r1_5min
  result[["r2_5min"]] = r2_5min
  result[["combined_5min"]] = combined_5min
  return(result)
}
# run loading script for each window length used
Waves50bp = LoadWaveRepData(r1File = paste(wave_dir, "/Rep1/", 50, "bpWindowing/", "cdk9_TCr1_WaveCalldat.RData", sep=""), 
                            r2File = paste(wave_dir, "/Rep2/", 50, "bpWindowing/", "cdk9_TCr2_WaveCalldat.RData", sep=""), 
                            combinedFile = paste(wave_dir, "/CombinedReps/", 50, "bpWindowing/", "cdk9_TC_WaveCalldat.RData", sep=""))
Waves100bp = LoadWaveRepData(r1File = paste(wave_dir, "/Rep1/", 100, "bpWindowing/", "cdk9_TCr1_WaveCalldat.RData", sep=""), 
                            r2File = paste(wave_dir, "/Rep2/", 100, "bpWindowing/", "cdk9_TCr2_WaveCalldat.RData", sep=""), 
                            combinedFile = paste(wave_dir, "/CombinedReps/", 100, "bpWindowing/", "cdk9_TC_WaveCalldat.RData", sep=""))
Waves250bp = LoadWaveRepData(r1File = paste(wave_dir, "/Rep1/", 250, "bpWindowing/", "cdk9_TCr1_WaveCalldat.RData", sep=""), 
                            r2File = paste(wave_dir, "/Rep2/", 250, "bpWindowing/", "cdk9_TCr2_WaveCalldat.RData", sep=""), 
                            combinedFile = paste(wave_dir, "/CombinedReps/", 250, "bpWindowing/", "cdk9_TC_WaveCalldat.RData", sep=""))
Waves500bp = LoadWaveRepData(r1File = paste(wave_dir, "/Rep1/", 500, "bpWindowing/", "cdk9_TCr1_WaveCalldat.RData", sep=""), 
                            r2File = paste(wave_dir, "/Rep2/", 500, "bpWindowing/", "cdk9_TCr2_WaveCalldat.RData", sep=""), 
                            combinedFile = paste(wave_dir, "/CombinedReps/", 500, "bpWindowing/", "cdk9_TC_WaveCalldat.RData", sep=""))

# get only genes which have TRUE for minOfMax and minOfAvg for both replicates at all time points. 
filterWaves = function(WaveDatList, GeneList, TSSoff = 500){
  result = list()
  allDat = cbind(data.frame(WaveDatList["r1_30sec"]), data.frame(WaveDatList["r2_30sec"]), data.frame(WaveDatList["combined_30sec"]),
                 data.frame(WaveDatList["r1_1min"]), data.frame(WaveDatList["r2_1min"]), data.frame(WaveDatList["combined_1min"]),
                 data.frame(WaveDatList["r1_2.5min"]), data.frame(WaveDatList["r2_2.5min"]), data.frame(WaveDatList["combined_2.5min"]),
                 data.frame(WaveDatList["r1_5min"]), data.frame(WaveDatList["r2_5min"]), data.frame(WaveDatList["combined_5min"]))
  # quality filter:
  ## All time points must have a successful wave call (minOfMax and minOfAvg are TRUE for combined reps ) 
  ## Advancing wave must move forward at subsequent time points
  qualFiltAll =  (allDat$combined_30sec.EndWave <= allDat$combined_1min.EndWave & 
                    allDat$combined_1min.EndWave <= allDat$combined_2.5min.EndWave & 
                    allDat$combined_2.5min.EndWave <= allDat$combined_5min.EndWave &
                    allDat$combined_30sec.minOfMax == TRUE & allDat$combined_30sec.minOfAvg == TRUE &
                    allDat$combined_1min.minOfMax == TRUE & allDat$combined_1min.minOfAvg == TRUE &
                    allDat$combined_2.5min.minOfMax == TRUE & allDat$combined_2.5min.minOfAvg == TRUE &
                    allDat$combined_5min.minOfMax == TRUE & allDat$combined_5min.minOfAvg == TRUE)
  #qualFilt30sec =  (#allDat$r1_30sec.minOfMax == TRUE & allDat$r1_30sec.minOfAvg == TRUE & 
                    #allDat$r2_30sec.minOfMax == TRUE & allDat$r2_30sec.minOfAvg == TRUE &
  #                  allDat$combined_30sec.minOfMax == TRUE & allDat$combined_30sec.minOfAvg == TRUE)
  #qualFilt1min =  (#allDat$r1_1min.minOfMax == TRUE & allDat$r1_1min.minOfAvg == TRUE & 
                      #allDat$r2_1min.minOfMax == TRUE & allDat$r2_1min.minOfAvg == TRUE &
  #                    allDat$combined_1min.minOfMax == TRUE & allDat$combined_1min.minOfAvg == TRUE)
  #qualFilt2.5min =  (#allDat$r1_2.5min.minOfMax == TRUE & allDat$r1_2.5min.minOfAvg == TRUE & 
                      #allDat$r2_2.5min.minOfMax == TRUE & allDat$r2_2.5min.minOfAvg == TRUE &
  #                    allDat$combined_2.5min.minOfMax == TRUE & allDat$combined_2.5min.minOfAvg == TRUE)
  #qualFilt5min =  (#allDat$r1_5min.minOfMax == TRUE & allDat$r1_5min.minOfAvg == TRUE & 
                    #  allDat$r2_5min.minOfMax == TRUE & allDat$r2_5min.minOfAvg == TRUE &
   #                   allDat$combined_5min.minOfMax == TRUE & allDat$combined_5min.minOfAvg == TRUE)
  Clean_combined_30sec_Dat = allDat[qualFiltAll, c("combined_30sec.ID", "combined_30sec.StartWave", "combined_30sec.EndWave", "combined_30sec.Rate", "combined_30sec.KLdiv")]
  Clean_combined_1min_Dat = allDat[qualFiltAll, c("combined_1min.ID", "combined_1min.StartWave", "combined_1min.EndWave", "combined_1min.Rate", "combined_1min.KLdiv")]
  Clean_combined_2.5min_Dat = allDat[qualFiltAll, c("combined_2.5min.ID", "combined_2.5min.StartWave", "combined_2.5min.EndWave", "combined_2.5min.Rate", "combined_2.5min.KLdiv")]
  Clean_combined_5min_Dat = allDat[qualFiltAll, c("combined_5min.ID", "combined_5min.StartWave", "combined_5min.EndWave", "combined_5min.Rate", "combined_5min.KLdiv")]
  returnWaveLocs <- function(Dat, genes, TSSoff){
    mergeDat = merge(Dat, genes, by.x = 1, by.y = 5) # note that we have reordered the genelist columns
    plusStrand = mergeDat$V6 == "+"
    minusStrand = mergeDat$V6 == "-"
    mergeDat$AdvWaveEnd = 0 ## placeholder so I can compute on plus and minus strand genes separately
    mergeDat$AdvWaveDist = 0 ## placeholder so I can compute on plus and minus strand genes separately
    mergeDat[plusStrand,]$AdvWaveEnd = (mergeDat[plusStrand,]$V2 + mergeDat[plusStrand, 2] + mergeDat[plusStrand, 3] - TSSoff)
    mergeDat[minusStrand,]$AdvWaveEnd = (mergeDat[minusStrand,]$V3 - mergeDat[minusStrand, 2] - mergeDat[minusStrand, 3] + TSSoff)
    # Get Wave boundary distances relative to TSS
    mergeDat[plusStrand,]$AdvWaveDist = mergeDat[plusStrand,]$AdvWaveEnd - mergeDat[plusStrand,]$V2
    mergeDat[minusStrand,]$AdvWaveDist = mergeDat[minusStrand,]$V3 - mergeDat[minusStrand,]$AdvWaveEnd
    res = mergeDat[,c(1,2,3,4,5,6,7,8,9,11,12)]
    return(res)
  }
  Clean_combined_30sec_Dat_exp = returnWaveLocs(Dat = Clean_combined_30sec_Dat, genes = GeneList, TSSoff)
  Clean_combined_1min_Dat_exp = returnWaveLocs(Dat = Clean_combined_1min_Dat, genes = GeneList, TSSoff)
  Clean_combined_2.5min_Dat_exp = returnWaveLocs(Dat = Clean_combined_2.5min_Dat, genes = GeneList, TSSoff)
  Clean_combined_5min_Dat_exp = returnWaveLocs(Dat = Clean_combined_5min_Dat, genes = GeneList, TSSoff)
  result[["30sec"]] = Clean_combined_30sec_Dat_exp[complete.cases(Clean_combined_30sec_Dat_exp),]
  result[["1min"]] = Clean_combined_1min_Dat_exp[complete.cases(Clean_combined_1min_Dat_exp),]
  result[["2.5min"]] = Clean_combined_2.5min_Dat_exp[complete.cases(Clean_combined_2.5min_Dat_exp),]
  result[["5min"]] = Clean_combined_5min_Dat_exp[complete.cases(Clean_combined_5min_Dat_exp),]
  return(result)
}
  
#Clean50bpWaves = filterWaves(Waves50bp, GeneList = GL, TSSoff = 500)

###################################################################################################################
# scripts for analysing cleaned up data
###################################################################################################################
## function for plotting histograms of wave front calls for each time point
## input is the CleanWaves output
require(ggplot2)
library(plyr) # for calculating the group means in function below
histWaveFronts <- function(FiltWaveDat, filename = "HistWaveDists.pdf", Bins = 30){
  dists_30sec = cbind(FiltWaveDat[["30sec"]][,11], 0.5)
  dists_1min = cbind(FiltWaveDat[["1min"]][,11], 1)
  dists_2.5min = cbind(FiltWaveDat[["2.5min"]][,11], 2.5)
  dists_5min = cbind(FiltWaveDat[["5min"]][,11], 5)
  result = data.frame(rbind(dists_30sec, dists_1min, dists_2.5min, dists_5min))
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
getPolIIrate <- function(FiltWaveDat, filename = "scatterWaveDistsVsTime.pdf"){
  dists_30sec = cbind(FiltWaveDat[["30sec"]][,11], 0.5)
  dists_1min = cbind(FiltWaveDat[["1min"]][,11], 1)
  dists_2.5min = cbind(FiltWaveDat[["2.5min"]][,11], 2.5)
  dists_5min = cbind(FiltWaveDat[["5min"]][,11], 5)
  result = data.frame(rbind(dists_30sec, dists_1min, dists_2.5min, dists_5min))
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
    ylab("Wave Distance (bp)")
  p1 = p + geom_text(data=eq,aes(x = 3, y = 4900,label=V1), parse = TRUE, inherit.aes=FALSE)
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
  dat_5min  = FiltWaveData[[4]]
  getWaveCoords <- function(Dat){ ## returns a list of dataframes.  one gives TSS to advancing Wave end, the other give Clear wave start to CPS, both formatted for geneLists
    plusStrand = Dat[,"V6"] == "+"
    minusStrand = Dat[,"V6"] == "-"
    AdvCoords = Dat
    AdvCoords$Advstart = 0 # set placeholder columns (col 14)
    AdvCoords$AdvEnd = 0 # (col 15)
    AdvCoords[plusStrand,]$Advstart = AdvCoords[plusStrand,]$V2
    AdvCoords[minusStrand,]$Advstart = AdvCoords[minusStrand,]$AdvWaveEnd ## switch order of TSS and Wave end on minus strand
    AdvCoords[plusStrand,]$AdvEnd = AdvCoords[plusStrand,]$AdvWaveEnd
    AdvCoords[minusStrand,]$AdvEnd = AdvCoords[minusStrand,]$V3 ## switch order of TSS and Wave end on minus strand
    AdvTable = AdvCoords[,c(6,12,13,1,11,9)]
    return(AdvTable)
  }
  waves_30sec = getWaveCoords(Dat = dat_30sec)
  waves_1min = getWaveCoords(Dat = dat_1min)
  waves_2.5min = getWaveCoords(Dat = dat_2.5min)
  waves_5min = getWaveCoords(Dat = dat_5min)
  write.table(waves_30sec, file = paste(fig_dir, filePrefix, "30secAdvWaves.bed", sep = ""), quote = F, col.names = T, row.names = F, sep = "\t")
  write.table(waves_1min, file = paste(fig_dir, filePrefix, "1minAdvWaves.bed", sep = ""), quote = F, col.names = T, row.names = F, sep = "\t")
  write.table(waves_2.5min, file = paste(fig_dir, filePrefix, "2.5minAdvWaves.bed", sep = ""), quote = F, col.names = T, row.names = F, sep = "\t")
  write.table(waves_5min, file = paste(fig_dir, filePrefix, "5minAdvWaves.bed", sep = ""), quote = F, col.names = T, row.names = F, sep = "\t")
}


###################################################################################################################
## Run the processing and plotting functions for waves called using different window sizes: 
# 50 bp window calls
dir.create(paste(fig_dir,"figs_50bp_window_calls/", sep = ""))
Clean50bpWaves = filterWaves(Waves50bp, GeneList = GL, TSSoff = 500)
hist_50bpWaves = histWaveFronts(Clean50bpWaves, filename = "figs_50bp_window_calls/HistWaveDists_50bpWin.pdf")
scatter_50bpWaves = getPolIIrate(Clean50bpWaves, filename = "figs_50bp_window_calls/scatterWaveDistsVsTime_50bpWin.pdf")
WaveCoordinates(genes = filteredGL, Clean50bpWaves, filePrefix = "WaveLocations_50bpWin/") # writes the bed files with "high confidence" wave coords for all time points 
# 100 bp window calls
dir.create(paste(fig_dir,"figs_100bp_window_calls/", sep = ""))
Clean100bpWaves = filterWaves(Waves100bp, GeneList = GL, TSSoff = 500)
hist_100bpWaves = histWaveFronts(Clean100bpWaves, filename = "figs_100bp_window_calls/HistWaveDists_100bpWin.pdf")
scatter_100bpWaves = getPolIIrate(Clean100bpWaves, filename = "figs_100bp_window_calls/scatterWaveDistsVsTime_100bpWin.pdf")
WaveCoordinates(genes = filteredGL, Clean100bpWaves, filePrefix = "WaveLocations_100bpWin/")
# 250 bp window calls
dir.create(paste(fig_dir,"figs_250bp_window_calls/", sep = ""))
Clean250bpWaves = filterWaves(Waves250bp, GeneList = GL, TSSoff = 500)
hist_250bpWaves = histWaveFronts(Clean250bpWaves, filename = "figs_250bp_window_calls/HistWaveDists_250bpWin.pdf", Bins = 20)
scatter_250bpWaves = getPolIIrate(Clean250bpWaves, filename = "figs_250bp_window_calls/scatterWaveDistsVsTime_250bpWin.pdf")
WaveCoordinates(genes = filteredGL, Clean250bpWaves, filePrefix = "WaveLocations_250bpWin/") ## Every wave call using this window has a 250 bp offset, whereas with others, 500 bp (upstreamDist used) is usually called
# 500 bp window calls
dir.create(paste(fig_dir,"figs_500bp_window_calls/", sep = ""))
Clean500bpWaves = filterWaves(Waves500bp, GeneList = GL, TSSoff = 500)
hist_500bpWaves = histWaveFronts(Clean500bpWaves, filename = "figs_500bp_window_calls/HistWaveDists_500bpWin.pdf", Bins = 10)
scatter_500bpWaves = getPolIIrate(Clean500bpWaves, filename = "figs_500bp_window_calls/scatterWaveDistsVsTime_500bpWin.pdf")
WaveCoordinates(genes = filteredGL, Clean500bpWaves, filePrefix = "WaveLocations_500bpWin/")


