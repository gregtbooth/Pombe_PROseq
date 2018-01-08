## function is to be used to calculate a confidence interval for rate (slope) from wave location files.  
## script designed to look in specified director and grab files for each time point.
### should be able to perform function on various wave call parameters by just changing "basepath" to a different wavelocation folder
fig_dir = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/figures/cdk9asTC_waves/10-05-17/RateVariance/"
dir.create(fig_dir)
basepath_Adv = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/figures/cdk9asTC_waves/10-05-17/fixedApproxDist_2kb/WaveLocations_50bpWin/"
basepath_clear ="/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/figures/cdk9asTC_waves/10-05-17/BothWaves/TryTSmooth_6Kbgenes/fixedApprox_3kb/5_TSmooth/WaveLocations_50bpWin/"
  
formatTimePointDists = function(filePath, Advancing = T){
  datList = list()
  if (Advancing){
    datList[["waveDat_0.5"]] = read.table(paste(filePath, "30secAdvWaves.bed", sep = ""), head = T)
    datList[["waveDat_1"]] = read.table(paste(filePath, "1minAdvWaves.bed", sep = ""), head = T)
    datList[["waveDat_2.5"]] = read.table(paste(filePath, "2.5minAdvWaves.bed", sep = ""), head = T)
    datList[["waveDat_5"]] = read.table(paste(filePath, "5minAdvWaves.bed", sep = ""), head = T)
  }
  else{
    datList[["waveDat_0.5"]] = read.table(paste(filePath, "30secClearWaves.bed", sep = ""), head = T)
    datList[["waveDat_1"]] = read.table(paste(filePath, "1minClearWaves.bed", sep = ""), head = T)
    datList[["waveDat_2.5"]] = read.table(paste(filePath, "2.5minClearWaves.bed", sep = ""), head = T)
  }
  TimeTable = data.frame(datList[["waveDat_0.5"]][,4])
  for(dat in 1:length(datList)){
  TimeTable = cbind(TimeTable, datList[[dat]][,5])
  }
  if(dim(TimeTable)[2] == 4){colnames(TimeTable) <- c("gene","30secDist", "1minDist", "2.5minDist")}
  else{colnames(TimeTable) <- c("gene","30secDist", "1minDist", "2.5minDist", "5minDist")}
  return(TimeTable)
}

permuteAvg = function(data, nPermute = 10){
  N = dim(data)[1]
  sampleFrac = 0.1
  M = as.integer(round(N * sampleFrac, 0))
  result = matrix(nrow = nPermute, ncol = (dim(data)[2]-1))
  for (i in 1:nPermute) {
    idx <- sample(N, size = M, replace = T) # gives M random row numbers.
    result[i, ] = colSums(data[idx, c(2:dim(data)[2])])/M
  }
return(result)
}
# SlopeVar calculates a slope for each gene (Standard Slope) and for averages of bootstrapped sets of genes (permuted)
SlopeVar = function(RawRes, permuteRes, Advancing = T){
  rawDat = RawRes[,2:dim(RawRes)[2]] # get rid of first col
  if (Advancing){timepoints = c(0.5, 1, 2.5, 5)}
  else{timepoints = c(0.5, 1, 2.5)}
  calcSlope = function(row){lm(row ~ timepoints)$coefficients[[2]]} # grabs the slope
  bootstrappedRes = apply(X = permuteRes, MARGIN = 1, FUN = calcSlope)
  StandardRes = apply(X = rawDat, MARGIN = 1, FUN = calcSlope)
  bootstrappedSlopes = cbind(permuteRes, bootstrappedRes)
  StandardSlopes = cbind(RawRes, StandardRes)
  result = list("standard" = StandardSlopes, "bootStrapped" = bootstrappedSlopes)
  return(result)
}

SlopeCI = function(slopeVar_res, qLow = 0.05, qHigh = 0.95){
  #quartiles from geneBYgene variance
  lowerQ = quantile(slopeVar_res[["standard"]][,dim(slopeVar_res[["standard"]])[2]], qLow)
  med = quantile(slopeVar_res[["standard"]][,dim(slopeVar_res[["standard"]])[2]], 0.5)
  upperQ =quantile(slopeVar_res[["standard"]][,dim(slopeVar_res[["standard"]])[2]], qHigh)
  res = list("lower" = lowerQ, "median" = med,"upper" = upperQ)
  #quartiles from permuted variance
  lowerQ_perm = quantile(slopeVar_res[["bootStrapped"]][,dim(slopeVar_res[["bootStrapped"]])[2]], qLow)
  med_perm = quantile(slopeVar_res[["bootStrapped"]][,dim(slopeVar_res[["bootStrapped"]])[2]], 0.5)
  upperQ_perm =quantile(slopeVar_res[["bootStrapped"]][,dim(slopeVar_res[["bootStrapped"]])[2]], qHigh)
  res_perm = list("lower" = lowerQ_perm, "median" = med_perm,"upper" = upperQ_perm)
  result = list("genebygeneQuantiles" = res, "bootstrappedQuantiles" = res_perm)
  return(result)
}

# box plot of distributions of rate measurments for advancing and clearing waves
require(ggplot2)
rateBoxPlots = function(adv_dat, clr_dat, filename = "testBoxPlot.pdf", bootstrapOnly = T, labs = c("Bootstrapped")){
  advSlopes_bootSt = cbind(adv_dat[[2]][,dim(adv_dat[[2]])[2]], 1, 1)
  advSlopes_gbg = cbind(adv_dat[[1]][,dim(adv_dat[[1]])[2]], 1, 2)
  clrSlopes_bootSt = cbind(clr_dat[[2]][,dim(clr_dat[[2]])[2]], 2, 1)
  clrSlopes_gbg = cbind(clr_dat[[1]][,dim(clr_dat[[1]])[2]], 2, 2)
  if (bootstrapOnly){reformatRes = rbind(advSlopes_bootSt, clrSlopes_bootSt)}
  else{reformatRes = rbind(advSlopes_bootSt, advSlopes_gbg, clrSlopes_bootSt, clrSlopes_gbg)}
  colnames(reformatRes)<- c("Slope", "Adv_clr", "gbg_bootst")
  p <- ggplot(data = data.frame(reformatRes), aes(x = factor(Adv_clr, labels = c("Adv. Rate", "Clr. Rate")), y = Slope, 
                                    fill = factor(gbg_bootst, labels = labs))) + 
    geom_boxplot() + 
    xlab("Wave") + 
    ylab("Rate (bp/min)") +
    labs(fill = "") +# gets rid of legend name
    scale_fill_brewer(palette = "Accent")
  ggsave(filename = paste(fig_dir, filename, sep = ""), plot = p, width = 4, height = 4)
  return(reformatRes) 
}

## Advancing Wave
Adv_formatted = formatTimePointDists(filePath = basepath_Adv, Advancing = T)
Adv_Permuted = permuteAvg(Adv_formatted, nPermute = 1000)
Adv_slopes = SlopeVar(RawRes = Adv_formatted, permuteRes = Adv_Permuted)
Adv_SlopeCI = SlopeCI(Adv_slopes, qLow = 0.05, qHigh = 0.95)

## Clearing Wave
Clr_formatted = formatTimePointDists(filePath = basepath_clear, Advancing = F)
Clr_Permuted = permuteAvg(Clr_formatted, nPermute = 1000)
Clr_slopes = SlopeVar(RawRes = Clr_formatted, permuteRes = Clr_Permuted, Advancing = F)
Clr_SlopeCI = SlopeCI(Clr_slopes, qLow = 0.05, qHigh = 0.95)

# Make BoxPlot of sampled rate estimates:
results = rateBoxPlots(adv_dat = Adv_slopes, filename = "boxPlot_WaveRateVariance_bothMeasures.pdf",
                       clr_dat = Clr_slopes, bootstrapOnly = F, labs = c("Bootstrapped", "gene by gene"))
results2 = rateBoxPlots(adv_dat = Adv_slopes, filename = "boxPlot_WaveRateVariance_bootstrapped.pdf", 
                        clr_dat = Clr_slopes, bootstrapOnly = T, labs = c("Bootstrapped"))


# print confidence interval values and t-test for clearing vs Advancing wave rates to file
Ttest_GBG.res = t.test(x = Adv_slopes[[1]]$StandardRes, y = Clr_slopes[[1]]$StandardRes)
Ttest_BOOT.res = t.test(x = Adv_slopes[[2]][,5], y = Clr_slopes[[2]][,4])
sink(paste(fig_dir, "RateCI_values.txt", sep = ""))
cat("confidence intervals for advancing wave rates are below\n")
Adv_SlopeCI
cat("\n")
cat("confidence intervals for advancing wave rates are below\n")
Clr_SlopeCI
cat("t-test results for difference in rates for advancing and clearing waves\n Gene By Gene\n")
Ttest_GBG.res
cat("t-test results for difference in rates for advancing and clearing waves\n Bootstrapped\n")
Ttest_BOOT.res
sink()


