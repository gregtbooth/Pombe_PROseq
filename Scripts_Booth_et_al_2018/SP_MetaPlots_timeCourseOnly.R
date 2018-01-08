source("/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/scripts/GB_Functions.R")
library(gplots)
bwpath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/alignment/bw/"
infopath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/pombeinfo/"
fig_dir = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/figures/MetaPlotComparison/10-02-17/TimeCourse/"
dir.create(fig_dir)
SPmappability = load.bigWig(paste(infopath, "SPombeMappable_SacCerNonOverLap_expanded_ChrFix.bw", sep = ""))
filteredGL = read.table(paste(infopath, "SP_CompleteFilteredGenes.bed", sep = ""))
SPGL_1kbsep = read.table(file = paste(infopath, "SP_PROcapObservedTSS_1000bpSep.bed", sep = ""), stringsAsFactors = F)
SPGL_1kbsep_filt = merge(filteredGL, SPGL_1kbsep, by.x = 4, by.y = 4)[, c(2,3,4,1,5,6)]
ObsTSS_NOoverlap = read.table(file = paste(infopath, "SP_PROcapObservedTSS_noOverlap.bed", sep = ""))
# make an index for filtering genes to a certain length. 
indx = (filteredGL[,3]-filteredGL[,2] > 6000)
cat("dimensions of the filtered gene list are: ", dim(filteredGL[indx,]))

#####################################################################################################################
## prepare wig table
## load NormFactors: 
## Note: for the treated mcs6-as sample I am only considering the second biological replicate rather than combined reps, due to artifacts of the 1st sample.
NFs = read.table(file = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/countData/combined_spikeNormFactors.txt", head = T)

TimeCourse_wigset = rbind(c("ePROseq_Pombe_CDK9as_0min_DMSO_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_0min_DMSO_COMBINED_pombe_minus.bw", "Cdk9as_0min_combined", NFs[26,'Spikereads']/100000),
                          c("ePROseq_Pombe_CDK9as_30sec_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_30sec_3MBPP1_10uM_COMBINED_pombe_minus.bw", "Cdk9as_30sec_combined", NFs[27,'Spikereads']/100000),
                          c("ePROseq_Pombe_CDK9as_1min_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_1min_3MBPP1_10uM_COMBINED_pombe_minus.bw", "Cdk9as_1min_combined", NFs[28,'Spikereads']/100000),
                          c("ePROseq_Pombe_CDK9as_2min30sec_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_2min30sec_3MBPP1_10uM_COMBINED_pombe_minus.bw", "Cdk9as_2min30sec_combined", NFs[29,'Spikereads']/100000),
                          c("ePROseq_Pombe_CDK9as_5min_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_5min_3MBPP1_10uM_COMBINED_pombe_minus.bw", "Cdk9as_5min_combined", NFs[30,'Spikereads']/100000),
                          c("ePROseq_Pombe_CDK9as_7min30sec_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_7min30sec_3MBPP1_10uM_COMBINED_pombe_minus.bw", "Cdk9as_7min30sec_combined", NFs[31,'Spikereads']/100000),
                          c("ePROseq_Pombe_CDK9as_10min_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_10min_3MBPP1_10uM_COMBINED_pombe_minus.bw", "Cdk9as_10min_combined", NFs[32,'Spikereads']/100000),
                          c("ePROseq_Pombe_CDK9as_20min_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_20min_3MBPP1_10uM_COMBINED_pombe_minus.bw", "Cdk9as_20min_combined", NFs[33,'Spikereads']/100000))
#####################################################################################################################
load.wigset <- function(wigset, wigset_row) {
  file = wigset[wigset_row, 1]
  wig.p = NULL
  if (file != "")
    wig.p = load.bigWig(paste(bwpath, file, sep=''))
  file = wigset[wigset_row, 2]
  wig.m = NULL
  if (file != "")
    wig.m = load.bigWig(paste(bwpath, file, sep=''))
  return(list(wig.p, wig.m, wigset[wigset_row, 3], wigset[wigset_row, 4]))
}

require(lattice)
my.panel.bands <-
  function(x, y, upper, lower,
           fill, col,
           subscripts, ..., font, fontface)
  {
    upper <- upper[subscripts]
    lower <- lower[subscripts]
    panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
                  col = fill, border = FALSE,
                  ...)
  }

lattice_meta_TC.proSeq = function(filename= "test.pdf", df, ylim = c(0, 10), main = NULL, 
                               xlab = "Distance to TSS (bp)",  ylab = "Median PRO-seq intensity",
                               minuteLabs = c("0 min", "30 sec", "1 min", "2.5 min","5 min", "7.5 min", "10 min", "20 min")){
  pdf(paste(fig_dir, filename), width = 8, height = 6)
  #result <- xyplot(mean ~ x | factor(background), data = df,
  result <- xyplot(mean ~ x | factor(background), data = df,
                   #group = factor(treatment, labels = c("DMSO", "3MB-PP1")), 
                   group = factor(sample, labels = minuteLabs), 
                   scales = list(tck=c(1,0),alternating = c(1,1),x=list(relation='free',axs='i'),
                                 y=list(relation='free',axs='i')),
                   key=list(corner=c(0.98,0.85),padding.text=3,
                            #text=list(c("DMSO", "3MB-PP1"),col=c("black", "black"), cex=0.8, font=2),
                            text=list(minuteLabs, col=c("black", "black"), cex=0.8, font=2),
                            rectangles=list(col=c(rgb(0,0.5,0,0.9), rgb(0.1,0,0.7,0.9), rgb(0.2,0,0.6,0.9), rgb(0.3,0,0.5,0.9), rgb(0.4,0,0.4,0.9), rgb(0.5,0,0.2,0.9), rgb(0.6,0,0.1,0.9), rgb(0.7,0,0,0.9)), size=2.7, height=0.8, border='white')),
                   type = 'l', 
                   ylim = ylim,
                   col = c(rgb(0,0.5,0,0.9), rgb(0.1,0,0.7,0.9), rgb(0.2,0,0.6,0.9), rgb(0.3,0,0.5,0.9), rgb(0.4,0,0.4,0.9), rgb(0.5,0,0.2,0.9), rgb(0.6,0,0.1,0.9), rgb(0.7,0,0,0.9)),  ## just add or take away colors to match number of samples in each panel
                   fill = c(rgb(0,0.5,0,0.1), rgb(0.1,0,0.7,0.1), rgb(0.2,0,0.6,0.1), rgb(0.3,0,0.5,0.1), rgb(0.4,0,0.4,0.1), rgb(0.5,0,0.2,0.1), rgb(0.6,0,0.1,0.1), rgb(0.7,0,0,0.1)),
                   ylab = ylab,
                   xlab = xlab,
                   main = main,
                   aspect=1,
                   lwd=3,
                   par.settings = list(strip.background=list(col="lightgrey"),par.xlab.text=list(cex=1.1,font=2),par.ylab.text=list(cex=1.1,font=2),axis.text=list(cex=1)),
                   par.strip.text=list(cex=1.2,font=2),
                   upper = df$upper,
                   lower = df$lower,
                   panel = function(x, y, ...){
                     panel.grid(h=-1, v=-1)
                     panel.superpose(x, y, panel.groups = 'my.panel.bands', ...)
                     panel.xyplot(x, y, ...)
                   })
  print(result)
  dev.off()
}

lattice_Scaledmeta_TC.proSeq = function(filename= "test.pdf", df, ylim = c(0, 10), main = NULL, 
                                     xlab = "Distance to gene boundaries (bp)",  ylab = "Median PRO-seq intensity", 
                                     minuteLabs = c("0 min", "30 sec", "1 min", "2.5 min","5 min", "7.5 min", "10 min", "20 min")){
  pdf(paste(fig_dir, filename), width = 8, height = 6)
  #result <- xyplot(mean ~ x | factor(background), data = df,
  result <- xyplot(mean ~ x | factor(background), data = df, 
                   #group = factor(treatment, labels = c("DMSO", "3MB-PP1")), 
                   group = factor(sample, labels = minuteLabs),
                   scales = list(tck=c(1,0),alternating = c(1,1),
                                 x=list(relation='free',axs='i', labels = c('-1000', 'TSS', "+300", "-300", "CPS", "1000"), at = c(1, 101, 131, 191, 221, 320)),
                                 y=list(relation='free',axs='i')),
                   key=list(corner=c(0.98,0.85),padding.text=3,
                            #text=list(c("DMSO", "3MB-PP1"),col=c("black", "black"), cex=0.8, font=2),
                            text=list(minuteLabs,col=c("black", "black"), cex=0.8, font=2),
                            rectangles=list(col=c(rgb(0,0.5,0,0.9), rgb(0.1,0,0.7,0.9), rgb(0.2,0,0.6,0.9), rgb(0.3,0,0.5,0.9), rgb(0.4,0,0.4,0.9), rgb(0.5,0,0.2,0.9), rgb(0.6,0,0.1,0.9), rgb(0.7,0,0,0.9)), size=2.7, height=0.8, border='white')),
                   type = 'l', 
                   ylim = ylim,
                   col = c(rgb(0,0.5,0,0.9), rgb(0.1,0,0.7,0.9), rgb(0.2,0,0.6,0.9), rgb(0.3,0,0.5,0.9), rgb(0.4,0,0.4,0.9), rgb(0.5,0,0.2,0.9), rgb(0.6,0,0.1,0.9), rgb(0.7,0,0,0.9)),  ## just add or take away colors to match number of samples in each panel
                   fill = c(rgb(0,0.5,0,0.1), rgb(0.1,0,0.7,0.1), rgb(0.2,0,0.6,0.1), rgb(0.3,0,0.5,0.1), rgb(0.4,0,0.4,0.1), rgb(0.5,0,0.2,0.1), rgb(0.6,0,0.1,0.1), rgb(0.7,0,0,0.1)),
                   ylab = ylab,
                   xlab = xlab,
                   main = main,
                   aspect=0.5,
                   lwd=3,
                   par.settings = list(strip.background=list(col="lightgrey"),par.xlab.text=list(cex=1.1,font=2),par.ylab.text=list(cex=1.1,font=2),axis.text=list(cex=1)),
                   par.strip.text=list(cex=1.2,font=2),
                   upper = df$upper,
                   lower = df$lower,
                   panel = function(x, y, ...){
                     panel.grid(h=-1, v=-1)
                     panel.superpose(x, y, panel.groups = 'my.panel.bands', ...)
                     panel.xyplot(x, y, ...)
                   }) + latticeExtra::layer(panel.abline(v = c(101,131,191,221), lty = 2, lwd = 2, col = 'gray'))
  print(result)
  dev.off()
}
#####################################################################################################################
## Scaled Meta plot functions 

## formatting Scaled meta plots for lattice
FormatScaledMetaData_TC_lattice <- function(wigset, bed){
  N = dim(wigset)[1]
  df = data.frame(matrix(ncol = 7, nrow = 0), stringsAsFactors = F)
  sampleNames = c()
  BG_Indx = 1
  for (i in 1:N){
    if (i == ceiling(N/2)){cat("* 50% complete ... \n")}
    cat("* loading", i, "\n")
    wigs = load.wigset(wigset, i) #this should give all the bigwig files
    cat("* generating and combining all sample scaled meta-plot data ...\n")
    meta= meta.subsample.scaled(bed, bigWig.plus = wigs[[1]], bigWig.minus = wigs[[2]], 1000, 300, 10, do.sum = T)
    metaNorm = meta.normalize(result = meta, scaleFactor = 1/as.numeric(wigs[[4]]))
    sampleNames = c(sampleNames, wigs[[3]])
    sampleVal = i ### index for later retrieval of sample name from sampleNames list.
    if (grepl("_0min_", wigs[[3]])){ #grepl returns TRUE or false based on the presence or absence of "DMSO" in sample name
      treatment = 0
    }
    else{treatment = 1}
    xAxis = seq(from = 1, to = 320)
    sample_df = cbind(xAxis, metaNorm[[4]], metaNorm[[3]], metaNorm[[2]], sampleVal, treatment, BG_Indx)
    df = rbind.data.frame(df, sample_df)
  }
  colnames(df) = c("x", "mean", "lower", "upper", "sample", 'treatment', 'background')
  # lastly I will replace values in the 'sample', 'treatment', and 'background' columns with actual names
  df$treatment[df$treatment == 0] <- "DMSO"
  df$treatment[df$treatment == 1] <- "3-MB_PP1"
  df$background[df$background == 1] <- "CDK9_as"
  return(df)
}

FormatMetaData_TC_lattice <- function(wigset, bed, halfWindow, step){
  N = dim(wigset)[1]
  df = data.frame(matrix(ncol = 7, nrow = 0), stringsAsFactors = F)
  sampleNames = c()
  BG_Indx = 1
  for (i in 1:N){
    if (i == ceiling(N/2)){cat("* 50% complete ... \n")}
    cat("* loading", i, "\n")
    wigs = load.wigset(wigset, i) #this should give all the bigwig files
    cat("* generating and combining all sample scaled meta-plot data ...\n")
    meta= meta.subsample(bed, bigWig.plus = wigs[[1]], bigWig.minus = wigs[[2]], 
                         step = step, at.TSS = T, halfWindow = halfWindow, do.sum = T)
    metaNorm = meta.normalize(result = meta, scaleFactor = 1/as.numeric(wigs[[4]]))
    sampleVal = i ### index for later retrieval of sample name from sampleNames list.
    xAxis = seq(from = -halfWindow, to = halfWindow, by = step)
    if (grepl("_0min_", wigs[[3]])){ #grepl returns TRUE or false based on the presence or absence of "DMSO" in sample name
      treatment = 0
    }
    else{treatment = 1}
    xAxis = seq(from = -halfWindow, to = halfWindow, by = step)
    sample_df = cbind(xAxis, metaNorm[[4]], metaNorm[[3]], metaNorm[[2]], sampleVal, treatment, BG_Indx)
    df = rbind.data.frame(df, sample_df)
  }
  colnames(df) = c("x", "mean", "lower", "upper", "sample", 'treatment', 'background')
  # lastly I will replace values in the 'sample', 'treatment', and 'background' columns with actual names
  df$treatment[df$treatment == 0] <- "DMSO"
  df$treatment[df$treatment == 1] <- "3-MB_PP1"
  df$background[df$background == 1] <- "CDK9_as"
  return(df)
}


####################################################################################################################
# Kb sep genes:
#### Scaled gene plots
TimeCourseScaledMetaData = FormatScaledMetaData_TC_lattice(TimeCourse_wigset, bed = bed)
lattice_Scaledmeta_TC.proSeq(filename = "ePROseqAllTimePoints.pdf",df = TimeCourseScaledMetaData, main = NULL, ylim = c(0,20))
# make plots of each time point vs 0' 
lattice_Scaledmeta_TC.proSeq(filename = "ePROseq0vs30secTimePoints.pdf",df = TimeCourseScaledMetaData[TimeCourseScaledMetaData$sample == 1 | TimeCourseScaledMetaData$sample == 2,],
                             main = NULL, ylim = c(0,20), minuteLabs = c("0 min", "30 sec"))
lattice_Scaledmeta_TC.proSeq(filename = "ePROseq0vs1minTimePoints.pdf",df = TimeCourseScaledMetaData[TimeCourseScaledMetaData$sample == 1 | TimeCourseScaledMetaData$sample == 3,],
                             main = NULL, ylim = c(0,20), minuteLabs = c("0 min", "1 min"))
lattice_Scaledmeta_TC.proSeq(filename = "ePROseq0vs2min30secTimePoints.pdf",df = TimeCourseScaledMetaData[TimeCourseScaledMetaData$sample == 1 | TimeCourseScaledMetaData$sample == 4,],
                             main = NULL, ylim = c(0,20), minuteLabs = c("0 min", "2.5 min"))
lattice_Scaledmeta_TC.proSeq(filename = "ePROseq0vs5minTimePoints.pdf",df = TimeCourseScaledMetaData[TimeCourseScaledMetaData$sample == 1 | TimeCourseScaledMetaData$sample == 5,],
                             main = NULL, ylim = c(0,20), minuteLabs = c("0 min", "5 min"))
lattice_Scaledmeta_TC.proSeq(filename = "ePROseq0vs7min30secTimePoints.pdf",df = TimeCourseScaledMetaData[TimeCourseScaledMetaData$sample == 1 | TimeCourseScaledMetaData$sample == 6,],
                             main = NULL, ylim = c(0,20), minuteLabs = c("0 min", "7.5 min"))
lattice_Scaledmeta_TC.proSeq(filename = "ePROseq0vs10minTimePoints.pdf",df = TimeCourseScaledMetaData[TimeCourseScaledMetaData$sample == 1 | TimeCourseScaledMetaData$sample == 7,],
                             main = NULL, ylim = c(0,20), minuteLabs = c("0 min", "10 min"))
lattice_Scaledmeta_TC.proSeq(filename = "ePROseq0vs20minTimePoints.pdf",df = TimeCourseScaledMetaData[TimeCourseScaledMetaData$sample == 1 | TimeCourseScaledMetaData$sample == 8,],
                             main = NULL, ylim = c(0,20), minuteLabs = c("0 min", "20 min"))
#### unscaled 3KB to TSS
#indx = (bed[,3]-bed[,2] > 3000)
cat("dimensions of the filtered gene list are: ", dim(filteredGL[indx,]))
TimeCourseMetaData = FormatMetaData_TC_lattice(TimeCourse_wigset, bed = filteredGL[indx,], halfWindow = 6000, step = 300)
xAxisIDX = TimeCourseMetaData$x >= -900
TimeCourseMetaData_xLim = TimeCourseMetaData[xAxisIDX,]
lattice_meta_TC.proSeq(filename = "ePROseqAllTimePoints_longGenes6KB_300bpsmooth.pdf",df = TimeCourseMetaData, main = NULL, ylim = c(0,600))

# make plots of each time point vs 0' 
lattice_meta_TC.proSeq(filename = "ePROseq0vs30secTimePoints_longGenes6KB_300bpsmooth.pdf",df = TimeCourseMetaData_xLim[TimeCourseMetaData_xLim$sample == 1 | TimeCourseMetaData_xLim$sample == 2,],
                       main = NULL, ylim = c(0,600), minuteLabs = c("0 min", "30 sec"))
lattice_meta_TC.proSeq(filename = "ePROseq0vs1minTimePoints_longGenes6KB_300bpsmooth.pdf",df = TimeCourseMetaData_xLim[TimeCourseMetaData_xLim$sample == 1 | TimeCourseMetaData_xLim$sample == 3,],
                       main = NULL, ylim = c(0,600), minuteLabs = c("0 min", "1 min"))
lattice_meta_TC.proSeq(filename = "ePROseq0vs2min30secTimePoints_longGenes6KB_300bpsmooth.pdf",df = TimeCourseMetaData_xLim[TimeCourseMetaData_xLim$sample == 1 | TimeCourseMetaData_xLim$sample == 4,],
                       main = NULL, ylim = c(0,600), minuteLabs = c("0 min", "2.5 min"))
lattice_meta_TC.proSeq(filename = "ePROseq0vs5minTimePoints_longGenes6KB_300bpsmooth.pdf",df = TimeCourseMetaData_xLim[TimeCourseMetaData_xLim$sample == 1 | TimeCourseMetaData_xLim$sample == 5,],
                       main = NULL, ylim = c(0,600), minuteLabs = c("0 min", "5 min"))
lattice_meta_TC.proSeq(filename = "ePROseq0vs7min30secTimePoints_longGenes6KB_300bpsmooth.pdf",df = TimeCourseMetaData_xLim[TimeCourseMetaData_xLim$sample == 1 | TimeCourseMetaData_xLim$sample == 6,],
                       main = NULL, ylim = c(0,600), minuteLabs = c("0 min", "7.5 min"))
lattice_meta_TC.proSeq(filename = "ePROseq0vs10minTimePoints_longGenes6KB_300bpsmooth.pdf",df = TimeCourseMetaData_xLim[TimeCourseMetaData_xLim$sample == 1 | TimeCourseMetaData_xLim$sample == 7,],
                       main = NULL, ylim = c(0,600), minuteLabs = c("0 min", "10 min"))
lattice_meta_TC.proSeq(filename = "ePROseq0vs20minTimePoints_longGenes6KB_300bpsmooth.pdf",df = TimeCourseMetaData_xLim[TimeCourseMetaData_xLim$sample == 1 | TimeCourseMetaData_xLim$sample == 8,],
                       main = NULL, ylim = c(0,600), minuteLabs = c("0 min", "20 min"))





