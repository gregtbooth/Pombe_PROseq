source("/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/scripts/GB_Functions.R")
library(gplots)
bwpath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/alignment/bw/"
infopath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/pombeinfo/"
fig_dir = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/figures/Dis2Analysis/MetaPlots/01-29-18/"
dir.create(fig_dir)
SPmappability = load.bigWig(paste(infopath, "SPombeMappable_SacCerNonOverLap_expanded_ChrFix.bw", sep = ""))
filteredGL = read.table(paste(infopath, "SP_CompleteFilteredGenes.bed", sep = ""))
SPGL_1kbsep = read.table(file = paste(infopath, "SP_PROcapObservedTSS_1000bpSep.bed", sep = ""), stringsAsFactors = F)
SPGL_1kbsep_filt = merge(filteredGL, SPGL_1kbsep, by.x = 4, by.y = 4)[, c(2,3,4,1,5,6)]
ObsTSS_NOoverlap = read.table(file = paste(infopath, "SP_PROcapObservedTSS_noOverlap.bed", sep = ""))
Introns = read.table(file = paste(infopath, "/features/SPintrons/s_pombe_Introns.bed", sep = "")) 
pausedGenes = read.table(file = paste(infopath, "/SP_PausedGene_Filtered_bothReps.bed", sep = "")) #defined  in genome research paper
NotpausedGenes = read.table(file = paste(infopath, "SP_nonPausedGene_Filtered_bothReps.bed", sep = "")) #defined  in genome research paper
p1nucl = read.table(file = paste(infopath, "/features/SP_NucleosomePositions/SP_firstNuclcoords.bed", sep = ""))   
p2nucl = read.table(file = paste(infopath, "/features/SP_NucleosomePositions/SP_secondNuclcoords.bed", sep = ""))
p3nucl = read.table(file = paste(infopath, "/features/SP_NucleosomePositions/SP_thirdNuclcoords.bed", sep = ""))
p4nucl = read.table(file = paste(infopath, "/features/SP_NucleosomePositions/SP_fourthNuclcoords.bed", sep = ""))
gbnucls = rbind(p2nucl, p3nucl, p4nucl)
# make an index for filtering genes to a certain length. 
indx = (filteredGL[,3]-filteredGL[,2] > 6000)
cat("dimensions of the filtered gene list are: ", dim(filteredGL[indx,]))

#bed = ObsTSS_NOoverlap[indx,]
bed = SPGL_1kbsep_filt
cat("number of genes being used =", length(bed[,1]), "\n")
bed3p = bed[, c(1,3,2,4,5,6)]
Ngenes = length(bed[,1])
#####################################################################################################################
## prepare wig table
## load NormFactors: 
## Note: for the treated mcs6-as sample I am only considering the second biological replicate rather than combined reps, due to artifacts of the 1st sample.
NFs = read.table(file = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/countData/combined_spikeNormFactors.txt", head = T)
repNFs = read.table(file = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/countData/replicate_spikeNormFactors.txt", head = T)
wigset = rbind(c("ePROseq_Pombe_Fisher_WT_COMBINED_pombe_plus.bw", "ePROseq_Pombe_Fisher_WT_COMBINED_pombe_minus.bw", "Fisher_WT_combined", NFs[36,'Spikereads']/100000),
               c("ePROseq_Pombe_Fisher_dis2_11_COMBINED_pombe_plus.bw", "ePROseq_Pombe_Fisher_dis2_11_COMBINED_pombe_minus.bw", "Fisher_dis2_11_combined", NFs[37,'Spikereads']/100000),
               c("ePROseq_Pombe_Fisher_dis2Del_COMBINED_pombe_plus.bw", "ePROseq_Pombe_Fisher_dis2Del_COMBINED_pombe_minus.bw", "Fisher_dis2Del_combined", NFs[38,'Spikereads']/100000),
               c("ePROseq_Pombe_Fisher_dis2_T316A_COMBINED_pombe_plus.bw", "ePROseq_Pombe_Fisher_dis2_T316A_COMBINED_pombe_minus.bw", "Fisher_dis2_T316A_combined", NFs[39,'Spikereads']/100000),
               c("ePROseq_Pombe_Fisher_dis2_T316D_COMBINED_pombe_plus.bw", "ePROseq_Pombe_Fisher_dis2_T316D_COMBINED_pombe_minus.bw", "Fisher_dis2_T316D_combined", NFs[40,'Spikereads']/100000))
          
#####################################################################################################################

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
individualScaledMetas = function(wigset, bed){
  N = dim(wigset)[1]
  #pi.res = vector(mode="list", length=N)
  for (i in 1:N){
    cat("* loading", i, "\n")
    wigs = load.wigset(wigset, i) #this should give all the bigwig files
    cat("* computing meta-plot data ...\n")
    meta= meta.subsample.scaled(bed, bigWig.plus = wigs[[1]], bigWig.minus = wigs[[2]], 1000, 300, 10, do.sum = T)
    metaNorm = meta.normalize(result = meta, scaleFactor = 1/as.numeric(wigs[[4]]))
    cat("* plotting data ...\n")
    pdf(file = paste(fig_dir, sprintf("%s_KbSepGenes_scaled.pdf", wigs[[3]]), sep = ""), width = 8, height = 5)
    meta.plot.scaled.GROseq(metaNorm, metaNorm, 10, ylab = "Normalized Median Counts", main = sprintf("%s", wigs[[3]]), bothStrands = F)
    abline(v=c(-600,-300,300, 600), lwd=2, lty=2)
    dev.off()
    cat("* unloading.\n")
    unload.wigset(wigs)
  }
}

# Plot individual metaplots
individualScaledMetas(wigset = wigset, bed = bed)

#######################################################################################################

lattice_meta_dis2.proSeq = function(filename= "test.pdf", df, ylim = c(0, 10), main = NULL, 
                               xlab = "Distance to TSS (bp)",  ylab = "Median PRO-seq intensity",
                               strainLabs = c("WT", "dis2-11", "dis2_KO", "dis2_T316A","dis2_T316D")){
  pdf(paste(fig_dir, filename), width = 8, height = 6)
  #result <- xyplot(mean ~ x | factor(background), data = df,
  result <- xyplot(mean ~ x, data = df,
                   #group = factor(treatment, labels = c("DMSO", "3MB-PP1")), 
                   group = factor(sample, labels = strainLabs), 
                   scales = list(tck=c(1,0),alternating = c(1,1),x=list(relation='free',axs='i'),
                                 y=list(relation='free',axs='i')),
                   key=list(corner=c(0.98,0.85),padding.text=3,
                            #text=list(c("DMSO", "3MB-PP1"),col=c("black", "black"), cex=0.8, font=2),
                            text=list(strainLabs, col=c("black", "black"), cex=0.8, font=2),
                            rectangles=list(col=c(rgb(0,0.5,0,0.9), rgb(0.1,0,0.7,0.9), rgb(0.2,0,0.6,0.9), rgb(0.3,0,0.5,0.9), rgb(0.4,0,0.4,0.9)), size=2.7, height=0.8, border='white')),
                   type = 'l', 
                   ylim = ylim,
                   col = c(rgb(0,0.5,0,0.9), rgb(0.1,0,0.7,0.9), rgb(0.2,0,0.6,0.9), rgb(0.3,0,0.5,0.9), rgb(0.4,0,0.4,0.9)),  ## just add or take away colors to match number of samples in each panel
                   fill = c(rgb(0,0.5,0,0.1), rgb(0.1,0,0.7,0.1), rgb(0.2,0,0.6,0.1), rgb(0.3,0,0.5,0.1), rgb(0.4,0,0.4,0.1)),
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

lattice_Scaledmeta_dis2.proSeq = function(filename= "test.pdf", df, ylim = c(0, 10), main = NULL, 
                                     xlab = "Distance to gene boundaries (bp)",  ylab = "Median PRO-seq intensity", 
                                     strainLabs = c("WT", "dis2-11", "dis2_KO", "dis2_T316A","dis2_T316D")){
  pdf(paste(fig_dir, filename), width = 8, height = 6)
  #result <- xyplot(mean ~ x | factor(background), data = df,
  result <- xyplot(mean ~ x, data = df, 
                   #group = factor(treatment, labels = c("DMSO", "3MB-PP1")), 
                   group = factor(sample, labels = strainLabs),
                   scales = list(tck=c(1,0),alternating = c(1,1),
                                 x=list(relation='free',axs='i', labels = c('-1000', 'TSS', "+300", "-300", "CPS", "1000"), at = c(1, 101, 131, 191, 221, 320)),
                                 y=list(relation='free',axs='i')),
                   key=list(corner=c(0.98,0.85),padding.text=3,
                            #text=list(c("DMSO", "3MB-PP1"),col=c("black", "black"), cex=0.8, font=2),
                            text=list(strainLabs,col=c("black", "black"), cex=0.8, font=2),
                            rectangles=list(col=c(rgb(0,0.5,0,0.9), rgb(0.1,0,0.7,0.9), rgb(0.2,0,0.6,0.9), rgb(0.3,0,0.5,0.9), rgb(0.4,0,0.4,0.9)), size=2.7, height=0.8, border='white')),
                   type = 'l', 
                   ylim = ylim,
                   col = c(rgb(0,0.5,0,0.9), rgb(0.1,0,0.7,0.9), rgb(0.2,0,0.6,0.9), rgb(0.3,0,0.5,0.9), rgb(0.4,0,0.4,0.9)),  ## just add or take away colors to match number of samples in each panel
                   fill = c(rgb(0,0.5,0,0.1), rgb(0.1,0,0.7,0.1), rgb(0.2,0,0.6,0.1), rgb(0.3,0,0.5,0.1), rgb(0.4,0,0.4,0.1)),
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
FormatScaledMetaData_dis2_lattice <- function(wigset, bed){
  N = dim(wigset)[1]
  df = data.frame(matrix(ncol = 7, nrow = 0), stringsAsFactors = F)
  for (i in 1:N){
    if (i == ceiling(N/2)){cat("* 50% complete ... \n")}
    cat("* loading", i, "\n")
    wigs = load.wigset(wigset, i) #this should give all the bigwig files
    cat("* generating and combining all sample scaled meta-plot data ...\n")
    meta= meta.subsample.scaled(bed, bigWig.plus = wigs[[1]], bigWig.minus = wigs[[2]], 1000, 300, 10, do.sum = T)
    metaNorm = meta.normalize(result = meta, scaleFactor = 1/as.numeric(wigs[[4]]))
    sampleVal = i ### index for later retrieval of sample name from sampleNames list.
    xAxis = seq(from = 1, to = 320)
    sample_df = cbind(xAxis, metaNorm[[4]], metaNorm[[3]], metaNorm[[2]], sampleVal)
    df = rbind.data.frame(df, sample_df)
  }
  colnames(df) = c("x", "mean", "lower", "upper", "sample")
  # lastly I will replace values in the 'sample', 'treatment', and 'background' columns with actual names
  return(df)
}

FormatMetaData_dis2_lattice <- function(wigset, bed, halfWindow, step){
  N = dim(wigset)[1]
  df = data.frame(matrix(ncol = 7, nrow = 0), stringsAsFactors = F)
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
    sample_df = cbind(xAxis, metaNorm[[4]], metaNorm[[3]], metaNorm[[2]], sampleVal)
    df = rbind.data.frame(df, sample_df)
  }
  colnames(df) = c("x", "mean", "lower", "upper", "sample")
  return(df)
}


####################################################################################################################
# Kb sep genes:
#### Scaled gene plots
dis2ScaledMetaData = FormatScaledMetaData_dis2_lattice(wigset, bed = bed)
lattice_Scaledmeta_dis2.proSeq(filename = "All_Dis2_strains.pdf",df = dis2ScaledMetaData, main = NULL, ylim = c(0,20))
# make plots of each time point vs 0' 
lattice_Scaledmeta_dis2.proSeq(filename = "Dis2_WT_vs_dis2_11.pdf",df = dis2ScaledMetaData[dis2ScaledMetaData$sample == 1 | dis2ScaledMetaData$sample == 2,],
                             main = NULL, ylim = c(0,20), strainLabs = c("Dis2", "dis2-11"))
lattice_Scaledmeta_dis2.proSeq(filename = "Dis2_WT_vs_dis2Del.pdf",df = dis2ScaledMetaData[dis2ScaledMetaData$sample == 1 | dis2ScaledMetaData$sample == 3,],
                             main = NULL, ylim = c(0,20), strainLabs = c("Dis2", "dis2Del"))
lattice_Scaledmeta_dis2.proSeq(filename = "Dis2_WT_vs_dis2_T316A.pdf",df = dis2ScaledMetaData[dis2ScaledMetaData$sample == 1 | dis2ScaledMetaData$sample == 4,],
                               main = NULL, ylim = c(0,20), strainLabs = c("Dis2", "dis2_T316A"))
lattice_Scaledmeta_dis2.proSeq(filename = "Dis2_WT_vs_dis2_T316D.pdf",df = dis2ScaledMetaData[dis2ScaledMetaData$sample == 1 | dis2ScaledMetaData$sample == 5,],
                               main = NULL, ylim = c(0,20), strainLabs = c("Dis2", "dis2_T316D"))
####################################################################################################################
# Kb sep genes:
#### CPS centered plots
dis2CPSmeta = FormatMetaData_dis2_lattice(wigset, bed = bed[, c(1,3,2,4,5,6)], halfWindow = 1000, step = 10)
lattice_meta_dis2.proSeq(filename = "All_Dis2_strains_CPS.pdf",df = dis2CPSmeta, main = NULL, ylim = c(0,15))
# make plots of each time point vs 0' 
lattice_meta_dis2.proSeq(filename = "Dis2_WT_vs_dis2_11_CPS_trunc.pdf",df = dis2CPSmeta[dis2CPSmeta$sample == 1 & dis2CPSmeta$x >-260 | dis2CPSmeta$sample == 2 & dis2CPSmeta$x >-260,],
                               main = NULL, ylim = c(0,15), strainLabs = c("Dis2", "dis2-11"))
lattice_meta_dis2.proSeq(filename = "Dis2_WT_vs_dis2Del_CPS_trunc.pdf",df = dis2CPSmeta[dis2CPSmeta$sample == 1 & dis2CPSmeta$x >-260 | dis2CPSmeta$sample == 3 & dis2CPSmeta$x >-260,],
                               main = NULL, ylim = c(0,15), strainLabs = c("Dis2", "dis2Del"))
lattice_meta_dis2.proSeq(filename = "Dis2_WT_vs_dis2_T316A_CPS_trunc.pdf",df = dis2CPSmeta[dis2CPSmeta$sample == 1 & dis2CPSmeta$x >-260 | dis2CPSmeta$sample == 4 & dis2CPSmeta$x >-260,],
                               main = NULL, ylim = c(0,15), strainLabs = c("Dis2", "dis2_T316A"))
lattice_meta_dis2.proSeq(filename = "Dis2_WT_vs_dis2_T316D_CPS_trunc.pdf",df = dis2CPSmeta[dis2CPSmeta$sample == 1 & dis2CPSmeta$x >-260 | dis2CPSmeta$sample == 5 & dis2CPSmeta$x >-260,],
                               main = NULL, ylim = c(0,15), strainLabs = c("Dis2", "dis2_T316D"))




