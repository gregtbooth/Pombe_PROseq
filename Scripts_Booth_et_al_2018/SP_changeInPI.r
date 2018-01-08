### This Script contains the analysis intended for looking for any changes, or lack there of, in promoter proximal pausing after the inhibition of each kinase.
source("/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/scripts/GB_Functions.R")
library(bigWig)
fig_dir = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/figures/PI/10-02-17/"
dir.create(fig_dir)
bwpath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/alignment/bw/"
GEO_bwpath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/GEO_datasets/"
infopath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/pombeinfo/"
countpath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/countData/"
DEseq_Table_path = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/countData/AllGenesObsTSS/DESeqOutput/"
ObsTSS_Filt = read.table(file = paste(infopath, "SP_CompleteFilteredGenes.bed", sep = ""), stringsAsFactors = F)
PausedGL_Filt = read.table(file = paste(infopath, "SP_PausedGene_Filtered_bothReps.bed", sep = ""), stringsAsFactors = F)
NotPausedGL_Filt = read.table(file = paste(infopath, "SP_nonPausedGene_Filtered_bothReps.bed", sep = ""), stringsAsFactors = F)

NFlist = read.table(file = paste(countpath, "combined_spikeNormFactors.txt", sep = ""), head = T)
WT_DMSO = read.table(file = paste(countpath, "AllGenesObsTSS/WT_DMSO_combinedPIcountData.txt", sep = ""), head = T)
WT_MB3PP1 = read.table(file = paste(countpath, "AllGenesObsTSS/WT_3MBPP1_combinedPIcountData.txt", sep = ""), head = T)
MCS6_DMSO = read.table(file = paste(countpath, "AllGenesObsTSS/mcs6as_DMSO_combinedPIcountData.txt", sep = ""), head = T)
MCS6_MB3PP1 = read.table(file = paste(countpath, "AllGenesObsTSS/mcs6as_3MBPP1r2PIcountData.txt", sep = ""), head = T) # only using 2nd replicate because of issues with rep 1
CDK9_DMSO = read.table(file = paste(countpath, "AllGenesObsTSS/Cdk9as_DMSO_combinedPIcountData.txt", sep = ""), head = T)
CDK9_MB3PP1 = read.table(file = paste(countpath, "AllGenesObsTSS/Cdk9as_3MBPP1_combinedPIcountData.txt", sep = ""), head = T)

# function for returning change in PI from untreated to treated
PIchange_quantile.genes = function(bed, WT_df, quantiles = 4){
  bed = merge(x = bed, y = WT_df[,c(18,13)], by.x = 'V4', by.y = 'row.names')
  bed = bed[is.finite(bed[,8]) & bed[,8] > 0,][,c(2,3,4,1,5,6,7,8)]
  quants = quantile(bed[,8], prob = seq(0,1, length = quantiles+1))
  result = list()
  for (i in 1: quantiles){
    quantbed = bed[bed[,8] > quants[i] & bed[,8] < quants[i+1], ]
    result[[i]] = quantbed
  }
  return(result)
}

# function returns a list of gene lists based on their PI decile (based on WT_DMSO Pausing indeces)
PI_quantile.genes = function(bed, WT_df, quantiles = 4){
  bed = merge(x = bed, y = WT_df[,c(18,13,8,10)], by.x = 'V4', by.y = 'row.names')
  bed = bed[is.finite(bed[,8]) & bed[,8] > 0,][,c(2,3,4,1,5,6,7,8,9,10)]
  quants = quantile(bed[,8], prob = seq(0,1, length = quantiles+1))
  result = list()
  for (i in 1: quantiles){
    quantbed = bed[bed[,8] > quants[i] & bed[,8] < quants[i+1], ]
    result[[i]] = quantbed
  }
  return(result)
}

PI_quantile_GL = PI_quantile.genes(ObsTSS_Filt, WT_DMSO, quantiles = 4)

#####################################################################################################################
## prepare wig table
## load NormFactors: 
## Note: for the treated mcs6-as sample I am only considering the second biological replicate rather than combined reps, due to artifacts of the 1st sample.
NFs = read.table(file = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/countData/combined_spikeNormFactors.txt", head = T)
repNFs = read.table(file = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/countData/replicate_spikeNormFactors.txt", head = T)
wigset = rbind(c("5993_7157_26977_HJ72CBGXX_pombe_WT_COMBINED_DMSO_CAGATC_R1.fastq_pombe_plus.bw", "5993_7157_26977_HJ72CBGXX_pombe_WT_COMBINED_DMSO_CAGATC_R1.fastq_pombe_minus.bw", "WT_DMSO_combined", NFs[1,'Spikereads']/100000),
               c("5993_7157_26978_HJ72CBGXX_pombe_WT_COMBINED_3MB-PPI_ACTTGA_R1.fastq_pombe_plus.bw", "5993_7157_26978_HJ72CBGXX_pombe_WT_COMBINED_3MB-PPI_ACTTGA_R1.fastq_pombe_minus.bw", "WT_3MBPP1_combined", NFs[2,'Spikereads']/100000),
               c("5993_7157_26979_HJ72CBGXX_pombe_mcs6as_COMBINED_DMSO_GATCAG_R1.fastq_pombe_plus.bw", "5993_7157_26979_HJ72CBGXX_pombe_mcs6as_COMBINED_DMSO_GATCAG_R1.fastq_pombe_minus.bw", "mcs6as_DMSO_combined", NFs[3,'Spikereads']/100000),
               c("5993_7157_26980_HJ72CBGXX_pombe_mcs6as_rep2_3MB-PPI_TAGCTT_R1.fastq_pombe_plus.bw", "5993_7157_26980_HJ72CBGXX_pombe_mcs6as_rep2_3MB-PPI_TAGCTT_R1.fastq_pombe_minus.bw", "mcs6as_3MBPP1_rep2only", repNFs[8,'Spikereads']/100000),
               c("5993_7157_26981_HJ72CBGXX_pombe_CDK9as_COMBINED_DMSO_GGCTAC_R1.fastq_pombe_plus.bw", "5993_7157_26981_HJ72CBGXX_pombe_CDK9as_COMBINED_DMSO_GGCTAC_R1.fastq_pombe_minus.bw", "Cdk9as_DMSO_combined", NFs[5,'Spikereads']/100000),
               c("5993_7157_26982_HJ72CBGXX_pombe_CDK9as_COMBINED_3MB-PPI_CTTGTA_R1.fastq_pombe_plus.bw", "5993_7157_26982_HJ72CBGXX_pombe_CDK9as_COMBINED_3MB-PPI_CTTGTA_R1.fastq_pombe_minus.bw", "Cdk9as_3MBPP1_combined", NFs[6,'Spikereads']/100000),
               c("6713_7157_31629_HWV7YBGXX_pombe_Isk1as_DMSO_COMBINED_CAGATC_R1.fastq_pombe_plus.bw", "6713_7157_31629_HWV7YBGXX_pombe_Isk1as_DMSO_COMBINED_CAGATC_R1.fastq_pombe_minus.bw", "Isk1as_DMSO_combined", NFs[7,'Spikereads']/100000),
               c("6713_7157_31630_HWV7YBGXX_pombe_Isk1as_3MB-PP1_COMBINED_ACTTGA_R1.fastq_pombe_plus.bw", "6713_7157_31630_HWV7YBGXX_pombe_Isk1as_3MB-PP1_COMBINED_ACTTGA_R1.fastq_pombe_minus.bw", "Isk1as_3MBPP1_combined", NFs[8,'Spikereads']/100000),
               c("6713_7157_31625_HWV7YBGXX_pombe_Mcs6as_CDK9as_DMSO_COMBINED_TTAGGC_R1.fastq_pombe_plus.bw", "6713_7157_31625_HWV7YBGXX_pombe_Mcs6as_CDK9as_DMSO_COMBINED_TTAGGC_R1.fastq_pombe_minus.bw", "Mcs6as_Cdk9as_DMSO_combined", NFs[9,'Spikereads']/100000),
               c("6713_7157_31626_HWV7YBGXX_pombe_Mcs6as_CDK9as_3MB-PP1_COMBINED_TGACCA_R1.fastq_pombe_plus.bw", "6713_7157_31626_HWV7YBGXX_pombe_Mcs6as_CDK9as_3MB-PP1_COMBINED_TGACCA_R1.fastq_pombe_minus.bw", "Mcs6as_Cdk9as_3MBPP1_combined",NFs[10,'Spikereads']/100000))
## Load various ChIP and MNase datasets for assessing correlations with pausing (04-27-16)
GEO_wigset = rbind(c("S.pombeMnase_GSE49575_RAW/GSM1201969_MNase_WT_100Uml_N2_90_200bp_gcnorm19_KDE_bandw20sampled10_chrNumeric.bw", "S.pombeMnase_GSE49575_RAW/GSM1201969_MNase_WT_100Uml_N2_90_200bp_gcnorm19_KDE_bandw20sampled10_chrNumeric.bw", "pombe_MNase", 1), 
                   c("pombe_WinstonPolIIchipSeq/WT_H2B/GSM1201984_ChIP-seq_WT_H2B_II_over_WT_Input_II_FixChr.bw", "pombe_WinstonPolIIchipSeq/WT_H2B/GSM1201984_ChIP-seq_WT_H2B_II_over_WT_Input_II_FixChr.bw", "H2B",1),
                   c("pombe_WinstonPolIIchipSeq/WT_H3/GSM1201990_ChIP-seq_WT_H3_II_over_WT_Input_II_FixChr.bw", "pombe_WinstonPolIIchipSeq/WT_H3/GSM1201990_ChIP-seq_WT_H3_II_over_WT_Input_II_FixChr.bw", "H3", 1), 
                   c("pombe_WinstonPolIIchipSeq/WT_H3K4me3/GSM1201988_ChIP-seq_WT_H3K4me3_II_over_WT_Input_II_FixChr.bw", "pombe_WinstonPolIIchipSeq/WT_H3K4me3/GSM1201988_ChIP-seq_WT_H3K4me3_II_over_WT_Input_II_FixChr.bw", "H3K4me3",1),
                   c("pombe_WinstonPolIIchipSeq/WT_H3K36me3/GSM1201986_ChIP-seq_WT_H3K36me3_II_over_WT_Input_II_FixChr.bw", "pombe_WinstonPolIIchipSeq/WT_H3K36me3/GSM1201986_ChIP-seq_WT_H3K36me3_II_over_WT_Input_II_FixChr.bw", "H3K36me3", 1))
#####################################################################################################################
## meta plots with lattice
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

lattice_meta.proSeq = function(filename= "test.pdf", df, ylim = c(0, 10), main = NULL, 
                               xlab = "Distance to TSS (bp)",  ylab = "Median PRO-seq intensity", width = 5, height = 5){
  pdf(paste(fig_dir, filename), width = width, height = height)
  result <- xyplot(mean ~ x | factor(decile, labels = c("PI quartile 1", "PI quartile 2", "PI quartile 3", "PI quartile 4")) + factor(background), data = df, 
                   group = factor(treatment, labels = c("DMSO", "3MB-PP1")), 
                   scales = list(tck=c(1,0),alternating = c(1,1),x=list(relation='free',axs='i'),
                                 y=list(relation='free',axs='i')),
                   key=list(corner=c(0.98,0.85),padding.text=3,
                            text=list(c("DMSO", "3MB-PP1"),col=c("black", "black"), cex=0.8, font=2),
                            rectangles=list(col=c(rgb(0,0,0.5), rgb(0.5,0,0)), size=2.7, height=0.8, border='white')),
                   type = 'l', 
                   ylim = ylim,
                   col = c(rgb(0.5,0,0,0.9), rgb(0,0,0.5,0.9), rgb(0,0.5,0,0.9), rgb(0,0,0,0.9)),  ## just add or take away colors to match number of samples in each panel
                   fill = c(rgb(0.5,0,0,0.3), rgb(0,0,0.5,0.3), rgb(0,0.5,0,0.3), rgb(0,0,0,0.3)),
                   ylab = ylab,
                   xlab = xlab,
                   main = main,
                   aspect=1,
                   lwd=2.5,
                   par.settings = list(strip.background=list(col="lightgrey"),par.xlab.text=list(cex=1.1,font=2),par.ylab.text=list(cex=1.1,font=2),axis.text=list(cex=1)),
                   par.strip.text=list(cex=1,font=2),
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

lattice_meta.GEOdata = function(filename= "test.pdf", df, ylim = c(0, 10), main = NULL, 
                               xlab = "Distance to TSS (bp)",  ylab = "Median PRO-seq intensity", width = 5, height = 5){
  pdf(paste(fig_dir, filename), width = width, height = height)
  result <- xyplot(mean ~ x | factor(sample), data = df, 
                   group = factor(decile, labels = c("PI quartile 1", "PI quartile 2", "PI quartile 3", "PI quartile 4")), 
                   scales = list(tck=c(1,0),alternating = c(1,1),x=list(relation='free',axs='i'),
                                 y=list(relation='free',axs='i')),
                   key=list(corner=c(0.98,0.85),padding.text=3,
                            text=list(c("PI quartile 1", "PI quartile 2", "PI quartile 3", "PI quartile 4"),col=c("black", "black","black", "black"), cex=0.8, font=2),
                            rectangles=list(col=c(rgb(0,0,0.5), rgb(0.5,0,0), rgb(0,0.5,0), rgb(0,0,0)), size=2.7, height=0.8, border='white')),
                   type = 'l', 
                   ylim = ylim,
                   col = c(rgb(0.5,0,0,0.9), rgb(0,0,0.5,0.9), rgb(0,0.5,0,0.9), rgb(0,0,0,0.9)),  ## just add or take away colors to match number of samples in each panel
                   fill = c(rgb(0.5,0,0,0.3), rgb(0,0,0.5,0.3), rgb(0,0.5,0,0.3), rgb(0,0,0,0.3)),
                   ylab = ylab,
                   xlab = xlab,
                   main = main,
                   aspect=1,
                   lwd=2.5,
                   par.settings = list(strip.background=list(col="lightgrey"),par.xlab.text=list(cex=1.1,font=2),par.ylab.text=list(cex=1.1,font=2),axis.text=list(cex=1)),
                   par.strip.text=list(cex=1,font=2),
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
#####################################################################################################################

## Scaled Meta plot functions 
load.wigset <- function(wigset, wigset_row, bwpath = bwpath) {
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

## formatting meta plots for lattice (genes are also split into quantiles)
FormatMetaData_lattice <- function(wigset, bed_quantiles, halfWindow, step, bwpath = bwpath){
  N = dim(wigset)[1]
  df = data.frame(matrix(ncol = 8, nrow = 0), stringsAsFactors = F)
  sampleNames = c()
  BG_Indx = 0 # counter for tracking genotype, which changes every other sample.
  for (i in 1:N){
    #if (i == ceiling(N/2)){cat("* 50% complete ... \n")}
    cat("* loading", i, "\n")
    wigs = load.wigset(wigset, i, bwpath = bwpath) #this should give all the bigwig files
    cat("* generating and combining all sample meta-plot data ...\n")
    alldec_df = data.frame(matrix(ncol = 4, nrow = 0), stringsAsFactors = F)
    for (ii in 1:length(bed_quantiles)){ ## this loop is needed because we have a list of gene lists based on PI quantiles.
      meta= meta.subsample(bed_quantiles[[ii]], bigWig.plus = wigs[[1]], bigWig.minus = wigs[[2]], 
                           step = step, at.TSS = T, halfWindow = halfWindow, do.sum = T)
      metaNorm = meta.normalize(result = meta, scaleFactor = 1/as.numeric(wigs[[4]]))
      dec_df = cbind(metaNorm[[4]], metaNorm[[3]], metaNorm[[2]], ii) ## make a df with a column noting the quantile
      alldec_df = rbind.data.frame(alldec_df, dec_df)
      }
    sampleNames = c(sampleNames, wigs[[3]])
    sampleVal = i ### index for later retrieval of sample name from sampleNames list.
    if ((-1)**i == -1){BG_Indx = BG_Indx + 1}
    else{BG_Indx = BG_Indx}
    if (grepl("DMSO", wigs[[3]])){ #grepl returns TRUE or false based on the presence or absence of "DMSO" in sample name
      treatment = 0
    }
    else{treatment = 1}
    xAxis = seq(from = -1*(halfWindow), to= halfWindow-step, by = step)
    sample_df = cbind(xAxis, alldec_df, sampleVal, treatment, BG_Indx)
    df = rbind.data.frame(df, sample_df)
  }
  colnames(df) = c("x", "mean", "lower", "upper", 'decile', "sample", 'treatment', 'background')
  # lastly I will replace values in the 'sample', 'treatment', and 'background' columns with actual names
  for (ii in 1:length(sampleNames)){
    df$sample[df$sample == ii] <- sampleNames[ii]
  }
  df$treatment[df$treatment == 0] <- "DMSO"
  df$treatment[df$treatment == 1] <- "3-MB_PP1"
  df$background[df$background == 1] <- "WT"
  df$background[df$background == 2] <- "Mcs6as"
  df$background[df$background == 3] <- "CDK9as"
  df$background[df$background == 4] <- "Isk1as"
  df$background[df$background == 5] <- "Mcs6as_CDK9as"
  return(df)
}
#####################################################################################################################
MetaData = FormatMetaData_lattice(wigset = wigset, bed_quantiles = PI_quantile_GL, halfWindow = 1000, step = 10, bwpath = bwpath)
lattice_meta.proSeq(filename = "TreatedVsUntreatedByPauseIndexQuartiles.pdf", df = MetaData, ylim = c(0,75), width = 10, height = 15)
#look at nucleosomes for genes by pause index quantile. 
MetaData_GEOsamples = FormatMetaData_lattice(wigset = GEO_wigset, bed_quantiles = PI_quantile_GL, halfWindow = 1000, step = 10, bwpath = GEO_bwpath)
lattice_meta.GEOdata(filename = "MetaData_ChIP-MNase_PIquartile.pdf", df = MetaData_GEOsamples, ylim = c(0,6), width = 12, height = 8, ylab = "Median ChIP-seq coverage")

#####################################################################################################################
# Prepare heatmaps of log fold change (treated / untreated) around TSS sorted by WT_DMSO PI. 
# the function below (getHeat.Data)  prepares a list of heatData matrices
getHeat.Data = function(bed, wigset, WT_PI_df, halfwindow = 1000, step = 10){
  N = dim(wigset)[1]
  result = list()
  for (i in 1:N){
    cat("* loading", i, "\n")
    wigs = load.wigset(wigset, i, bwpath = bwpath) #this should give all the bigwig files
    cat("* generating heat data for each sample...\n")
    heatData = collect.many(bed, wigs[[1]], wigs[[2]], halfWindow = halfwindow, step = step, at.TSS=T, do.sum=T)
    heatData = cbind(heatData, WT_PI_df[,13]) ### this column will be used for sorting genes (based on PI, etc.) 
    heatData_sorted = heatData[order(heatData[,dim(heatData)[2]], decreasing = T),]
    heatData_norm = (heatData_sorted[, 1:dim(heatData_sorted)[2]-1])/as.numeric(wigs[[4]])
    heatData_norm[heatData_norm == 0] <- 0.1 ## prevents dividing by zero when taking ratios later on
    result[[i]] = heatData_norm
  }
return(result)
}
AllHeatData = getHeat.Data(ObsTSS_Filt, wigset, WT_DMSO)
WT_treatedFC_heat = log2(AllHeatData[[2]]/AllHeatData[[1]])
MCS6_treatedFC_heat = log2(AllHeatData[[4]]/AllHeatData[[3]])
CDK9_treatedFC_heat = log2(AllHeatData[[6]]/AllHeatData[[5]])
Isk1_treatedFC_heat = log2(AllHeatData[[8]]/AllHeatData[[7]])
Mcs6_CDK9_treatedFC_heat = log2(AllHeatData[[10]]/AllHeatData[[9]])
#####################################################################################################################
# Draw heatmaps 
require(gplots)
my_palette = colorRampPalette(c('lightblue4','lightblue3','lightblue1','white','pink','salmon1','tomato4'))(299)
br = c(seq(-10,-3,length=50), seq(-2.9,-1,length=50), seq(-0.9,0,length=50), seq(0.1,1,length=50), seq(1.1,3,length=50), seq(3.1,10,length=50))

DrawHeat = function(heatData, filename = "heat.pdf", br=br){
  pdf(file = paste(fig_dir, filename, sep = ""), width = 10, height = 10)
  heatmap.2(as.matrix(heatData), Rowv = NA, Colv = NA,  
            dendrogram=c("none"), breaks = br, labRow = F, labCol= F, 
            col =  my_palette, symkey = F, symm=F, trace = c("none"), 
            key.title = "Fold Change", key.xlab = "log2(treated / untreated)", 
            key.ylab = "",
            lmat = rbind(c(2,0),c(3,4),c(0,1)), 
            lwid = c(0.4,4),
            lhei = c(0.4,1,4),
            #key.par = list(cex = 1),
            density.info = "none", useRaster = T)
  dev.off()
}

DrawHeat(WT_treatedFC_heat, filename = "WT_TreatedFC_PIsort_heat.pdf", br=br)
DrawHeat(MCS6_treatedFC_heat, filename = "MCS6as_TreatedFC_PIsort_heat.pdf", br=br)
DrawHeat(CDK9_treatedFC_heat, filename = "CDK9as_TreatedFC_PIsort_heat.pdf", br=br)
DrawHeat(Isk1_treatedFC_heat, filename = "Isk1as_TreatedFC_PIsort_heat.pdf", br=br)
DrawHeat(Mcs6_CDK9_treatedFC_heat, filename = "Mcs6as_CDK9as_TreatedFC_PIsort_heat.pdf", br=br)

###########################################################################################################################
# Looking for correlation between pausing and expression level (use boxplots to assess shifts)
require(ggplot2)
WTPIdeciles = PI_quantile.genes(bed = ObsTSS_Filt, WT_df = WT_DMSO, quantiles = 10)


boxplot.quantiles <- function(QuantileList, filename = "PI_quantile_PRdensity_boxplot.pdf"){
  res = matrix(ncol = 3)
  colnames(res) <- c('PRdensity', "GBdensity", "quantile")
  for (i in 1:length(QuantileList)){
    df = as.matrix(cbind(QuantileList[[i]][,c(9,10)], i))
    res = rbind(res, df)
  }
  result = as.data.frame(res[c(2:dim(res)[1]),])  ## remove first row since it's just NAs
  result$quantile <- as.factor(result$quantile)
  p <- ggplot(result, aes(x = quantile, y = log(PRdensity))) + geom_boxplot()
  ggsave(filename = paste(fig_dir, filename, sep = ""), plot = p, width = 10, height = 10)
}

boxplot.quantiles(WTPIdeciles)

###########################################################################################################################
# added 03-28-17
###########################################################################################################################
require(ggplot2)
# Prepare CDFs of changes in promoter proximal counts for all filtered genes, paused genes and not-paused genes. 
## Load DEseq2 results to get the accurate fold change calculation. 
WT_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_WT_gb_res.txt", sep = ""), head = T)
WT_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_WT_pr_res.txt", sep = ""), head = T)
mcs6_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_mcs6_gb_res.txt", sep = ""), head = T)
mcs6_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_mcs6_pr_res.txt", sep = ""), head = T)
CDK9_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_gb_res.txt", sep = ""), head = T)
CDK9_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_pr_res.txt", sep = ""), head = T)
Isk1_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Isk1_gb_res.txt", sep = ""), head = T)
Isk1_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Isk1_pr_res.txt", sep = ""), head = T)
mcs6_CDK9_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_mcs6_CDK9_gb_res.txt", sep = ""), head = T)
mcs6_CDK9_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_mcs6_CDK9_pr_res.txt", sep = ""), head = T)
resList = list(WT_pr_res, mcs6_pr_res, CDK9_pr_res, Isk1_pr_res, mcs6_CDK9_pr_res)


## function below just returns DESeq caluculated fold change values for a list of DESeq comparisons and filters to include only genes within "genelist"
CompileFC_dat = function(DESeq_res_list, geneList = ObsTSS_Filt){
  FC_dat = matrix(ncol = 2, nrow = 0)
  for (sample in 1:length(DESeq_res_list)){
    FC = merge(DESeq_res_list[sample], geneList, by.x = 0, by.y = 4)[,c(1,3)]
    FC_res_merge = cbind(FC[,2], sample)
    #rownames(FC_res_merge) = FC[,1]
    FC_dat = rbind(FC_dat, FC_res_merge)
  }
  FC_dat = as.data.frame(FC_dat)
  FC_dat[FC_dat[,2]==1, 2] = "WT"
  FC_dat[FC_dat[,2]==2, 2] = "Mcs6as"
  FC_dat[FC_dat[,2]==3, 2] = "Cdk9as"
  FC_dat[FC_dat[,2]==4, 2] = "Lsk1as"
  FC_dat[FC_dat[,2]==5, 2] = "Mcs6as_Cdk9as"
  FC_dat[,2] = as.factor(FC_dat[,2])
  return(FC_dat)
}

## Function for plotting CDFs using ggplot 
plotCDFs = function(filename = "test.pdf",  df, samples = factor(c("WT","Cdk9as")), ylab = "fraction of genes", xlab = "log2 FC promoter"){
  indx = (df[,2] %in% samples)
  res = ggplot(df[indx,], aes(x = V1, colour = sample)) + stat_ecdf() +
    labs(x = xlab, y = ylab)
  ggsave(filename = paste(fig_dir, filename, sep = ""), plot = res, width = 5, height = 5)
}

# get fold change cdf data for specific gene sets 
FiltGenes_cdf_FC_dat = CompileFC_dat(resList, geneList = ObsTSS_Filt)
pausedGenes_cdf_FC_dat = CompileFC_dat(resList, geneList = PausedGL_Filt)
NotpausedGenes_cdf_FC_dat = CompileFC_dat(resList, geneList = NotPausedGL_Filt)

# plot cdfs All samples together
plotCDFs(filename = "Allsamples_AllGenesFilt_prFC_cdf.pdf",  df = FiltGenes_cdf_FC_dat, samples = factor(c("WT","Cdk9as","Lsk1as", "Mcs6as", "Mcs6as_Cdk9as")), ylab = "fraction of genes", xlab = "log2 FC promoter")
plotCDFs(filename = "Allsamples_PausedGenesFilt_prFC_cdf.pdf",  df = pausedGenes_cdf_FC_dat, samples = factor(c("WT","Cdk9as","Lsk1as", "Mcs6as", "Mcs6as_Cdk9as")), ylab = "fraction of genes", xlab = "log2 FC promoter")
plotCDFs(filename = "Allsamples_NotPausedGenesFilt_prFC_cdf.pdf",  df = NotpausedGenes_cdf_FC_dat, samples = factor(c("WT","Cdk9as","Lsk1as", "Mcs6as", "Mcs6as_Cdk9as")), ylab = "fraction of genes", xlab = "log2 FC promoter")

# plot cdfs All samples together
plotCDFs(filename = "WTvsCdk9as_AllGenesFilt_prFC_cdf.pdf",  df = FiltGenes_cdf_FC_dat, samples = factor(c("WT","Cdk9as")), ylab = "fraction of genes", xlab = "log2 FC promoter")
plotCDFs(filename = "WTvsCdk9as_PausedGenesFilt_prFC_cdf.pdf",  df = pausedGenes_cdf_FC_dat, samples = factor(c("WT","Cdk9as")), ylab = "fraction of genes", xlab = "log2 FC promoter")
plotCDFs(filename = "WTvsCdk9as_NotPausedGenesFilt_prFC_cdf.pdf",  df = NotpausedGenes_cdf_FC_dat, samples = factor(c("WT","Cdk9as")), ylab = "fraction of genes", xlab = "log2 FC promoter")
###########################################################################################################################
## added 06-19-17
# Get CDFs for Pause Index for each strain before and after treatment:
plotCDFs_PI = function(DMSOdat, PP1dat, geneList, filename = "test.pdf", ylab = "fraction of genes", xlab = "log10(Pause Index)"){
  Ddat = merge(DMSOdat, geneList, by.x = 0, by.y = 4) ## filter out unused genes
  Pdat = merge(PP1dat, geneList, by.x = 0, by.y = 4)
  DMSO_PI = cbind(log10(Ddat[,14]), 1)
  colnames(DMSO_PI) = c("PI", "treatment") 
  PP1_PI = cbind(log10(Pdat[,14]), 2)
  colnames(PP1_PI) = c("PI", "treatment") 
  dat = data.frame(rbind(DMSO_PI, PP1_PI))
  dat = as.data.frame(dat)
  dat[dat[,2]==1, 2] = "DMSO"
  dat[dat[,2]==2, 2] = "3MBPP1"
  dat[,2] = as.factor(dat[,2])
  #indx = dat[is.finite(dat[,1]),]
  res = ggplot(dat, aes(x = PI, colour = treatment)) + stat_ecdf() +
    labs(x = xlab, y = ylab) + xlim(-2,2)
  ggsave(filename = paste(fig_dir, filename, sep = ""), plot = res, width = 5, height = 5)
  return(res)
}

## run script
WTPIcdf = plotCDFs_PI(DMSOdat = WT_DMSO, PP1dat = WT_MB3PP1, geneList = ObsTSS_Filt, filename = "WT_PI_DMSOvs3MBPP1_cdf_1.pdf", ylab = "fraction of genes", xlab = "log10(Pause Index)")
CDK9asPIcdf = plotCDFs_PI(DMSOdat = CDK9_DMSO, PP1dat = CDK9_MB3PP1, geneList = ObsTSS_Filt, filename = "Cdk9as_PI_DMSOvs3MBPP1_cdf_1.pdf")
Mcs6asPIcdf = plotCDFs_PI(DMSOdat = MCS6_DMSO, PP1dat = MCS6_MB3PP1, geneList = ObsTSS_Filt, filename = "MCS6as_PI_DMSOvs3MBPP1_cdf_1.pdf")


###########################################################################################################################
# Looking for correlation between pausing and change in GB (use boxplots to assess shifts)
WTPIquartiles = PI_quantile.genes(bed = ObsTSS_Filt, WT_df = WT_DMSO, quantiles = 4)
WTPIdeciles = PI_quantile.genes(bed = ObsTSS_Filt, WT_df = WT_DMSO, quantiles = 10)

boxplot.DESeq_vs_PIquantiles <- function(QuantileList, DESeqRes, filename = "PI_quantile_GBfoldchange_boxplot.pdf"){
  res = matrix(ncol = 4)
  colnames(res) <- c('PRdensity', "GBdensity","gbFoldChange", "quantile")
  for (i in 1:length(QuantileList)){
    combineDat = merge(QuantileList[i], DESeqRes, by.x = "V4", by.y = 0)
    df = as.matrix(cbind(combineDat[,c(9,10,12)], i)) # prDensity, gbDensity, gbFoldChange (fromDEseq), quantile
    res = rbind(res, df)
  }
  result = as.data.frame(res[c(2:dim(res)[1]),])  ## remove first row since it's just NAs
  result$quantile <- as.factor(result$quantile)
  p <- ggplot(result, aes(x = quantile, y = gbFoldChange)) + geom_boxplot()
  ggsave(filename = paste(fig_dir, filename, sep = ""), plot = p, width = 5, height = 5)
  return(res)
}
## compare relationship between PI and gb fold change for each mutant. (deciles and quantiles)
Cdk9res_4 = boxplot.DESeq_vs_PIquantiles(QuantileList = WTPIquartiles, DESeqRes = CDK9_gb_res, filename = "PI_quantile_Cdk9asGBfoldchange_boxplot.pdf")
Cdk9res_10 = boxplot.DESeq_vs_PIquantiles(QuantileList = WTPIdeciles, DESeqRes = CDK9_gb_res, filename = "PI_decile_Cdk9asGBfoldchange_boxplot.pdf")
Lsk1res_4 = boxplot.DESeq_vs_PIquantiles(QuantileList = WTPIquartiles, DESeqRes = Isk1_gb_res, filename = "PI_quantile_Lsk1asGBfoldchange_boxplot.pdf")
Lsk1res_10 = boxplot.DESeq_vs_PIquantiles(QuantileList = WTPIdeciles, DESeqRes = Isk1_gb_res, filename = "PI_decile_Lsk1asGBfoldchange_boxplot.pdf")
Mcs6res_4 = boxplot.DESeq_vs_PIquantiles(QuantileList = WTPIquartiles, DESeqRes = mcs6_gb_res, filename = "PI_quantile_Mcs6asGBfoldchange_boxplot.pdf")
Mcs6res_10 = boxplot.DESeq_vs_PIquantiles(QuantileList = WTPIdeciles, DESeqRes = mcs6_gb_res, filename = "PI_decile_Mcs6asGBfoldchange_boxplot.pdf")
Mcs6Cdk9res_4 = boxplot.DESeq_vs_PIquantiles(QuantileList = WTPIquartiles, DESeqRes = mcs6_CDK9_gb_res, filename = "PI_quantile_Mcs6asCdk9asGBfoldchange_boxplot.pdf")
Mcs6Cdk9res_10 = boxplot.DESeq_vs_PIquantiles(QuantileList = WTPIdeciles, DESeqRes = mcs6_CDK9_gb_res, filename = "PI_decile_Mcs6asCdk9asGBfoldchange_boxplot.pdf")
## Get actual Correlations and generate a scatter plot for PI (WT DMSO) and fold change after treatment (each strain)
## make single data frame with PI in first col and each strain's fold change in subsequent columns (will make individual plots relative to PI column)
# Note: join_all should allow to merge multiple dataframes simultaneously
library(plyr)
PI = data.frame(WT_DMSO[,12])
PI$rn = row.names(WT_DMSO)
WT = data.frame(WT_gb_res[,2])
WT$rn = row.names(WT_gb_res)
lsk1 = data.frame(Isk1_gb_res[,2])
lsk1$rn = row.names(Isk1_gb_res)
mcs6 = data.frame(mcs6_gb_res[,2])
mcs6$rn = row.names(mcs6_gb_res)
cdk9 = data.frame(CDK9_gb_res[,2])
cdk9$rn = row.names(CDK9_gb_res)
mcs6_cdk9 = data.frame(mcs6_CDK9_gb_res[,2])
mcs6_cdk9$rn = row.names(mcs6_CDK9_gb_res)
FCvsPI_df = join_all(list(PI, WT, lsk1, mcs6, cdk9, mcs6_cdk9), by = "rn", type = "full")
FCvsPI_df_filt = merge(FCvsPI_df, ObsTSS_Filt, by.x = "rn", by.y = "V4")[,c(1:7)]

## plotting function for comparison between strain FC (treated/untreated) and WT PI
x = log2(FCvsPI_df_filt[,2])
indx = is.finite(x)
library(hexbin)
ScatterPlot.cor = function(x = FCvsPI_df_filt[indx, 2], y = FCvsPI_df_filt[indx, 3], 
                           filename = "test_scatter.pdf", main = NULL, xlab = "log2(WT PI)", ylab = "log2(gb FC)"){
  Ngenes = length(x)
  pdf(file = paste(fig_dir, filename, sep = ""), width = 7, height = 7)
  result = hexbinplot(y ~ log2(x),
             aspect = 1, xbins = 35, 
             xlim = c(-10,10),
             ylim = c(-2.5,2.5),
             main= list(label = main, fontsize = 25),  #fontsize mustbe specified here to be displayed (cex won't work)
             xlab = list(label = xlab, fontsize = 20), ## lists must be used to set the various parameters for each label
             ylab = list(label = ylab, fontsize =20), type = c('r', 'g'), lwd = 3,
             scales = list(cex = 1.5)) + latticeExtra::layer(panel.key(c(sprintf("Rho = %f", cor.test(x, y, method = 'spearman')[[4]]), 
                                                                         sprintf("cor_Pval = %f", cor.test(x, y, method = 'spearman')[[3]]), sprintf("N = %i", length(x))), corner = c(0,.98), lines = FALSE, points = FALSE, cex = 1.5), packets = 1)
  print(result)
  dev.off()
}
# make plot for each comparsion
ScatterPlot.cor(x = FCvsPI_df_filt[indx, 2], y = FCvsPI_df_filt[indx, 3], filename = "WTgbFC_vs_WTPI_scatter.pdf", main = NULL, xlab = "log2(WT PI)", ylab = "log2(gb FC)")
ScatterPlot.cor(x = FCvsPI_df_filt[indx, 2], y = FCvsPI_df_filt[indx, 4], filename = "Lsk1gbFC_vs_WTPI_scatter.pdf", main = NULL, xlab = "log2(WT PI)", ylab = "log2(gb FC)")
ScatterPlot.cor(x = FCvsPI_df_filt[indx, 2], y = FCvsPI_df_filt[indx, 5], filename = "Mcs6gbFC_vs_WTPI_scatter.pdf", main = NULL, xlab = "log2(WT PI)", ylab = "log2(gb FC)")
ScatterPlot.cor(x = FCvsPI_df_filt[indx, 2], y = FCvsPI_df_filt[indx, 6], filename = "Cdk9gbFC_vs_WTPI_scatter.pdf", main = NULL, xlab = "log2(WT PI)", ylab = "log2(gb FC)")
ScatterPlot.cor(x = FCvsPI_df_filt[indx, 2], y = FCvsPI_df_filt[indx, 7], filename = "Mcs6Cdk9gbFC_vs_WTPI_scatter.pdf", main = NULL, xlab = "log2(WT PI)", ylab = "log2(gb FC)")

cor.test(log2(FCvsPI_df_filt[indx, 2]), FCvsPI_df_filt[indx, 5], method = 'spearman')









