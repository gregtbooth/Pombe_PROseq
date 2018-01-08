source("/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/scripts/GB_Functions.R")
bwpath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/alignment/bw/"
infopath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/pombeinfo/"
fig_dir = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/figures/PostCPS_Distributions/10-02-17/"
dir.create(fig_dir)
SPmappability = load.bigWig(paste(infopath, "SPombeMappable_SacCerNonOverLap_expanded_ChrFix.bw", sep = ""))
filteredGL = read.table(paste(infopath, "SP_CompleteFilteredGenes.bed", sep = ""))
SPGL_1kbsep = read.table(file = paste(infopath, "SP_PROcapObservedTSS_1000bpSep.bed", sep = ""), stringsAsFactors = F)
#
usedbed = merge(SPGL_1kbsep, filteredGL, by.x  ="V4", by.y = "V4")[,c(2,3,4,1,5,6)]
#####################################################################################################################
## prepare wig table
## load NormFactors: 
## Note: for the treated mcs6-as sample I am only considering the second biological replicate rather than combined reps, due to artifacts of the 1st sample.
NFs = read.table(file = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/countData/combined_spikeNormFactors.txt", head = T)
repNFs = read.table(file = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/countData/replicate_spikeNormFactors.txt", head = T)
TimeCourse_wigset = rbind(c("ePROseq_Pombe_CDK9as_0min_DMSO_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_0min_DMSO_COMBINED_pombe_minus.bw", "Cdk9as_0min_combined", NFs[26,'Spikereads']/100000),
                          c("ePROseq_Pombe_CDK9as_30sec_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_30sec_3MBPP1_10uM_COMBINED_pombe_minus.bw", "Cdk9as_30sec_combined", NFs[27,'Spikereads']/100000),
                          c("ePROseq_Pombe_CDK9as_1min_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_1min_3MBPP1_10uM_COMBINED_pombe_minus.bw", "Cdk9as_1min_combined", NFs[28,'Spikereads']/100000),
                          c("ePROseq_Pombe_CDK9as_2min30sec_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_2min30sec_3MBPP1_10uM_COMBINED_pombe_minus.bw", "Cdk9as_2min30sec_combined", NFs[29,'Spikereads']/100000),
                          c("ePROseq_Pombe_CDK9as_5min_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_5min_3MBPP1_10uM_COMBINED_pombe_minus.bw", "Cdk9as_5min_combined", NFs[30,'Spikereads']/100000),
                          c("ePROseq_Pombe_CDK9as_7min30sec_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_7min30sec_3MBPP1_10uM_COMBINED_pombe_minus.bw", "Cdk9as_7min30sec_combined", NFs[31,'Spikereads']/100000),
                          c("ePROseq_Pombe_CDK9as_10min_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_10min_3MBPP1_10uM_COMBINED_pombe_minus.bw", "Cdk9as_10min_combined", NFs[32,'Spikereads']/100000),
                          c("ePROseq_Pombe_CDK9as_20min_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_20min_3MBPP1_10uM_COMBINED_pombe_minus.bw", "Cdk9as_20min_combined", NFs[33,'Spikereads']/100000))
# Also look at original 5 minute treatment data (WT and Cdk9as) for comparison
Orig_wigset = rbind(c("5993_7157_26977_HJ72CBGXX_pombe_WT_COMBINED_DMSO_CAGATC_R1.fastq_pombe_plus.bw", "5993_7157_26977_HJ72CBGXX_pombe_WT_COMBINED_DMSO_CAGATC_R1.fastq_pombe_minus.bw", "WT_DMSO_combined", NFs[1,'Spikereads']/100000),
               c("5993_7157_26978_HJ72CBGXX_pombe_WT_COMBINED_3MB-PPI_ACTTGA_R1.fastq_pombe_plus.bw", "5993_7157_26978_HJ72CBGXX_pombe_WT_COMBINED_3MB-PPI_ACTTGA_R1.fastq_pombe_minus.bw", "WT_3MBPP1_combined", NFs[2,'Spikereads']/100000),
               c("5993_7157_26981_HJ72CBGXX_pombe_CDK9as_COMBINED_DMSO_GGCTAC_R1.fastq_pombe_plus.bw", "5993_7157_26981_HJ72CBGXX_pombe_CDK9as_COMBINED_DMSO_GGCTAC_R1.fastq_pombe_minus.bw", "Cdk9as_DMSO_combined", NFs[5,'Spikereads']/100000),
               c("5993_7157_26982_HJ72CBGXX_pombe_CDK9as_COMBINED_3MB-PPI_CTTGTA_R1.fastq_pombe_plus.bw", "5993_7157_26982_HJ72CBGXX_pombe_CDK9as_COMBINED_3MB-PPI_CTTGTA_R1.fastq_pombe_minus.bw", "Cdk9as_3MBPP1_combined", NFs[6,'Spikereads']/100000))
               
#####################################################################################################################
### Functions for Meta plots:
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
#
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
                                  minuteLabs = c("0 min", "30 sec", "1 min", "2.5 min", "5 min", "7.5 min", "10 min", "20 min")){
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

#####################################################################################################################
# Make a metaplot around CPS of each time point compared to the 0 minute. 
TimeCourseMetaData = FormatMetaData_TC_lattice(TimeCourse_wigset, bed = usedbed[, c(1,3,2,4,5,6)], halfWindow = 1000, step = 10)
lattice_meta_TC.proSeq(filename = "ePROseqAllTimePoints_1kb_aroundCPS.pdf",df = TimeCourseMetaData, main = NULL, ylim = c(0,15))

# make plots of each time point vs 0' 
lattice_meta_TC.proSeq(filename = "ePROseq0vs30secTimePoints_1kb_aroundCPS.pdf",df = TimeCourseMetaData[TimeCourseMetaData$sample == 1 | TimeCourseMetaData$sample == 2,],
                       main = NULL, ylim = c(0,15), minuteLabs = c("0 min", "0.5 min"))
lattice_meta_TC.proSeq(filename = "ePROseq0vs1minTimePoints_1kb_aroundCPS.pdf",df = TimeCourseMetaData[TimeCourseMetaData$sample == 1 | TimeCourseMetaData$sample == 3,],
                       main = NULL, ylim = c(0,15), minuteLabs = c("0 min", "1 min"))
lattice_meta_TC.proSeq(filename = "ePROseq0vs2.5minTimePoints_1kb_aroundCPS.pdf",df = TimeCourseMetaData[TimeCourseMetaData$sample == 1 | TimeCourseMetaData$sample == 4,],
                       main = NULL, ylim = c(0,15), minuteLabs = c("0 min", "2.5 min"))
lattice_meta_TC.proSeq(filename = "ePROseq0vs5minTimePoints_1kb_aroundCPS.pdf",df = TimeCourseMetaData[TimeCourseMetaData$sample == 1 | TimeCourseMetaData$sample == 5,],
                       main = NULL, ylim = c(0,15), minuteLabs = c("0 min", "5 min"))
lattice_meta_TC.proSeq(filename = "ePROseq0vs7.5minTimePoints_1kb_aroundCPS.pdf",df = TimeCourseMetaData[TimeCourseMetaData$sample == 1 | TimeCourseMetaData$sample == 6,],
                       main = NULL, ylim = c(0,15), minuteLabs = c("0 min", "7.5 min"))
lattice_meta_TC.proSeq(filename = "ePROseq0vs10minTimePoints_1kb_aroundCPS.pdf",df = TimeCourseMetaData[TimeCourseMetaData$sample == 1 | TimeCourseMetaData$sample == 7,],
                       main = NULL, ylim = c(0,15), minuteLabs = c("0 min", "10 min"))
lattice_meta_TC.proSeq(filename = "ePROseq0vs20minTimePoints_1kb_aroundCPS.pdf",df = TimeCourseMetaData[TimeCourseMetaData$sample == 1 | TimeCourseMetaData$sample == 8,],
                       main = NULL, ylim = c(0,15), minuteLabs = c("0 min", "20 min"))

# Consider Original data. 
OrigMetaData = FormatMetaData_TC_lattice(Orig_wigset, bed = usedbed[, c(1,3,2,4,5,6)], halfWindow = 1000, step = 10)
lattice_meta_TC.proSeq(filename = "Original_WT_Cdk9asData_1kb_aroundCPS.pdf", df = OrigMetaData, main = NULL, ylim = c(0,15), minuteLabs = c("WT_DMSO", "WT_3MBPP1", "Cdk9as_DMSO", "Cdk9as_3MBPP1"))

# make plots of each time point vs 0' 
lattice_meta_TC.proSeq(filename = "Original_WTDMSO_WT3MBPP1_1kb_aroundCPS.pdf",df = OrigMetaData[OrigMetaData$sample == 1 | OrigMetaData$sample == 2,],
                       main = NULL, ylim = c(0,15), minuteLabs = c("WT_DMSO", "WT_3MBPP1"))
lattice_meta_TC.proSeq(filename = "Original_WTDMSO_Cdk9DMSO_1kb_aroundCPS.pdf",df = OrigMetaData[OrigMetaData$sample == 1 | OrigMetaData$sample == 3,],
                       main = NULL, ylim = c(0,15), minuteLabs = c("WT_DMSO", "Cdk9as_DMSO"))
lattice_meta_TC.proSeq(filename = "Original_WTDMSO_Cdk93MBPP1_1kb_aroundCPS.pdf",df = OrigMetaData[OrigMetaData$sample == 1 | OrigMetaData$sample == 4,],
                       main = NULL, ylim = c(0,15), minuteLabs = c("WT_DMSO", "Cdk9as_3MBPP1"))
lattice_meta_TC.proSeq(filename = "Original_Cdk9DMSO_Cdk93MBPP1_1kb_aroundCPS.pdf",df = OrigMetaData[OrigMetaData$sample == 3 | OrigMetaData$sample == 4,],
                       main = NULL, ylim = c(0,15), minuteLabs = c("Cdk9as_DMSO", "Cdk9as_3MBPP1"))




#####################################################################################################################
# Make fold change heat maps around CPS of each time point compared to the 0 minute. 
# sort genes by length, even though we are looking at the gene ends
## expect longer genes to have unchanged post-CPS distributions based on transcription rates (only really at early time points).
lenSortGL = orderGenesbyLen(usedbed, StoL = F) ## note because we have to transpose the matrix for plotting, sorting must be opposite what you want on plot (top to bottom)
## running function in GBfunctions to prepare and format heat data. for each comparison
halfMin_data = compareHeats(lenSortGL[,c(1,3,2,4,5,6)], TimeCourse_wigset[1,], TimeCourse_wigset[2,], step = 10, maxLength = 1000)
oneMin_data = compareHeats(lenSortGL[,c(1,3,2,4,5,6)], TimeCourse_wigset[1,], TimeCourse_wigset[3,], step = 10, maxLength = 1000)
two.5Min_data = compareHeats(lenSortGL[,c(1,3,2,4,5,6)], TimeCourse_wigset[1,], TimeCourse_wigset[4,], step = 10, maxLength = 1000)
fiveMin_data = compareHeats(lenSortGL[,c(1,3,2,4,5,6)], TimeCourse_wigset[1,], TimeCourse_wigset[5,], step = 10,maxLength = 1000)
seven.5Min_data = compareHeats(lenSortGL[,c(1,3,2,4,5,6)], TimeCourse_wigset[1,], TimeCourse_wigset[6,], step = 10, maxLength = 1000)
tenMin_data = compareHeats(lenSortGL[,c(1,3,2,4,5,6)], TimeCourse_wigset[1,], TimeCourse_wigset[7,], step = 10, maxLength = 1000)
twentyMin_data = compareHeats(lenSortGL[,c(1,3,2,4,5,6)], TimeCourse_wigset[1,], TimeCourse_wigset[8,], step = 10, maxLength = 1000)
# run plotting function from GBfunctions for each above time point comparison
dir.create(paste(fig_dir,"/lengthSort/", sep = ""))
Plot3Heat(halfMin_data, filename = "/lengthSort/ePROseq_Cdk9as_TC_0vs0.5min_kbsepGenes_Heatmapst_postCPS.pdf")
Plot3Heat(oneMin_data, filename = "/lengthSort/ePROseq_Cdk9as_TC_0vs1min_kbsepGenes_Heatmapst_postCPS.pdf")
Plot3Heat(two.5Min_data, filename = "/lengthSort/ePROseq_Cdk9as_TC_0vs2.5min_kbsepGenes_Heatmapst_postCPS.pdf")
Plot3Heat(fiveMin_data, filename = "/lengthSort/ePROseq_Cdk9as_TC_0vs5min_kbsepGenes_Heatmaps_postCPS.pdf")
Plot3Heat(seven.5Min_data, filename = "/lengthSort/ePROseq_Cdk9as_TC_0vs7.5min_kbsepGenes_Heatmapst_postCPS.pdf")
Plot3Heat(tenMin_data, filename = "/lengthSort/ePROseq_Cdk9as_TC_0vs10min_kbsepGenes_Heatmaps_postCPS.pdf")
Plot3Heat(twentyMin_data, filename = "/lengthSort/ePROseq_Cdk9as_TC_0vs20min_kbsepGenes_Heatmapst_postCPS.pdf")


# Run the same set of functions, but for the original data (only comparing Cdk9as DMSO vs 3MBPP1 5' treatment)
Orig_data = compareHeats(lenSortGL[,c(1,3,2,4,5,6)], Orig_wigset[3,], Orig_wigset[4,], step = 10, maxLength = 1000)
# run plotting function from GBfunctions for each above time point comparison
Plot3Heat(Orig_data, filename = "/lengthSort/Original_Cdk9DMSO_Cdk93MBPP1_Heatmapst_postCPS.pdf")









