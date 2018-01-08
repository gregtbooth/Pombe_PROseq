source("/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/scripts/GB_Functions.R")
library(gplots)
bwpath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/alignment/bw/"
infopath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/pombeinfo/"
fig_dir = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/figures/MetaPlotComparison/12-20-17/"
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
wigset = rbind(c("5993_7157_26977_HJ72CBGXX_pombe_WT_COMBINED_DMSO_CAGATC_R1.fastq_pombe_plus.bw", "5993_7157_26977_HJ72CBGXX_pombe_WT_COMBINED_DMSO_CAGATC_R1.fastq_pombe_minus.bw", "WT_DMSO_combined", NFs[1,'Spikereads']/100000),
               c("5993_7157_26978_HJ72CBGXX_pombe_WT_COMBINED_3MB-PPI_ACTTGA_R1.fastq_pombe_plus.bw", "5993_7157_26978_HJ72CBGXX_pombe_WT_COMBINED_3MB-PPI_ACTTGA_R1.fastq_pombe_minus.bw", "WT_3MBPP1_combined", NFs[2,'Spikereads']/100000),
               c("5993_7157_26979_HJ72CBGXX_pombe_mcs6as_COMBINED_DMSO_GATCAG_R1.fastq_pombe_plus.bw", "5993_7157_26979_HJ72CBGXX_pombe_mcs6as_COMBINED_DMSO_GATCAG_R1.fastq_pombe_minus.bw", "mcs6as_DMSO_combined", NFs[3,'Spikereads']/100000),
               c("5993_7157_26980_HJ72CBGXX_pombe_mcs6as_rep2_3MB-PPI_TAGCTT_R1.fastq_pombe_plus.bw", "5993_7157_26980_HJ72CBGXX_pombe_mcs6as_rep2_3MB-PPI_TAGCTT_R1.fastq_pombe_minus.bw", "mcs6as_3MBPP1_rep2only", repNFs[8,'Spikereads']/100000),
               c("5993_7157_26981_HJ72CBGXX_pombe_CDK9as_COMBINED_DMSO_GGCTAC_R1.fastq_pombe_plus.bw", "5993_7157_26981_HJ72CBGXX_pombe_CDK9as_COMBINED_DMSO_GGCTAC_R1.fastq_pombe_minus.bw", "Cdk9as_DMSO_combined", NFs[5,'Spikereads']/100000),
               c("5993_7157_26982_HJ72CBGXX_pombe_CDK9as_COMBINED_3MB-PPI_CTTGTA_R1.fastq_pombe_plus.bw", "5993_7157_26982_HJ72CBGXX_pombe_CDK9as_COMBINED_3MB-PPI_CTTGTA_R1.fastq_pombe_minus.bw", "Cdk9as_3MBPP1_combined", NFs[6,'Spikereads']/100000),
               c("6713_7157_31629_HWV7YBGXX_pombe_Isk1as_DMSO_COMBINED_CAGATC_R1.fastq_pombe_plus.bw", "6713_7157_31629_HWV7YBGXX_pombe_Isk1as_DMSO_COMBINED_CAGATC_R1.fastq_pombe_minus.bw", "Isk1as_DMSO_combined", NFs[7,'Spikereads']/100000),
               c("6713_7157_31630_HWV7YBGXX_pombe_Isk1as_3MB-PP1_COMBINED_ACTTGA_R1.fastq_pombe_plus.bw", "6713_7157_31630_HWV7YBGXX_pombe_Isk1as_3MB-PP1_COMBINED_ACTTGA_R1.fastq_pombe_minus.bw", "Isk1as_3MBPP1_combined", NFs[8,'Spikereads']/100000),
               c("6713_7157_31625_HWV7YBGXX_pombe_Mcs6as_CDK9as_DMSO_COMBINED_TTAGGC_R1.fastq_pombe_plus.bw", "6713_7157_31625_HWV7YBGXX_pombe_Mcs6as_CDK9as_DMSO_COMBINED_TTAGGC_R1.fastq_pombe_minus.bw", "Mcs6as_Cdk9as_DMSO_combined", NFs[9,'Spikereads']/100000),
               c("6713_7157_31626_HWV7YBGXX_pombe_Mcs6as_CDK9as_3MB-PP1_COMBINED_TGACCA_R1.fastq_pombe_plus.bw", "6713_7157_31626_HWV7YBGXX_pombe_Mcs6as_CDK9as_3MB-PP1_COMBINED_TGACCA_R1.fastq_pombe_minus.bw", "Mcs6as_Cdk9as_3MBPP1_combined", NFs[10,'Spikereads']/100000),
               c("7772_7157_43149_HVLYCBGXY_pombe_18C_Cdk9as_DMSO_COMBINED_GGCTAC_R1_pombe_plus.bw", "7772_7157_43149_HVLYCBGXY_pombe_18C_Cdk9as_DMSO_COMBINED_GGCTAC_R1_pombe_minus.bw", "Cdk9as_DMSO_18C_combined", NFs[11,'Spikereads']/100000),
               c("7772_7157_43150_HVLYCBGXY_pombe_18C_Cdk9as_3MBPP1_COMBINED_CTTGTA_R1_pombe_plus.bw", "7772_7157_43150_HVLYCBGXY_pombe_18C_Cdk9as_3MBPP1_COMBINED_CTTGTA_R1_pombe_minus.bw", "Cdk9as_3MBPP1_18C_combined", NFs[12,'Spikereads']/100000),
               c("7772_7157_43151_HVLYCBGXY_pombe_18C_dis2ts_DMSO_COMBINED_ATTCCT_R1_pombe_plus.bw", "7772_7157_43151_HVLYCBGXY_pombe_18C_dis2ts_DMSO_COMBINED_ATTCCT_R1_pombe_minus.bw", "dis2ts_DMSO_18C_combined", NFs[13,'Spikereads']/100000),
               c("7772_7157_43152_HVLYCBGXY_pombe_18C_dis2ts_3MBPP1_COMBINED_CAAAAG_R1_pombe_plus.bw", "7772_7157_43152_HVLYCBGXY_pombe_18C_dis2ts_3MBPP1_COMBINED_CAAAAG_R1_pombe_minus.bw", "dis2ts_3MBPP1_18C_combined", NFs[14,'Spikereads']/100000),
               c("7772_7157_43153_HVLYCBGXY_pombe_30C_dis2ts_DMSO_COMBINED_CAACTA_R1_pombe_plus.bw", "7772_7157_43153_HVLYCBGXY_pombe_30C_dis2ts_DMSO_COMBINED_CAACTA_R1_pombe_minus.bw", "dis2ts_DMSO_30C_combined", NFs[15,'Spikereads']/100000),
               c("7772_7157_43154_HVLYCBGXY_pombe_30C_dis2ts_3MBPP1_COMBINED_CACCGG_R1_pombe_plus.bw", "7772_7157_43154_HVLYCBGXY_pombe_30C_dis2ts_3MBPP1_COMBINED_CACCGG_R1_pombe_minus.bw", "dis2ts_3MBPP1_30C_combined", NFs[16,'Spikereads']/100000),
               c("7772_7157_43155_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_DMSO_COMBINED_CACGAT_R1_pombe_plus.bw", "7772_7157_43155_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_DMSO_COMBINED_CACGAT_R1_pombe_minus.bw", "CDK9as_dis2ts_DMSO_18C_combined", NFs[17,'Spikereads']/100000),
               c("7772_7157_43156_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_3MBPP1_COMBINED_CACTCA_R1_pombe_plus.bw", "7772_7157_43156_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_3MBPP1_COMBINED_CACTCA_R1_pombe_minus.bw", "CDK9as_dis2ts_3MBPP1_18C_combined", NFs[18,'Spikereads']/100000),
               c("7772_7157_43157_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_DMSO_COMBINED_CAGGCG_R1_pombe_plus.bw", "7772_7157_43157_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_DMSO_COMBINED_CAGGCG_R1_pombe_minus.bw", "CDK9as_dis2ts_DMSO_30C_combined", NFs[19,'Spikereads']/100000),
               c("7772_7157_43158_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_3MBPP1_COMBINED_CATGGC_R1_pombe_plus.bw", "7772_7157_43158_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_3MBPP1_COMBINED_CATGGC_R1_pombe_minus.bw", "CDK9as_dis2ts_3MBPP1_30C_combined", NFs[20,'Spikereads']/100000),
               c("ePROseq_Pombe_Lsk1asCDK9as_DMSO_COMBINED_pombe_plus.bw", "ePROseq_Pombe_Lsk1asCDK9as_DMSO_COMBINED_pombe_minus.bw", "Lsk1as_Cdk9as_DMSO_combined", NFs[34,'Spikereads']/100000),
               c("ePROseq_Pombe_Lsk1asCDK9as_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_Lsk1asCDK9as_3MBPP1_10uM_COMBINED_pombe_minus.bw", "Lsk1as_Cdk9as_3MBPP1_combined", NFs[35,'Spikereads']/100000),)
          
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
                               xlab = "Distance to TSS (bp)",  ylab = "Median PRO-seq intensity"){
  pdf(paste(fig_dir, filename), width = 30, height = 10)
  result <- xyplot(mean ~ x | factor(background), data = df,
  #result <- xyplot(mean ~ x | factor(strain) + factor(treatment), data = df, 
                   group = factor(treatment, labels = c("DMSO", "3MB-PP1")), 
                   #group = factor(temp, labels = c("30 deg C", "18 deg C")), 
                   scales = list(tck=c(1,0),alternating = c(1,1),x=list(relation='free',axs='i'),
                                 y=list(relation='free',axs='i')),
                   key=list(corner=c(0.98,0.85),padding.text=3,
                            text=list(c("DMSO", "3MB-PP1"),col=c("black", "black"), cex=0.8, font=2),
                            #text=list(c("30 deg C", "18 deg C"),col=c("black", "black"), cex=0.8, font=2),
                            rectangles=list(col=c(rgb(0,0,0.5), rgb(0.5,0,0)), size=2.7, height=0.8, border='white')),
                   type = 'l', 
                   ylim = ylim,
                   col = c(rgb(0.5,0,0,0.9), rgb(0,0,0.5,0.9)),# rgb(0,0.5,0,0.9)),  ## just add or take away colors to match number of samples in each panel
                   fill = c(rgb(0.5,0,0,0.3), rgb(0,0,0.5,0.3)),# rgb(0,0.5,0,0.3)),
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

lattice_Scaledmeta.proSeq = function(filename= "test.pdf", df, ylim = c(0, 10), main = NULL, 
                               xlab = "Distance to gene boundaries (bp)",  ylab = "Median PRO-seq intensity"){
  pdf(paste(fig_dir, filename), width = 35, height = 10)
  result <- xyplot(mean ~ x | factor(background), data = df,
  #result <- xyplot(mean ~ x | factor(strain) + factor(treatment), data = df, 
                   group = factor(treatment, labels = c("DMSO", "3MB-PP1")), 
                   #group = factor(temp, labels = c("30 deg C", "18 deg C")),
                   scales = list(tck=c(1,0),alternating = c(1,1),
                                 x=list(relation='free',axs='i', labels = c('-1000', 'TSS', "+300", "-300", "CPS", "1000"), at = c(1, 101, 131, 191, 221, 320)),
                                 y=list(relation='free',axs='i')),
                   key=list(corner=c(0.98,0.85),padding.text=3,
                            text=list(c("DMSO", "3MB-PP1"),col=c("black", "black"), cex=0.8, font=2),
                            rectangles=list(col=c(rgb(0,0,0.5), rgb(0.5,0,0)), size=2.7, height=0.8, border='white')),
                            #text=list(c("Cdk9as", "Mcs6as", "WT"),col=c("black", "black", "black"), cex=0.8, font=2),
                            #rectangles=list(col=c(rgb(0.5,0,0), rgb(0,0,0.5), rgb(0,0.5,0)), size=2.7, height=0.8, border='white')),
                   type = 'l', 
                   ylim = ylim,
                   col = c(rgb(0.5,0,0,0.9), rgb(0,0,0.5,0.9)),# rgb(0,0.5,0,0.9)),  ## just add or take away colors to match number of samples in each panel
                   fill = c(rgb(0.5,0,0,0.3), rgb(0,0,0.5,0.3)),# rgb(0,0.5,0,0.3)),
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

## formatting Scaled meta plots for lattice
FormatScaledMetaData_lattice <- function(wigset, bed){
  N = dim(wigset)[1]
  df = data.frame(matrix(ncol = 9, nrow = 0), stringsAsFactors = F)
  sampleNames = c()
  BG_Indx = 0 # counter for tracking genotype, which changes every other sample.
  for (i in 1:N){
    if (i == ceiling(N/2)){cat("* 50% complete ... \n")}
    cat("* loading", i, "\n")
    wigs = load.wigset(wigset, i) #this should give all the bigwig files
    cat("* generating and combining all sample scaled meta-plot data ...\n")
    meta= meta.subsample.scaled(bed, bigWig.plus = wigs[[1]], bigWig.minus = wigs[[2]], 1000, 300, 10, do.sum = T)
    metaNorm = meta.normalize(result = meta, scaleFactor = 1/as.numeric(wigs[[4]]))
    sampleNames = c(sampleNames, wigs[[3]])
    sampleVal = i ### index for later retrieval of sample name from sampleNames list.
    if ((-1)**i == -1){BG_Indx = BG_Indx + 1}
    else{BG_Indx = BG_Indx}
    if (grepl("DMSO", wigs[[3]])){ #grepl returns TRUE or false based on the presence or absence of "DMSO" in sample name
      treatment = 0
    }
    else{treatment = 1}
    # this section is to distinguish between strains for the purpose of later separating based on temperature
    cdk9 = c("CDK9", "Cdk9")
    mcs6 = c("mcs6", "Mcs6")
    if(grepl(paste(cdk9, collapse = "|"), wigs[[3]]) & grepl("Mcs6as_Cdk9as", wigs[[3]])) {strain = 4}
    else if(grepl(paste(cdk9, collapse = "|"), wigs[[3]]) & grepl("CDK9as_dis2ts", wigs[[3]])) {strain = 7}
    else if(grepl(paste(cdk9, collapse = "|"), wigs[[3]]) & grepl("Lsk1as_Cdk9as", wigs[[3]])) {strain = 8}
    else  if(grepl(paste(cdk9, collapse = "|"), wigs[[3]])) {strain = 2}
    if (grepl(paste(mcs6, collapse = "|"), wigs[[3]]) & grepl("Mcs6as_Cdk9as", wigs[[3]])) {strain = 4}
    else if (grepl(paste(mcs6, collapse = "|"), wigs[[3]])) {strain = 3}
    if(grepl("dis2", wigs[[3]]) & grepl("CDK9as_dis2ts", wigs[[3]])) {strain = 7}
    else if(grepl("dis2",  wigs[[3]])) {strain = 6}
    if(grepl("WT", wigs[[3]])) {strain = 1}
    if(grepl("Isk1", wigs[[3]])) {strain = 5}
    if(grepl("18C", wigs[[3]])){temp = 18}
    else{temp = 30}
    print(strain)
    print(temp)
    xAxis = seq(from = 1, to = 320)
    sample_df = cbind(xAxis, metaNorm[[4]], metaNorm[[3]], metaNorm[[2]], sampleVal, treatment, BG_Indx, strain, temp)
    df = rbind.data.frame(df, sample_df)
  }
  colnames(df) = c("x", "mean", "lower", "upper", "sample", 'treatment', 'background', "strain", "temp")
  # lastly I will replace values in the 'sample', 'treatment', and 'background' columns with actual names
  for (ii in 1:length(sampleNames)){
    df$sample[df$sample == ii] <- sampleNames[ii]
  }
  df$treatment[df$treatment == 0] <- "DMSO"
  df$treatment[df$treatment == 1] <- "3-MB_PP1"
  df$background[df$background == 1] <- "WT"
  df$background[df$background == 2] <- "Mcs6_as"
  df$background[df$background == 3] <- "CDK9_as"
  df$background[df$background == 4] <- "Lsk1_as"
  df$background[df$background == 5] <- "Mcs6as_CDK9_as"
  df$background[df$background == 6] <- "CDK9_as_18C"
  df$background[df$background == 7] <- "dis2_ts_18C"
  df$background[df$background == 8] <- "dis2_ts_30C"
  df$background[df$background == 9] <- "CDK9as_dis2ts_18C"
  df$background[df$background == 10] <- "CDK9as_dis2ts_30C"
  df$background[df$background == 11] <- "Lsk1as_Cdk9as"
  df$strain[df$strain == 1] <- "WT"
  df$strain[df$strain == 2] <- "CDK9_as"
  df$strain[df$strain == 3] <- "Mcs6_as"
  df$strain[df$strain == 4] <- "Mcs6as_CDK9_as"
  df$strain[df$strain == 5] <- "Lsk1_as"
  df$strain[df$strain == 6] <- "dis2_ts"
  df$strain[df$strain == 7] <- "CDK9as_dis2ts"
  df$strain[df$strain == 8] <- "Lsk1as_Cdk9as"
  return(df)
}
#####################################################################################################################
#plot scaled meta-plots with lattice. 
scaledMetaData = FormatScaledMetaData_lattice(wigset, bed)
lattice_Scaledmeta.proSeq(filename = "Experiments_byDrugTreatment_ScaledMetas_kbSepGenes.pdf", df = scaledMetaData, ylim = c(0, 30), main = "scaled Gene Length", 
                          xlab = "Distance to Gene boundaries (bp)", ylab = "Median PRO-seq intensity")
## Plot scaled meta plots comparing just WT, Mcs6as, and Cdk9as (original samples) before and after adding 3MBPP1.
## this plot is used in Extended Data Fig2. 
## Need to first adjust the grouping in the lattice_Scaledmeta.proSeq function. 
scaledMetaData_ED2 = FormatScaledMetaData_lattice(wigset[c(1:6),], bed)
lattice_Scaledmeta.proSeq(filename = "ED2_OrigExperiments_byDrugTreatment_ScaledMetas_kbSepGenes.pdf", df = scaledMetaData_ED2, ylim = c(0, 25), main = "scaled Gene Length", 
                          xlab = "Distance to Gene boundaries (bp)", ylab = "Median PRO-seq intensity")


#####################################################################################################################
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

#####################################################################################################################
## formatting meta plots for lattice (not scaled Metaplots)
FormatMetaData_lattice <- function(wigset, bed, halfWindow, step){
  N = dim(wigset)[1]
  df = data.frame(matrix(ncol = 7, nrow = 0), stringsAsFactors = F)
  sampleNames = c()
  BG_Indx = 0 # counter for tracking genotype, which changes every other sample.
  for (i in 1:N){
    if (i == ceiling(N/2)){cat("* 50% complete ... \n")}
    cat("* loading", i, "\n")
    wigs = load.wigset(wigset, i) #this should give all the bigwig files
    cat("* generating and combining all sample meta-plot data ...\n")
    meta= meta.subsample(bed, bigWig.plus = wigs[[1]], bigWig.minus = wigs[[2]], 
                         step = step, at.TSS = T, halfWindow = halfWindow, do.sum = T)
    metaNorm = meta.normalize(result = meta, scaleFactor = 1/as.numeric(wigs[[4]]))
    sampleNames = c(sampleNames, wigs[[3]])
    sampleVal = i ### index for later retrieval of sample name from sampleNames list.
    if ((-1)**i == -1){BG_Indx = BG_Indx + 1}
    else{BG_Indx = BG_Indx}
    if (grepl("DMSO", wigs[[3]])){ #grepl returns TRUE or false based on the presence or absence of "DMSO" in sample name
      treatment = 0
    }
    else{treatment = 1}
    # this section is to distinguish between strains for the purpose of later separating based on temperature
    cdk9 = c("CDK9", "Cdk9")
    mcs6 = c("mcs6", "Mcs6")
    if(grepl(paste(cdk9, collapse = "|"), wigs[[3]]) & grepl("Mcs6as_Cdk9as", wigs[[3]])) {strain = 4}
    else if(grepl(paste(cdk9, collapse = "|"), wigs[[3]]) & grepl("CDK9as_dis2ts", wigs[[3]])) {strain = 7}
    else if(grepl(paste(cdk9, collapse = "|"), wigs[[3]]) & grepl("Lsk1as_Cdk9as", wigs[[3]])) {strain = 8}
    else  if(grepl(paste(cdk9, collapse = "|"), wigs[[3]])) {strain = 2}
    if (grepl(paste(mcs6, collapse = "|"), wigs[[3]]) & grepl("Mcs6as_Cdk9as", wigs[[3]])) {strain = 4}
    else if (grepl(paste(mcs6, collapse = "|"), wigs[[3]])) {strain = 3}
    if(grepl("dis2", wigs[[3]]) & grepl("CDK9as_dis2ts", wigs[[3]])) {strain = 7}
    else if(grepl("dis2",  wigs[[3]])) {strain = 6}
    if(grepl("WT", wigs[[3]])) {strain = 1}
    if(grepl("Isk1", wigs[[3]])) {strain = 5}
    if(grepl("18C", wigs[[3]])){temp = 18}
    else{temp = 30}
    print(strain)
    print(temp)
    xAxis = seq(from = -halfWindow, to = halfWindow, by = step)
    sample_df = cbind(xAxis, metaNorm[[4]], metaNorm[[3]], metaNorm[[2]], sampleVal, treatment, BG_Indx, strain, temp)
    df = rbind.data.frame(df, sample_df)
  }
  colnames(df) = c("x", "mean", "lower", "upper", "sample", 'treatment', 'background', 'strain', 'temp')
  # lastly I will replace values in the 'sample', 'treatment', and 'background' columns with actual names
  for (ii in 1:length(sampleNames)){
    df$sample[df$sample == ii] <- sampleNames[ii]
  }
  df$treatment[df$treatment == 0] <- "DMSO"
  df$treatment[df$treatment == 1] <- "3-MB_PP1"
  df$background[df$background == 1] <- "WT"
  df$background[df$background == 2] <- "Mcs6_as"
  df$background[df$background == 3] <- "CDK9_as"
  df$background[df$background == 4] <- "Isk1_as"
  df$background[df$background == 5] <- "Mcs6as_CDK9_as"
  df$background[df$background == 6] <- "CDK9_as_18C"
  df$background[df$background == 7] <- "dis2_ts_18C"
  df$background[df$background == 8] <- "dis2_ts_30C"
  df$background[df$background == 9] <- "CDK9as_dis2ts_18C"
  df$background[df$background == 10] <- "CDK9as_dis2ts_30C"
  df$background[df$background == 11] <- "Lsk1as_Cdk9as"
  df$strain[df$strain == 1] <- "WT"
  df$strain[df$strain == 2] <- "CDK9_as"
  df$strain[df$strain == 3] <- "Mcs6_as"
  df$strain[df$strain == 4] <- "Mcs6as_CDK9_as"
  df$strain[df$strain == 5] <- "Lsk1_as"
  df$strain[df$strain == 6] <- "dis2_ts"
  df$strain[df$strain == 7] <- "CDK9as_dis2ts"
  df$strain[df$strain == 8] <- "Lsk1as_Cdk9as"
  return(df)
}
#####################################################################################################################

### Plots around TSS and CPS (for KBsep genes)
# TSS
TSSmetaData = FormatMetaData_lattice(wigset, bed, halfWindow = 1000, step = 10)
lattice_meta.proSeq(filename = "SP_TSS_experiments_byDrugTreatment_Meta_expY.pdf", df = TSSmetaData, ylim = c(0, 30), main = 'Distance to TSS', 
                    xlab = "Distance to TSS (bp)",  ylab = "Median PRO-seq intensity")
# CPS
CPSmetaData = FormatMetaData_lattice(wigset, bed[, c(1,3,2,4,5,6)], halfWindow = 1000, step = 10)
lattice_meta.proSeq(filename = "SP_CPS_experiments_byDrugTreatment_Meta_expY.pdf", df = CPSmetaData, ylim = c(0, 30), main = 'Distance to CPS', 
                    xlab = "Distance to CPS (bp)",  ylab = "Median PRO-seq intensity")
### Plots around GB nucleosomes
GBnuclmetaData = FormatMetaData_lattice(wigset, gbnucls, halfWindow = 80, step = 5)
lattice_meta.proSeq(filename = "SP_PROaroundGBnuclDyad_Allexperiments_Meta.pdf", df = GBnuclmetaData, ylim = c(4, 15), main = 'Distance to Gene-body Nucleosome Dyad', 
                    xlab = "Distance to Dyad (bp)",  ylab = "Median PRO-seq intensity")

### Plots around Splice Sites
# 3pSS
SS3PmetaData = FormatMetaData_lattice(wigset, Introns[, c(1,3,2,4,5,6)], halfWindow = 80, step = 5)
lattice_meta.proSeq(filename = "SP_3pSS_Allexperiments_Meta.pdf", df = SS3PmetaData, ylim = c(0, 15), main = 'Distance to 3p Splice Site', 
                    xlab = "Distance to 3p SS (bp)",  ylab = "Median PRO-seq intensity")
# 5pSS
SS5PmetaData = FormatMetaData_lattice(wigset, Introns, halfWindow = 80, step = 5)
lattice_meta.proSeq(filename = "SP_5pSS_Allexperiments_Meta.pdf", df = SS5PmetaData, ylim = c(0, 15), main = 'Distance to 5p Splice Site', 
                    xlab = "Distance to 5p SS (bp)",  ylab = "Median PRO-seq intensity")

#####################################################################################################################
#####################################################################################################################
# Functions Specifically for Time course data
#####################################################################################################################
#####################################################################################################################
TimeCourse_wigset = rbind(c("ePROseq_Pombe_CDK9as_0min_DMSO_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_0min_DMSO_COMBINED_pombe_minus.bw", "Cdk9as_0min_combined", NFs[26,'Spikereads']/100000),
                          c("ePROseq_Pombe_CDK9as_30sec_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_30sec_3MBPP1_10uM_COMBINED_pombe_minus.bw", "Cdk9as_30sec_combined", NFs[27,'Spikereads']/100000),
                          c("ePROseq_Pombe_CDK9as_1min_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_1min_3MBPP1_10uM_COMBINED_pombe_minus.bw", "Cdk9as_1min_combined", NFs[28,'Spikereads']/100000),
                          c("ePROseq_Pombe_CDK9as_2min30sec_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_2min30sec_3MBPP1_10uM_COMBINED_pombe_minus.bw", "Cdk9as_2min30sec_combined", NFs[29,'Spikereads']/100000),
                          c("ePROseq_Pombe_CDK9as_5min_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_5min_3MBPP1_10uM_COMBINED_pombe_minus.bw", "Cdk9as_5min_combined", NFs[30,'Spikereads']/100000),
                          c("ePROseq_Pombe_CDK9as_7min30sec_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_7min30sec_3MBPP1_10uM_COMBINED_pombe_minus.bw", "Cdk9as_7min30sec_combined", NFs[31,'Spikereads']/100000),
                          c("ePROseq_Pombe_CDK9as_10min_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_10min_3MBPP1_10uM_COMBINED_pombe_minus.bw", "Cdk9as_10min_combined", NFs[32,'Spikereads']/100000),
                          c("ePROseq_Pombe_CDK9as_20min_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_20min_3MBPP1_10uM_COMBINED_pombe_minus.bw", "Cdk9as_20min_combined", NFs[33,'Spikereads']/100000))
#####################################################################################################################

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
lattice_meta_TC.proSeq(filename = "ePROseqAllTimePoints_longGenes6KB_300bpsmooth.pdf",df = TimeCourseMetaData, main = NULL, ylim = c(0,600))

# make plots of each time point vs 0' 
lattice_meta_TC.proSeq(filename = "ePROseq0vs30secTimePoints_longGenes6KB_300bpsmooth.pdf",df = TimeCourseMetaData[TimeCourseMetaData$sample == 1 | TimeCourseMetaData$sample == 2,],
                       main = NULL, ylim = c(0,600), minuteLabs = c("0 min", "30 sec"))
lattice_meta_TC.proSeq(filename = "ePROseq0vs1minTimePoints_longGenes6KB_300bpsmooth.pdf",df = TimeCourseMetaData[TimeCourseMetaData$sample == 1 | TimeCourseMetaData$sample == 3,],
                       main = NULL, ylim = c(0,600), minuteLabs = c("0 min", "1 min"))
lattice_meta_TC.proSeq(filename = "ePROseq0vs2min30secTimePoints_longGenes6KB_300bpsmooth.pdf",df = TimeCourseMetaData[TimeCourseMetaData$sample == 1 | TimeCourseMetaData$sample == 4,],
                       main = NULL, ylim = c(0,600), minuteLabs = c("0 min", "2.5 min"))
lattice_meta_TC.proSeq(filename = "ePROseq0vs5minTimePoints_longGenes6KB_300bpsmooth.pdf",df = TimeCourseMetaData[TimeCourseMetaData$sample == 1 | TimeCourseMetaData$sample == 5,],
                       main = NULL, ylim = c(0,600), minuteLabs = c("0 min", "5 min"))
lattice_meta_TC.proSeq(filename = "ePROseq0vs7min30secTimePoints_longGenes6KB_300bpsmooth.pdf",df = TimeCourseMetaData[TimeCourseMetaData$sample == 1 | TimeCourseMetaData$sample == 6,],
                       main = NULL, ylim = c(0,600), minuteLabs = c("0 min", "7.5 min"))
lattice_meta_TC.proSeq(filename = "ePROseq0vs10minTimePoints_longGenes6KB_300bpsmooth.pdf",df = TimeCourseMetaData[TimeCourseMetaData$sample == 1 | TimeCourseMetaData$sample == 7,],
                       main = NULL, ylim = c(0,600), minuteLabs = c("0 min", "10 min"))
lattice_meta_TC.proSeq(filename = "ePROseq0vs20minTimePoints_longGenes6KB_300bpsmooth.pdf",df = TimeCourseMetaData[TimeCourseMetaData$sample == 1 | TimeCourseMetaData$sample == 8,],
                       main = NULL, ylim = c(0,600), minuteLabs = c("0 min", "20 min"))





