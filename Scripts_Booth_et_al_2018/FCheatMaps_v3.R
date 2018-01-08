source("/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/scripts/GB_Functions.R")
library(gplots)
library(ggplot2)
library(lattice)
bwpath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/alignment/bw/"
infopath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/pombeinfo/"
fig_dir = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/figures/FC_heatMaps/10-02-17/"
dir.create(fig_dir)
SPmappability = load.bigWig(paste(infopath, "SPombeMappable_SacCerNonOverLap_expanded_ChrFix.bw", sep = ""))
SPGL_1kbsep = read.table(file = paste(infopath, "SP_PROcapObservedTSS_1000bpSep.bed", sep = ""), stringsAsFactors = F)
ObsTSS_NOoverlap = read.table(file = paste(infopath, "SP_PROcapObservedTSS_noOverlap.bed", sep = ""))
filteredGL = read.table(paste(infopath, "SP_CompleteFilteredGenes.bed", sep = ""))
UntreatedCountDat = read.table(file = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/countData/AllGenesObsTSS/ePRO_Cdk9as_0min_combinedPIcountData.txt")
Min20CountDat = read.table(file = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/countData/AllGenesObsTSS/ePRO_Cdk9as_20min_combinedPIcountData.txt")
###################################################################################################################
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
               c("ePROseq_Pombe_Lsk1asCDK9as_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_Lsk1asCDK9as_3MBPP1_10uM_COMBINED_pombe_minus.bw", "Lsk1as_Cdk9as_3MBPP1_combined", NFs[35,'Spikereads']/100000))

TimeCourse_wigset = rbind(c("ePROseq_Pombe_CDK9as_0min_DMSO_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_0min_DMSO_COMBINED_pombe_minus.bw", "Cdk9as_0min_combined", NFs[26,'Spikereads']/100000),
                          c("ePROseq_Pombe_CDK9as_30sec_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_30sec_3MBPP1_10uM_COMBINED_pombe_minus.bw", "Cdk9as_30sec_combined", NFs[27,'Spikereads']/100000),
                          c("ePROseq_Pombe_CDK9as_1min_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_1min_3MBPP1_10uM_COMBINED_pombe_minus.bw", "Cdk9as_1min_combined", NFs[28,'Spikereads']/100000),
                          c("ePROseq_Pombe_CDK9as_2min30sec_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_2min30sec_3MBPP1_10uM_COMBINED_pombe_minus.bw", "Cdk9as_2min30sec_combined", NFs[29,'Spikereads']/100000),
                          c("ePROseq_Pombe_CDK9as_5min_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_5min_3MBPP1_10uM_COMBINED_pombe_minus.bw", "Cdk9as_5min_combined", NFs[30,'Spikereads']/100000),
                          c("ePROseq_Pombe_CDK9as_7min30sec_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_7min30sec_3MBPP1_10uM_COMBINED_pombe_minus.bw", "Cdk9as_7min30sec_combined", NFs[31,'Spikereads']/100000),
                          c("ePROseq_Pombe_CDK9as_10min_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_10min_3MBPP1_10uM_COMBINED_pombe_minus.bw", "Cdk9as_10min_combined", NFs[32,'Spikereads']/100000),
                          c("ePROseq_Pombe_CDK9as_20min_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_20min_3MBPP1_10uM_COMBINED_pombe_minus.bw", "Cdk9as_20min_combined", NFs[33,'Spikereads']/100000))
###################################################################################################################
### Metaplotting function for comparing with corresponding heatmaps (only comparing 2 conditions)
lattice_meta.proSeq = function(filename= "test.pdf", df, ylim = c(0, 10), main = NULL, 
                               xlab = "Distance to TSS (bp)",  ylab = "Median PRO-seq intensity"){
  pdf(paste(fig_dir, filename, sep =""), width = 10, height = 5)
  result <- xyplot(mean ~ x, data = df, 
                   group = factor(treatment, labels = c("DMSO", "FP")),
                   scales = list(tck=c(1,0),alternating = c(1,1),x=list(relation='free',axs='i'),
                                 y=list(relation='free',axs='i')),
                   key=list(corner=c(0.98,0.85),padding.text=3,
                            text=list(c("DMSO", "FP"),col=c("black", "black"), cex=0.8, font=2),
                            rectangles=list(col=c(rgb(0,0,0.5), rgb(0.5,0,0)), size=2.7, height=0.8, border='white')),
                   type = 'l', 
                   ylim = ylim,
                   col = c(rgb(0,0,0.5,0.75), rgb(0.5,0,0,0.75), rgb(0,0.5,0,0.75), rgb(0.5,0,0.5,0.75)),  ## just add or take away colors to match number of samples in each panel
                   fill = c(rgb(0,0,0.5,0.3), rgb(0.5,0,0,0.3), rgb(0,0.5,0,0.3), rgb(0.5,0,0.5,0.3)),
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
                   })
  print(result)
  dev.off()
}
#####################################################################################################################
## Want to be able to generate a figure with 3 heat maps
## the first and second will be absolute counts per bin for each sample. 
##The third will be log2 fold change per bin between the two samples. 
## heatmaps will be plotted from -250 to + maxLength.
# 06-02-17  The boundaries of each gene will now be colored black. 

compareHeats = function(genes, BWset1, BWset2, maxLength = 4000, step = 10, BWpath = bwpath){
  winNum = maxLength/step
  WindowGeneEnd = 250/step + ceiling((genes[,3]- genes[,2])/step) # this defines the CPS of each gene
  BW1pl = load.bigWig(paste(BWpath, BWset1[[1]], sep = ""))
  BW1mn = load.bigWig(paste(BWpath,BWset1[[2]], sep = ""))
  BW2pl = load.bigWig(paste(BWpath,BWset2[[1]], sep = ""))
  BW2mn = load.bigWig(paste(BWpath,BWset2[[2]], sep = ""))
  Untreated = collect.many(genes, BW1pl, BW1mn, halfWindow = maxLength, step = step, at.TSS = T, do.sum = T)[,c((winNum - (250/step)+1):(2*winNum))]
  Un_normed = Untreated*(1/as.numeric(BWset1[[4]]))
  Un_normed[Un_normed == 0] <- 0.1 # permits calculation of log2 FC for all bins 
  row.names(Un_normed) <- genes[,4]
  Treated = collect.many(genes, BW2pl, BW2mn, halfWindow = maxLength, step = step, at.TSS = T, do.sum = T)[,c((winNum - (250/step)+1):(2*winNum))]
  Tr_normed = Treated*(1/as.numeric(BWset2[[4]]))
  Tr_normed[Tr_normed == 0] <- 0.1
  row.names(Tr_normed) <- genes[,4]
  FC = log2(Tr_normed/Un_normed)
  row.names(FC) <- genes[,4]
  # for each gene change the value of the cells for the start and end of each gene (change to 10000)
  for (i in 1:length(WindowGeneEnd)){
    if (WindowGeneEnd[i] < (winNum + (250/step))){
      Un_normed[i, WindowGeneEnd[i]] = 10000
      Tr_normed[i, WindowGeneEnd[i]] = 10000
      FC[i, WindowGeneEnd[i]] = 10000
    }
  }
  Un_normed[,(250/step)] = 10000
  Tr_normed[,(250/step)] = 10000
  FC[,(250/step)] = 10000
  ## also prepare corresponding metaplot data for comparisons
  Un_meta = meta.subsample(genes, BW1pl, BW1mn, step = step, at.TSS = T, halfWindow = maxLength, do.sum = T)
  Un_metaNorm = meta.normalize(result = Un_meta, scaleFactor = 1/as.numeric(BWset1[[4]]))
  Un_metaDat = cbind(seq(-240, maxLength, step), Un_metaNorm[[4]][c((winNum - (250/step)+1):(2*winNum))], 
                     Un_metaNorm[[3]][c((winNum - (250/step)+1):(2*winNum))], Un_metaNorm[[2]][c((winNum - (250/step)+1):(2*winNum))], 1, 1)
  Tr_meta = meta.subsample(genes, BW2pl, BW2mn, step = step, at.TSS = T, halfWindow = maxLength, do.sum = T)
  Tr_metaNorm = meta.normalize(result = Tr_meta, scaleFactor = 1/as.numeric(BWset2[[4]]))
  Tr_metaDat = cbind(seq(-240, maxLength, step), Tr_metaNorm[[4]][c((winNum - (250/step)+1):(2*winNum))], 
                     Tr_metaNorm[[3]][c((winNum - (250/step)+1):(2*winNum))], Tr_metaNorm[[2]][c((winNum - (250/step)+1):(2*winNum))], 1, 2)
  metaData = rbind(Un_metaDat, Tr_metaDat)
  colnames(metaData) <- c("x", "mean", "lower", "upper", "sample", "treatment")
  #result = list(as.data.frame(Un_normed), as.data.frame(Tr_normed), as.data.frame(FC))
  result = list(Un_normed, Tr_normed, FC, metaData)
  return(result)
}

#Function below will resort the 3 heat maps generated by "compareHeats" based on the geneList provided
sortHeats = function(genes, heatList){
  sorted = list()
  for (i in 1:length(heatList)){
    heat = merge(genes, heatList[[i]], by.x = "V4", by.y = 0)
    heat = heat[,c(7:dim(heat)[2])]
    row.names(heat) = genes[,4]
    sorted[[i]] = heat
  }  
  return(sorted)
}
#####################################################################################################################
# Sort gene list by diff criteria: 
## sort filtered genes by length (short to long)
orderGenesbyLen = function(geneList, StoL = F){
  len = geneList[,3]-geneList[,2]
  geneList$length = len
  if (StoL){
    result = geneList[order(geneList$length),]
  }
  else{
    result = geneList[order(-geneList$length),]
  }
  return(result)
}

orderGenesbyPI = function(geneList, CountDat = UntreatedCountDat, HighToLow = F){
  GL = merge(geneList, CountDat, by.x = "V4", by.y = 0)[,c(2,3,4,1,5,6,18)] ## col 18 is "scaled" pause index
  if (HighToLow){
  result = GL[order(GL[,7]),]
  }
  else{
  result = GL[order(-GL[,7]),]
  }
  return(result)
}

orderGenesbyActivity = function(geneList, CountDat = UntreatedCountDat, HighToLow = T){
  GL = merge(geneList, CountDat, by.x = "V4", by.y = 0)[,c(2,3,4,1,5,6,15)] ## col 15 is gene body density (i.e. "activity" )
  if (HighToLow){
    result = GL[order(GL[,7]),]
  }
  else{
    result = GL[order(-GL[,7]),]
  }
  return(result)
  }

#####################################################################################################################
# Plotting parameters #
library(RColorBrewer)
#### function for printing all three plots (untreated, treated, foldChange) on one fig
Plot3Heat = function(compareHeatRes, filename = "/test_LevPlot.jpeg", allFC = F, distToTSS = 4000, step = 10){
  Col_breaksRaw = c(seq(0, 3, length.out = 100), 10000) ## this gives 101 total values for mapping colors to
  Col_breaksFC = c(seq(-10, 10, length.out = 100), 10000) 
  col.raw <- colorRampPalette(c("white", "#EECD86",	"#E18942",	"#B95835",	"#3D3242"))(100) #set 100 values to within this color range. 
  col.raw <- c(col.raw, 'black') # add black as the last value in the color mapping list. 
  col.fc <- colorRampPalette(c('blue', 'white', 'red'))(100)
  col.fc <- c(col.fc, 'black')
  ## use distToTSS and step to determine x axis ticks and labels
  if (distToTSS >= 10000){
    tcks = c(26)
    labs = c("TSS")
    for (i in 1: (distToTSS/5000)){
      tcks = c(tcks, i*(5000/step)+26)
      labs = c(labs, sprintf("%sKb", i*5))
    }
  }
  else {
    tcks = c(0, 26*(10/step), 126*(10/step), 226*(10/step), 326*(10/step))
    labs = c("-250", "TSS", "1Kb", "2Kb", "3Kb")
  }
  #generate each heat map
  if (distToTSS > 10000){
    if (allFC){
      Unplot = levelplot(t(compareHeatRes[[1]]), ylab = "genes log2(Fold Change)", xlab = "Dist. from TSS (bp)", main = "LacZ fold change (FP/DMSO)",
                         at = Col_breaksFC, col.regions = col.fc, colorkey = list(title = "LacZ (Fold Change)", at = seq(-10, 10, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks, labels = labs)), aspect = 0.33, useRaster = T)
      Trplot = levelplot(t(compareHeatRes[[2]]), ylab = "genes log2(Fold Change)", xlab = "Dist. from TSS (bp)", main = "KD fold change (FP/DMSO)",
                         at = Col_breaksFC, col.regions = col.fc, colorkey = list(title = "KD (Fold Change)", at = seq(-10, 10, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks, labels = labs)), aspect = 0.33, useRaster = T)
      FCplot = levelplot(t(compareHeatRes[[3]]), ylab = "genes (FC KD / FC LacZ)", xlab = "Dist. from TSS (bp)", main = "fold change in fold change",
                         at = Col_breaksFC, col.regions = col.fc, colorkey = list(title = "FC_KD / FC_LacZ", at = seq(-10, 10, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks, labels = labs)), aspect = 0.33, useRaster = T)
    }
    else{
      Unplot = levelplot(t(log10(compareHeatRes[[1]])), ylab = "genes log10(reads)", xlab = "Dist. from TSS (bp)", main = "PROseq DMSO",
                         at = Col_breaksRaw, col.regions = col.raw, colorkey = list(title = "log10(reads)", at = seq(0, 3, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks, labels = labs)), aspect = 0.33, useRaster = T)
      Trplot = levelplot(t(log10(compareHeatRes[[2]])), ylab = "genes log10(reads)", xlab = "Dist. from TSS (bp)", main = "PROseq FP (300nM)",
                         at = Col_breaksRaw, col.regions = col.raw, colorkey = list(title = "log10(reads)", at = seq(0, 3, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks, labels = labs)), aspect = 0.33, useRaster = T)
      FCplot = levelplot(t(compareHeatRes[[3]]), ylab = "genes log2(fold change)", xlab = "Dist. from TSS (bp)", main = "fold change (FP/DMSO)",
                         at = Col_breaksFC, col.regions = col.fc, colorkey = list(title = "log2(FC)", at = seq(-10, 10, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks , labels = labs)), aspect = 0.33, useRaster = T)
    }
    ##plotting to jpeg
    #jpeg(file = paste(fig_dir, filename, sep = ''), width = 7, height = 10, units = "in", res = 300)
    pdf(file = paste(fig_dir, filename, sep = ''), width = 7, height = 10)
    plot.new()
    print(Unplot,  split = c(1,1,1,3), more = T)
    print(Trplot,  split = c(1,2,1,3), more = T)
    print(FCplot,  split = c(1,3,1,3), more = F)
    dev.off()
  }
  else{
    if (allFC){
      Unplot = levelplot(t(compareHeatRes[[1]]), ylab = "genes log2(Fold Change)", xlab = "Dist. from TSS (bp)", main = "WT fold change (FP/DMSO)",
                         at = Col_breaksFC, col.regions = col.fc, colorkey = list(title = "LacZ (Fold Change)", at = seq(-10, 10, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks, labels = labs)), aspect = 3, useRaster = T)
      Trplot = levelplot(t(compareHeatRes[[2]]), ylab = "genes log2(Fold Change)", xlab = "Dist. from TSS (bp)", main = "Mut fold change (Treated/Untreated)",
                         at = Col_breaksFC, col.regions = col.fc, colorkey = list(title = "KD (Fold Change)", at = seq(-10, 10, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks, labels = labs)), aspect = 3, useRaster = T)
      FCplot = levelplot(t(compareHeatRes[[3]]), ylab = "genes (FC KD / FC LacZ)", xlab = "Dist. from TSS (bp)", main = "fold change in fold change",
                         at = Col_breaksFC, col.regions = col.fc, colorkey = list(title = "FC_KD / FC_LacZ", at = seq(-10, 10, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks, labels = labs)), aspect = 3, useRaster = T)
    }
    else{
      Unplot = levelplot(t(log10(compareHeatRes[[1]])), ylab = "genes log10(reads)", xlab = "Dist. from TSS (bp)", main = "PROseq Untreated",
                         at = Col_breaksRaw, col.regions = col.raw, colorkey = list(title = "log10(reads)", at = seq(0, 3, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks, labels = labs)), aspect = 3, useRaster = T)
      Trplot = levelplot(t(log10(compareHeatRes[[2]])), ylab = "genes log10(reads)", xlab = "Dist. from TSS (bp)", main = "PROseq Treated",
                         at = Col_breaksRaw, col.regions = col.raw, colorkey = list(title = "log10(reads)", at = seq(0, 3, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks, labels = labs)), aspect = 3, useRaster = T)
      FCplot = levelplot(t(compareHeatRes[[3]]), ylab = "genes log2(fold change)", xlab = "Dist. from TSS (bp)", main = "fold change (Treated/Untreated)",
                         at = Col_breaksFC, col.regions = col.fc, colorkey = list(title = "log2(FC)", at = seq(-10, 10, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks , labels = labs)), aspect = 3, useRaster = T)
    }
    ##plotting to jpeg
    #jpeg(file = paste(fig_dir, filename, sep = ''), width = 10, height = 7, units = "in", res = 300)
    pdf(file = paste(fig_dir, filename, sep = ''), width = 7, height = 10)
    plot.new()
    print(Unplot,  split = c(1,1,3,1), more = T)
    print(Trplot,  split = c(2,1,3,1), more = T)
    print(FCplot,  split = c(3,1,3,1), more = F)
    dev.off()
  }
}

#############################################################################################################
# Heatmaps of changes within Knockdown samples (between treatments)
#############################################################################################################
## sort geneLists for running obove scripts 
lenSortGL = orderGenesbyLen(filteredGL, StoL = F) ## note because we have to transpose the matrix for plotting, sorting must be opposite what you want on plot (top to bottom)
PIsortGL = orderGenesbyPI(filteredGL, CountDat = UntreatedCountDat, HighToLow = F)
# restrict genes to specific length
filteredGL$length = filteredGL[,3] - filteredGL[,2]
long6Kbgenes = filteredGL[filteredGL$length>6000,]
ActivitySortLongGL = orderGenesbyActivity(long6Kbgenes[,c(1:6)], HighToLow = T)
FinalActivitySortLongGL = orderGenesbyActivity(CountDat = Min20CountDat, long6Kbgenes[,c(1:6)], HighToLow = T)
# run heatmatrix generating script separately for each mutant (treated vs untreated)
WT_data = compareHeats(lenSortGL, wigset[1,], wigset[2,])
Cdk9as_data = compareHeats(lenSortGL, wigset[5,], wigset[6,])
Mcs6as_data = compareHeats(lenSortGL, wigset[3,], wigset[4,])
Lsk1as_data = compareHeats(lenSortGL, wigset[7,], wigset[8,])
Mcs6as_Cdk9as_data = compareHeats(lenSortGL, wigset[9,], wigset[10,])
Lsk1as_Cdk9as_data = compareHeats(lenSortGL, wigset[21,], wigset[22,])
# run plotting function above for each AS mutant comparison
dir.create(paste(fig_dir,"/lengthSort/", sep = ""))
Plot3Heat(WT_data, filename = "/lengthSort/WT_Heatmaps_LenSort_raster.pdf")
Plot3Heat(Cdk9as_data, filename = "/lengthSort/Cdk9as_Heatmaps_LenSort_raster.pdf")
Plot3Heat(Mcs6as_data, filename = "/lengthSort/Mcs6as_Heatmaps_LenSort_raster.pdf")
Plot3Heat(Lsk1as_data, filename = "/lengthSort/Lsk1as_Heatmaps_LenSort_raster.pdf")
Plot3Heat(Lsk1as_Cdk9as_data, filename = "/lengthSort/Lsk1as_Cdk9as_Heatmaps_LenSort_raster.pdf")

# run heatmatrix generating script separately for time point sample (against the 0') in the Cdk9as inhibition time course with length sorted genes
halfMin_data_len = compareHeats(lenSortGL, TimeCourse_wigset[1,], TimeCourse_wigset[2,], step = 10, maxLength = 4000)
oneMin_data_len = compareHeats(lenSortGL, TimeCourse_wigset[1,], TimeCourse_wigset[3,], step = 10, maxLength = 4000)
two.5Min_data_len = compareHeats(lenSortGL, TimeCourse_wigset[1,], TimeCourse_wigset[4,], step = 10, maxLength = 4000)
fiveMin_data_len = compareHeats(lenSortGL, TimeCourse_wigset[1,], TimeCourse_wigset[5,], step = 10,maxLength = 4000)
seven.5Min_data_len = compareHeats(lenSortGL, TimeCourse_wigset[1,], TimeCourse_wigset[6,], step = 10, maxLength = 4000)
tenMin_data_len = compareHeats(lenSortGL, TimeCourse_wigset[1,], TimeCourse_wigset[7,], step = 10, maxLength = 4000)
twentyMin_data_len = compareHeats(lenSortGL, TimeCourse_wigset[1,], TimeCourse_wigset[8,], step = 10, maxLength = 4000)
# run plotting function above for each time point comparison
Plot3Heat(halfMin_data_len, filename = "/lengthSort/ePROseq_Cdk9as_TC_0vs30sec_Heatmaps_LenSort_raster.pdf")
Plot3Heat(oneMin_data_len, filename = "/lengthSort/ePROseq_Cdk9as_TC_0vs1min_Heatmaps_LenSort_raster.pdf")
Plot3Heat(two.5Min_data_len, filename = "/lengthSort/ePROseq_Cdk9as_TC_0vs2.5min_Heatmaps_LenSort_raster.pdf")
Plot3Heat(fiveMin_data_len, filename = "/lengthSort/ePROseq_Cdk9as_TC_0vs5min_Heatmaps_LenSort_raster.pdf")
Plot3Heat(seven.5Min_data_len, filename = "/lengthSort/ePROseq_Cdk9as_TC_0vs7.5min_Heatmaps_LenSort_raster.pdf")
Plot3Heat(tenMin_data_len, filename = "/lengthSort/ePROseq_Cdk9as_TC_0vs10min_Heatmaps_LenSort_raster.pdf")
Plot3Heat(twentyMin_data_len, filename = "/lengthSort/ePROseq_Cdk9as_TC_0vs20min_Heatmaps_LenSort_raster.pdf")

# run heatmatrix generating script separately for time point sample (against the 0') in the Cdk9as inhibition time course 6kb genes sorted by activity
halfMin_data = compareHeats(FinalActivitySortLongGL, TimeCourse_wigset[1,], TimeCourse_wigset[2,], step = 250, maxLength = 6000)
oneMin_data = compareHeats(FinalActivitySortLongGL, TimeCourse_wigset[1,], TimeCourse_wigset[3,], step = 250, maxLength = 6000)
two.5Min_data = compareHeats(FinalActivitySortLongGL, TimeCourse_wigset[1,], TimeCourse_wigset[4,], step = 250, maxLength = 6000)
fiveMin_data = compareHeats(FinalActivitySortLongGL, TimeCourse_wigset[1,], TimeCourse_wigset[5,], step = 250,maxLength = 6000)
seven.5Min_data = compareHeats(FinalActivitySortLongGL, TimeCourse_wigset[1,], TimeCourse_wigset[6,], step = 250, maxLength = 6000)
tenMin_data = compareHeats(FinalActivitySortLongGL, TimeCourse_wigset[1,], TimeCourse_wigset[7,], step = 250, maxLength = 6000)
twentyMin_data = compareHeats(FinalActivitySortLongGL, TimeCourse_wigset[1,], TimeCourse_wigset[8,], step = 250, maxLength = 6000)
# run plotting function above for each time point comparison
dir.create(paste(fig_dir,"/activitySort/", sep = ""))
Plot3Heat(halfMin_data, filename = "/activitySort/ePROseq_Cdk9as_TC_0vs30sec_Heatmaps_20minActivitySort_6kb_genes_smooth250bp_raster.pdf")
Plot3Heat(oneMin_data, filename = "/activitySort/ePROseq_Cdk9as_TC_0vs1min_Heatmaps_20minActivitySort_6kb_genes_smooth250bp_raster.pdf")
Plot3Heat(two.5Min_data, filename = "/activitySort/ePROseq_Cdk9as_TC_0vs2.5min_Heatmaps_20minActivitySort_6kb_genes_smooth250bp_raster.pdf")
Plot3Heat(fiveMin_data, filename = "/activitySort/ePROseq_Cdk9as_TC_0vs5min_Heatmaps_20minActivitySort_6kb_genes_smooth250bp_raster.pdf")
Plot3Heat(seven.5Min_data, filename = "/activitySort/ePROseq_Cdk9as_TC_0vs7.5min_Heatmaps_20minActivitySort_6kb_genes_smooth250bp_raster.pdf")
Plot3Heat(tenMin_data, filename = "/activitySort/ePROseq_Cdk9as_TC_0vs10min_Heatmaps_20minActivitySort_6kb_genes_smooth250bp_raster.pdf")
Plot3Heat(twentyMin_data, filename = "/activitySort/ePROseq_Cdk9as_TC_0vs20min_Heatmaps_20minActivitySort_6kb_genes_smooth250bp_raster.pdf")
# run heatmatrix generating script to compare Cdk9as Dis2ts strain at 30C and 18C etc. (treated vs untreated)
Dis2ts_DMSO_data = compareHeats(lenSortGL, wigset[15,], wigset[13,])
Dis2ts_3MB_data = compareHeats(lenSortGL, wigset[16,], wigset[14,])
Cdk9as_Dis2ts_DMSO_data = compareHeats(lenSortGL, wigset[19,], wigset[17,])
Cdk9as_Dis2ts_3MB_data = compareHeats(lenSortGL, wigset[20,], wigset[18,])
# run plotting function for each above comparison
Plot3Heat(Dis2ts_DMSO_data, filename = "/lengthSort/Dis2ts_DMSO_18v30C_Heatmaps_LenSort_raster.pdf")
Plot3Heat(Dis2ts_3MB_data, filename = "/lengthSort/Dis2ts_3MBPP1_18v30C_Heatmaps_LenSort_raster.pdf")
Plot3Heat(Cdk9as_Dis2ts_DMSO_data, filename = "/lengthSort/Cdk9as_Dis2ts_DMSO_18v30C_Heatmaps_LenSort_raster.pdf")
Plot3Heat(Cdk9as_Dis2ts_3MB_data, filename = "/lengthSort/Cdk9as_Dis2ts_3MBPP1_18v30C_Heatmaps_LenSort_raster.pdf")

