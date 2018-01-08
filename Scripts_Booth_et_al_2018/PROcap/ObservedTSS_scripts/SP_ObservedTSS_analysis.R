source("/Volumes/SEAGATE_EXP/Greg_YeastData/YeastAnalysisScripts/GB_Functions.R")
library(bigWig)
library(lattice)
library(grid)
dir.create("/Volumes/SEAGATE_EXP/Greg_YeastData/PRO_Cap/Pombe/Analysis/observedTSS/figures/SP_CapAnalysis_05-19-16/")
directory = "/Volumes/SEAGATE_EXP/Greg_YeastData/PRO_Cap/Pombe/Analysis/observedTSS/figures/SP_CapAnalysis_05-19-16/"

fullGL=read.table("/Volumes/SEAGATE_EXP/Greg_YeastData/GeneLists_final/Pombe_08_23-15/pombe.ASM294v1.16.cleaned_sorted.bed",head=F)
SPGL_1kbsep = read.table("/Volumes/SEAGATE_EXP/Greg_YeastData/GeneLists_final/Pombe_08_23-15/observedTSS/SP_PROcapObservedTSS_1000bpSep.bed", stringsAsFactors = F)
ObsTSS_fullData = read.table(file = "/Volumes/SEAGATE_EXP/Greg_YeastData/GeneLists_final/Pombe_05-19-16/pombe.ASM294v1.16.cleaned.bed_PROcapObservedTSS.bed", sep = "\t")
ObsTSS_Filtered = read.table(file = "/Volumes/SEAGATE_EXP/Greg_YeastData/GeneLists_final/Pombe_05-19-16/SP_CompleteFilteredGenes.bed", sep = "\t") # genes are active and not influenced by upstream transcription. 

ObsTSS_data_filtered = ObsTSS_fullData[ObsTSS_fullData[,4] %in% ObsTSS_Filtered[,4] ,]
MatchedFullgList = merge(fullGL, ObsTSS_Filtered, by.x = "V4", by.y = "V4")[,c(2,3,4,1,5,6)]

bed = ObsTSS_data_filtered
cat("number of genes being used =", length(bed[,1]), "\n")
#write.table(MatchedFullgList[, c(1,2,3,4,5,6)], "/Volumes/SEAGATE_EXP/Greg_YeastData/GeneLists_final/Pombe_08_23-15/SP_CompleteFilteredGenes_AnnotatedTSS.bed", quote = F, sep = "\t", row.names =F, col.names =F)
################################################################################################################
### distributions of observed TSS features
get_TSSdif = function(ObservedGL, AnnotatedGL){
  N = dim(ObservedGL)[1]
  result = vector(mode="integer", length=N)
  for (i in 1:N){
    if (ObservedGL[i,6] == "+"){
      result[i] = ObservedGL[i,2] - AnnotatedGL[i,2]
    }
    else{
      result[i] = AnnotatedGL[i,3] - ObservedGL[i,3]
    }
  }
  return(cbind(result))
}

pdf(file = paste(directory, "ObservedTSS_countHist.pdf", sep = ""), width = 5.5, height = 5.5)
hist(log10(ObsTSS_data_filtered[,5]), breaks = 50, xlab = "log10(counts(RPM) - background)", main = "Obs. TSS Counts", cex.lab = 1.5, cex.main = 2, cex.axis = 1.5, col = "blue4")
dev.off()

Dists = get_TSSdif(ObsTSS_Filtered, MatchedFullgList)
pdf(file = paste(directory,"ObservedTSS_DistanceHist.pdf", sep = ""), width = 5.5, height = 5.5)
hist(Dists, breaks = 50, xlab = "distance bp", main = "Distance From Annotation", cex.lab = 1.5, cex.main = 2, cex.axis = 1.5, col = rgb(0.25, 0, 1, 0.5), xlim = c(-250,250))
dev.off()

#pdf(file = "/Volumes/SEAGATE_EXP/Greg_YeastData/PRO_Cap/Pombe/Analysis/observedTSS/figures/ObservedTSS_DispersionHist.pdf", width = 7, height = 7)
#hist(ObservedTSS$TSSdispersion, breaks = 50, xlab = "dispersion (stDev)", main = "Observed Signal Dispersion", cex.lab = 1.5, cex.main = 2, cex.axis = 1.5, col = "orange4")
#dev.off()

################################################################################################################
## SP_PROcap +TAP data 

SPcapPbw = load.bigWig("/Volumes/SEAGATE_EXP/Greg_YeastData/PRO_Cap/Pombe/bigWig/3870_7157_12385_C53ARACXX_pombe972h-ProCap_ATCACG_R1.fastq_plus.bw")
SPcapMbw = load.bigWig("/Volumes/SEAGATE_EXP/Greg_YeastData/PRO_Cap/Pombe/bigWig/3870_7157_12385_C53ARACXX_pombe972h-ProCap_ATCACG_R1.fastq_minus.bw")
SPcap_RPM = 1/(total.reads(SPcapPbw,SPcapMbw)/1000000)

#scaledbins = collect.many.scaled(SPgList, Pbw, Mbw, 1000, 300, 10, do.sum = T)
ProCapplus = meta.subsample(ObsTSS_Filtered, SPcapPbw, SPcapMbw, 50, 1, do.sum = T, at.TSS=T)
ProCapminus = meta.subsample(ObsTSS_Filtered, SPcapMbw, SPcapPbw, 50, 1, do.sum = T, at.TSS=T)
ProCapplusRPM = meta.normalize(ProCapplus, SPcap_RPM)
ProCapminusRPM = meta.normalize(ProCapminus, SPcap_RPM)

pdf(file = paste(directory, "metaplot_ObservedTSS.pdf", sep = ""), width = 5, height = 5)
meta.plot.GROseq(ProCapplusRPM, ProCapminusRPM, 1, main="SP PROcap Observed TSS",  
                 ylab = "Avg. read count (RPM)", cex.main=2, cex.lab=1.5, cex.axis=1.5, 
                 bothStrands = F, col1 = 'sienna4')
abline(v=c(0), lwd=2, lty=2)
dev.off()

#Annotated TSS
ProCapplusAnn = meta.subsample(MatchedFullgList, SPcapPbw, SPcapMbw, 50, 1, do.sum = T, at.TSS=T)
ProCapminusAnn = meta.subsample(MatchedFullgList, SPcapMbw, SPcapPbw, 50, 1, do.sum = T, at.TSS=T)
ProCapplusRPMAnn = meta.normalize(ProCapplusAnn, SPcap_RPM)
ProCapminusRPMAnn = meta.normalize(ProCapminusAnn, SPcap_RPM)

pdf(file = paste(directory, "metaplot_AnnotatedTSS.pdf", sep = ""), width = 5, height = 5)
meta.plot.GROseq(ProCapplusRPMAnn, ProCapminusRPMAnn, 1, main="SP PROcap Annotated TSS",  
                 ylab = "Avg. read count (RPM)", cex.main=2, cex.lab=1.5, cex.axis=1.5, 
                 bothStrands = F, col1 = 'sienna4')
abline(v=c(0), lwd=2, lty=2)
dev.off()

## SP_PROcap -TAP data 

SPcapPbw_mTAP = load.bigWig("/Volumes/SEAGATE_EXP/Greg_YeastData/PRO_Cap/Pombe/bigWig/3870_7157_12386_C53ARACXX_pombe972h-ProCap-TAP_CGATGT_R1.fastq_plus.bw")
SPcapMbw_mTAP = load.bigWig("/Volumes/SEAGATE_EXP/Greg_YeastData/PRO_Cap/Pombe/bigWig/3870_7157_12386_C53ARACXX_pombe972h-ProCap-TAP_CGATGT_R1.fastq_minus.bw")
SPcap_RPM_mTAP = 1/(total.reads(SPcapPbw_mTAP,SPcapMbw_mTAP)/1000000)

#scaledbins = collect.many.scaled(SPgList, Pbw, Mbw, 1000, 300, 10, do.sum = T)
ProCapplus_mTAP = meta.subsample(ObsTSS_Filtered, SPcapPbw_mTAP, SPcapMbw_mTAP, 50, 1, do.sum = T, at.TSS=T)
ProCapminus_mTAP = meta.subsample(ObsTSS_Filtered, SPcapMbw_mTAP, SPcapPbw_mTAP, 50, 1, do.sum = T, at.TSS=T)
ProCapplusRPM_mTAP = meta.normalize(ProCapplus_mTAP, SPcap_RPM_mTAP)
ProCapminusRPM_mTAP = meta.normalize(ProCapminus_mTAP, SPcap_RPM_mTAP)

pdf(file = paste(directory,"metaplot_kb_ObservedTSS_VsBackground.pdf", sep = ""), width = 5, height = 5)
meta.plot.GROseq(ProCapplusRPM, ProCapminusRPM, 1, main="SP PROcap +/- TAP Obs. TSS",  
                 ylab = "Avg. read count (RPM)", cex.main=2, cex.lab=1.5, cex.axis=1.5, 
                 bothStrands = F, col1 = 'sienna4')
meta.overlay(ProCapplusRPM_mTAP, step = 1, strand = "+", col1 = 'blue4', col2 = rgb(0,0,1,0.2))
#meta.overlay(ProCapminusRPM_mTAP, step = 1, strand = "-")
abline(v=c(0), lwd=2, lty=2)
legend('topleft', c("+TAP", "-TAP"), text.col = c('sienna4', 'blue4'), cex = 2)
dev.off()

#Annotated TSS
ProCapplus_mTAP_Ann = meta.subsample(MatchedFullgList, SPcapPbw_mTAP, SPcapMbw_mTAP, 50, 1, do.sum = T, at.TSS=T)
ProCapminus_mTAP_Ann = meta.subsample(MatchedFullgList, SPcapMbw_mTAP, SPcapPbw_mTAP, 50, 1, do.sum = T, at.TSS=T)
ProCapplusRPM_mTAP_Ann = meta.normalize(ProCapplus_mTAP_Ann, SPcap_RPM_mTAP)
ProCapminusRPM_mTAP_Ann = meta.normalize(ProCapminus_mTAP_Ann, SPcap_RPM_mTAP)

pdf(file = paste(directory,"metaplot_AnnotatedTSS_VsBackground.pdf", sep = ""), width = 5, height = 5)
meta.plot.GROseq(ProCapplusRPMAnn, ProCapminusRPMAnn, 1, main="SP PROcap +/- TAP Ann. TSS",  
                 ylab = "Avg. read count (RPM)", cex.main=2, cex.lab=1.5, cex.axis=1.5, 
                 bothStrands = F, col1 = 'sienna4')
meta.overlay(ProCapplusRPM_mTAP_Ann, step = 1, strand = "+", col1 = 'blue4', col2 = rgb(0,0,1,0.2))
#meta.overlay(ProCapminusRPM_mTAP_Ann, step = 1, strand = "-")
abline(v=c(0), lwd=2, lty=2)
legend('topleft', c("+TAP", "-TAP"), text.col = c('sienna4', 'blue4'), cex = 2)
dev.off()


################################################################################################################
#Annotated TSS heatmap (Use background subtracted data)
SPcapPbw_mBackground = load.bigWig("/Volumes/SEAGATE_EXP/Greg_YeastData/PRO_Cap/Pombe/bigWig/pombe972h-ProCap_ATCACG_R1_normed_plus_BackgroundSubtracted.bw")
SPcapMbw_mBackground = load.bigWig("/Volumes/SEAGATE_EXP/Greg_YeastData/PRO_Cap/Pombe/bigWig/pombe972h-ProCap_ATCACG_R1_normed_minus_BackgroundSubtracted.bw")

SPWTheatDataAnn = collect.many(MatchedFullgList, SPcapPbw_mBackground, SPcapMbw_mBackground, 250, 1, at.TSS=T, do.sum=T)
#sort by centerbin
#SPWTheatDataSortedAnn = SPWTheatDataAnn[order(SPWTheatDataAnn[,(ncol(SPWTheatDataAnn)/2)], decreasing = T),]
# Sort by distance from Annotation to Obs
AnnToObsDist = get_TSSdif(ObservedGL = ObsTSS_data_filtered, AnnotatedGL = MatchedFullgList)
SPWTheatDataAnn = cbind(SPWTheatDataAnn, AnnToObsDist)
SPWTheatDataSortedAnn = SPWTheatDataAnn[order(SPWTheatDataAnn[,501], decreasing = F),]
library(gplots)
#my_palette = colorRampPalette(c('lightblue4','lightblue3','lightblue1','white','pink','salmon1','tomato4'))(299)
my_palette = colorRampPalette(c('white','lightblue1','lightblue2','lightblue3', 'lightblue4', 'midnightblue'))(299)
br = c(seq(0,0.09,length=50), seq(0.1,0.19,length=50), seq(0.2,0.49,length=50), seq(0.5,0.99,length=50), seq(1,4.99,length=50), seq(5,10,length=50))
br1 = c(seq(0,1,length=50), seq(2,3,length=50), seq(4,5,length=50), seq(6,12,length=50), seq(13,18,length=50), seq(19,20,length=50))

jpeg(file = paste(directory,"heatmap_AnnDistToObsSort_250bp_AnnotatedTSS.jpeg", sep = ""), units = "in", width = 10, height = 10, res = 300)
#pdf(file = paste(directory,"heatmap_AnnDistToObsSort_250bp_AnnotatedTSS.pdf", sep = ""), width = 10, height = 10)
heatmap.2(as.matrix(SPWTheatDataSortedAnn[,1:500]), Rowv = NA, Colv = NA,  dendrogram=c("none"), 
          breaks = br, labRow = F, labCol= F, col =  my_palette, symkey = F, symm=F, 
          trace = c("none"), key.title = "", key.xlab = "", key.ylab = "", 
          density.info = "none", useRaster = T)
dev.off()
## make heatmap centered on observed TSS
#indx1 <- (bed[,8] != 'NaN')
SPWTheatData = collect.many(ObsTSS_Filtered, SPcapPbw_mBackground, SPcapMbw_mBackground, 250, 1, at.TSS=T, do.sum=T)
#sort by centerbin
#SPWTheatDataSorted = SPWTheatData[order(SPWTheatData[,(ncol(SPWTheatData)/2)], decreasing = T),]
#sort by increasing dispersion
#SPWTheatData1 = cbind(SPWTheatData, bed[indx1,]$TSSdispersion)
#SPWTheatDataSorted_disp = SPWTheatData1[order(SPWTheatData1[,201], decreasing = F), 1:200)]
#sort by distance from Annotation to Obs
SPWTheatData = cbind(SPWTheatData, AnnToObsDist)
SPWTheatDataSorted = SPWTheatData[order(SPWTheatData[,501], decreasing = F),]
jpeg(file = paste(directory,"heatmap_AnnDistToObsSort_250bp_ObservedTSS.jpeg", sep = ""), units = "in", width = 10, height = 10, res = 300)
#pdf(file = paste(directory, "heatmap_AnnDistToObsSort_250bp_ObservedTSS.pdf", sep = ""), width = 10, height = 10)
heatmap.2(as.matrix(SPWTheatDataSorted[,1:500]), Rowv = NA, Colv = NA,  dendrogram=c("none"), 
          breaks = br, labRow = F, labCol= F, col =  my_palette, symkey = F, symm=F, trace = c("none"), 
          key.title = "", key.xlab = "", key.ylab = "", density.info = "none", useRaster = T)
dev.off()

################################################################################################################
# S. pombe WT combined PRO-seq replicates:

SPPbw = load.bigWig("/Volumes/SEAGATE_EXP/Greg_YeastData/PRO_Seq/Pombe/bigWig/3870_7157_12395_C53ARACXX_SP_WT_972h-COMBINED_REPS_GGCTAC_R1.fastq_sorted.bed_plus.bw")
SPMbw = load.bigWig("/Volumes/SEAGATE_EXP/Greg_YeastData/PRO_Seq/Pombe/bigWig/3870_7157_12395_C53ARACXX_SP_WT_972h-COMBINED_REPS_GGCTAC_R1.fastq_sorted.bed_minus.bw")
SP_WT_RPM = 1/(total.reads(SPPbw,SPMbw)/1000000)

#scaledbins = collect.many.scaled(SPgList, Pbw, Mbw, 1000, 300, 10, do.sum = T)
WTplus = meta.subsample(ObsTSS_Filtered, SPPbw, SPMbw, 500, 10, do.sum = T, at.TSS=T)
WTminus = meta.subsample(ObsTSS_Filtered, SPMbw, SPPbw, 500, 10, do.sum = T, at.TSS=T)
WTplusRPM = meta.normalize(WTplus, SP_WT_RPM)
WTminusRPM = meta.normalize(WTminus, SP_WT_RPM)

pdf(file = paste(directory, "PROseq_ObservedTSS_metaplot.pdf", sep = ""), width = 5, height = 5)
meta.plot.GROseq(WTplus, WTminus, 10, main="SP PROseq Observed TSS",  cex.main=2, cex.lab=1.5, cex.axis=1.5, bothStrands = F)
abline(v=0, lwd=2, lty=2)
dev.off()

### overlay PROcap and PROseq metaplots
pdf(file = paste(directory,"PROcapVsSeq_kb_ObservedTSS_meta.pdf", sep = ""), width = 10, height = 10)
meta.plot.GROseq(ProCapplusRPM, ProCapminusRPM, 10, main="SP PROseqVsCap Observed TSS",  ylab = "Avg. read count (RPM)", cex.main=2, cex.lab=1.5, cex.axis=1.5, bothStrands = F)
meta.overlay(WTplusRPM, step = 10, strand = "+")
abline(v=c(0), lwd=2, lty=2)
legend('topleft', c("PRO-Cap", "PRO-seq"), text.col= c("red", "green"), cex = 1.5)
dev.off()

#scaledbins = collect.many.scaled(SPgList, Pbw, Mbw, 1000, 300, 10, do.sum = T)
WTplusAnn = meta.subsample(MatchedFullgList, SPPbw, SPMbw, 500, 10, do.sum = T, at.TSS=T)
WTminusAnn = meta.subsample(MatchedFullgList, SPMbw, SPPbw, 500, 10, do.sum = T, at.TSS=T)
WTplusRPMAnn = meta.normalize(WTplusAnn, SP_WT_RPM)
WTminusRPMAnn = meta.normalize(WTminusAnn, SP_WT_RPM)

pdf(file = paste(directory, "PROseq_AnnoatedTSS_metaplot.pdf", sep = ""), width = 7, height = 7)
meta.plot.GROseq(WTplusAnn, WTminusAnn, 10, main="S.pombe PROseq Annotated TSS",  cex.main=2, cex.lab=1.5, cex.axis=1.5, bothStrands = F)
abline(v=0, lwd=2, lty=2)
dev.off()

#Annotated TSS
pdf(file = paste(directory,"PROcapVsSeq_AnnotatedTSS_meta.pdf", sep = ""), width = 10, height = 10)
meta.plot.GROseq(WTplusRPMAnn, WTminusRPMAnn, 10, main="SP PROseqVsCap Annotated TSS",  ylab = "Avg. read count (RPM)", col1 = 'green', cex.main=2, cex.lab=1.5, cex.axis=1.5, bothStrands = F)
meta.overlay(ProCapplusRPMAnn, step = 10, strand = "+", col1 = 'red', col2 = rgb(1,0,0,0.25))
abline(v=c(0), lwd=2, lty=2)
legend('topleft', c("PRO-Cap", "PRO-seq"), text.col= c("red", "green"), cex = 1.5)
dev.off()

### Calculate fraction of Maximum (See Fig. 1; Core, Martins, et. al. Nat. Gen 2014)
## overlay fraction of maximum meta plots for pro seq and pro-cap.  This should normalize them to improve visualization.
WTseqplusObs = meta.subsample(bed, SPPbw, SPMbw, 300, 5, do.sum = T, at.TSS=T)
WTcapplusObs = meta.subsample(bed, SPcapPbw, SPcapMbw, 300, 5, do.sum = T, at.TSS=T)
WTseqplusAnn = meta.subsample(MatchedFullgList, SPPbw, SPMbw, 300, 5, do.sum = T, at.TSS=T)
WTcapplusAnn = meta.subsample(MatchedFullgList, SPcapPbw, SPcapMbw, 300, 5, do.sum = T, at.TSS=T)
WT_PROseq_FracMax_Obs_pl = WTseqplusObs[[4]]/(max(WTseqplusObs[[4]]))
WT_PROcap_FracMax_Obs_pl = WTcapplusObs[[4]]/(max(WTcapplusObs[[4]]))
WT_PROseq_FracMax_Ann_pl = WTseqplusAnn[[4]]/(max(WTseqplusAnn[[4]]))
WT_PROcap_FracMax_Ann_pl = WTcapplusAnn[[4]]/(max(WTcapplusAnn[[4]]))
base = seq(-295, 300,5) 
FracMaxTable = cbind(base, WT_PROseq_FracMax_Obs_pl, WT_PROcap_FracMax_Obs_pl, WT_PROseq_FracMax_Ann_pl, WT_PROcap_FracMax_Ann_pl)
# FracMax around Observed TSS
pdf(file = paste(directory, "SP_PROcapVsSeq_300bp_ObservedTSS_FractionOfMaximum_meta.pdf", sep = ""), width = 5, height = 5)
plot(FracMaxTable[,1], FracMaxTable[,3], type = 'l', lwd = 4, col = 'sienna4',
     main="SP Observed TSS",  ylab = "Fraction of Maximum", xlab = "Distance to TSS",
     cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(FracMaxTable[,1], FracMaxTable[,2], type = 'l', lwd = 4, col = 'blue4')
abline(v=c(0), lwd=2, lty=2)
legend('topleft', c("PRO-Cap", "PRO-seq"), text.col= c("sienna4", "blue4"), cex = 1.5)
dev.off()

# FracMax around annotated TSS
pdf(file = paste(directory, "SP_PROcapVsSeq_300bp_AnnotatedTSS_FractionOfMaximum_meta.pdf", sep = ""), width = 5, height = 5)
plot(FracMaxTable[,1], FracMaxTable[,5], type = 'l', lwd = 4, col = 'sienna4',
     main="SP Annotated TSS",  ylab = "Fraction of Maximum", xlab = "Distance to TSS",
     cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(FracMaxTable[,1], FracMaxTable[,4], type = 'l', lwd = 4, col = 'blue4')
abline(v=c(0), lwd=2, lty=2)
legend('topleft', c("PRO-Cap", "PRO-seq"), text.col= c("sienna4", "blue4"), cex = 1.5)
dev.off()

################################################################################################################

## SP_MNase-seq data 

MnPbw = load.bigWig("/Volumes/SEAGATE_EXP/GEO_datasets/S.pombeMnase_GSE49575_RAW/GSM1201970_MNase_WT_150Uml_N1_90_200bp_gcnorm19_KDE_bandw20sampled10_chrNumeric.bw")
MnMbw = load.bigWig("/Volumes/SEAGATE_EXP/GEO_datasets/S.pombeMnase_GSE49575_RAW/GSM1201970_MNase_WT_150Uml_N1_90_200bp_gcnorm19_KDE_bandw20sampled10_chrNumeric.bw")
SP_MN_scale = 30

#scaledbins = collect.many.scaled(SPgList, Pbw, Mbw, 1000, 300, 10, do.sum = T)
Mnplus = meta.subsample(ObsTSS_Filtered, MnPbw, MnMbw, 1000, 10, do.sum = T, at.TSS=T)
Mnminus = meta.subsample(ObsTSS_Filtered, MnMbw, MnPbw, 1000, 10, do.sum = T, at.TSS=T)
MnplusScaled = meta.normalize(Mnplus, SP_MN_scale)
MnminusScaled = meta.normalize(Mnminus, SP_MN_scale)
MnplusAnn = meta.subsample(MatchedFullgList, MnPbw, MnMbw, 1000, 10, do.sum = T, at.TSS=T)
MnminusAnn = meta.subsample(MatchedFullgList, MnMbw, MnPbw, 1000, 10, do.sum = T, at.TSS=T)
MnplusScaledAnn = meta.normalize(MnplusAnn, SP_MN_scale)
MnminusScaledAnn = meta.normalize(MnminusAnn, SP_MN_scale)

pdf(file = paste(directory,"MNase_kb_ObservedvsAnnotatedTSS_meta.pdf", sep = "") , width = 5, height = 5)
meta.plot.GROseq(Mnplus, Mnminus, 10, main="SP Mnase around TSS", cex.main=2, cex.lab=1.5, 
                 cex.axis=1.5, bothStrands = F, col1 = "sienna4")
meta.overlay(MnplusAnn, step = 10, strand = "+", col1 = 'blue4', col2 = rgb(0,0,1,0.2))
legend("topleft", c("Annotated TSS", "Observed TSS"), text.col = c("sienna4", "blue4"), cex = 1.5)
abline(v=0, lwd=2, lty=2)
dev.off()

################################################################################################################
################################################################################################################
# Check Correlations between PRO-seq GB density and PRO-cap signal (proving that PRO-cap is a good filter for "gene activity")

WT_PROseqCountData = read.table(file = "/Volumes/SEAGATE_EXP/Greg_YeastData/YeastAnalysisScripts/SC_WT_Spt4/ObservedTSS_based/countData/SC_WTcombined_CountData.txt", head = T, stringsAsFactors = F)
cor.test(ObsTSS_data_filtered[,5], WT_PROseqCountData$gb_densityRPKM, method = "spearman")

plot(log10(ObsTSS_data_filtered[,5]), log10(WT_PROseqCountData$gb_densityRPKM))


################################################################################################################

## Split up Genes based on "focused" (1st quartile) or "dispersed" (4th quartile) signal around observed TSS. 
indx <- (ObservedTSS[,8] != "NaN")
SummaryDisp = summary(ObservedTSS[indx,]$TSSdispersion)
focusedTSS = ObservedTSS[indx,][ObservedTSS[indx,]$TSSdispersion < SummaryDisp[2],]
dispersedTSS = ObservedTSS[indx,][ObservedTSS[indx,]$TSSdispersion > SummaryDisp[5],]

boxplot(log10(focusedTSS$X50mer_sumCounts), log10(dispersedTSS$X50mer_sumCounts), col = c("red2", "blue2"), names = c('focused', "dispersed"), ylab = "log10(TSS counts)", cex.lab = 1.5, cex.axis = 1.5)
### Need to match groups based on total reads in 50mers around TSS
summaryFoc = summary(focusedTSS$X50mer_sumCounts)
summaryDisp = summary(dispersedTSS$X50mer_sumCounts)

MatchedFocTSS = focusedTSS[focusedTSS$X50mer_sumCounts > 131 & focusedTSS$X50mer_sumCounts < 318, ]
MatchedDispTSS = dispersedTSS[dispersedTSS$X50mer_sumCounts > 131 & dispersedTSS$X50mer_sumCounts < 318, ]


### metaplot these two groups. 
focusedPlus = meta.subsample(MatchedFocTSS, SPcapPbw, SPcapMbw, 100, 5, do.sum = T, at.TSS=T)
focusedMinus = meta.subsample(MatchedFocTSS, SPcapMbw, SPcapPbw, 100, 5, do.sum = T, at.TSS=T)
dispersedPlus = meta.subsample(MatchedDispTSS, SPcapPbw, SPcapMbw, 100, 5, do.sum = T, at.TSS=T)
dispersedMinus = meta.subsample(MatchedDispTSS, SPcapMbw, SPcapPbw, 100, 5, do.sum = T, at.TSS=T)

pdf(file = "/Volumes/SEAGATE_EXP/Greg_YeastData/PRO_Cap/Pombe/Analysis/observedTSS/figures/focusedVsDispersed_observedTSS_Meta.pdf", width = 7, height = 7)
meta.plot.GROseq(focusedPlus, focusedMinus, 5, main="SP PROcap Focused vs Disp. Obs. TSS",  ylab = "Avg. read count", cex.main=2, cex.lab=1.5, cex.axis=1.5, bothStrands = T)
meta.overlay(dispersedPlus, step = 5, strand = "+")
meta.overlay(dispersedMinus, step = 5, strand = "-")
abline(v=c(0), lwd=2, lty=2)
legend('topleft', c("focused N = 1139", "dispersed N = 1139"), text.col= c("red", "green"), cex = 1.5)
dev.off()

#### metaplot of other Data around these groups 

PROseqfocusedPlus = meta.subsample(focusedTSS, SPPbw, SPMbw, 250, 5, do.sum = T, at.TSS=T)
PROseqfocusedMinus = meta.subsample(focusedTSS, SPMbw, SPPbw, 250, 5, do.sum = T, at.TSS=T)
PROseqdispersedPlus = meta.subsample(dispersedTSS, SPPbw, SPMbw, 250, 5, do.sum = T, at.TSS=T)
PROseqdispersedMinus = meta.subsample(dispersedTSS, SPMbw, SPPbw, 250, 5, do.sum = T, at.TSS=T)

pdf(file = "/Volumes/SEAGATE_EXP/Greg_YeastData/PRO_Cap/Pombe/Analysis/observedTSS/figures/PROseqfocusedVsDispersed_observedTSS_Meta.pdf", width = 7, height = 7)
meta.plot.GROseq(PROseqfocusedPlus, PROseqfocusedMinus, 5, main="SP PROseq Focused vs Disp. Obs. TSS",  ylab = "Avg. read count", cex.main=2, cex.lab=1.5, cex.axis=1.5, bothStrands = T)
meta.overlay(PROseqdispersedPlus, step = 5, strand = "+")
meta.overlay(PROseqdispersedMinus, step = 5, strand = "-")
abline(v=c(0), lwd=2, lty=2)
legend('topleft', c("focused N = 1139", "dispersed N = 1139"), text.col= c("red", "green"), cex = 1.5)
dev.off()


MNasefocusedPlus = meta.subsample(focusedTSS, MnPbw, MnMbw, 250, 1, do.sum = T, at.TSS=T)
MNasefocusedMinus = meta.subsample(focusedTSS, MnMbw, MnPbw, 250, 1, do.sum = T, at.TSS=T)
MNasedispersedPlus = meta.subsample(dispersedTSS, MnPbw, MnMbw, 250, 1, do.sum = T, at.TSS=T)
MNasedispersedMinus = meta.subsample(dispersedTSS, MnMbw, MnPbw, 250, 1, do.sum = T, at.TSS=T)

pdf(file = "/Volumes/SEAGATE_EXP/Greg_YeastData/PRO_Cap/Pombe/Analysis/observedTSS/figures/MNasefocusedVsDispersed_observedTSS_Meta.pdf", width = 7, height = 7)
meta.plot.GROseq(MNasefocusedPlus, MNasefocusedMinus, 5, main="SP MNase Focused vs Disp. Obs. TSS",  ylab = "Avg. read count", cex.main=2, cex.lab=1.5, cex.axis=1.5, bothStrands = F)
meta.overlay(MNasedispersedPlus, step = 5, strand = "+")
meta.overlay(MNasedispersedMinus, step = 5, strand = "-")
abline(v=c(0), lwd=2, lty=2)
legend('topleft', c("focused N = 1139", "dispersed N = 1139"), text.col= c("red", "green"), cex = 1.5)
dev.off()

