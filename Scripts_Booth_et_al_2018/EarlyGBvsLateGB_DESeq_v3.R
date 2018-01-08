# 06-19-17 edit:  added another heat map showing significant changes in Promoter Proximal regions.
source("/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/scripts/GB_Functions.R")
library(gplots)
library(ggplot2)
library(lattice)
library(RColorBrewer)
dir.create("/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/figures/EarlyGBvsLateGB_DESeq/10-02-17/")
bwpath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/alignment/bw/"
fig_dir = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/figures/EarlyGBvsLateGB_DESeq/10-02-17/"
infopath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/pombeinfo/"
countpath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/countData/"
ObsTSS_Filt = read.table(file = paste(infopath, "SP_CompleteFilteredGenes.bed", sep = ""))
ObsTSS_Filt_long = ObsTSS_Filt[(ObsTSS_Filt[,3] - ObsTSS_Filt[,2])>800,] ### restrict to genes longer than 800 bp (no overlap of early and late gb)
summary(ObsTSS_Filt_long[,3] - ObsTSS_Filt_long[,2]) # print out the quartiles of lengths
SPmappability = load.bigWig(file = paste(infopath, "SPombeMappable_SacCerNonOverLap_expanded_ChrFix.bw", sep = ""))
AllGenes = read.table(paste(infopath, "pombe.ASM294v1.16.cleaned_sorted_PROcapObservedTSS.bed", sep = ""))[,c(1:6)]
NoOverlapGenes = read.table(paste(infopath, "SP_PROcapObservedTSS_noOverlap.bed", sep = ""))
NFlist = read.table(file = paste(countpath, "replicate_spikeNormFactors.txt", sep = ""), head = T)
CombinedNFs = read.table(file = paste(countpath, "combined_spikeNormFactors.txt", sep = ""), head = T)
#indx = (NoOverlapGenes[,3]-NoOverlapGenes[,2] > 3000) ## filter for using only genes longer than 3 kb
wigset = rbind(c("5993_7157_26971_HJ72CBGXX_pombe_WT_rep1_DMSO_ATCACG_R1.fastq_pombe_plus.bw", "5993_7157_26971_HJ72CBGXX_pombe_WT_rep1_DMSO_ATCACG_R1.fastq_pombe_minus.bw", "WT_DMSOr1"),
               c("5993_7157_26977_HJ72CBGXX_pombe_WT_rep2_DMSO_CAGATC_R1.fastq_pombe_plus.bw", "5993_7157_26977_HJ72CBGXX_pombe_WT_rep2_DMSO_CAGATC_R1.fastq_pombe_minus.bw", "WT_DMSOr2"),
               c("5993_7157_26972_HJ72CBGXX_pombe_WT_rep1_3MB-PPI_CGATGT_R1.fastq_pombe_plus.bw", "5993_7157_26972_HJ72CBGXX_pombe_WT_rep1_3MB-PPI_CGATGT_R1.fastq_pombe_minus.bw", "WT_3MBPP1r1"),
               c("5993_7157_26978_HJ72CBGXX_pombe_WT_rep2_3MB-PPI_ACTTGA_R1.fastq_pombe_plus.bw", "5993_7157_26978_HJ72CBGXX_pombe_WT_rep2_3MB-PPI_ACTTGA_R1.fastq_pombe_minus.bw", "WT_3MBPP1r2"),
               c("5993_7157_26973_HJ72CBGXX_pombe_mcs6as_rep1_DMSO_TTAGGC_R1.fastq_pombe_plus.bw", "5993_7157_26973_HJ72CBGXX_pombe_mcs6as_rep1_DMSO_TTAGGC_R1.fastq_pombe_minus.bw", "mcs6as_DMSOr1"),
               c("5993_7157_26979_HJ72CBGXX_pombe_mcs6as_rep2_DMSO_GATCAG_R1.fastq_pombe_plus.bw", "5993_7157_26979_HJ72CBGXX_pombe_mcs6as_rep2_DMSO_GATCAG_R1.fastq_pombe_minus.bw", "mcs6as_DMSOr2"),
               c("5993_7157_26974_HJ72CBGXX_pombe_mcs6as_rep1_3MB-PPI_TGACCA_R1.fastq_pombe_plus.bw", "5993_7157_26974_HJ72CBGXX_pombe_mcs6as_rep1_3MB-PPI_TGACCA_R1.fastq_pombe_minus.bw", "mcs6as_3MBPP1r1"),
               c("5993_7157_26980_HJ72CBGXX_pombe_mcs6as_rep2_3MB-PPI_TAGCTT_R1.fastq_pombe_plus.bw", "5993_7157_26980_HJ72CBGXX_pombe_mcs6as_rep2_3MB-PPI_TAGCTT_R1.fastq_pombe_minus.bw", "mcs6as_3MBPP1r2"),
               c("5993_7157_26975_HJ72CBGXX_pombe_CDK9as_rep1_DMSO_ACAGTG_R1.fastq_pombe_plus.bw", "5993_7157_26975_HJ72CBGXX_pombe_CDK9as_rep1_DMSO_ACAGTG_R1.fastq_pombe_minus.bw", "CDK9as_DMSOr1"),
               c("5993_7157_26981_HJ72CBGXX_pombe_CDK9as_rep2_DMSO_GGCTAC_R1.fastq_pombe_plus.bw", "5993_7157_26981_HJ72CBGXX_pombe_CDK9as_rep2_DMSO_GGCTAC_R1.fastq_pombe_minus.bw", "CDK9as_DMSOr2"),
               c("5993_7157_26976_HJ72CBGXX_pombe_CDK9as_rep1_3MB-PPI_GCCAAT_R1.fastq_pombe_plus.bw", "5993_7157_26976_HJ72CBGXX_pombe_CDK9as_rep1_3MB-PPI_GCCAAT_R1.fastq_pombe_minus.bw", "CDK9as_3MBPP1r1"),
               c("5993_7157_26982_HJ72CBGXX_pombe_CDK9as_rep2_3MB-PPI_CTTGTA_R1.fastq_pombe_plus.bw", "5993_7157_26982_HJ72CBGXX_pombe_CDK9as_rep2_3MB-PPI_CTTGTA_R1.fastq_pombe_minus.bw", "CDK9as_3MBPP1r2"),
               c("6713_7157_31627_HWV7YBGXX_pombe_Isk1as_DMSO_rep1_ACAGTG_R1.fastq_pombe_plus.bw", "6713_7157_31627_HWV7YBGXX_pombe_Isk1as_DMSO_rep1_ACAGTG_R1.fastq_pombe_minus.bw", "Isk1as_DMSOr1"),
               c("6713_7157_31629_HWV7YBGXX_pombe_Isk1as_DMSO_rep2_CAGATC_R1.fastq_pombe_plus.bw", "6713_7157_31629_HWV7YBGXX_pombe_Isk1as_DMSO_rep2_CAGATC_R1.fastq_pombe_minus.bw", "Isk1as_DMSOr2"),
               c("6713_7157_31628_HWV7YBGXX_pombe_Isk1as_3MB-PP1_rep1_GCCAAT_R1.fastq_pombe_plus.bw", "6713_7157_31628_HWV7YBGXX_pombe_Isk1as_3MB-PP1_rep1_GCCAAT_R1.fastq_pombe_minus.bw", "Isk1as_3MBPP1r1"),
               c("6713_7157_31630_HWV7YBGXX_pombe_Isk1as_3MB-PP1_rep2_ACTTGA_R1.fastq_pombe_plus.bw", "6713_7157_31630_HWV7YBGXX_pombe_Isk1as_3MB-PP1_rep2_ACTTGA_R1.fastq_pombe_minus.bw", "Isk1as_3MBPP1r2"),
               c("6713_7157_31623_HWV7YBGXX_pombe_Mcs6as_CDK9as_DMSO_rep1_ATCACG_R1.fastq_pombe_plus.bw", "6713_7157_31623_HWV7YBGXX_pombe_Mcs6as_CDK9as_DMSO_rep1_ATCACG_R1.fastq_pombe_minus.bw", "Mcs6as_Cdk9as_DMSOr1"),
               c("6713_7157_31625_HWV7YBGXX_pombe_Mcs6as_CDK9as_DMSO_rep2_TTAGGC_R1.fastq_pombe_plus.bw", "6713_7157_31625_HWV7YBGXX_pombe_Mcs6as_CDK9as_DMSO_rep2_TTAGGC_R1.fastq_pombe_minus.bw", "Mcs6as_Cdk9as_DMSOr2"),
               c("6713_7157_31624_HWV7YBGXX_pombe_Mcs6as_CDK9as_3MB-PP1_rep1_CGATGT_R1.fastq_pombe_plus.bw", "6713_7157_31624_HWV7YBGXX_pombe_Mcs6as_CDK9as_3MB-PP1_rep1_CGATGT_R1.fastq_pombe_minus.bw", "Mcs6as_Cdk9as_3MBPP1r1"),
               c("6713_7157_31626_HWV7YBGXX_pombe_Mcs6as_CDK9as_3MB-PP1_rep2_TGACCA_R1.fastq_pombe_plus.bw", "6713_7157_31626_HWV7YBGXX_pombe_Mcs6as_CDK9as_3MB-PP1_rep2_TGACCA_R1.fastq_pombe_minus.bw", "Mcs6as_Cdk9as_3MBPP1r2"),
               c("ePROseq_Pombe_Lsk1asCDK9as_DMSO_rep1_pombe_plus.bw", "ePROseq_Pombe_Lsk1asCDK9as_DMSO_rep1_pombe_minus.bw", "ePRO_Lsk1as_Cdk9as_DMSOr1"),
               c("ePROseq_Pombe_Lsk1asCDK9as_DMSO_rep2_pombe_plus.bw", "ePROseq_Pombe_Lsk1asCDK9as_DMSO_rep2_pombe_minus.bw", "ePRO_Lsk1as_Cdk9as_DMSOr2"),
               c("ePROseq_Pombe_Lsk1asCDK9as_3MBPP1_10uM_rep1_pombe_plus.bw", "ePROseq_Pombe_Lsk1asCDK9as_3MBPP1_10uM_rep1_pombe_minus.bw", "ePRO_Lsk1as_Cdk9as_3MBPP1r1"),
               c("ePROseq_Pombe_Lsk1asCDK9as_3MBPP1_10uM_rep2_pombe_plus.bw", "ePROseq_Pombe_Lsk1asCDK9as_3MBPP1_10uM_rep2_pombe_minus.bw", "ePRO_Lsk1as_Cdk9as_3MBPP1r2"))

TCwigset = rbind(c("ePROseq_Pombe_CDK9as_0min_DMSO_rep1_pombe_plus.bw", "ePROseq_Pombe_CDK9as_0min_DMSO_rep1_pombe_minus.bw", "ePRO_Cdk9as_0min_r1"),
                 c("ePROseq_Pombe_CDK9as_0min_DMSO_rep2_pombe_plus.bw", "ePROseq_Pombe_CDK9as_0min_DMSO_rep2_pombe_minus.bw", "ePRO_Cdk9as_0min_r2"),
                 c("ePROseq_Pombe_CDK9as_30sec_3MBPP1_10uM_rep1_pombe_plus.bw", "ePROseq_Pombe_CDK9as_30sec_3MBPP1_10uM_rep1_pombe_minus.bw", "ePRO_Cdk9as_30sec_r1"),
                 c("ePROseq_Pombe_CDK9as_30sec_3MBPP1_10uM_rep2_pombe_plus.bw", "ePROseq_Pombe_CDK9as_30sec_3MBPP1_10uM_rep2_pombe_minus.bw", "ePRO_Cdk9as_30sec_r2"),
                 c("ePROseq_Pombe_CDK9as_1min_3MBPP1_10uM_rep1_pombe_plus.bw", "ePROseq_Pombe_CDK9as_1min_3MBPP1_10uM_rep1_pombe_minus.bw", "ePRO_Cdk9as_1min_r1"),
                 c("ePROseq_Pombe_CDK9as_1min_3MBPP1_10uM_rep2_pombe_plus.bw", "ePROseq_Pombe_CDK9as_1min_3MBPP1_10uM_rep2_pombe_minus.bw", "ePRO_Cdk9as_1min_r2"),
                 c("ePROseq_Pombe_CDK9as_2min30sec_3MBPP1_10uM_rep1_pombe_plus.bw", "ePROseq_Pombe_CDK9as_2min30sec_3MBPP1_10uM_rep1_pombe_minus.bw", "ePRO_Cdk9as_2min30sec_r1"),
                 c("ePROseq_Pombe_CDK9as_2min30sec_3MBPP1_10uM_rep2_pombe_plus.bw", "ePROseq_Pombe_CDK9as_2min30sec_3MBPP1_10uM_rep2_pombe_minus.bw", "ePRO_Cdk9as_2min30sec_r2"),
                 c("ePROseq_Pombe_CDK9as_5min_3MBPP1_10uM_rep1_pombe_plus.bw", "ePROseq_Pombe_CDK9as_5min_3MBPP1_10uM_rep1_pombe_minus.bw", "ePRO_Cdk9as_5min_r1"),
                 c("ePROseq_Pombe_CDK9as_5min_3MBPP1_10uM_rep2_pombe_plus.bw", "ePROseq_Pombe_CDK9as_5min_3MBPP1_10uM_rep2_pombe_minus.bw", "ePRO_Cdk9as_5min_r2"),
                 c("ePROseq_Pombe_CDK9as_7min30sec_3MBPP1_10uM_rep1_pombe_plus.bw", "ePROseq_Pombe_CDK9as_7min30sec_3MBPP1_10uM_rep1_pombe_minus.bw", "ePRO_Cdk9as_7min30sec_r1"),
                 c("ePROseq_Pombe_CDK9as_7min30sec_3MBPP1_10uM_rep2_pombe_plus.bw", "ePROseq_Pombe_CDK9as_7min30sec_3MBPP1_10uM_rep2_pombe_minus.bw", "ePRO_Cdk9as_7min30sec_r2"),
                 c("ePROseq_Pombe_CDK9as_10min_3MBPP1_10uM_rep1_pombe_plus.bw", "ePROseq_Pombe_CDK9as_10min_3MBPP1_10uM_rep1_pombe_minus.bw", "ePRO_Cdk9as_10min_r1"),
                 c("ePROseq_Pombe_CDK9as_10min_3MBPP1_10uM_rep2_pombe_plus.bw", "ePROseq_Pombe_CDK9as_10min_3MBPP1_10uM_rep2_pombe_minus.bw", "ePRO_Cdk9as_10min_r2"),
                 c("ePROseq_Pombe_CDK9as_20min_3MBPP1_10uM_rep1_pombe_plus.bw", "ePROseq_Pombe_CDK9as_20min_3MBPP1_10uM_rep1_pombe_minus.bw", "ePRO_Cdk9as_20min_r1"),
                 c("ePROseq_Pombe_CDK9as_20min_3MBPP1_10uM_rep2_pombe_plus.bw", "ePROseq_Pombe_CDK9as_20min_3MBPP1_10uM_rep2_pombe_minus.bw", "ePRO_Cdk9as_20min_r2"))
###########################################################################################################################################
###########################################################################################################################################

load.wigset <- function(wigset, wigset_row) {
  file = wigset[wigset_row, 1]
  wig.p = NULL
  if (file != "")
    wig.p = load.bigWig(paste(bwpath, file, sep=''))
  file = wigset[wigset_row, 2]
  wig.m = NULL
  if (file != "")
    wig.m = load.bigWig(paste(bwpath, file, sep=''))
  return(list(wig.p, wig.m, wigset[wigset_row, 3]))
}

# function to split up genes (note that the genes must be long enough to be broken into desired segments)
SepGeneRegions = function(bed){
  pr = matrix(nrow = 0, ncol = 2)
  earlyGB = matrix(nrow = 0, ncol = 2)
  lateGB = matrix(nrow = 0, ncol = 2)
  plus = bed[,6] == "+"
  N = dim(bed)[1]
  for (i in 1:N){
    start = bed[i,2]
    end = bed[i,3]
    if (plus[i]){
      pr = rbind(pr, c(start, start+100))
      earlyGB = rbind(earlyGB, c(start+150, start+450)) #do not include the pause region
      lateGB = rbind(lateGB, c(end-300, end))
    }
    else{
      pr = rbind(pr, c(end-100, end))
      earlyGB = rbind(earlyGB, c(end-450, end-150)) #do not include the pause region
      lateGB = rbind(lateGB, c(start, start+300))
    }
  }
  pr_bed = cbind(bed, pr)[,c(1,7,8,4,5,6)]
  early_bed = cbind(bed, earlyGB)[,c(1,7,8,4,5,6)]
  late_bed = cbind(bed, lateGB)[,c(1,7,8,4,5,6)]
  result = list(pr_bed, early_bed, late_bed)
  return(result)
}
# run the above script to get a list of early and late gene lists
EarlyLateGB_list = SepGeneRegions(ObsTSS_Filt_long)

# function for returning early and late region counts/densities for each sample
## outputs a list of dataframes, each containing the gene name, early count, early map, early density, late count, late map and late density
returnEarlyLateData = function(wigset = wigset, map = SPmappability, genes = EarlyLateGB_list){
  N = dim(wigset)[1]
  resultList = list()
  for (i in 1:N){
    cat("* loading", i, "\n")
    wigs = load.wigset(wigset, i) #this should give all the bigwig files
    cat("* computing ...\n")
    PRgenes = genes[[1]][,c(1,6,2,3,4)]
    Egenes = genes[[2]][,c(1,6,2,3,4)]
    Lgenes = genes[[3]][,c(1,6,2,3,4)]
    PRregion.counts = getCounts(wig.p= wigs[[1]], wig.m = wigs[[2]], PRgenes)
    PRregion.map = getCounts.mappability(map, PRgenes)
    PRregion.density = PRregion.counts[,1] / PRregion.map[,1]
    Eregion.counts = getCounts(wig.p= wigs[[1]], wig.m = wigs[[2]], Egenes)
    Eregion.map = getCounts.mappability(map, Egenes)
    Eregion.density = Eregion.counts[,1] / Eregion.map[,1]
    Lregion.counts = getCounts(wig.p= wigs[[1]], wig.m = wigs[[2]], Lgenes)
    Lregion.map = getCounts.mappability(map, Lgenes)
    Lregion.density = Lregion.counts[,1] / Lregion.map[,1]
    result = cbind(PRregion.counts, PRregion.map, PRregion.density, Eregion.counts, Eregion.map, Eregion.density, Lregion.counts, Lregion.map, Lregion.density)
    rownames(result) <- Egenes[,5]
    colnames(result) <- c("PrCounts", "PrMap", "PrDense", "EarlyCounts", "EarlyMap", "EarlyDense", "LateCount", "LateMap", "LateDense")
    resultList[[wigs[[3]]]] = result
    cat("* unloading.\n")
    unload.wigset(wigs)
  }
  return(resultList)
}
## run the above script to get a list of early and late region count dataframes for each sample
EarlyLateData = returnEarlyLateData(wigset = wigset, map = SPmappability, genes = EarlyLateGB_list)
EarlyLateDataTC = returnEarlyLateData(wigset = TCwigset, map = SPmappability, genes = EarlyLateGB_list)
### Format Early counts and Late counts separately for DESeq analyisis (i.e. cbind all same strain samples)
### using list names to combine strain info
FormatForDESeq = function(countData_list){
  result = list()
  WTres = list()
  Mcs6res = list()
  CDK9res = list()
  Lsk1res = list()
  Mcs6CDK9res = list()
  Lsk1CDK9res = list()
  N = dim(countData_list[[1]])[1]
  LateData = matrix(nrow = N, ncol = 0)
  EarlyData = matrix(nrow = N, ncol = 0)
  PrData = matrix(nrow = N, ncol = 0)
  for (i in grep(pattern = "WT", names(countData_list))){
    PcolData = countData_list[[i]][,1]
    EcolData = countData_list[[i]][,4]
    LcolData = countData_list[[i]][,7]
    PrData = cbind(PrData, PcolData)
    EarlyData = cbind(EarlyData, EcolData)
    LateData = cbind(LateData, LcolData)
  }
  rownames(PrData) <- rownames(countData_list[[1]]) 
  colnames(PrData) <- c("DMSO1", "DMSO2", "3MBPP1_1", "3MBPP1_2")
  rownames(EarlyData) <- rownames(countData_list[[1]]) 
  colnames(EarlyData) <- c("DMSO1", "DMSO2", "3MBPP1_1", "3MBPP1_2")
  rownames(LateData) <- rownames(countData_list[[1]]) 
  colnames(LateData) <- c("DMSO1", "DMSO2", "3MBPP1_1", "3MBPP1_2")
  WTres[["pr"]] = PrData
  WTres[["early"]] =  EarlyData
  WTres[["late"]] =  LateData
  EarlyData = matrix(nrow = N, ncol = 0)
  LateData = matrix(nrow = N, ncol = 0)
  PrData = matrix(nrow = N, ncol = 0)
  for (i in grep(pattern = "mcs6as", names(countData_list))){
    PcolData = countData_list[[i]][,1]
    EcolData = countData_list[[i]][,4]
    LcolData = countData_list[[i]][,7]
    PrData = cbind(PrData, PcolData)
    EarlyData = cbind(EarlyData, EcolData)
    LateData = cbind(LateData, LcolData)
  }
  EarlyData = EarlyData[,c(1,2,4)] ## Not using rep 1 of the treated samples
  LateData = LateData[,c(1,2,4)] ## Not using rep 1 of the treated samples
  PrData = PrData[,c(1,2,4)]## Not using rep 1 of the treated samples
  rownames(PrData) <- rownames(countData_list[[1]]) 
  colnames(PrData) <- c("DMSO1", "DMSO2",  "3MBPP1_2")
  rownames(EarlyData) <- rownames(countData_list[[1]]) 
  colnames(EarlyData) <- c("DMSO1", "DMSO2", "3MBPP1_2")
  rownames(LateData) <- rownames(countData_list[[1]]) 
  colnames(LateData) <- c("DMSO1", "DMSO2", "3MBPP1_2")
  Mcs6res[["pr"]] = PrData
  Mcs6res[["early"]] =  EarlyData
  Mcs6res[["late"]] =  LateData
  EarlyData = matrix(nrow = N, ncol = 0)
  LateData = matrix(nrow = N, ncol = 0)
  PrData = matrix(nrow = N, ncol = 0)
  for (i in grep(pattern = "CDK9as", names(countData_list))){
    PcolData = countData_list[[i]][,1]
    EcolData = countData_list[[i]][,4]
    LcolData = countData_list[[i]][,7]
    PrData = cbind(PrData, PcolData)
    EarlyData = cbind(EarlyData, EcolData)
    LateData = cbind(LateData, LcolData)
  }
  rownames(PrData) <- rownames(countData_list[[1]]) 
  colnames(PrData) <- c("DMSO1", "DMSO2", "3MBPP1_1", "3MBPP1_2")
  rownames(EarlyData) <- rownames(countData_list[[1]]) 
  colnames(EarlyData) <- c("DMSO1", "DMSO2", "3MBPP1_1", "3MBPP1_2")
  rownames(LateData) <- rownames(countData_list[[1]]) 
  colnames(LateData) <- c("DMSO1", "DMSO2", "3MBPP1_1", "3MBPP1_2")
  CDK9res[["pr"]] = PrData
  CDK9res[["early"]] =  EarlyData
  CDK9res[["late"]] =  LateData
  EarlyData = matrix(nrow = N, ncol = 0)
  LateData = matrix(nrow = N, ncol = 0)
  PrData = matrix(nrow = N, ncol = 0)
  for (i in grep(pattern = "Isk1as", names(countData_list))){
    PcolData = countData_list[[i]][,1]
    EcolData = countData_list[[i]][,4]
    LcolData = countData_list[[i]][,7]
    PrData = cbind(PrData, PcolData)
    EarlyData = cbind(EarlyData, EcolData)
    LateData = cbind(LateData, LcolData)
  }
  rownames(PrData) <- rownames(countData_list[[1]]) 
  colnames(PrData) <- c("DMSO1", "DMSO2", "3MBPP1_1", "3MBPP1_2")
  rownames(EarlyData) <- rownames(countData_list[[1]]) 
  colnames(EarlyData) <- c("DMSO1", "DMSO2", "3MBPP1_1", "3MBPP1_2")
  rownames(LateData) <- rownames(countData_list[[1]]) 
  colnames(LateData) <- c("DMSO1", "DMSO2", "3MBPP1_1", "3MBPP1_2")
  Lsk1res[["pr"]] = PrData
  Lsk1res[["early"]] =  EarlyData
  Lsk1res[["late"]] =  LateData
  EarlyData = matrix(nrow = N, ncol = 0)
  LateData = matrix(nrow = N, ncol = 0)
  PrData = matrix(nrow = N, ncol = 0)
  for (i in grep(pattern = "Mcs6as_Cdk9as", names(countData_list))){
    PcolData = countData_list[[i]][,1]
    EcolData = countData_list[[i]][,4]
    LcolData = countData_list[[i]][,7]
    PrData = cbind(PrData, PcolData)
    EarlyData = cbind(EarlyData, EcolData)
    LateData = cbind(LateData, LcolData)
  }
  rownames(PrData) <- rownames(countData_list[[1]]) 
  colnames(PrData) <- c("DMSO1", "DMSO2", "3MBPP1_1", "3MBPP1_2")
  rownames(EarlyData) <- rownames(countData_list[[1]]) 
  colnames(EarlyData) <- c("DMSO1", "DMSO2", "3MBPP1_1", "3MBPP1_2")
  rownames(LateData) <- rownames(countData_list[[1]]) 
  colnames(LateData) <- c("DMSO1", "DMSO2", "3MBPP1_1", "3MBPP1_2")
  Mcs6CDK9res[["pr"]] = PrData
  Mcs6CDK9res[["early"]] =  EarlyData
  Mcs6CDK9res[["late"]] =  LateData
  EarlyData = matrix(nrow = N, ncol = 0)
  LateData = matrix(nrow = N, ncol = 0)
  PrData = matrix(nrow = N, ncol = 0)
  for (i in grep(pattern = "Lsk1as_Cdk9as", names(countData_list))){
    PcolData = countData_list[[i]][,1]
    EcolData = countData_list[[i]][,4]
    LcolData = countData_list[[i]][,7]
    PrData = cbind(PrData, PcolData)
    EarlyData = cbind(EarlyData, EcolData)
    LateData = cbind(LateData, LcolData)
  }
  rownames(PrData) <- rownames(countData_list[[1]]) 
  colnames(PrData) <- c("DMSO1", "DMSO2", "3MBPP1_1", "3MBPP1_2")
  rownames(EarlyData) <- rownames(countData_list[[1]]) 
  colnames(EarlyData) <- c("DMSO1", "DMSO2", "3MBPP1_1", "3MBPP1_2")
  rownames(LateData) <- rownames(countData_list[[1]]) 
  colnames(LateData) <- c("DMSO1", "DMSO2", "3MBPP1_1", "3MBPP1_2")
  Lsk1CDK9res[["pr"]] = PrData
  Lsk1CDK9res[["early"]] =  EarlyData
  Lsk1CDK9res[["late"]] =  LateData
  result[["WT"]] = WTres
  result[["Mcs6"]] = Mcs6res
  result[["CDK9"]] = CDK9res
  result[["Lsk1"]] = Lsk1res
  result[["Mcs6CDK9"]] = Mcs6CDK9res
  result[["Lsk1CDK9"]] = Lsk1CDK9res
  return(result)
}

## Run Formatting Script
FormattedELData = FormatForDESeq(EarlyLateData)

###########################################################################################################################################
# DESeq2 Analysis
###########################################################################################################################################
# Normalization factors (size factors)
NFs = list()
NFs[["WT"]] = NFlist[1:4,"Spikereads"]/100000
NFs[["Mcs6"]] = NFlist[c(5,6,8), "Spikereads"]/100000 # omitting rep 1 of the treated sample
NFs[["CDK9"]] = NFlist[9:12,"Spikereads"]/100000
NFs[["Lsk1"]] = NFlist[13:16,"Spikereads"]/100000
NFs[["Mcs6CDK9"]] = NFlist[17:20,"Spikereads"]/100000
NFs[["Lsk1CDK9"]] = NFlist[67:70,"Spikereads"]/100000

require(DESeq2)
DEseq2.analyze = function(df, sizeFactors, D_reps = 2, M_reps = 2){
  coldata = data.frame(row.names=colnames(df),condition=as.factor(c(rep("DMSO",D_reps),rep("3MB-PP1",M_reps))))
  df = df[complete.cases(df),]
  dds.df <- DESeqDataSetFromMatrix(
    countData = df,
    colData = coldata,
    design = ~ condition)
  dds.df$condition <- relevel(dds.df$condition, "DMSO")
  sizeFactors(dds.df) = c(sizeFactors)
  DESeq.df = DESeq(dds.df)
  res.df = results(DESeq.df, independentFiltering = F) # indipendent filtering seems to create a lot of NAs in the Padj column
  return(res.df)
}

# function for applying DEseq2.analyze to my formated list of dataframes
BulkAnalyze = function(Data, Normfactors, strains = c("WT", "Mcs6", "CDK9", "Lsk1", "Mcs6CDK9", "Lsk1CDK9")){
  result = list()
  for (strain in strains){
    if (strain == "Mcs6"){ ## need this loop because there is one less sample for this experiment
      prRes = DEseq2.analyze(df = Data[[strain]][["pr"]], sizeFactors = Normfactors[[strain]], D_reps = 2, M_reps = 1)
      earlyRes = DEseq2.analyze(df = Data[[strain]][["early"]], sizeFactors = Normfactors[[strain]], D_reps = 2, M_reps = 1)
      lateRes = DEseq2.analyze(df = Data[[strain]][["late"]], sizeFactors = Normfactors[[strain]], D_reps = 2, M_reps = 1)
      Res = cbind(prRes, earlyRes, lateRes)
      #rownames(Res) <- rownames(Data[[strain]][["early"]])
      result[[strain]] = Res
    }
    else{
      prRes = DEseq2.analyze(df = Data[[strain]][["pr"]], sizeFactors = Normfactors[[strain]])
      earlyRes = DEseq2.analyze(df = Data[[strain]][["early"]], sizeFactors = Normfactors[[strain]])
      lateRes = DEseq2.analyze(df = Data[[strain]][["late"]], sizeFactors = Normfactors[[strain]])
      Res = cbind(prRes, earlyRes, lateRes)
      #rownames(Res) <- rownames(Data[[strain]][["early"]])
      result[[strain]] = Res
    }
  }
return(result)
}

DESeqResults = BulkAnalyze(Data = FormattedELData, Normfactors = NFs, strains = c("WT", "Mcs6", "CDK9", "Lsk1", "Mcs6CDK9", "Lsk1CDK9"))


##########################################################################################################
# 03-29-2017 New version of plotting heat maps 
# levelplotting function for plotting DESeq calculated Fold change and Pvalues for early and late gene body regions for all genes (sorted by length to match other heatmaps) 
##########################################################################################################
## prepare a matrix of columns for early and late fold changes for each strain. 
FCdata = cbind(DESeqResults[["WT"]][,c(2,6,8,12,14,18)], DESeqResults[["Lsk1"]][,c(2,6,8,12,14,18)], DESeqResults[["Mcs6"]][,c(2,6,8,12,14,18)], 
               DESeqResults[["CDK9"]][,c(2,6,8,12,14,18)], DESeqResults[["Mcs6CDK9"]][,c(2,6,8,12,14,18)], DESeqResults[["Lsk1CDK9"]][,c(2,6,8,12,14,18)])

### function for sorting heatmaps
sortData = function(genes, df, StoL = F){
  sorted = list()
  res = merge(genes, df, by.x = "V4", by.y = 0)
  if (StoL){
    res = res[order(res$length),]
  }
  else{
    res = res[order(-res$length),]
  }
  result = res[,c(7:dim(res)[2])]
  row.names(result) = res[,1]
  return(result)
}
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
#Sort FC data by gene length
LenSortedGL = orderGenesbyLen(ObsTSS_Filt)
sortedFCdata = sortData(genes = LenSortedGL, FCdata)
# segregate data by strain (each will have early FC, early pVal, late FC, late pVal)
WTdat_sorted = sortedFCdata[,c(2,3,4,5,6,7)]
Lsk1dat_sorted = sortedFCdata[,c(8,9,10,11,12,13)]
Mcs6dat_sorted = sortedFCdata[,c(14,15,16,17,18,19)]
Cdk9dat_sorted = sortedFCdata[,c(20,21,22,23,24,25)]
Mcs6Cdk9dat_sorted = sortedFCdata[,c(26,27,28,29,30,31)]

# levelplotting function for plotting DESeq calculated Fold change and Pvalues for early and late gene body regions for all genes (sorted by length to match other heatmaps) 
plotEarlyLate_FC.Heat = function(HeatRes, filename = "/test_LevPlot.jpeg", simple = T){
  SigSimplify <- function(dat){
    dat = data.frame(dat)
    dat = dat[complete.cases(dat),]
    dat$DirectionP = 0
    dat$DirectionE = 0
    dat$DirectionL = 0
    P_SigUp = (dat[,1] > 0 & dat[,2] < 0.01)
    P_SigDown = (dat[,1] < 0 & dat[,2] < 0.01)
    E_SigUp = (dat[,3] > 0 & dat[,4] < 0.01)
    E_SigDown = (dat[,3] < 0 & dat[,4] < 0.01)
    L_SigUp = (dat[,5] > 0 & dat[,6] < 0.01)
    L_SigDown = (dat[,5] < 0 & dat[,6] < 0.01)
    if(dim(dat[P_SigUp,])[1] > 0){dat[P_SigUp,]$DirectionP = 1}
    if(dim(dat[P_SigDown,])[1] > 0){dat[P_SigDown,]$DirectionP = -1}
    if(dim(dat[E_SigUp,])[1] > 0){dat[E_SigUp,]$DirectionE = 1}
    if(dim(dat[E_SigDown,])[1] > 0){dat[E_SigDown,]$DirectionE = -1}
    if(dim(dat[L_SigUp,])[1] > 0){dat[L_SigUp,]$DirectionL = 1}
    if(dim(dat[L_SigDown,])[1] > 0){dat[L_SigDown,]$DirectionL = -1}
    return(dat)
  }
  #generate each heat map
  if (simple){ ## plots single heat map for early and late regions (blue = sig up, red = sig down)
    Dat = SigSimplify(HeatRes)
    Col_breaks = c(-1.01, -0.99, 0, .99, 1.01)
    col <- colorRampPalette(c('blue', 'white', 'white','red'))
    simpleDirectionPR = levelplot(t(Dat[,7]), ylab = "PR Sig. Up or Down", xlab = NULL, main = "Pr FC",
                                     at = Col_breaks, col.regions = col, colorkey = list(title = "P<0.01"),
                                     scales = list(y = list(draw = F), x = list(draw = F)), aspect = 10)#, useRaster = T)
    simpleDirectionEarly = levelplot(t(Dat[,8]), ylab = "Early GB Sig. Up or Down", xlab = NULL, main = "Early FC",
                            at = Col_breaks, col.regions = col, colorkey = list(title = "P<0.01"),
                            scales = list(y = list(draw = F), x = list(draw = F)), aspect = 10)#, useRaster = T)
    simpleDirectionLate = levelplot(t(Dat[,9]), ylab = "Late GB Sig. Up or Down", xlab = NULL, main = "Late FC",
                           at = Col_breaks, col.regions = col, colorkey = list(title = "P<0.01"), 
                           scales = list(y = list(draw = F), x = list(draw = F)), aspect = 10)#, useRaster = T)
    #jpeg(file = paste(fig_dir, filename, sep = ''), width = 5, height = 10, units = "in")
    pdf(file = paste(fig_dir, filename, sep = ''), width = 5, height = 10)
    plot.new()
    print(simpleDirectionPR,  split = c(1,1,3,1), more = T)
    print(simpleDirectionEarly,  split = c(2,1,3,1), more = T)
    print(simpleDirectionLate,  split = c(3,1,3,1), more = T)
    dev.off()
  }
  else{
    Col_breaksPval = c(0, 0.01)
    Col_breaksFC = seq(-10, 10, length.out = 100)
    col.pVal <- colorRampPalette(c("black", "white"))
    col.fc <- colorRampPalette(c('blue', 'white', 'red'))
    FCplotEarly = levelplot(t(HeatRes[,1]), ylab = "Early GB log2(3MBPP1 / DMSO)", xlab = NULL, main = "Early FC",
                          at = Col_breaksFC, col.regions = col.fc, colorkey = list(title = "log2(3MBPP1 / DMSO)"),
                          scales = list(y = list(draw = F), x = list(draw = F)), aspect = 10)
    pValplotEarly = levelplot(t(HeatRes[,2]),ylab = "Early GB adjusted pVal", xlab = NULL, main = "Early pVal",
                            at = Col_breaksPval, col.regions = col.pVal, colorkey = list(title = "P<0.01"),
                            scales = list(y = list(draw = F), x = list(draw = F)), aspect = 10)
    FCplotLate = levelplot(t(HeatRes[,3]), ylab = "Late GB log2(3MBPP1 / DMSO)", xlab = NULL, main = "Late FC",
                         at = Col_breaksFC, col.regions = col.fc, colorkey = list(title = "log2(3MBPP1 / DMSO)"), 
                         scales = list(y = list(draw = F), x = list(draw = F)), aspect = 10)
    pValplotLate = levelplot(t(HeatRes[,4]), ylab = "Late GB adjusted pVal", xlab = NULL, main = "Late pVal",
                           at = Col_breaksPval, col.regions = col.pVal, colorkey = list(title = "P<0.01"),
                           scales = list(y = list(draw = F), x = list(draw = F)), aspect = 10)
    ##plotting to jpeg
    jpeg(file = paste(fig_dir, filename, sep = ''), width = 10, height = 10, units = "in", res = 300)
    plot.new()
    print(FCplotEarly,  split = c(1,1,4,1), more = T)
    print(pValplotEarly,  split = c(2,1,4,1), more = T)
    print(FCplotLate,  split = c(3,1,4,1), more = T)
    print(pValplotLate,  split = c(4,1,4,1), more = F)
    dev.off()
  }
}
### Prepare early vs late FC heatmaps for each strain. 
plotEarlyLate_FC.Heat(HeatRes = WTdat_sorted, filename = "/WTearlyVsLateGB_FCheatMaps.pdf", simple = T)
plotEarlyLate_FC.Heat(HeatRes = Lsk1dat_sorted, filename = "/Lsk1earlyVsLateGB_FCheatMaps.pdf", simple = T)
plotEarlyLate_FC.Heat(HeatRes = Mcs6dat_sorted, filename = "/Mcs6earlyVsLateGB_FCheatMaps.pdf", simple = T)
plotEarlyLate_FC.Heat(HeatRes = Cdk9dat_sorted, filename = "/Cdk9earlyVsLateGB_FCheatMaps.pdf", simple = T)
plotEarlyLate_FC.Heat(HeatRes = Mcs6Cdk9dat_sorted, filename = "/Mcs6Cdk9earlyVsLateGB_FCheatMaps.pdf", simple = T)


