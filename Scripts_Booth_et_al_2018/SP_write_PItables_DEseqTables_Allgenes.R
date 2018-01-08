source("/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/scripts/GB_Functions.R")
library(lattice)
library(grid)

bwpath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/alignment/bw/"
dir.create("/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/countData/AllGenesObsTSS/")
count_outpath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/countData/AllGenesObsTSS/"

SPmappability = load.bigWig("/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/pombeInfo/SPombeMappable_SacCerNonOverLap_expanded_ChrFix.bw")
AllGenesObs = read.table(file = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/pombeInfo/pombe.ASM294v1.16.cleaned_sorted_PROcapObservedTSS.bed", sep = "\t")
SPGL = AllGenesObs[,c(1,6,2,3,4)]
############
#using above functions to generate count and PI data tables for all samples 
# load all individual sample bigWig files
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

unload.wigset <- function(set) {
  if (!is.null(set[[1]]))
    unload.bigWig(set[[1]])
  if (!is.null(set[[2]]))
    unload.bigWig(set[[2]])
}
calc.pause.index = function(pbw, mbw, mappability = SPmappability, genes = SPGL, name){ 
  PIdata = pause.index(pbw, mbw, mappability, genes)
  result = calcPI.sig(PIdata)
  colnames(result) <- c(paste(name,"_pr_counts", sep=""), paste(name,"gb_counts", sep=""), paste(name,"term_counts", sep=""), paste(name,"pr_map", sep=""), 
                        paste(name,"gb_map", sep=""), paste(name,"term_map", sep=""), paste(name,"pr_density", sep=""), paste(name,"pr_densityRPKM", sep=""), 
                        paste(name,"gb_density", sep=""), paste(name,"gb_densityRPKM", sep=""), paste(name,"pause_index_CI_low", sep=""), 
                        paste(name,"pause_index_CI_high", sep=""), paste(name,"pause_index_scaled", sep=""), paste(name,"pause_index_CI_low_scaled", sep=""), 
                        paste(name,"pause_index_CI_high_scaled", sep=""), paste(name,"_ExpPR", sep=""), paste(name,"_ExpGB", sep=""), paste(name,"Fexact_Pval", sep=""))
  return(result)
}


wigset = rbind(c("5993_7157_26971_HJ72CBGXX_pombe_WT_rep1_DMSO_ATCACG_R1.fastq_pombe_plus.bw", "5993_7157_26971_HJ72CBGXX_pombe_WT_rep1_DMSO_ATCACG_R1.fastq_pombe_minus.bw", "WT_DMSOr1"),
               c("5993_7157_26977_HJ72CBGXX_pombe_WT_rep2_DMSO_CAGATC_R1.fastq_pombe_plus.bw", "5993_7157_26977_HJ72CBGXX_pombe_WT_rep2_DMSO_CAGATC_R1.fastq_pombe_minus.bw", "WT_DMSOr2"),
               c("5993_7157_26972_HJ72CBGXX_pombe_WT_rep1_3MB-PPI_CGATGT_R1.fastq_pombe_plus.bw", "5993_7157_26972_HJ72CBGXX_pombe_WT_rep1_3MB-PPI_CGATGT_R1.fastq_pombe_minus.bw", "WT_3MBPP1r1"),
               c("5993_7157_26978_HJ72CBGXX_pombe_WT_rep2_3MB-PPI_ACTTGA_R1.fastq_pombe_plus.bw", "5993_7157_26978_HJ72CBGXX_pombe_WT_rep2_3MB-PPI_ACTTGA_R1.fastq_pombe_minus.bw", "WT_3MBPP1r2"),
               c("5993_7157_26973_HJ72CBGXX_pombe_mcs6as_rep1_DMSO_TTAGGC_R1.fastq_pombe_plus.bw", "5993_7157_26973_HJ72CBGXX_pombe_mcs6as_rep1_DMSO_TTAGGC_R1.fastq_pombe_minus.bw", "mcs6as_DMSOr1"),
               c("5993_7157_26979_HJ72CBGXX_pombe_mcs6as_rep2_DMSO_GATCAG_R1.fastq_pombe_plus.bw", "5993_7157_26979_HJ72CBGXX_pombe_mcs6as_rep2_DMSO_GATCAG_R1.fastq_pombe_minus.bw", "mcs6as_DMSOr2"),
               c("5993_7157_26974_HJ72CBGXX_pombe_mcs6as_rep1_3MB-PPI_TGACCA_R1.fastq_pombe_plus.bw", "5993_7157_26974_HJ72CBGXX_pombe_mcs6as_rep1_3MB-PPI_TGACCA_R1.fastq_pombe_minus.bw", "mcs6as_3MBPP1r1"),
               c("5993_7157_26980_HJ72CBGXX_pombe_mcs6as_rep2_3MB-PPI_TAGCTT_R1.fastq_pombe_plus.bw", "5993_7157_26980_HJ72CBGXX_pombe_mcs6as_rep2_3MB-PPI_TAGCTT_R1.fastq_pombe_minus.bw", "mcs6as_3MBPP1r2"),
               c("5993_7157_26975_HJ72CBGXX_pombe_CDK9as_rep1_DMSO_ACAGTG_R1.fastq_pombe_plus.bw", "5993_7157_26975_HJ72CBGXX_pombe_CDK9as_rep1_DMSO_ACAGTG_R1.fastq_pombe_minus.bw", "Cdk9as_DMSOr1"),
               c("5993_7157_26981_HJ72CBGXX_pombe_CDK9as_rep2_DMSO_GGCTAC_R1.fastq_pombe_plus.bw", "5993_7157_26981_HJ72CBGXX_pombe_CDK9as_rep2_DMSO_GGCTAC_R1.fastq_pombe_minus.bw", "Cdk9as_DMSOr2"),
               c("5993_7157_26976_HJ72CBGXX_pombe_CDK9as_rep1_3MB-PPI_GCCAAT_R1.fastq_pombe_plus.bw", "5993_7157_26976_HJ72CBGXX_pombe_CDK9as_rep1_3MB-PPI_GCCAAT_R1.fastq_pombe_minus.bw", "Cdk9as_3MBPP1r1"),
               c("5993_7157_26982_HJ72CBGXX_pombe_CDK9as_rep2_3MB-PPI_CTTGTA_R1.fastq_pombe_plus.bw", "5993_7157_26982_HJ72CBGXX_pombe_CDK9as_rep2_3MB-PPI_CTTGTA_R1.fastq_pombe_minus.bw", "Cdk9as_3MBPP1r2"),
               c("6713_7157_31627_HWV7YBGXX_pombe_Isk1as_DMSO_rep1_ACAGTG_R1.fastq_pombe_plus.bw", "6713_7157_31627_HWV7YBGXX_pombe_Isk1as_DMSO_rep1_ACAGTG_R1.fastq_pombe_minus.bw", "Isk1as_DMSOr1"),
               c("6713_7157_31629_HWV7YBGXX_pombe_Isk1as_DMSO_rep2_CAGATC_R1.fastq_pombe_plus.bw", "6713_7157_31629_HWV7YBGXX_pombe_Isk1as_DMSO_rep2_CAGATC_R1.fastq_pombe_minus.bw", "Isk1as_DMSOr2"),
               c("6713_7157_31628_HWV7YBGXX_pombe_Isk1as_3MB-PP1_rep1_GCCAAT_R1.fastq_pombe_plus.bw", "6713_7157_31628_HWV7YBGXX_pombe_Isk1as_3MB-PP1_rep1_GCCAAT_R1.fastq_pombe_minus.bw", "Isk1as_3MBPP1r1"),
               c("6713_7157_31630_HWV7YBGXX_pombe_Isk1as_3MB-PP1_rep2_ACTTGA_R1.fastq_pombe_plus.bw", "6713_7157_31630_HWV7YBGXX_pombe_Isk1as_3MB-PP1_rep2_ACTTGA_R1.fastq_pombe_minus.bw", "Isk1as_3MBPP1r2"),
               c("6713_7157_31623_HWV7YBGXX_pombe_Mcs6as_CDK9as_DMSO_rep1_ATCACG_R1.fastq_pombe_plus.bw", "6713_7157_31623_HWV7YBGXX_pombe_Mcs6as_CDK9as_DMSO_rep1_ATCACG_R1.fastq_pombe_minus.bw", "Mcs6as_Cdk9as_DMSOr1"),
               c("6713_7157_31625_HWV7YBGXX_pombe_Mcs6as_CDK9as_DMSO_rep2_TTAGGC_R1.fastq_pombe_plus.bw", "6713_7157_31625_HWV7YBGXX_pombe_Mcs6as_CDK9as_DMSO_rep2_TTAGGC_R1.fastq_pombe_minus.bw", "Mcs6as_Cdk9as_DMSOr2"),
               c("6713_7157_31624_HWV7YBGXX_pombe_Mcs6as_CDK9as_3MB-PP1_rep1_CGATGT_R1.fastq_pombe_plus.bw", "6713_7157_31624_HWV7YBGXX_pombe_Mcs6as_CDK9as_3MB-PP1_rep1_CGATGT_R1.fastq_pombe_minus.bw", "Mcs6as_Cdk9as_3MBPP1r1"),
               c("6713_7157_31626_HWV7YBGXX_pombe_Mcs6as_CDK9as_3MB-PP1_rep2_TGACCA_R1.fastq_pombe_plus.bw", "6713_7157_31626_HWV7YBGXX_pombe_Mcs6as_CDK9as_3MB-PP1_rep2_TGACCA_R1.fastq_pombe_minus.bw", "Mcs6as_Cdk9as_3MBPP1r2"),
               c("7772_7157_43139_HVLYCBGXY_pombe_18C_Cdk9as_DMSO_rep1_ATCACG_R1_pombe_plus.bw", "7772_7157_43139_HVLYCBGXY_pombe_18C_Cdk9as_DMSO_rep1_ATCACG_R1_pombe_minus.bw", "Cdk9as_DMSO_18Cr1"),
               c("7772_7157_43149_HVLYCBGXY_pombe_18C_Cdk9as_DMSO_rep2_GGCTAC_R1_pombe_plus.bw", "7772_7157_43149_HVLYCBGXY_pombe_18C_Cdk9as_DMSO_rep2_GGCTAC_R1_pombe_minus.bw", "Cdk9as_DMSO_18Cr2"),
               c("7772_7157_43140_HVLYCBGXY_pombe_18C_Cdk9as_3MBPP1_rep1_CGATGT_R1_pombe_plus.bw", "7772_7157_43140_HVLYCBGXY_pombe_18C_Cdk9as_3MBPP1_rep1_CGATGT_R1_pombe_minus.bw", "Cdk9as_3MBPP1_18Cr1"),
               c("7772_7157_43150_HVLYCBGXY_pombe_18C_Cdk9as_3MBPP1_rep2_CTTGTA_R1_pombe_plus.bw", "7772_7157_43150_HVLYCBGXY_pombe_18C_Cdk9as_3MBPP1_rep2_CTTGTA_R1_pombe_minus.bw", "Cdk9as_3MBPP1_18Cr2"),
               c("7772_7157_43141_HVLYCBGXY_pombe_18C_dis2ts_DMSO_rep1_TTAGGC_R1_pombe_plus.bw", "7772_7157_43141_HVLYCBGXY_pombe_18C_dis2ts_DMSO_rep1_TTAGGC_R1_pombe_minus.bw", "dis2ts_DMSO_18Cr1"),
               c("7772_7157_43151_HVLYCBGXY_pombe_18C_dis2ts_DMSO_rep2_ATTCCT_R1_pombe_plus.bw", "7772_7157_43151_HVLYCBGXY_pombe_18C_dis2ts_DMSO_rep2_ATTCCT_R1_pombe_minus.bw", "dis2ts_DMSO_18Cr2"),
               c("7772_7157_43142_HVLYCBGXY_pombe_18C_dis2ts_3MBPP1_rep1_TGACCA_R1_pombe_plus.bw", "7772_7157_43142_HVLYCBGXY_pombe_18C_dis2ts_3MBPP1_rep1_TGACCA_R1_pombe_minus.bw", "dis2ts_3MBPP1_18Cr1"),
               c("7772_7157_43152_HVLYCBGXY_pombe_18C_dis2ts_3MBPP1_rep2_CAAAAG_R1_pombe_plus.bw", "7772_7157_43152_HVLYCBGXY_pombe_18C_dis2ts_3MBPP1_rep2_CAAAAG_R1_pombe_minus.bw", "dis2ts_3MBPP1_18Cr2"),
               c("7772_7157_43143_HVLYCBGXY_pombe_30C_dis2ts_DMSO_rep1_ACAGTG_R1_pombe_plus.bw", "7772_7157_43143_HVLYCBGXY_pombe_30C_dis2ts_DMSO_rep1_ACAGTG_R1_pombe_minus.bw", "dis2ts_DMSO_30Cr1"),
               c("7772_7157_43153_HVLYCBGXY_pombe_30C_dis2ts_DMSO_rep2_CAACTA_R1_pombe_plus.bw", "7772_7157_43153_HVLYCBGXY_pombe_30C_dis2ts_DMSO_rep2_CAACTA_R1_pombe_minus.bw", "dis2ts_DMSO_30Cr2"),
               c("7772_7157_43144_HVLYCBGXY_pombe_30C_dis2ts_3MBPP1_rep1_GCCAAT_R1_pombe_plus.bw", "7772_7157_43144_HVLYCBGXY_pombe_30C_dis2ts_3MBPP1_rep1_GCCAAT_R1_pombe_minus.bw", "dis2ts_3MBPP1_30Cr1"),
               c("7772_7157_43154_HVLYCBGXY_pombe_30C_dis2ts_3MBPP1_rep2_CACCGG_R1_pombe_plus.bw", "7772_7157_43154_HVLYCBGXY_pombe_30C_dis2ts_3MBPP1_rep2_CACCGG_R1_pombe_minus.bw", "dis2ts_3MBPP1_30Cr2"),
               c("7772_7157_43145_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_DMSO_rep1_CAGATC_R1_pombe_plus.bw", "7772_7157_43145_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_DMSO_rep1_CAGATC_R1_pombe_minus.bw", "CDK9as_dis2ts_DMSO_18Cr1"),
               c("7772_7157_43155_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_DMSO_rep2_CACGAT_R1_pombe_plus.bw", "7772_7157_43155_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_DMSO_rep2_CACGAT_R1_pombe_minus.bw", "CDK9as_dis2ts_DMSO_18Cr2"),
               c("7772_7157_43146_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_3MBPP1_rep1_ACTTGA_R1_pombe_plus.bw", "7772_7157_43146_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_3MBPP1_rep1_ACTTGA_R1_pombe_minus.bw", "CDK9as_dis2ts_3MBPP1_18Cr1"),
               c("7772_7157_43156_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_3MBPP1_rep2_CACTCA_R1_pombe_plus.bw", "7772_7157_43156_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_3MBPP1_rep2_CACTCA_R1_pombe_minus.bw", "CDK9as_dis2ts_3MBPP1_18Cr2"),
               c("7772_7157_43147_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_DMSO_rep1_GATCAG_R1_pombe_plus.bw", "7772_7157_43147_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_DMSO_rep1_GATCAG_R1_pombe_minus.bw", "CDK9as_dis2ts_DMSO_30Cr1"),
               c("7772_7157_43157_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_DMSO_rep2_CAGGCG_R1_pombe_plus.bw", "7772_7157_43157_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_DMSO_rep2_CAGGCG_R1_pombe_minus.bw", "CDK9as_dis2ts_DMSO_30Cr2"),
               c("7772_7157_43148_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_3MBPP1_rep1_TAGCTT_R1_pombe_plus.bw", "7772_7157_43148_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_3MBPP1_rep1_TAGCTT_R1_pombe_minus.bw", "CDK9as_dis2ts_3MBPP1_30Cr1"),
               c("7772_7157_43158_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_3MBPP1_rep2_CATGGC_R1_pombe_plus.bw", "7772_7157_43158_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_3MBPP1_rep2_CATGGC_R1_pombe_minus.bw", "CDK9as_dis2ts_3MBPP1_30Cr2"),
               c("ePROseq_Pombe_CDK9as_0min_DMSO_rep1_pombe_plus.bw", "ePROseq_Pombe_CDK9as_0min_DMSO_rep1_pombe_minus.bw", "ePRO_Cdk9as_0min_r1"),
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
               c("ePROseq_Pombe_CDK9as_20min_3MBPP1_10uM_rep2_pombe_plus.bw", "ePROseq_Pombe_CDK9as_20min_3MBPP1_10uM_rep2_pombe_minus.bw", "ePRO_Cdk9as_20min_r2"),
               c("ePROseq_Pombe_Lsk1asCDK9as_DMSO_rep1_pombe_plus.bw", "ePROseq_Pombe_Lsk1asCDK9as_DMSO_rep1_pombe_minus.bw", "ePRO_Lsk1as_Cdk9as_DMSOr1"),
               c("ePROseq_Pombe_Lsk1asCDK9as_DMSO_rep2_pombe_plus.bw", "ePROseq_Pombe_Lsk1asCDK9as_DMSO_rep2_pombe_minus.bw", "ePRO_Lsk1as_Cdk9as_DMSOr2"),
               c("ePROseq_Pombe_Lsk1asCDK9as_3MBPP1_10uM_rep1_pombe_plus.bw", "ePROseq_Pombe_Lsk1asCDK9as_3MBPP1_10uM_rep1_pombe_minus.bw", "ePRO_Lsk1as_Cdk9as_3MBPP1r1"),
               c("ePROseq_Pombe_Lsk1asCDK9as_3MBPP1_10uM_rep2_pombe_plus.bw", "ePROseq_Pombe_Lsk1asCDK9as_3MBPP1_10uM_rep2_pombe_minus.bw", "ePRO_Lsk1as_Cdk9as_3MBPP1r2"),
               c("ePROseq_Pombe_Booth_Cdk9as_spt5_WT7_0min_DMSO_rep1_pombe_plus.bw", "ePROseq_Pombe_Booth_Cdk9as_spt5_WT7_0min_DMSO_rep1_pombe_minus.bw", "Cdk9as_spt5_WT7_0min_DMSO_r1"),
               c("ePROseq_Pombe_Booth_Cdk9as_spt5_WT7_0min_DMSO_rep2_pombe_plus.bw", "ePROseq_Pombe_Booth_Cdk9as_spt5_WT7_0min_DMSO_rep2_pombe_minus.bw", "Cdk9as_spt5_WT7_0min_DMSO_r2"),
               c("ePROseq_Pombe_Booth_Cdk9as_spt5_WT7_1min_3MBPP1_rep1_pombe_plus.bw", "ePROseq_Pombe_Booth_Cdk9as_spt5_WT7_1min_3MBPP1_rep1_pombe_minus.bw", "Cdk9as_spt5_WT7_1min_3MBPP1_r1"),
               c("ePROseq_Pombe_Booth_Cdk9as_spt5_WT7_1min_3MBPP1_rep2_pombe_plus.bw", "ePROseq_Pombe_Booth_Cdk9as_spt5_WT7_1min_3MBPP1_rep2_pombe_minus.bw", "Cdk9as_spt5_WT7_1min_3MBPP1_r2"),
               c("ePROseq_Pombe_Booth_Cdk9as_spt5_WT7_5min_3MBPP1_rep1_pombe_plus.bw", "ePROseq_Pombe_Booth_Cdk9as_spt5_WT7_5min_3MBPP1_rep1_pombe_minus.bw", "Cdk9as_spt5_WT7_5min_3MBPP1_r1"),
               c("ePROseq_Pombe_Booth_Cdk9as_spt5_WT7_5min_3MBPP1_rep2_pombe_plus.bw", "ePROseq_Pombe_Booth_Cdk9as_spt5_WT7_5min_3MBPP1_rep2_pombe_minus.bw", "Cdk9as_spt5_WT7_5min_3MBPP1_r2"),
               c("ePROseq_Pombe_Booth_Cdk9as_spt5_T1A_0min_DMSO_rep1_pombe_plus.bw", "ePROseq_Pombe_Booth_Cdk9as_spt5_T1A_0min_DMSO_rep1_pombe_minus.bw", "Cdk9as_spt5_T1A_0min_DMSO_r1"),
               c("ePROseq_Pombe_Booth_Cdk9as_spt5_T1A_0min_DMSO_rep2_pombe_plus.bw", "ePROseq_Pombe_Booth_Cdk9as_spt5_T1A_0min_DMSO_rep2_pombe_minus.bw", "Cdk9as_spt5_T1A_0min_DMSO_r2"),
               c("ePROseq_Pombe_Booth_Cdk9as_spt5_T1A_1min_3MBPP1_rep1_pombe_plus.bw", "ePROseq_Pombe_Booth_Cdk9as_spt5_T1A_1min_3MBPP1_rep1_pombe_minus.bw", "Cdk9as_spt5_T1A_1min_3MBPP1_r1"),
               c("ePROseq_Pombe_Booth_Cdk9as_spt5_T1A_1min_3MBPP1_rep2_pombe_plus.bw", "ePROseq_Pombe_Booth_Cdk9as_spt5_T1A_1min_3MBPP1_rep2_pombe_minus.bw", "Cdk9as_spt5_T1A_1min_3MBPP1_r2"),
               c("ePROseq_Pombe_Booth_Cdk9as_spt5_T1A_5min_3MBPP1_rep1_pombe_plus.bw", "ePROseq_Pombe_Booth_Cdk9as_spt5_T1A_5min_3MBPP1_rep1_pombe_minus.bw", "Cdk9as_spt5_T1A_5min_3MBPP1_r1"),
               c("ePROseq_Pombe_Booth_Cdk9as_spt5_T1A_5min_3MBPP1_rep2_pombe_plus.bw", "ePROseq_Pombe_Booth_Cdk9as_spt5_T1A_5min_3MBPP1_rep2_pombe_minus.bw", "Cdk9as_spt5_T1A_5min_3MBPP1_r2"),
               c("ePROseq_Pombe_Booth_Cdk9as_spt5_T1E_0min_DMSO_rep1_pombe_plus.bw", "ePROseq_Pombe_Booth_Cdk9as_spt5_T1E_0min_DMSO_rep1_pombe_minus.bw", "Cdk9as_spt5_T1E_0min_DMSO_r1"),
               c("ePROseq_Pombe_Booth_Cdk9as_spt5_T1E_0min_DMSO_rep2_pombe_plus.bw", "ePROseq_Pombe_Booth_Cdk9as_spt5_T1E_0min_DMSO_rep2_pombe_minus.bw", "Cdk9as_spt5_T1E_0min_DMSO_r2"),
               c("ePROseq_Pombe_Booth_Cdk9as_spt5_T1E_1min_3MBPP1_rep1_pombe_plus.bw", "ePROseq_Pombe_Booth_Cdk9as_spt5_T1E_1min_3MBPP1_rep1_pombe_minus.bw", "Cdk9as_spt5_T1E_1min_3MBPP1_r1"),
               c("ePROseq_Pombe_Booth_Cdk9as_spt5_T1E_1min_3MBPP1_rep2_pombe_plus.bw", "ePROseq_Pombe_Booth_Cdk9as_spt5_T1E_1min_3MBPP1_rep2_pombe_minus.bw", "Cdk9as_spt5_T1E_1min_3MBPP1_r2"),
               c("ePROseq_Pombe_Booth_Cdk9as_spt5_T1E_5min_3MBPP1_rep1_pombe_plus.bw", "ePROseq_Pombe_Booth_Cdk9as_spt5_T1E_5min_3MBPP1_rep1_pombe_minus.bw", "Cdk9as_spt5_T1E_5min_3MBPP1_r1"),
               c("ePROseq_Pombe_Booth_Cdk9as_spt5_T1E_5min_3MBPP1_rep2_pombe_plus.bw", "ePROseq_Pombe_Booth_Cdk9as_spt5_T1E_5min_3MBPP1_rep2_pombe_minus.bw", "Cdk9as_spt5_T1E_5min_3MBPP1_r2"),
               c("ePROseq_Pombe_Fisher_WT_rep1_pombe_plus.bw", "ePROseq_Pombe_Fisher_WT_rep1_pombe_minus.bw", "Fisher_WT_r1"),
               c("ePROseq_Pombe_Fisher_WT_rep2_pombe_plus.bw", "ePROseq_Pombe_Fisher_WT_rep2_pombe_minus.bw", "Fisher_WT_r2"),
               c("ePROseq_Pombe_Fisher_dis2_11_rep1_pombe_plus.bw", "ePROseq_Pombe_Fisher_dis2_11_rep1_pombe_minus.bw", "Fisher_dis2_11_r1"),
               c("ePROseq_Pombe_Fisher_dis2_11_rep2_pombe_plus.bw", "ePROseq_Pombe_Fisher_dis2_11_rep2_pombe_minus.bw", "Fisher_dis2_11_r2"),
               c("ePROseq_Pombe_Fisher_dis2Del_rep1_pombe_plus.bw", "ePROseq_Pombe_Fisher_dis2Del_rep1_pombe_minus.bw", "Fisher_dis2Del_r1"),
               c("ePROseq_Pombe_Fisher_dis2Del_rep2_pombe_plus.bw", "ePROseq_Pombe_Fisher_dis2Del_rep2_pombe_minus.bw", "Fisher_dis2Del_r2"),
               c("ePROseq_Pombe_Fisher_dis2_T316A_rep1_pombe_plus.bw", "ePROseq_Pombe_Fisher_dis2_T316A_rep1_pombe_minus.bw", "Fisher_dis2_T316A_r1"),
               c("ePROseq_Pombe_Fisher_dis2_T316A_rep2_pombe_plus.bw", "ePROseq_Pombe_Fisher_dis2_T316A_rep2_pombe_minus.bw", "Fisher_dis2_T316A_r2"),
               c("ePROseq_Pombe_Fisher_dis2_T316D_rep1_pombe_plus.bw", "ePROseq_Pombe_Fisher_dis2_T316D_rep1_pombe_minus.bw", "Fisher_dis2_T316D_r1"),
               c("ePROseq_Pombe_Fisher_dis2_T316D_rep2_pombe_plus.bw", "ePROseq_Pombe_Fisher_dis2_T316D_rep2_pombe_minus.bw", "Fisher_dis2_T316D_r2"),
               c("5993_7157_26977_HJ72CBGXX_pombe_WT_COMBINED_DMSO_CAGATC_R1.fastq_pombe_plus.bw", "5993_7157_26977_HJ72CBGXX_pombe_WT_COMBINED_DMSO_CAGATC_R1.fastq_pombe_minus.bw", "WT_DMSO_combined"),
               c("5993_7157_26978_HJ72CBGXX_pombe_WT_COMBINED_3MB-PPI_ACTTGA_R1.fastq_pombe_plus.bw", "5993_7157_26978_HJ72CBGXX_pombe_WT_COMBINED_3MB-PPI_ACTTGA_R1.fastq_pombe_minus.bw", "WT_3MBPP1_combined"),
               c("5993_7157_26979_HJ72CBGXX_pombe_mcs6as_COMBINED_DMSO_GATCAG_R1.fastq_pombe_plus.bw", "5993_7157_26979_HJ72CBGXX_pombe_mcs6as_COMBINED_DMSO_GATCAG_R1.fastq_pombe_minus.bw", "mcs6as_DMSO_combined"),
               c("5993_7157_26980_HJ72CBGXX_pombe_mcs6as_COMBINED_3MB-PPI_TAGCTT_R1.fastq_pombe_plus.bw", "5993_7157_26980_HJ72CBGXX_pombe_mcs6as_COMBINED_3MB-PPI_TAGCTT_R1.fastq_pombe_minus.bw", "mcs6as_3MBPP1_combined"),
               c("5993_7157_26981_HJ72CBGXX_pombe_CDK9as_COMBINED_DMSO_GGCTAC_R1.fastq_pombe_plus.bw", "5993_7157_26981_HJ72CBGXX_pombe_CDK9as_COMBINED_DMSO_GGCTAC_R1.fastq_pombe_minus.bw", "Cdk9as_DMSO_combined"),
               c("5993_7157_26982_HJ72CBGXX_pombe_CDK9as_COMBINED_3MB-PPI_CTTGTA_R1.fastq_pombe_plus.bw", "5993_7157_26982_HJ72CBGXX_pombe_CDK9as_COMBINED_3MB-PPI_CTTGTA_R1.fastq_pombe_minus.bw", "Cdk9as_3MBPP1_combined"),
               c("6713_7157_31629_HWV7YBGXX_pombe_Isk1as_DMSO_COMBINED_CAGATC_R1.fastq_pombe_plus.bw", "6713_7157_31629_HWV7YBGXX_pombe_Isk1as_DMSO_COMBINED_CAGATC_R1.fastq_pombe_minus.bw", "Isk1as_DMSO_combined"),
               c("6713_7157_31630_HWV7YBGXX_pombe_Isk1as_3MB-PP1_COMBINED_ACTTGA_R1.fastq_pombe_plus.bw", "6713_7157_31630_HWV7YBGXX_pombe_Isk1as_3MB-PP1_COMBINED_ACTTGA_R1.fastq_pombe_minus.bw", "Isk1as_3MBPP1_combined"),
               c("6713_7157_31625_HWV7YBGXX_pombe_Mcs6as_CDK9as_DMSO_COMBINED_TTAGGC_R1.fastq_pombe_plus.bw", "6713_7157_31625_HWV7YBGXX_pombe_Mcs6as_CDK9as_DMSO_COMBINED_TTAGGC_R1.fastq_pombe_minus.bw", "Mcs6as_Cdk9as_DMSO_combined"),
               c("6713_7157_31626_HWV7YBGXX_pombe_Mcs6as_CDK9as_3MB-PP1_COMBINED_TGACCA_R1.fastq_pombe_plus.bw", "6713_7157_31626_HWV7YBGXX_pombe_Mcs6as_CDK9as_3MB-PP1_COMBINED_TGACCA_R1.fastq_pombe_minus.bw", "Mcs6as_Cdk9as_3MBPP1_combined"),
               c("7772_7157_43149_HVLYCBGXY_pombe_18C_Cdk9as_DMSO_COMBINED_GGCTAC_R1_pombe_plus.bw", "7772_7157_43149_HVLYCBGXY_pombe_18C_Cdk9as_DMSO_COMBINED_GGCTAC_R1_pombe_minus.bw", "Cdk9as_DMSO_18C_combined"),
               c("7772_7157_43150_HVLYCBGXY_pombe_18C_Cdk9as_3MBPP1_COMBINED_CTTGTA_R1_pombe_plus.bw", "7772_7157_43150_HVLYCBGXY_pombe_18C_Cdk9as_3MBPP1_COMBINED_CTTGTA_R1_pombe_minus.bw", "Cdk9as_3MBPP1_18C_combined"),
               c("7772_7157_43151_HVLYCBGXY_pombe_18C_dis2ts_DMSO_COMBINED_ATTCCT_R1_pombe_plus.bw", "7772_7157_43151_HVLYCBGXY_pombe_18C_dis2ts_DMSO_COMBINED_ATTCCT_R1_pombe_minus.bw", "dis2ts_DMSO_18C_combined"),
               c("7772_7157_43152_HVLYCBGXY_pombe_18C_dis2ts_3MBPP1_COMBINED_CAAAAG_R1_pombe_plus.bw", "7772_7157_43152_HVLYCBGXY_pombe_18C_dis2ts_3MBPP1_COMBINED_CAAAAG_R1_pombe_minus.bw", "dis2ts_3MBPP1_18C_combined"),
               c("7772_7157_43153_HVLYCBGXY_pombe_30C_dis2ts_DMSO_COMBINED_CAACTA_R1_pombe_plus.bw", "7772_7157_43153_HVLYCBGXY_pombe_30C_dis2ts_DMSO_COMBINED_CAACTA_R1_pombe_minus.bw", "dis2ts_DMSO_30C_combined"),
               c("7772_7157_43154_HVLYCBGXY_pombe_30C_dis2ts_3MBPP1_COMBINED_CACCGG_R1_pombe_plus.bw", "7772_7157_43154_HVLYCBGXY_pombe_30C_dis2ts_3MBPP1_COMBINED_CACCGG_R1_pombe_minus.bw", "dis2ts_3MBPP1_30C_combined"),
               c("7772_7157_43155_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_DMSO_COMBINED_CACGAT_R1_pombe_plus.bw", "7772_7157_43155_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_DMSO_COMBINED_CACGAT_R1_pombe_minus.bw", "CDK9as_dis2ts_DMSO_18C_combined"),
               c("7772_7157_43156_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_3MBPP1_COMBINED_CACTCA_R1_pombe_plus.bw", "7772_7157_43156_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_3MBPP1_COMBINED_CACTCA_R1_pombe_minus.bw", "CDK9as_dis2ts_3MBPP1_18C_combined"),
               c("7772_7157_43157_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_DMSO_COMBINED_CAGGCG_R1_pombe_plus.bw", "7772_7157_43157_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_DMSO_COMBINED_CAGGCG_R1_pombe_minus.bw", "CDK9as_dis2ts_DMSO_30C_combined"),
               c("7772_7157_43158_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_3MBPP1_COMBINED_CATGGC_R1_pombe_plus.bw", "7772_7157_43158_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_3MBPP1_COMBINED_CATGGC_R1_pombe_minus.bw", "CDK9as_dis2ts_3MBPP1_30C_combined"),
               c("ePROseq_Pombe_CDK9as_0min_DMSO_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_0min_DMSO_COMBINED_pombe_minus.bw", "ePRO_Cdk9as_0min_combined"),
               c("ePROseq_Pombe_CDK9as_30sec_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_30sec_3MBPP1_10uM_COMBINED_pombe_minus.bw", "ePRO_Cdk9as_30sec_combined"),
               c("ePROseq_Pombe_CDK9as_1min_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_1min_3MBPP1_10uM_COMBINED_pombe_minus.bw", "ePRO_Cdk9as_1min_combined"),
               c("ePROseq_Pombe_CDK9as_2min30sec_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_2min30sec_3MBPP1_10uM_COMBINED_pombe_minus.bw", "ePRO_Cdk9as_2min30sec_combined"),
               c("ePROseq_Pombe_CDK9as_5min_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_5min_3MBPP1_10uM_COMBINED_pombe_minus.bw", "ePRO_Cdk9as_5min_combined"),
               c("ePROseq_Pombe_CDK9as_7min30sec_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_7min30sec_3MBPP1_10uM_COMBINED_pombe_minus.bw", "ePRO_Cdk9as_7min30sec_combined"),
               c("ePROseq_Pombe_CDK9as_10min_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_10min_3MBPP1_10uM_COMBINED_pombe_minus.bw", "ePRO_Cdk9as_10min_combined"),
               c("ePROseq_Pombe_CDK9as_20min_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_CDK9as_20min_3MBPP1_10uM_COMBINED_pombe_minus.bw", "ePRO_Cdk9as_20min_combined"),
               c("ePROseq_Pombe_Lsk1asCDK9as_DMSO_COMBINED_pombe_plus.bw", "ePROseq_Pombe_Lsk1asCDK9as_DMSO_COMBINED_pombe_minus.bw", "ePRO_Lsk1as_Cdk9as_DMSO_combined"),
               c("ePROseq_Pombe_Lsk1asCDK9as_3MBPP1_10uM_COMBINED_pombe_plus.bw", "ePROseq_Pombe_Lsk1asCDK9as_3MBPP1_10uM_COMBINED_pombe_minus.bw", "ePRO_Lsk1as_Cdk9as_3MBPP1_combined"),
               c("ePROseq_Pombe_Booth_Cdk9as_spt5_WT7_0min_DMSO_COMBINED_pombe_plus.bw", "ePROseq_Pombe_Booth_Cdk9as_spt5_WT7_0min_DMSO_COMBINED_pombe_minus.bw", "Cdk9as_spt5_WT7_0min_DMSO_combined"),
               c("ePROseq_Pombe_Booth_Cdk9as_spt5_WT7_1min_3MBPP1_COMBINED_pombe_plus.bw", "ePROseq_Pombe_Booth_Cdk9as_spt5_WT7_1min_3MBPP1_COMBINED_pombe_minus.bw", "Cdk9as_spt5_WT7_1min_3MBPP1_combined"),
               c("ePROseq_Pombe_Booth_Cdk9as_spt5_WT7_5min_3MBPP1_COMBINED_pombe_plus.bw", "ePROseq_Pombe_Booth_Cdk9as_spt5_WT7_5min_3MBPP1_COMBINED_pombe_minus.bw", "Cdk9as_spt5_WT7_5min_3MBPP1_combined"),
               c("ePROseq_Pombe_Booth_Cdk9as_spt5_T1A_0min_DMSO_COMBINED_pombe_plus.bw", "ePROseq_Pombe_Booth_Cdk9as_spt5_T1A_0min_DMSO_COMBINED_pombe_minus.bw", "Cdk9as_spt5_T1A_0min_DMSO_combined"),
               c("ePROseq_Pombe_Booth_Cdk9as_spt5_T1A_1min_3MBPP1_COMBINED_pombe_plus.bw", "ePROseq_Pombe_Booth_Cdk9as_spt5_T1A_1min_3MBPP1_COMBINED_pombe_minus.bw", "Cdk9as_spt5_T1A_1min_3MBPP1_combined"),
               c("ePROseq_Pombe_Booth_Cdk9as_spt5_T1A_5min_3MBPP1_COMBINED_pombe_plus.bw", "ePROseq_Pombe_Booth_Cdk9as_spt5_T1A_5min_3MBPP1_COMBINED_pombe_minus.bw", "Cdk9as_spt5_T1A_5min_3MBPP1_combined"),
               c("ePROseq_Pombe_Booth_Cdk9as_spt5_T1E_0min_DMSO_COMBINED_pombe_plus.bw", "ePROseq_Pombe_Booth_Cdk9as_spt5_T1E_0min_DMSO_COMBINED_pombe_minus.bw", "Cdk9as_spt5_T1E_0min_DMSO_combined"),
               c("ePROseq_Pombe_Booth_Cdk9as_spt5_T1E_1min_3MBPP1_COMBINED_pombe_plus.bw", "ePROseq_Pombe_Booth_Cdk9as_spt5_T1E_1min_3MBPP1_COMBINED_pombe_minus.bw", "Cdk9as_spt5_T1E_1min_3MBPP1_combined"),
               c("ePROseq_Pombe_Booth_Cdk9as_spt5_T1E_5min_3MBPP1_COMBINED_pombe_plus.bw", "ePROseq_Pombe_Booth_Cdk9as_spt5_T1E_5min_3MBPP1_COMBINED_pombe_minus.bw", "Cdk9as_spt5_T1E_5min_3MBPP1_combined"),
               c("ePROseq_Pombe_Fisher_WT_COMBINED_pombe_plus.bw", "ePROseq_Pombe_Fisher_WT_COMBINED_pombe_minus.bw", "Fisher_WT_combined"),
               c("ePROseq_Pombe_Fisher_dis2_11_COMBINED_pombe_plus.bw", "ePROseq_Pombe_Fisher_dis2_11_COMBINED_pombe_minus.bw", "Fisher_dis2_11_combined"),
               c("ePROseq_Pombe_Fisher_dis2Del_COMBINED_pombe_plus.bw", "ePROseq_Pombe_Fisher_dis2Del_COMBINED_pombe_minus.bw", "Fisher_dis2Del_combined"),
               c("ePROseq_Pombe_Fisher_dis2_T316A_COMBINED_pombe_plus.bw", "ePROseq_Pombe_Fisher_dis2_T316A_COMBINED_pombe_minus.bw", "Fisher_dis2_T316A_combined"),
               c("ePROseq_Pombe_Fisher_dis2_T316D_COMBINED_pombe_plus.bw", "ePROseq_Pombe_Fisher_dis2_T316D_COMBINED_pombe_minus.bw", "Fisher_dis2_T316D_combined"))

write_PIfiles = function(wigset = wigset){
  N = dim(wigset)[1]
  #pi.res = vector(mode="list", length=N)
  for (i in 1:N){
    cat("* loading", i, "\n")
    wigs = load.wigset(wigset, i) #this should give all the bigwig files
    cat("* computing ...\n")
    PI_results = calc.pause.index(pbw = wigs[[1]], mbw = wigs[[2]], name = wigs[[3]])
    write.table(PI_results, file = paste(count_outpath, wigs[[3]], "PIcountData.txt",sep = ""), sep = "\t", quote = F, row.names = T)
    cat("* unloading.\n")
    unload.wigset(wigs)
  }
}

write_PIfiles(wigset)

############################################################################################
# Run DESeq2 on all gene data and save dataframes containing logFC

library(DESeq2)
dir.create("/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/countData/AllGenesObsTSS/DESeqOutput/")
DEseq_outpath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/countData/AllGenesObsTSS/DESeqOutput/"
infopath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/pombeinfo/"
countpath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/countData/AllGenesObsTSS/"

# 1st experiment
WT_D1 = read.table(file = paste(countpath, "WT_DMSOr1PIcountData.txt", sep = ""), head = T)
WT_D2 = read.table(file = paste(countpath, "WT_DMSOr2PIcountData.txt", sep = ""), head = T)
WT_M1 = read.table(file = paste(countpath, "WT_3MBPP1r1PIcountData.txt", sep = ""), head = T)
WT_M2 = read.table(file = paste(countpath, "WT_3MBPP1r2PIcountData.txt", sep = ""), head = T)
MCS6_D1 = read.table(file = paste(countpath, "mcs6as_DMSOr1PIcountData.txt", sep = ""), head = T)
MCS6_D2 = read.table(file = paste(countpath, "mcs6as_DMSOr2PIcountData.txt", sep = ""), head = T)
MCS6_M1 = read.table(file = paste(countpath, "mcs6as_3MBPP1r1PIcountData.txt", sep = ""), head = T)
MCS6_M2 = read.table(file = paste(countpath, "mcs6as_3MBPP1r2PIcountData.txt", sep = ""), head = T)
CDK9_D1 = read.table(file = paste(countpath, "Cdk9as_DMSOr1PIcountData.txt", sep = ""), head = T)
CDK9_D2 = read.table(file = paste(countpath, "Cdk9as_DMSOr2PIcountData.txt", sep = ""), head = T)
CDK9_M1 = read.table(file = paste(countpath, "Cdk9as_3MBPP1r1PIcountData.txt", sep = ""), head = T)
CDK9_M2 = read.table(file = paste(countpath, "Cdk9as_3MBPP1r2PIcountData.txt", sep = ""), head = T)
# 2nd experiment
Isk1as_D1 = read.table(file = paste(countpath, "Isk1as_DMSOr1PIcountData.txt", sep = ""), head = T)
Isk1as_D2 = read.table(file = paste(countpath, "Isk1as_DMSOr2PIcountData.txt", sep = ""), head = T)
Isk1as_M1 = read.table(file = paste(countpath, "Isk1as_3MBPP1r1PIcountData.txt", sep = ""), head = T)
Isk1as_M2 = read.table(file = paste(countpath, "Isk1as_3MBPP1r2PIcountData.txt", sep = ""), head = T)
mcs6as_Cdk9as_D1 = read.table(file = paste(countpath, "Mcs6as_Cdk9as_DMSOr1PIcountData.txt", sep = ""), head = T)
mcs6as_Cdk9as_D2 = read.table(file = paste(countpath, "Mcs6as_Cdk9as_DMSOr2PIcountData.txt", sep = ""), head = T)
mcs6as_Cdk9as_M1 = read.table(file = paste(countpath, "Mcs6as_Cdk9as_3MBPP1r1PIcountData.txt", sep = ""), head = T)
mcs6as_Cdk9as_M2 = read.table(file = paste(countpath, "Mcs6as_Cdk9as_3MBPP1r2PIcountData.txt", sep = ""), head = T)
# 3rd experiment
CDK9_18C_D1 = read.table(file = paste(countpath, "Cdk9as_DMSO_18Cr1PIcountData.txt", sep = ""), head = T)
CDK9_18C_D2 = read.table(file = paste(countpath, "Cdk9as_DMSO_18Cr2PIcountData.txt", sep = ""), head = T)
CDK9_18C_M1 = read.table(file = paste(countpath, "Cdk9as_3MBPP1_18Cr1PIcountData.txt", sep = ""), head = T)
CDK9_18C_M2 = read.table(file = paste(countpath, "Cdk9as_3MBPP1_18Cr2PIcountData.txt", sep = ""), head = T)
Dis2ts_18_D1 = read.table(file = paste(countpath, "dis2ts_DMSO_18Cr1PIcountData.txt", sep = ""), head = T)
Dis2ts_18_D2 = read.table(file = paste(countpath, "dis2ts_DMSO_18Cr2PIcountData.txt", sep = ""), head = T)
Dis2ts_30_D1 = read.table(file = paste(countpath, "dis2ts_DMSO_30Cr1PIcountData.txt", sep = ""), head = T)
Dis2ts_30_D2 = read.table(file = paste(countpath, "dis2ts_DMSO_30Cr2PIcountData.txt", sep = ""), head = T)
Dis2ts_18_M1 = read.table(file = paste(countpath, "dis2ts_3MBPP1_18Cr1PIcountData.txt", sep = ""), head = T)
Dis2ts_18_M2 = read.table(file = paste(countpath, "dis2ts_3MBPP1_18Cr2PIcountData.txt", sep = ""), head = T)
Dis2ts_30_M1 = read.table(file = paste(countpath, "dis2ts_3MBPP1_30Cr1PIcountData.txt", sep = ""), head = T)
Dis2ts_30_M2 = read.table(file = paste(countpath, "dis2ts_3MBPP1_30Cr2PIcountData.txt", sep = ""), head = T)
Cdk9as_dis2ts_18_D1 = read.table(file = paste(countpath, "CDK9as_dis2ts_DMSO_18Cr1PIcountData.txt", sep = ""), head = T)
Cdk9as_dis2ts_18_D2 = read.table(file = paste(countpath, "CDK9as_dis2ts_DMSO_18Cr2PIcountData.txt", sep = ""), head = T)
Cdk9as_dis2ts_30_D1 = read.table(file = paste(countpath, "CDK9as_dis2ts_DMSO_30Cr1PIcountData.txt", sep = ""), head = T)
Cdk9as_dis2ts_30_D2 = read.table(file = paste(countpath, "CDK9as_dis2ts_DMSO_30Cr2PIcountData.txt", sep = ""), head = T)
Cdk9as_dis2ts_18_M1 = read.table(file = paste(countpath, "CDK9as_dis2ts_3MBPP1_18Cr1PIcountData.txt", sep = ""), head = T)
Cdk9as_dis2ts_18_M2 = read.table(file = paste(countpath, "CDK9as_dis2ts_3MBPP1_18Cr2PIcountData.txt", sep = ""), head = T)
Cdk9as_dis2ts_30_M1 = read.table(file = paste(countpath, "CDK9as_dis2ts_3MBPP1_30Cr1PIcountData.txt", sep = ""), head = T)
Cdk9as_dis2ts_30_M2 = read.table(file = paste(countpath, "CDK9as_dis2ts_3MBPP1_30Cr2PIcountData.txt", sep = ""), head = T)
# 4th experiment: Cdk9as time course (Full timecourse with new data)
ePRO_Cdk9as_0min_r1 = read.table(file = paste(countpath, "ePRO_Cdk9as_0min_r1PIcountData.txt", sep = ""), head = T)
ePRO_Cdk9as_0min_r2 = read.table(file = paste(countpath, "ePRO_Cdk9as_0min_r2PIcountData.txt", sep = ""), head = T)
ePRO_Cdk9as_30sec_r1 = read.table(file = paste(countpath, "ePRO_Cdk9as_30sec_r1PIcountData.txt", sep = ""), head = T)
ePRO_Cdk9as_30sec_r2 = read.table(file = paste(countpath, "ePRO_Cdk9as_30sec_r2PIcountData.txt", sep = ""), head = T)
ePRO_Cdk9as_1min_r1 = read.table(file = paste(countpath, "ePRO_Cdk9as_1min_r1PIcountData.txt", sep = ""), head = T)
ePRO_Cdk9as_1min_r2 = read.table(file = paste(countpath, "ePRO_Cdk9as_1min_r2PIcountData.txt", sep = ""), head = T)
ePRO_Cdk9as_2min30sec_r1 = read.table(file = paste(countpath, "ePRO_Cdk9as_2min30sec_r1PIcountData.txt", sep = ""), head = T)
ePRO_Cdk9as_2min30sec_r2 = read.table(file = paste(countpath, "ePRO_Cdk9as_2min30sec_r2PIcountData.txt", sep = ""), head = T)
ePRO_Cdk9as_5min_r1 = read.table(file = paste(countpath, "ePRO_Cdk9as_5min_r1PIcountData.txt", sep = ""), head = T)
ePRO_Cdk9as_5min_r2 = read.table(file = paste(countpath, "ePRO_Cdk9as_5min_r2PIcountData.txt", sep = ""), head = T)
ePRO_Cdk9as_7min30sec_r1 = read.table(file = paste(countpath, "ePRO_Cdk9as_7min30sec_r1PIcountData.txt", sep = ""), head = T)
ePRO_Cdk9as_7min30sec_r2 = read.table(file = paste(countpath, "ePRO_Cdk9as_7min30sec_r2PIcountData.txt", sep = ""), head = T)
ePRO_Cdk9as_10min_r1 = read.table(file = paste(countpath, "ePRO_Cdk9as_10min_r1PIcountData.txt", sep = ""), head = T)
ePRO_Cdk9as_10min_r2 = read.table(file = paste(countpath, "ePRO_Cdk9as_10min_r2PIcountData.txt", sep = ""), head = T)
ePRO_Cdk9as_20min_r1 = read.table(file = paste(countpath, "ePRO_Cdk9as_20min_r1PIcountData.txt", sep = ""), head = T)
ePRO_Cdk9as_20min_r2 = read.table(file = paste(countpath, "ePRO_Cdk9as_20min_r2PIcountData.txt", sep = ""), head = T)
# Lsk1as_Cdk9as  double mutant
Lsk1as_Cdk9as_D1 = read.table(file = paste(countpath, "ePRO_Lsk1as_Cdk9as_DMSOr1PIcountData.txt", sep = ""), head = T)
Lsk1as_Cdk9as_D2 = read.table(file = paste(countpath, "ePRO_Lsk1as_Cdk9as_DMSOr2PIcountData.txt", sep = ""), head = T)
Lsk1as_Cdk9as_M1 = read.table(file = paste(countpath, "ePRO_Lsk1as_Cdk9as_3MBPP1r1PIcountData.txt", sep = ""), head = T)
Lsk1as_Cdk9as_M2 = read.table(file = paste(countpath, "ePRO_Lsk1as_Cdk9as_3MBPP1r2PIcountData.txt", sep = ""), head = T)
# 5th experiment: Cdk9as Time course with different Spt5 CTR mutants
Cdk9as_spt5_WT7_0min_DMSO_r1 = read.table(file = paste(countpath, "Cdk9as_spt5_WT7_0min_DMSO_r1PIcountData.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_WT7_0min_DMSO_r2 = read.table(file = paste(countpath, "Cdk9as_spt5_WT7_0min_DMSO_r2PIcountData.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_WT7_1min_3MBPP1_r1 = read.table(file = paste(countpath, "Cdk9as_spt5_WT7_1min_3MBPP1_r1PIcountData.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_WT7_1min_3MBPP1_r2 = read.table(file = paste(countpath, "Cdk9as_spt5_WT7_1min_3MBPP1_r2PIcountData.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_WT7_5min_3MBPP1_r1 = read.table(file = paste(countpath, "Cdk9as_spt5_WT7_5min_3MBPP1_r1PIcountData.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_WT7_5min_3MBPP1_r2 = read.table(file = paste(countpath, "Cdk9as_spt5_WT7_5min_3MBPP1_r2PIcountData.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_T1A_0min_DMSO_r1 = read.table(file = paste(countpath, "Cdk9as_spt5_T1A_0min_DMSO_r1PIcountData.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_T1A_0min_DMSO_r2 = read.table(file = paste(countpath, "Cdk9as_spt5_T1A_0min_DMSO_r2PIcountData.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_T1A_1min_3MBPP1_r1 = read.table(file = paste(countpath, "Cdk9as_spt5_T1A_1min_3MBPP1_r1PIcountData.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_T1A_1min_3MBPP1_r2 = read.table(file = paste(countpath, "Cdk9as_spt5_T1A_1min_3MBPP1_r2PIcountData.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_T1A_5min_3MBPP1_r1 = read.table(file = paste(countpath, "Cdk9as_spt5_T1A_5min_3MBPP1_r1PIcountData.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_T1A_5min_3MBPP1_r2 = read.table(file = paste(countpath, "Cdk9as_spt5_T1A_5min_3MBPP1_r2PIcountData.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_T1E_0min_DMSO_r1 = read.table(file = paste(countpath, "Cdk9as_spt5_T1E_0min_DMSO_r1PIcountData.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_T1E_0min_DMSO_r2 = read.table(file = paste(countpath, "Cdk9as_spt5_T1E_0min_DMSO_r2PIcountData.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_T1E_1min_3MBPP1_r1 = read.table(file = paste(countpath, "Cdk9as_spt5_T1E_1min_3MBPP1_r1PIcountData.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_T1E_1min_3MBPP1_r2 = read.table(file = paste(countpath, "Cdk9as_spt5_T1E_1min_3MBPP1_r2PIcountData.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_T1E_5min_3MBPP1_r1 = read.table(file = paste(countpath, "Cdk9as_spt5_T1E_5min_3MBPP1_r1PIcountData.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_T1E_5min_3MBPP1_r2 = read.table(file = paste(countpath, "Cdk9as_spt5_T1E_5min_3MBPP1_r2PIcountData.txt", sep = ""), head = T, row.names=1)
# 6th experiment: Comparing variuos Dis2 mutants to WT (Parua, Fisher et al.)
Fisher_WT_r1 = read.table(file = paste(countpath, "Fisher_WT_r1PIcountData.txt", sep = ""), head = T, row.names=1)
Fisher_WT_r2 = read.table(file = paste(countpath, "Fisher_WT_r2PIcountData.txt", sep = ""), head = T, row.names=1)
Fisher_dis2_11_r1 = read.table(file = paste(countpath, "Fisher_dis2_11_r1PIcountData.txt", sep = ""), head = T, row.names=1)
Fisher_dis2_11_r2 = read.table(file = paste(countpath, "Fisher_dis2_11_r2PIcountData.txt", sep = ""), head = T, row.names=1)
Fisher_dis2Del_r1 = read.table(file = paste(countpath, "Fisher_dis2Del_r1PIcountData.txt", sep = ""), head = T, row.names=1)
Fisher_dis2Del_r2 = read.table(file = paste(countpath, "Fisher_dis2Del_r2PIcountData.txt", sep = ""), head = T, row.names=1)
Fisher_dis2_T316A_r1 = read.table(file = paste(countpath, "Fisher_dis2_T316A_r1PIcountData.txt", sep = ""), head = T, row.names=1)
Fisher_dis2_T316A_r2 = read.table(file = paste(countpath, "Fisher_dis2_T316A_r2PIcountData.txt", sep = ""), head = T, row.names=1)
Fisher_dis2_T316D_r1 = read.table(file = paste(countpath, "Fisher_dis2_T316D_r1PIcountData.txt", sep = ""), head = T, row.names=1)
Fisher_dis2_T316D_r2 = read.table(file = paste(countpath, "Fisher_dis2_T316D_r2PIcountData.txt", sep = ""), head = T, row.names=1)

df_List = list(WT_D1, WT_D2, WT_M1, WT_M2, MCS6_D1, MCS6_D2, MCS6_M1, MCS6_M2, CDK9_D1, CDK9_D2, CDK9_M1, CDK9_M2,
               Isk1as_D1, Isk1as_D2, Isk1as_M1, Isk1as_M2, mcs6as_Cdk9as_D1, mcs6as_Cdk9as_D2, mcs6as_Cdk9as_M1, mcs6as_Cdk9as_M2,
               CDK9_18C_D1, CDK9_18C_D2, CDK9_18C_M1, CDK9_18C_M2, Dis2ts_18_D1, Dis2ts_18_D2, Dis2ts_18_M1, Dis2ts_18_M2, 
               Dis2ts_30_D1, Dis2ts_30_D2, Dis2ts_30_M1, Dis2ts_30_M2, Cdk9as_dis2ts_18_D1, Cdk9as_dis2ts_18_D2, Cdk9as_dis2ts_18_M1, Cdk9as_dis2ts_18_M2,
               Cdk9as_dis2ts_30_D1, Cdk9as_dis2ts_30_D2, Cdk9as_dis2ts_30_M1, Cdk9as_dis2ts_30_M2, 
               ePRO_Cdk9as_0min_r1, ePRO_Cdk9as_0min_r2, ePRO_Cdk9as_30sec_r1,ePRO_Cdk9as_30sec_r2, ePRO_Cdk9as_1min_r1, ePRO_Cdk9as_1min_r2, 
               ePRO_Cdk9as_2min30sec_r1, ePRO_Cdk9as_2min30sec_r2, ePRO_Cdk9as_5min_r1, ePRO_Cdk9as_5min_r2, ePRO_Cdk9as_7min30sec_r1, ePRO_Cdk9as_7min30sec_r2,
               ePRO_Cdk9as_10min_r1, ePRO_Cdk9as_10min_r2, ePRO_Cdk9as_20min_r1, ePRO_Cdk9as_20min_r2, Lsk1as_Cdk9as_D1, Lsk1as_Cdk9as_D2, Lsk1as_Cdk9as_M1, Lsk1as_Cdk9as_M2,
               Cdk9as_spt5_WT7_0min_DMSO_r1, Cdk9as_spt5_WT7_0min_DMSO_r2, Cdk9as_spt5_WT7_1min_3MBPP1_r1, Cdk9as_spt5_WT7_1min_3MBPP1_r2, Cdk9as_spt5_WT7_5min_3MBPP1_r1,
               Cdk9as_spt5_WT7_5min_3MBPP1_r2, Cdk9as_spt5_T1A_0min_DMSO_r1, Cdk9as_spt5_T1A_0min_DMSO_r2, Cdk9as_spt5_T1A_1min_3MBPP1_r1, Cdk9as_spt5_T1A_1min_3MBPP1_r2, 
               Cdk9as_spt5_T1A_5min_3MBPP1_r1,Cdk9as_spt5_T1A_5min_3MBPP1_r2, Cdk9as_spt5_T1E_0min_DMSO_r1, Cdk9as_spt5_T1E_0min_DMSO_r2, Cdk9as_spt5_T1E_1min_3MBPP1_r1, 
               Cdk9as_spt5_T1E_1min_3MBPP1_r2, Cdk9as_spt5_T1E_5min_3MBPP1_r1, Cdk9as_spt5_T1E_5min_3MBPP1_r2, Fisher_WT_r1, Fisher_WT_r2, Fisher_dis2_11_r1, Fisher_dis2_11_r2, 
               Fisher_dis2Del_r1, Fisher_dis2Del_r2, Fisher_dis2_T316A_r1, Fisher_dis2_T316A_r2, Fisher_dis2_T316D_r1, Fisher_dis2_T316D_r2)
SampleList = c("WT_DMSOr1", "WT_DMSOr2", "WT_3MBPP1r1", "WT_3MBPP1r2", "MCS6_DMSOr1", "MCS6_DMSOr2", "MCS6_3MBPP1r1", "MCS6_3MBPP1r2", "CDK9_DMSOr1", "CDK9_DMSOr2", "CDK9_3MBPP1r1", "CDK9_3MBPP1r2",
               "Isk1as_DMSOr1", "Isk1as_DMSOr2", "Isk1as_3MBPP1r1", "Isk1as_3MBPP1r2", "mcs6as_Cdk9as_DMSOr1", "mcs6as_Cdk9as_DMSOr2", "mcs6as_Cdk9as_3MBPP1r1", "mcs6as_Cdk9as_3MBPP1r2",
               "CDK9_18C_DMSOr1", "CDK9_18C_DMSOr2", "CDK9_18C_3MBPP1r1", "CDK9_18C_3MBPP1r2", "Dis2ts_18_DMSOr1", "Dis2ts_18_DMSOr2", "Dis2ts_18_3MBPP1r1", "Dis2ts_18_3MBPP1r2", 
               "Dis2ts_30_DMSOr1", "Dis2ts_30_DMSOr2", "Dis2ts_30_3MBPP1r1", "Dis2ts_30_3MBPP1r2", "Cdk9as_dis2ts_18_DMSOr1", "Cdk9as_dis2ts_18_DMSOr2", "Cdk9as_dis2ts_18_3MBPP1r1", "Cdk9as_dis2ts_18_3MBPP1r2",
               "Cdk9as_dis2ts_30_DMSOr1", "Cdk9as_dis2ts_30_DMSOr2", "Cdk9as_dis2ts_30_3MBPP1r1", "Cdk9as_dis2ts_30_3MBPP1r2",
               "ePRO_Cdk9as_0min_r1", "ePRO_Cdk9as_0min_r2", "ePRO_Cdk9as_30sec_r1", "ePRO_Cdk9as_30sec_r2", "ePRO_Cdk9as_1min_r1", "ePRO_Cdk9as_1min_r2", 
               "ePRO_Cdk9as_2min30sec_r1", "ePRO_Cdk9as_2min30sec_r2", "ePRO_Cdk9as_5min_r1", "ePRO_Cdk9as_5min_r2", "ePRO_Cdk9as_7min30sec_r1", "ePRO_Cdk9as_7min30sec_r2",
               "ePRO_Cdk9as_10min_r1", "ePRO_Cdk9as_10min_r2", "ePRO_Cdk9as_20min_r1", "ePRO_Cdk9as_20min_r2", "Lsk1as_Cdk9as_DMSOr1", "Lsk1as_Cdk9as_DMSOr2", "Lsk1as_Cdk9as_3MBPP1r1", "Lsk1as_Cdk9as_3MBPP1r2",
               "Cdk9as_spt5_WT7_0min_DMSO_r1", "Cdk9as_spt5_WT7_0min_DMSO_r2", "Cdk9as_spt5_WT7_1min_3MBPP1_r1", "Cdk9as_spt5_WT7_1min_3MBPP1_r2", "Cdk9as_spt5_WT7_5min_3MBPP1_r1",
               "Cdk9as_spt5_WT7_5min_3MBPP1_r2", "Cdk9as_spt5_T1A_0min_DMSO_r1", "Cdk9as_spt5_T1A_0min_DMSO_r2", "Cdk9as_spt5_T1A_1min_3MBPP1_r1", "Cdk9as_spt5_T1A_1min_3MBPP1_r2", 
               "Cdk9as_spt5_T1A_5min_3MBPP1_r1", "Cdk9as_spt5_T1A_5min_3MBPP1_r2", "Cdk9as_spt5_T1E_0min_DMSO_r1", "Cdk9as_spt5_T1E_0min_DMSO_r2", "Cdk9as_spt5_T1E_1min_3MBPP1_r1", 
               "Cdk9as_spt5_T1E_1min_3MBPP1_r2", "Cdk9as_spt5_T1E_5min_3MBPP1_r1", "Cdk9as_spt5_T1E_5min_3MBPP1_r2", "Fisher_WT_r1", "Fisher_WT_r2", "Fisher_dis2_11_r1", "Fisher_dis2_11_r2", 
               "Fisher_dis2Del_r1", "Fisher_dis2Del_r2", "Fisher_dis2_T316A_r1", "Fisher_dis2_T316A_r2", "Fisher_dis2_T316D_r1", "Fisher_dis2_T316D_r2")
rewrite_CountTables = function(dflist, samplelist, geneList){
  mTable = merge(geneList[,c(1,4)], dflist[[1]][,c(2,5)], by.x = "V4", by.y = "row.names")
  for (i in 2:length(dflist)){
    mTable = merge(mTable, dflist[[i]][,c(2,5)], by.x = "V4", by.y = "row.names")
  }
  mTableFin = mTable[complete.cases(mTable),c(1,2)]
  for (i in 1:length(dflist)){
    reduceTable = merge(mTableFin, dflist[[i]], by.x = "V4", by.y = "row.names")[,c(1,seq(3,20))]
    write.table(reduceTable, file = paste(countpath, samplelist[[i]], "PIcountData_Fix.txt",sep = ""), sep = "\t", quote = F, row.names = F)
  }
}
rewrite_CountTables(df_List, samplelist = SampleList, AllGenesObs)
#########################################################################################################
# over-write count tables with fixed lists (all have same genes in same order)
# 1st experiment
WT_D1 = read.table(file = paste(countpath, "WT_DMSOr1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
WT_D2 = read.table(file = paste(countpath, "WT_DMSOr2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
WT_M1 = read.table(file = paste(countpath, "WT_3MBPP1r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
WT_M2 = read.table(file = paste(countpath, "WT_3MBPP1r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
MCS6_D1 = read.table(file = paste(countpath, "MCS6_DMSOr1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
MCS6_D2 = read.table(file = paste(countpath, "MCS6_DMSOr2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
MCS6_M1 = read.table(file = paste(countpath, "MCS6_3MBPP1r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
MCS6_M2 = read.table(file = paste(countpath, "MCS6_3MBPP1r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
CDK9_D1 = read.table(file = paste(countpath, "CDK9_DMSOr1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
CDK9_D2 = read.table(file = paste(countpath, "CDK9_DMSOr2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
CDK9_M1 = read.table(file = paste(countpath, "CDK9_3MBPP1r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
CDK9_M2 = read.table(file = paste(countpath, "CDK9_3MBPP1r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
# 2nd experiment
Isk1as_D1 = read.table(file = paste(countpath, "Isk1as_DMSOr1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Isk1as_D2 = read.table(file = paste(countpath, "Isk1as_DMSOr2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Isk1as_M1 = read.table(file = paste(countpath, "Isk1as_3MBPP1r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Isk1as_M2 = read.table(file = paste(countpath, "Isk1as_3MBPP1r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
mcs6as_Cdk9as_D1 = read.table(file = paste(countpath, "Mcs6as_Cdk9as_DMSOr1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
mcs6as_Cdk9as_D2 = read.table(file = paste(countpath, "Mcs6as_Cdk9as_DMSOr2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
mcs6as_Cdk9as_M1 = read.table(file = paste(countpath, "Mcs6as_Cdk9as_3MBPP1r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
mcs6as_Cdk9as_M2 = read.table(file = paste(countpath, "Mcs6as_Cdk9as_3MBPP1r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
# 3rd experiment
CDK9_18C_D1 = read.table(file = paste(countpath, "CDK9_18C_DMSOr1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
CDK9_18C_D2 = read.table(file = paste(countpath, "CDK9_18C_DMSOr2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
CDK9_18C_M1 = read.table(file = paste(countpath, "CDK9_18C_3MBPP1r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
CDK9_18C_M2 = read.table(file = paste(countpath, "CDK9_18C_3MBPP1r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Dis2ts_18_D1 = read.table(file = paste(countpath, "Dis2ts_18_DMSOr1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Dis2ts_18_D2 = read.table(file = paste(countpath, "Dis2ts_18_DMSOr2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Dis2ts_30_D1 = read.table(file = paste(countpath, "Dis2ts_30_DMSOr1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Dis2ts_30_D2 = read.table(file = paste(countpath, "Dis2ts_30_DMSOr2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Dis2ts_18_M1 = read.table(file = paste(countpath, "Dis2ts_18_3MBPP1r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Dis2ts_18_M2 = read.table(file = paste(countpath, "Dis2ts_18_3MBPP1r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Dis2ts_30_M1 = read.table(file = paste(countpath, "Dis2ts_30_3MBPP1r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Dis2ts_30_M2 = read.table(file = paste(countpath, "Dis2ts_30_3MBPP1r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Cdk9as_dis2ts_18_D1 = read.table(file = paste(countpath, "Cdk9as_dis2ts_18_DMSOr1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Cdk9as_dis2ts_18_D2 = read.table(file = paste(countpath, "Cdk9as_dis2ts_18_DMSOr2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Cdk9as_dis2ts_30_D1 = read.table(file = paste(countpath, "Cdk9as_dis2ts_30_DMSOr1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Cdk9as_dis2ts_30_D2 = read.table(file = paste(countpath, "Cdk9as_dis2ts_30_DMSOr2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Cdk9as_dis2ts_18_M1 = read.table(file = paste(countpath, "Cdk9as_dis2ts_18_3MBPP1r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Cdk9as_dis2ts_18_M2 = read.table(file = paste(countpath, "Cdk9as_dis2ts_18_3MBPP1r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Cdk9as_dis2ts_30_M1 = read.table(file = paste(countpath, "Cdk9as_dis2ts_30_3MBPP1r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Cdk9as_dis2ts_30_M2 = read.table(file = paste(countpath, "Cdk9as_dis2ts_30_3MBPP1r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
# 4th experiment: Cdk9as Time course (New data: Full time course)
ePRO_Cdk9as_0min_r1 = read.table(file = paste(countpath, "ePRO_Cdk9as_0min_r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
ePRO_Cdk9as_0min_r2 = read.table(file = paste(countpath, "ePRO_Cdk9as_0min_r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
ePRO_Cdk9as_30sec_r1 = read.table(file = paste(countpath, "ePRO_Cdk9as_30sec_r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
ePRO_Cdk9as_30sec_r2 = read.table(file = paste(countpath, "ePRO_Cdk9as_30sec_r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
ePRO_Cdk9as_1min_r1 = read.table(file = paste(countpath, "ePRO_Cdk9as_1min_r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
ePRO_Cdk9as_1min_r2 = read.table(file = paste(countpath, "ePRO_Cdk9as_1min_r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
ePRO_Cdk9as_2min30sec_r1 = read.table(file = paste(countpath, "ePRO_Cdk9as_2min30sec_r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
ePRO_Cdk9as_2min30sec_r2 = read.table(file = paste(countpath, "ePRO_Cdk9as_2min30sec_r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
ePRO_Cdk9as_5min_r1 = read.table(file = paste(countpath, "ePRO_Cdk9as_5min_r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
ePRO_Cdk9as_5min_r2 = read.table(file = paste(countpath, "ePRO_Cdk9as_5min_r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
ePRO_Cdk9as_7min30sec_r1 = read.table(file = paste(countpath, "ePRO_Cdk9as_7min30sec_r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
ePRO_Cdk9as_7min30sec_r2 = read.table(file = paste(countpath, "ePRO_Cdk9as_7min30sec_r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
ePRO_Cdk9as_10min_r1 = read.table(file = paste(countpath, "ePRO_Cdk9as_10min_r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
ePRO_Cdk9as_10min_r2 = read.table(file = paste(countpath, "ePRO_Cdk9as_10min_r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
ePRO_Cdk9as_20min_r1 = read.table(file = paste(countpath, "ePRO_Cdk9as_20min_r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
ePRO_Cdk9as_20min_r2 = read.table(file = paste(countpath, "ePRO_Cdk9as_20min_r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
# Lsk1as_Cdk9as  double mutant
Lsk1as_Cdk9as_D1 = read.table(file = paste(countpath, "Lsk1as_Cdk9as_DMSOr1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Lsk1as_Cdk9as_D2 = read.table(file = paste(countpath, "Lsk1as_Cdk9as_DMSOr2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Lsk1as_Cdk9as_M1 = read.table(file = paste(countpath, "Lsk1as_Cdk9as_3MBPP1r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Lsk1as_Cdk9as_M2 = read.table(file = paste(countpath, "Lsk1as_Cdk9as_3MBPP1r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
# 5th experiment: Cdk9as Time course with different Spt5 CTR mutants
Cdk9as_spt5_WT7_0min_DMSO_r1 = read.table(file = paste(countpath, "Cdk9as_spt5_WT7_0min_DMSO_r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_WT7_0min_DMSO_r2 = read.table(file = paste(countpath, "Cdk9as_spt5_WT7_0min_DMSO_r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_WT7_1min_3MBPP1_r1 = read.table(file = paste(countpath, "Cdk9as_spt5_WT7_1min_3MBPP1_r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_WT7_1min_3MBPP1_r2 = read.table(file = paste(countpath, "Cdk9as_spt5_WT7_1min_3MBPP1_r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_WT7_5min_3MBPP1_r1 = read.table(file = paste(countpath, "Cdk9as_spt5_WT7_5min_3MBPP1_r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_WT7_5min_3MBPP1_r2 = read.table(file = paste(countpath, "Cdk9as_spt5_WT7_5min_3MBPP1_r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_T1A_0min_DMSO_r1 = read.table(file = paste(countpath, "Cdk9as_spt5_T1A_0min_DMSO_r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_T1A_0min_DMSO_r2 = read.table(file = paste(countpath, "Cdk9as_spt5_T1A_0min_DMSO_r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_T1A_1min_3MBPP1_r1 = read.table(file = paste(countpath, "Cdk9as_spt5_T1A_1min_3MBPP1_r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_T1A_1min_3MBPP1_r2 = read.table(file = paste(countpath, "Cdk9as_spt5_T1A_1min_3MBPP1_r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_T1A_5min_3MBPP1_r1 = read.table(file = paste(countpath, "Cdk9as_spt5_T1A_5min_3MBPP1_r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_T1A_5min_3MBPP1_r2 = read.table(file = paste(countpath, "Cdk9as_spt5_T1A_5min_3MBPP1_r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_T1E_0min_DMSO_r1 = read.table(file = paste(countpath, "Cdk9as_spt5_T1E_0min_DMSO_r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_T1E_0min_DMSO_r2 = read.table(file = paste(countpath, "Cdk9as_spt5_T1E_0min_DMSO_r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_T1E_1min_3MBPP1_r1 = read.table(file = paste(countpath, "Cdk9as_spt5_T1E_1min_3MBPP1_r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_T1E_1min_3MBPP1_r2 = read.table(file = paste(countpath, "Cdk9as_spt5_T1E_1min_3MBPP1_r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_T1E_5min_3MBPP1_r1 = read.table(file = paste(countpath, "Cdk9as_spt5_T1E_5min_3MBPP1_r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Cdk9as_spt5_T1E_5min_3MBPP1_r2 = read.table(file = paste(countpath, "Cdk9as_spt5_T1E_5min_3MBPP1_r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
# 6th experiment: Comparing variuos Dis2 mutants to WT (Parua, Fisher et al.)
Fisher_WT_r1 = read.table(file = paste(countpath, "Fisher_WT_r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Fisher_WT_r2 = read.table(file = paste(countpath, "Fisher_WT_r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Fisher_dis2_11_r1 = read.table(file = paste(countpath, "Fisher_dis2_11_r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Fisher_dis2_11_r2 = read.table(file = paste(countpath, "Fisher_dis2_11_r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Fisher_dis2Del_r1 = read.table(file = paste(countpath, "Fisher_dis2Del_r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Fisher_dis2Del_r2 = read.table(file = paste(countpath, "Fisher_dis2Del_r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Fisher_dis2_T316A_r1 = read.table(file = paste(countpath, "Fisher_dis2_T316A_r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Fisher_dis2_T316A_r2 = read.table(file = paste(countpath, "Fisher_dis2_T316A_r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Fisher_dis2_T316D_r1 = read.table(file = paste(countpath, "Fisher_dis2_T316D_r1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Fisher_dis2_T316D_r2 = read.table(file = paste(countpath, "Fisher_dis2_T316D_r2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)

#########################################################################################################
## Prepare Tables with columns of raw counts for samples you want to compare (format:  untreated, untreated, treated, treated)
# WT 
WT_gb = cbind(WT_D1[,2], WT_D2[,2], WT_M1[,2], WT_M2[,2]) 
colnames(WT_gb) <- c("D1", "D2", "M1", "M2")
rownames(WT_gb) <- rownames(WT_D1)
WT_pr = cbind(WT_D1[,1], WT_D2[,1], WT_M1[,1], WT_M2[,1])
colnames(WT_pr) <- c("D1", "D2", "M1", "M2")
rownames(WT_pr) <- rownames(WT_D1)
WT_term = cbind(WT_D1[,3], WT_D2[,3], WT_M1[,3], WT_M2[,3])
colnames(WT_term) <- c("D1", "D2", "M1", "M2")
rownames(WT_term) <- rownames(WT_D1)
# MCS6  ### Exclude MCS6_M1.  It is a huge outlier in terms of spike-in reads per pombe reads
MCS6_gb = cbind(MCS6_D1[,2], MCS6_D2[,2], MCS6_M2[,2]) 
colnames(MCS6_gb) <- c("D1", "D2", "M1")
rownames(MCS6_gb) <- rownames(WT_D1)
MCS6_pr = cbind(MCS6_D1[,1], MCS6_D2[,1], MCS6_M2[,1]) 
colnames(MCS6_pr) <- c("D1", "D2", "M1")
rownames(MCS6_pr) <- rownames(WT_D1)
MCS6_term = cbind(MCS6_D1[,3], MCS6_D2[,3], MCS6_M2[,3]) 
colnames(MCS6_term) <- c("D1", "D2", "M1")
rownames(MCS6_term) <- rownames(WT_D1)
# cdk9
CDK9_gb = cbind(CDK9_D1[,2], CDK9_D2[,2], CDK9_M1[,2], CDK9_M2[,2]) 
colnames(CDK9_gb) <- c("D1", "D2", "M1", "M2")
rownames(CDK9_gb) <- rownames(WT_D1)
CDK9_pr = cbind(CDK9_D1[,1], CDK9_D2[,1], CDK9_M1[,1], CDK9_M2[,1]) 
colnames(CDK9_pr) <- c("D1", "D2", "M1", "M2")
rownames(CDK9_pr) <- rownames(WT_D1)
CDK9_term = cbind(CDK9_D1[,3], CDK9_D2[,3], CDK9_M1[,3], CDK9_M2[,3]) 
colnames(CDK9_term) <- c("D1", "D2", "M1", "M2")
rownames(CDK9_term) <- rownames(WT_D1)
############ 2nd experiment
# Isk1
Isk1_gb = cbind(Isk1as_D1[,2], Isk1as_D2[,2], Isk1as_M1[,2], Isk1as_M2[,2]) 
colnames(Isk1_gb) <- c("D1", "D2", "M1", "M2")
rownames(Isk1_gb) <- rownames(WT_D1)
Isk1_pr = cbind(Isk1as_D1[,1], Isk1as_D2[,1], Isk1as_M1[,1], Isk1as_M2[,1]) 
colnames(Isk1_pr) <- c("D1", "D2", "M1", "M2")
rownames(Isk1_pr) <- rownames(WT_D1)
Isk1_term = cbind(Isk1as_D1[,3], Isk1as_D2[,3], Isk1as_M1[,3], Isk1as_M2[,3]) 
colnames(Isk1_term) <- c("D1", "D2", "M1", "M2")
rownames(Isk1_term) <- rownames(WT_D1)
# Mcs6_CDK9 (double Mutant)
Mcs6_CDK9_gb = cbind(mcs6as_Cdk9as_D1[,2], mcs6as_Cdk9as_D2[,2], mcs6as_Cdk9as_M1[,2], mcs6as_Cdk9as_M2[,2]) 
colnames(Mcs6_CDK9_gb) <- c("D1", "D2", "M1", "M2")
rownames(Mcs6_CDK9_gb) <- rownames(WT_D1)
Mcs6_CDK9_pr = cbind(mcs6as_Cdk9as_D1[,1], mcs6as_Cdk9as_D2[,1], mcs6as_Cdk9as_M1[,1], mcs6as_Cdk9as_M2[,1]) 
colnames(Mcs6_CDK9_pr) <- c("D1", "D2", "M1", "M2")
rownames(Mcs6_CDK9_pr) <- rownames(WT_D1)
Mcs6_CDK9_term = cbind(mcs6as_Cdk9as_D1[,3], mcs6as_Cdk9as_D2[,3], mcs6as_Cdk9as_M1[,3], mcs6as_Cdk9as_M2[,3]) 
colnames(Mcs6_CDK9_term) <- c("D1", "D2", "M1", "M2")
rownames(Mcs6_CDK9_term) <- rownames(WT_D1)

######### compare AS mutants to WT (untreated)
# WT vs Mcs6_CDK9, both untreated
WTvMcs6_gb = cbind(WT_D1[,2], WT_D2[,2], MCS6_D1[,2], MCS6_D2[,2]) 
colnames(WTvMcs6_gb) <- c("D1", "D2", "M1", "M2")
rownames(WTvMcs6_gb) <- rownames(WT_D1)
WTvMcs6_pr = cbind(WT_D1[,1], WT_D2[,1], MCS6_D1[,1], MCS6_D2[,1]) 
colnames(WTvMcs6_pr) <- c("D1", "D2", "M1", "M2")
rownames(WTvMcs6_pr) <- rownames(WT_D1)
WTvMcs6_term = cbind(WT_D1[,3], WT_D2[,3], MCS6_D1[,3], MCS6_D2[,3]) 
colnames(WTvMcs6_term) <- c("D1", "D2", "M1", "M2")
rownames(WTvMcs6_term) <- rownames(WT_D1)
# WT vs CDK9, both untreated
WTvCDK9_gb = cbind(WT_D1[,2], WT_D2[,2], CDK9_D1[,2], CDK9_D2[,2]) 
colnames(WTvCDK9_gb) <- c("D1", "D2", "M1", "M2")
rownames(WTvCDK9_gb) <- rownames(WT_D1)
WTvCDK9_pr = cbind(WT_D1[,1], WT_D2[,1], CDK9_D1[,1], CDK9_D2[,1]) 
colnames(WTvCDK9_pr) <- c("D1", "D2", "M1", "M2")
rownames(WTvCDK9_pr) <- rownames(WT_D1)
WTvCDK9_term = cbind(WT_D1[,3], WT_D2[,3], CDK9_D1[,3], CDK9_D2[,3]) 
colnames(WTvCDK9_term) <- c("D1", "D2", "M1", "M2")
rownames(WTvCDK9_term) <- rownames(WT_D1)
# WT vs Mcs6_CDK9, both untreated
WTvMcs6_CDK9_gb = cbind(WT_D1[,2], WT_D2[,2], mcs6as_Cdk9as_D1[,2], mcs6as_Cdk9as_D2[,2]) 
colnames(WTvMcs6_CDK9_gb) <- c("D1", "D2", "M1", "M2")
rownames(WTvMcs6_CDK9_gb) <- rownames(WT_D1)
WTvMcs6_CDK9_pr = cbind(WT_D1[,1], WT_D2[,1], mcs6as_Cdk9as_D1[,1], mcs6as_Cdk9as_D2[,1]) 
colnames(WTvMcs6_CDK9_pr) <- c("D1", "D2", "M1", "M2")
rownames(WTvMcs6_CDK9_pr) <- rownames(WT_D1)
WTvMcs6_CDK9_term = cbind(WT_D1[,3], WT_D2[,3], mcs6as_Cdk9as_D1[,3], mcs6as_Cdk9as_D2[,3]) 
colnames(WTvMcs6_CDK9_term) <- c("D1", "D2", "M1", "M2")
rownames(WTvMcs6_CDK9_term) <- rownames(WT_D1)
# WT vs Lsk1, both untreated
WTvIsk1_gb = cbind(WT_D1[,2], WT_D2[,2], Isk1as_D1[,2], Isk1as_D2[,2]) 
colnames(WTvIsk1_gb) <- c("D1", "D2", "M1", "M2")
rownames(WTvIsk1_gb) <- rownames(WT_D1)
WTvIsk1_pr = cbind(WT_D1[,1], WT_D2[,1], Isk1as_D1[,1], Isk1as_D2[,1]) 
colnames(WTvIsk1_pr) <- c("D1", "D2", "M1", "M2")
rownames(WTvIsk1_pr) <- rownames(WT_D1)
WTvIsk1_term = cbind(WT_D1[,3], WT_D2[,3], Isk1as_D1[,3], Isk1as_D2[,3]) 
colnames(WTvIsk1_term) <- c("D1", "D2", "M1", "M2")
rownames(WTvIsk1_term) <- rownames(WT_D1)

######### 3rd experiment
################### Within 18C treatmetns
# Cdk9as 18C 
Cdk9as_18C_gb = cbind(CDK9_18C_D1[,2], CDK9_18C_D2[,2], CDK9_18C_M1[,2], CDK9_18C_M2[,2]) 
colnames(Cdk9as_18C_gb) <- c("D1", "D2", "M1", "M2")
rownames(Cdk9as_18C_gb) <- rownames(WT_D1)
Cdk9as_18C_pr = cbind(CDK9_18C_D1[,1], CDK9_18C_D2[,1], CDK9_18C_M1[,1], CDK9_18C_M2[,1]) 
colnames(Cdk9as_18C_pr) <- c("D1", "D2", "M1", "M2")
rownames(Cdk9as_18C_pr) <- rownames(WT_D1)
Cdk9as_18C_term = cbind(CDK9_18C_D1[,3], CDK9_18C_D2[,3], CDK9_18C_M1[,3], CDK9_18C_M2[,3]) 
colnames(Cdk9as_18C_term) <- c("D1", "D2", "M1", "M2")
rownames(Cdk9as_18C_term) <- rownames(WT_D1)
# dis2ts 18C 
dis2ts_18C_gb = cbind(Dis2ts_18_D1[,2], Dis2ts_18_D2[,2], Dis2ts_18_M1[,2], Dis2ts_18_M2[,2]) 
colnames(dis2ts_18C_gb) <- c("D1", "D2", "M1", "M2")
rownames(dis2ts_18C_gb) <- rownames(WT_D1)
dis2ts_18C_pr = cbind(Dis2ts_18_D1[,1], Dis2ts_18_D2[,1], Dis2ts_18_M1[,1], Dis2ts_18_M2[,1]) 
colnames(dis2ts_18C_pr) <- c("D1", "D2", "M1", "M2")
rownames(dis2ts_18C_pr) <- rownames(WT_D1)
dis2ts_18C_term = cbind(Dis2ts_18_D1[,3], Dis2ts_18_D2[,3], Dis2ts_18_M1[,3], Dis2ts_18_M2[,3]) 
colnames(dis2ts_18C_term) <- c("D1", "D2", "M1", "M2")
rownames(dis2ts_18C_term) <- rownames(WT_D1)
# dis2ts Cdk9as 18C 
dis2ts_Cdk9as_18C_gb = cbind(Cdk9as_dis2ts_18_D1[,2], Cdk9as_dis2ts_18_D2[,2], Cdk9as_dis2ts_18_M1[,2], Cdk9as_dis2ts_18_M2[,2]) 
colnames(dis2ts_Cdk9as_18C_gb) <- c("D1", "D2", "M1", "M2")
rownames(dis2ts_Cdk9as_18C_gb) <- rownames(WT_D1)
dis2ts_Cdk9as_18C_pr = cbind(Cdk9as_dis2ts_18_D1[,1], Cdk9as_dis2ts_18_D2[,1], Cdk9as_dis2ts_18_M1[,1], Cdk9as_dis2ts_18_M2[,1]) 
colnames(dis2ts_Cdk9as_18C_pr) <- c("D1", "D2", "M1", "M2")
rownames(dis2ts_Cdk9as_18C_pr) <- rownames(WT_D1)
dis2ts_Cdk9as_18C_term = cbind(Cdk9as_dis2ts_18_D1[,3], Cdk9as_dis2ts_18_D2[,3], Cdk9as_dis2ts_18_M1[,3], Cdk9as_dis2ts_18_M2[,3]) 
colnames(dis2ts_Cdk9as_18C_term) <- c("D1", "D2", "M1", "M2")
rownames(dis2ts_Cdk9as_18C_term) <- rownames(WT_D1)
################### Within 30C treatmetns
# dis2ts 30C 
dis2ts_30C_gb = cbind(Dis2ts_30_D1[,2], Dis2ts_30_D2[,2], Dis2ts_30_M1[,2], Dis2ts_30_M2[,2]) 
colnames(dis2ts_30C_gb) <- c("D1", "D2", "M1", "M2")
rownames(dis2ts_30C_gb) <- rownames(WT_D1)
dis2ts_30C_pr = cbind(Dis2ts_30_D1[,1], Dis2ts_30_D2[,1], Dis2ts_30_M1[,1], Dis2ts_30_M2[,1]) 
colnames(dis2ts_30C_pr) <- c("D1", "D2", "M1", "M2")
rownames(dis2ts_30C_pr) <- rownames(WT_D1)
dis2ts_30C_term = cbind(Dis2ts_30_D1[,3], Dis2ts_30_D2[,3], Dis2ts_30_M1[,3], Dis2ts_30_M2[,3]) 
colnames(dis2ts_30C_term) <- c("D1", "D2", "M1", "M2")
rownames(dis2ts_30C_term) <- rownames(WT_D1)
# dis2ts Cdk9as 30C 
dis2ts_Cdk9as_30C_gb = cbind(Cdk9as_dis2ts_30_D1[,2], Cdk9as_dis2ts_30_D2[,2], Cdk9as_dis2ts_30_M1[,2], Cdk9as_dis2ts_30_M2[,2]) 
colnames(dis2ts_Cdk9as_30C_gb) <- c("D1", "D2", "M1", "M2")
rownames(dis2ts_Cdk9as_30C_gb) <- rownames(WT_D1)
dis2ts_Cdk9as_30C_pr = cbind(Cdk9as_dis2ts_30_D1[,1], Cdk9as_dis2ts_30_D2[,1], Cdk9as_dis2ts_30_M1[,1], Cdk9as_dis2ts_30_M2[,1]) 
colnames(dis2ts_Cdk9as_30C_pr) <- c("D1", "D2", "M1", "M2")
rownames(dis2ts_Cdk9as_30C_pr) <- rownames(WT_D1)
dis2ts_Cdk9as_30C_term = cbind(Cdk9as_dis2ts_30_D1[,3], Cdk9as_dis2ts_30_D2[,3], Cdk9as_dis2ts_30_M1[,3], Cdk9as_dis2ts_30_M2[,3]) 
colnames(dis2ts_Cdk9as_30C_term) <- c("D1", "D2", "M1", "M2")
rownames(dis2ts_Cdk9as_30C_term) <- rownames(WT_D1)
################### Between 18C and 30C treatmetns
# dis2ts 30C vs 18C DMSO treated
dis2ts_30Cv18C_DMSO_gb = cbind(Dis2ts_30_D1[,2], Dis2ts_30_D2[,2], Dis2ts_18_D1[,2], Dis2ts_18_D2[,2]) 
colnames(dis2ts_30Cv18C_DMSO_gb) <- c("D1", "D2", "M1", "M2")
rownames(dis2ts_30Cv18C_DMSO_gb) <- rownames(WT_D1)
dis2ts_30Cv18C_DMSO_pr = cbind(Dis2ts_30_D1[,1], Dis2ts_30_D2[,1], Dis2ts_18_D1[,1], Dis2ts_18_D2[,1]) 
colnames(dis2ts_30Cv18C_DMSO_pr) <- c("D1", "D2", "M1", "M2")
rownames(dis2ts_30Cv18C_DMSO_pr) <- rownames(WT_D1)
dis2ts_30Cv18C_DMSO_term = cbind(Dis2ts_30_D1[,3], Dis2ts_30_D2[,3], Dis2ts_18_D1[,3], Dis2ts_18_D2[,3]) 
colnames(dis2ts_30Cv18C_DMSO_term) <- c("D1", "D2", "M1", "M2")
rownames(dis2ts_30Cv18C_DMSO_term) <- rownames(WT_D1)
# dis2ts 30C vs 18C 3MB-PP1 treated
dis2ts_30Cv18C_3MBPP1_gb = cbind(Dis2ts_30_M1[,2], Dis2ts_30_M2[,2], Dis2ts_18_M1[,2], Dis2ts_18_M2[,2]) 
colnames(dis2ts_30Cv18C_3MBPP1_gb) <- c("D1", "D2", "M1", "M2")
rownames(dis2ts_30Cv18C_3MBPP1_gb) <- rownames(WT_D1)
dis2ts_30Cv18C_3MBPP1_pr = cbind(Dis2ts_30_M1[,1], Dis2ts_30_M2[,1], Dis2ts_18_M1[,1], Dis2ts_18_M2[,1]) 
colnames(dis2ts_30Cv18C_3MBPP1_pr) <- c("D1", "D2", "M1", "M2")
rownames(dis2ts_30Cv18C_3MBPP1_pr) <- rownames(WT_D1)
dis2ts_30Cv18C_3MBPP1_term = cbind(Dis2ts_30_M1[,3], Dis2ts_30_M2[,3], Dis2ts_18_M1[,3], Dis2ts_18_M2[,3]) 
colnames(dis2ts_30Cv18C_3MBPP1_term) <- c("D1", "D2", "M1", "M2")
rownames(dis2ts_30Cv18C_3MBPP1_term) <- rownames(WT_D1)
# dis2ts Cdk9as 30C vs 18C DMSO treated
dis2ts_Cdk9as_30Cv18C_DMSO_gb = cbind(Cdk9as_dis2ts_30_D1[,2], Cdk9as_dis2ts_30_D2[,2], Cdk9as_dis2ts_18_D1[,2], Cdk9as_dis2ts_18_D2[,2]) 
colnames(dis2ts_Cdk9as_30Cv18C_DMSO_gb) <- c("D1", "D2", "M1", "M2")
rownames(dis2ts_Cdk9as_30Cv18C_DMSO_gb) <- rownames(WT_D1)
dis2ts_Cdk9as_30Cv18C_DMSO_pr = cbind(Cdk9as_dis2ts_30_D1[,1], Cdk9as_dis2ts_30_D2[,1], Cdk9as_dis2ts_18_D1[,1], Cdk9as_dis2ts_18_D2[,1]) 
colnames(dis2ts_Cdk9as_30Cv18C_DMSO_pr) <- c("D1", "D2", "M1", "M2")
rownames(dis2ts_Cdk9as_30Cv18C_DMSO_pr) <- rownames(WT_D1)
dis2ts_Cdk9as_30Cv18C_DMSO_term = cbind(Cdk9as_dis2ts_30_D1[,3], Cdk9as_dis2ts_30_D2[,3], Cdk9as_dis2ts_18_D1[,3], Cdk9as_dis2ts_18_D2[,3]) 
colnames(dis2ts_Cdk9as_30Cv18C_DMSO_term) <- c("D1", "D2", "M1", "M2")
rownames(dis2ts_Cdk9as_30Cv18C_DMSO_term) <- rownames(WT_D1)
# dis2ts Cdk9as 30C vs 18C DMSO treated
dis2ts_Cdk9as_30Cv18C_3MBPP1_gb = cbind(Cdk9as_dis2ts_30_M1[,2], Cdk9as_dis2ts_30_M2[,2], Cdk9as_dis2ts_18_M1[,2], Cdk9as_dis2ts_18_M2[,2]) 
colnames(dis2ts_Cdk9as_30Cv18C_3MBPP1_gb) <- c("D1", "D2", "M1", "M2")
rownames(dis2ts_Cdk9as_30Cv18C_3MBPP1_gb) <- rownames(WT_D1)
dis2ts_Cdk9as_30Cv18C_3MBPP1_pr = cbind(Cdk9as_dis2ts_30_M1[,1], Cdk9as_dis2ts_30_M2[,1], Cdk9as_dis2ts_18_M1[,1], Cdk9as_dis2ts_18_M2[,1]) 
colnames(dis2ts_Cdk9as_30Cv18C_3MBPP1_pr) <- c("D1", "D2", "M1", "M2")
rownames(dis2ts_Cdk9as_30Cv18C_3MBPP1_pr) <- rownames(WT_D1)
dis2ts_Cdk9as_30Cv18C_3MBPP1_term = cbind(Cdk9as_dis2ts_30_M1[,3], Cdk9as_dis2ts_30_M2[,3], Cdk9as_dis2ts_18_M1[,3], Cdk9as_dis2ts_18_M2[,3]) 
colnames(dis2ts_Cdk9as_30Cv18C_3MBPP1_term) <- c("D1", "D2", "M1", "M2")
rownames(dis2ts_Cdk9as_30Cv18C_3MBPP1_term) <- rownames(WT_D1)

#### 4th experiment: Time course of Cdk9as inhibition with 10 uM) #####
### Compare each time point/ treatment with the 0' timepoint 
# 30sec inhibition
Cdk9as_30secVs0min_gb = cbind(ePRO_Cdk9as_0min_r1[,2], ePRO_Cdk9as_0min_r2[,2], ePRO_Cdk9as_30sec_r1[,2], ePRO_Cdk9as_30sec_r2[,2]) 
colnames(Cdk9as_30secVs0min_gb) <- c("D1", "D2", "M1", "M2")
rownames(Cdk9as_30secVs0min_gb) <- rownames(WT_D1)
Cdk9as_30secVs0min_pr = cbind(ePRO_Cdk9as_0min_r1[,1], ePRO_Cdk9as_0min_r2[,1], ePRO_Cdk9as_30sec_r1[,1], ePRO_Cdk9as_30sec_r2[,1]) 
colnames(Cdk9as_30secVs0min_pr) <- c("D1", "D2", "M1", "M2")
rownames(Cdk9as_30secVs0min_pr) <- rownames(WT_D1)
Cdk9as_30secVs0min_term = cbind(ePRO_Cdk9as_0min_r1[,3], ePRO_Cdk9as_0min_r2[,3], ePRO_Cdk9as_30sec_r1[,3], ePRO_Cdk9as_30sec_r2[,3]) 
colnames(Cdk9as_30secVs0min_term) <- c("D1", "D2", "M1", "M2")
rownames(Cdk9as_30secVs0min_term) <- rownames(WT_D1)
# 1 minute inhibition
Cdk9as_1minVs0min_gb = cbind(ePRO_Cdk9as_0min_r1[,2], ePRO_Cdk9as_0min_r2[,2], ePRO_Cdk9as_1min_r1[,2], ePRO_Cdk9as_1min_r2[,2]) 
colnames(Cdk9as_1minVs0min_gb) <- c("D1", "D2", "M1", "M2")
rownames(Cdk9as_1minVs0min_gb) <- rownames(WT_D1)
Cdk9as_1minVs0min_pr = cbind(ePRO_Cdk9as_0min_r1[,1], ePRO_Cdk9as_0min_r2[,1], ePRO_Cdk9as_1min_r1[,1], ePRO_Cdk9as_1min_r2[,1]) 
colnames(Cdk9as_1minVs0min_pr) <- c("D1", "D2", "M1", "M2")
rownames(Cdk9as_1minVs0min_pr) <- rownames(WT_D1)
Cdk9as_1minVs0min_term = cbind(ePRO_Cdk9as_0min_r1[,3], ePRO_Cdk9as_0min_r2[,3], ePRO_Cdk9as_1min_r1[,3], ePRO_Cdk9as_1min_r2[,3]) 
colnames(Cdk9as_1minVs0min_term) <- c("D1", "D2", "M1", "M2")
rownames(Cdk9as_1minVs0min_term) <- rownames(WT_D1)
# 2.5 minute inhibition
Cdk9as_2.5minVs0min_gb = cbind(ePRO_Cdk9as_0min_r1[,2], ePRO_Cdk9as_0min_r2[,2], ePRO_Cdk9as_2min30sec_r1[,2], ePRO_Cdk9as_2min30sec_r2[,2]) 
colnames(Cdk9as_2.5minVs0min_gb) <- c("D1", "D2", "M1", "M2")
rownames(Cdk9as_2.5minVs0min_gb) <- rownames(WT_D1)
Cdk9as_2.5minVs0min_pr = cbind(ePRO_Cdk9as_0min_r1[,1], ePRO_Cdk9as_0min_r2[,1], ePRO_Cdk9as_2min30sec_r1[,1], ePRO_Cdk9as_2min30sec_r2[,1]) 
colnames(Cdk9as_2.5minVs0min_pr) <- c("D1", "D2", "M1", "M2")
rownames(Cdk9as_2.5minVs0min_pr) <- rownames(WT_D1)
Cdk9as_2.5minVs0min_term = cbind(ePRO_Cdk9as_0min_r1[,3], ePRO_Cdk9as_0min_r2[,3], ePRO_Cdk9as_2min30sec_r1[,3], ePRO_Cdk9as_2min30sec_r2[,3]) 
colnames(Cdk9as_2.5minVs0min_term) <- c("D1", "D2", "M1", "M2")
rownames(Cdk9as_2.5minVs0min_term) <- rownames(WT_D1)
# 5 minute inhibition
Cdk9as_5minVs0min_gb = cbind(ePRO_Cdk9as_0min_r1[,2], ePRO_Cdk9as_0min_r2[,2], ePRO_Cdk9as_5min_r1[,2], ePRO_Cdk9as_5min_r2[,2]) 
colnames(Cdk9as_5minVs0min_gb) <- c("D1", "D2", "M1", "M2")
rownames(Cdk9as_5minVs0min_gb) <- rownames(WT_D1)
Cdk9as_5minVs0min_pr = cbind(ePRO_Cdk9as_0min_r1[,1], ePRO_Cdk9as_0min_r2[,1], ePRO_Cdk9as_5min_r1[,1], ePRO_Cdk9as_5min_r2[,1]) 
colnames(Cdk9as_5minVs0min_pr) <- c("D1", "D2", "M1", "M2")
rownames(Cdk9as_5minVs0min_pr) <- rownames(WT_D1)
Cdk9as_5minVs0min_term = cbind(ePRO_Cdk9as_0min_r1[,3], ePRO_Cdk9as_0min_r2[,3], ePRO_Cdk9as_5min_r1[,3], ePRO_Cdk9as_5min_r2[,3]) 
colnames(Cdk9as_5minVs0min_term) <- c("D1", "D2", "M1", "M2")
rownames(Cdk9as_5minVs0min_term) <- rownames(WT_D1)
# 7.5 minute inhibition
Cdk9as_7.5minVs0min_gb = cbind(ePRO_Cdk9as_0min_r1[,2], ePRO_Cdk9as_0min_r2[,2], ePRO_Cdk9as_7min30sec_r1[,2], ePRO_Cdk9as_7min30sec_r2[,2]) 
colnames(Cdk9as_7.5minVs0min_gb) <- c("D1", "D2", "M1", "M2")
rownames(Cdk9as_7.5minVs0min_gb) <- rownames(WT_D1)
Cdk9as_7.5minVs0min_pr = cbind(ePRO_Cdk9as_0min_r1[,1], ePRO_Cdk9as_0min_r2[,1], ePRO_Cdk9as_7min30sec_r1[,1], ePRO_Cdk9as_7min30sec_r2[,1]) 
colnames(Cdk9as_7.5minVs0min_pr) <- c("D1", "D2", "M1", "M2")
rownames(Cdk9as_7.5minVs0min_pr) <- rownames(WT_D1)
Cdk9as_7.5minVs0min_term = cbind(ePRO_Cdk9as_0min_r1[,3], ePRO_Cdk9as_0min_r2[,3], ePRO_Cdk9as_7min30sec_r1[,3], ePRO_Cdk9as_7min30sec_r2[,3]) 
colnames(Cdk9as_7.5minVs0min_term) <- c("D1", "D2", "M1", "M2")
rownames(Cdk9as_7.5minVs0min_term) <- rownames(WT_D1)
# 10 minute inhibition
Cdk9as_10minVs0min_gb = cbind(ePRO_Cdk9as_0min_r1[,2], ePRO_Cdk9as_0min_r2[,2], ePRO_Cdk9as_10min_r1[,2], ePRO_Cdk9as_10min_r2[,2]) 
colnames(Cdk9as_10minVs0min_gb) <- c("D1", "D2", "M1", "M2")
rownames(Cdk9as_10minVs0min_gb) <- rownames(WT_D1)
Cdk9as_10minVs0min_pr = cbind(ePRO_Cdk9as_0min_r1[,1], ePRO_Cdk9as_0min_r2[,1], ePRO_Cdk9as_10min_r1[,1], ePRO_Cdk9as_10min_r2[,1]) 
colnames(Cdk9as_10minVs0min_pr) <- c("D1", "D2", "M1", "M2")
rownames(Cdk9as_10minVs0min_pr) <- rownames(WT_D1)
Cdk9as_10minVs0min_term = cbind(ePRO_Cdk9as_0min_r1[,3], ePRO_Cdk9as_0min_r2[,3], ePRO_Cdk9as_10min_r1[,3], ePRO_Cdk9as_10min_r2[,3]) 
colnames(Cdk9as_10minVs0min_term) <- c("D1", "D2", "M1", "M2")
rownames(Cdk9as_10minVs0min_term) <- rownames(WT_D1)
# 20 minute inhibition
Cdk9as_20minVs0min_gb = cbind(ePRO_Cdk9as_0min_r1[,2], ePRO_Cdk9as_0min_r2[,2], ePRO_Cdk9as_20min_r1[,2], ePRO_Cdk9as_20min_r2[,2]) 
colnames(Cdk9as_20minVs0min_gb) <- c("D1", "D2", "M1", "M2")
rownames(Cdk9as_20minVs0min_gb) <- rownames(WT_D1)
Cdk9as_20minVs0min_pr = cbind(ePRO_Cdk9as_0min_r1[,1], ePRO_Cdk9as_0min_r2[,1], ePRO_Cdk9as_20min_r1[,1], ePRO_Cdk9as_20min_r2[,1]) 
colnames(Cdk9as_20minVs0min_pr) <- c("D1", "D2", "M1", "M2")
rownames(Cdk9as_20minVs0min_pr) <- rownames(WT_D1)
Cdk9as_20minVs0min_term = cbind(ePRO_Cdk9as_0min_r1[,3], ePRO_Cdk9as_0min_r2[,3], ePRO_Cdk9as_20min_r1[,3], ePRO_Cdk9as_20min_r2[,3]) 
colnames(Cdk9as_20minVs0min_term) <- c("D1", "D2", "M1", "M2")
rownames(Cdk9as_20minVs0min_term) <- rownames(WT_D1)
# Lsk1as Cdk9as double mutant
Lsk1_CDK9_gb = cbind(Lsk1as_Cdk9as_D1[,2], Lsk1as_Cdk9as_D2[,2], Lsk1as_Cdk9as_M1[,2], Lsk1as_Cdk9as_M2[,2]) 
colnames(Lsk1_CDK9_gb) <- c("D1", "D2", "M1", "M2")
rownames(Lsk1_CDK9_gb) <- rownames(WT_D1)
Lsk1_CDK9_pr = cbind(Lsk1as_Cdk9as_D1[,1], Lsk1as_Cdk9as_D2[,1], Lsk1as_Cdk9as_M1[,1], Lsk1as_Cdk9as_M2[,1]) 
colnames(Lsk1_CDK9_pr) <- c("D1", "D2", "M1", "M2")
rownames(Lsk1_CDK9_pr) <- rownames(WT_D1)
Lsk1_CDK9_term = cbind(Lsk1as_Cdk9as_D1[,3], Lsk1as_Cdk9as_D2[,3], Lsk1as_Cdk9as_M1[,3], Lsk1as_Cdk9as_M2[,3]) 
colnames(Lsk1_CDK9_term) <- c("D1", "D2", "M1", "M2")
rownames(Lsk1_CDK9_term) <- rownames(WT_D1)

#### 5th experiment: Time course of Cdk9as inhibition with 10 uM in various spt5 CTR mutants #####
### For each strain, compare each time point/ treatment with the 0' timepoint 
# Cdk9as spt5WT7 
# 1 min
Spt5WT7_1minVs0min_gb = cbind(Cdk9as_spt5_WT7_0min_DMSO_r1[,2], Cdk9as_spt5_WT7_0min_DMSO_r2[,2], Cdk9as_spt5_WT7_1min_3MBPP1_r1[,2], Cdk9as_spt5_WT7_1min_3MBPP1_r2[,2]) 
colnames(Spt5WT7_1minVs0min_gb) <- c("D1", "D2", "M1", "M2")
rownames(Spt5WT7_1minVs0min_gb) <- rownames(WT_D1)
Spt5WT7_1minVs0min_pr = cbind(Cdk9as_spt5_WT7_0min_DMSO_r1[,1], Cdk9as_spt5_WT7_0min_DMSO_r2[,1], Cdk9as_spt5_WT7_1min_3MBPP1_r1[,1], Cdk9as_spt5_WT7_1min_3MBPP1_r2[,1]) 
colnames(Spt5WT7_1minVs0min_pr) <- c("D1", "D2", "M1", "M2")
rownames(Spt5WT7_1minVs0min_pr) <- rownames(WT_D1)
Spt5WT7_1minVs0min_term = cbind(Cdk9as_spt5_WT7_0min_DMSO_r1[,3], Cdk9as_spt5_WT7_0min_DMSO_r2[,3], Cdk9as_spt5_WT7_1min_3MBPP1_r1[,3], Cdk9as_spt5_WT7_1min_3MBPP1_r2[,3]) 
colnames(Spt5WT7_1minVs0min_term) <- c("D1", "D2", "M1", "M2")
rownames(Spt5WT7_1minVs0min_term) <- rownames(WT_D1)
# 5 min
Spt5WT7_5minVs0min_gb = cbind(Cdk9as_spt5_WT7_0min_DMSO_r1[,2], Cdk9as_spt5_WT7_0min_DMSO_r2[,2], Cdk9as_spt5_WT7_1min_3MBPP1_r1[,2], Cdk9as_spt5_WT7_1min_3MBPP1_r2[,2]) 
colnames(Spt5WT7_5minVs0min_gb) <- c("D1", "D2", "M1", "M2")
rownames(Spt5WT7_5minVs0min_gb) <- rownames(WT_D1)
Spt5WT7_5minVs0min_pr = cbind(Cdk9as_spt5_WT7_0min_DMSO_r1[,1], Cdk9as_spt5_WT7_0min_DMSO_r2[,1], Cdk9as_spt5_WT7_1min_3MBPP1_r1[,1], Cdk9as_spt5_WT7_1min_3MBPP1_r2[,1]) 
colnames(Spt5WT7_5minVs0min_pr) <- c("D1", "D2", "M1", "M2")
rownames(Spt5WT7_5minVs0min_pr) <- rownames(WT_D1)
Spt5WT7_5minVs0min_term = cbind(Cdk9as_spt5_WT7_0min_DMSO_r1[,3], Cdk9as_spt5_WT7_0min_DMSO_r2[,3], Cdk9as_spt5_WT7_1min_3MBPP1_r1[,3], Cdk9as_spt5_WT7_1min_3MBPP1_r2[,3]) 
colnames(Spt5WT7_5minVs0min_term) <- c("D1", "D2", "M1", "M2")
rownames(Spt5WT7_5minVs0min_term) <- rownames(WT_D1)
# Cdk9as spt5T1A 
# 1 min
Spt5T1A_1minVs0min_gb = cbind(Cdk9as_spt5_T1A_0min_DMSO_r1[,2], Cdk9as_spt5_T1A_0min_DMSO_r2[,2], Cdk9as_spt5_T1A_1min_3MBPP1_r1[,2], Cdk9as_spt5_T1A_1min_3MBPP1_r2[,2]) 
colnames(Spt5T1A_1minVs0min_gb) <- c("D1", "D2", "M1", "M2")
rownames(Spt5T1A_1minVs0min_gb) <- rownames(WT_D1)
Spt5T1A_1minVs0min_pr = cbind(Cdk9as_spt5_T1A_0min_DMSO_r1[,1], Cdk9as_spt5_T1A_0min_DMSO_r2[,1], Cdk9as_spt5_T1A_1min_3MBPP1_r1[,1], Cdk9as_spt5_T1A_1min_3MBPP1_r2[,1]) 
colnames(Spt5T1A_1minVs0min_pr) <- c("D1", "D2", "M1", "M2")
rownames(Spt5T1A_1minVs0min_pr) <- rownames(WT_D1)
Spt5T1A_1minVs0min_term = cbind(Cdk9as_spt5_T1A_0min_DMSO_r1[,3], Cdk9as_spt5_T1A_0min_DMSO_r2[,3], Cdk9as_spt5_T1A_1min_3MBPP1_r1[,3], Cdk9as_spt5_T1A_1min_3MBPP1_r2[,3]) 
colnames(Spt5T1A_1minVs0min_term) <- c("D1", "D2", "M1", "M2")
rownames(Spt5T1A_1minVs0min_term) <- rownames(WT_D1)
# 5 min
Spt5T1A_5minVs0min_gb = cbind(Cdk9as_spt5_T1A_0min_DMSO_r1[,2], Cdk9as_spt5_T1A_0min_DMSO_r2[,2], Cdk9as_spt5_T1A_5min_3MBPP1_r1[,2], Cdk9as_spt5_T1A_5min_3MBPP1_r2[,2]) 
colnames(Spt5T1A_5minVs0min_gb) <- c("D1", "D2", "M1", "M2")
rownames(Spt5T1A_5minVs0min_gb) <- rownames(WT_D1)
Spt5T1A_5minVs0min_pr = cbind(Cdk9as_spt5_T1A_0min_DMSO_r1[,1], Cdk9as_spt5_T1A_0min_DMSO_r2[,1], Cdk9as_spt5_T1A_5min_3MBPP1_r1[,1], Cdk9as_spt5_T1A_5min_3MBPP1_r2[,1]) 
colnames(Spt5T1A_5minVs0min_pr) <- c("D1", "D2", "M1", "M2")
rownames(Spt5T1A_5minVs0min_pr) <- rownames(WT_D1)
Spt5T1A_5minVs0min_term = cbind(Cdk9as_spt5_T1A_0min_DMSO_r1[,3], Cdk9as_spt5_T1A_0min_DMSO_r2[,3], Cdk9as_spt5_T1A_5min_3MBPP1_r1[,3], Cdk9as_spt5_T1A_5min_3MBPP1_r2[,3]) 
colnames(Spt5T1A_5minVs0min_term) <- c("D1", "D2", "M1", "M2")
rownames(Spt5T1A_5minVs0min_term) <- rownames(WT_D1)
# Cdk9as spt5T1E 
# 1 min
Spt5T1E_1minVs0min_gb = cbind(Cdk9as_spt5_T1E_0min_DMSO_r1[,2], Cdk9as_spt5_T1E_0min_DMSO_r2[,2], Cdk9as_spt5_T1E_1min_3MBPP1_r1[,2], Cdk9as_spt5_T1E_1min_3MBPP1_r2[,2]) 
colnames(Spt5T1E_1minVs0min_gb) <- c("D1", "D2", "M1", "M2")
rownames(Spt5T1E_1minVs0min_gb) <- rownames(WT_D1)
Spt5T1E_1minVs0min_pr = cbind(Cdk9as_spt5_T1E_0min_DMSO_r1[,1], Cdk9as_spt5_T1E_0min_DMSO_r2[,1], Cdk9as_spt5_T1E_1min_3MBPP1_r1[,1], Cdk9as_spt5_T1E_1min_3MBPP1_r2[,1]) 
colnames(Spt5T1E_1minVs0min_pr) <- c("D1", "D2", "M1", "M2")
rownames(Spt5T1E_1minVs0min_pr) <- rownames(WT_D1)
Spt5T1E_1minVs0min_term = cbind(Cdk9as_spt5_T1E_0min_DMSO_r1[,3], Cdk9as_spt5_T1E_0min_DMSO_r2[,3], Cdk9as_spt5_T1E_1min_3MBPP1_r1[,3], Cdk9as_spt5_T1E_1min_3MBPP1_r2[,3]) 
colnames(Spt5T1E_1minVs0min_term) <- c("D1", "D2", "M1", "M2")
rownames(Spt5T1E_1minVs0min_term) <- rownames(WT_D1)
# 5 min
Spt5T1E_5minVs0min_gb = cbind(Cdk9as_spt5_T1E_0min_DMSO_r1[,2], Cdk9as_spt5_T1E_0min_DMSO_r2[,2], Cdk9as_spt5_T1E_5min_3MBPP1_r1[,2], Cdk9as_spt5_T1E_5min_3MBPP1_r2[,2]) 
colnames(Spt5T1E_5minVs0min_gb) <- c("D1", "D2", "M1", "M2")
rownames(Spt5T1E_5minVs0min_gb) <- rownames(WT_D1)
Spt5T1E_5minVs0min_pr = cbind(Cdk9as_spt5_T1E_0min_DMSO_r1[,1], Cdk9as_spt5_T1E_0min_DMSO_r2[,1], Cdk9as_spt5_T1E_5min_3MBPP1_r1[,1], Cdk9as_spt5_T1E_5min_3MBPP1_r2[,1]) 
colnames(Spt5T1E_5minVs0min_pr) <- c("D1", "D2", "M1", "M2")
rownames(Spt5T1E_5minVs0min_pr) <- rownames(WT_D1)
Spt5T1E_5minVs0min_term = cbind(Cdk9as_spt5_T1E_0min_DMSO_r1[,3], Cdk9as_spt5_T1E_0min_DMSO_r2[,3], Cdk9as_spt5_T1E_5min_3MBPP1_r1[,3], Cdk9as_spt5_T1E_5min_3MBPP1_r2[,3]) 
colnames(Spt5T1E_5minVs0min_term) <- c("D1", "D2", "M1", "M2")
rownames(Spt5T1E_5minVs0min_term) <- rownames(WT_D1)
# Compare Spt5 mutants (all untreated)
# Spt5WT7 vs Spt5T1A
Spt5WT7_Spt5T1A_0min_gb = cbind(Cdk9as_spt5_WT7_0min_DMSO_r1[,2], Cdk9as_spt5_WT7_0min_DMSO_r2[,2], Cdk9as_spt5_T1A_0min_DMSO_r1[,2], Cdk9as_spt5_T1A_0min_DMSO_r2[,2]) 
colnames(Spt5WT7_Spt5T1A_0min_gb) <- c("D1", "D2", "M1", "M2")
rownames(Spt5WT7_Spt5T1A_0min_gb) <- rownames(WT_D1)
Spt5WT7_Spt5T1A_0min_pr = cbind(Cdk9as_spt5_WT7_0min_DMSO_r1[,1], Cdk9as_spt5_WT7_0min_DMSO_r2[,1], Cdk9as_spt5_T1A_0min_DMSO_r1[,1], Cdk9as_spt5_T1A_0min_DMSO_r2[,1]) 
colnames(Spt5WT7_Spt5T1A_0min_pr) <- c("D1", "D2", "M1", "M2")
rownames(Spt5WT7_Spt5T1A_0min_pr) <- rownames(WT_D1)
Spt5WT7_Spt5T1A_0min_term = cbind(Cdk9as_spt5_WT7_0min_DMSO_r1[,3], Cdk9as_spt5_WT7_0min_DMSO_r2[,3], Cdk9as_spt5_T1A_0min_DMSO_r1[,3], Cdk9as_spt5_T1A_0min_DMSO_r2[,3]) 
colnames(Spt5WT7_Spt5T1A_0min_term) <- c("D1", "D2", "M1", "M2")
rownames(Spt5WT7_Spt5T1A_0min_term) <- rownames(WT_D1)
# Spt5WT7 vs Spt5T1E
Spt5WT7_Spt5T1E_0min_gb = cbind(Cdk9as_spt5_WT7_0min_DMSO_r1[,2], Cdk9as_spt5_WT7_0min_DMSO_r2[,2], Cdk9as_spt5_T1E_0min_DMSO_r1[,2], Cdk9as_spt5_T1E_0min_DMSO_r2[,2]) 
colnames(Spt5WT7_Spt5T1E_0min_gb) <- c("D1", "D2", "M1", "M2")
rownames(Spt5WT7_Spt5T1E_0min_gb) <- rownames(WT_D1)
Spt5WT7_Spt5T1E_0min_pr = cbind(Cdk9as_spt5_WT7_0min_DMSO_r1[,1], Cdk9as_spt5_WT7_0min_DMSO_r2[,1], Cdk9as_spt5_T1E_0min_DMSO_r1[,1], Cdk9as_spt5_T1E_0min_DMSO_r2[,1]) 
colnames(Spt5WT7_Spt5T1E_0min_pr) <- c("D1", "D2", "M1", "M2")
rownames(Spt5WT7_Spt5T1E_0min_pr) <- rownames(WT_D1)
Spt5WT7_Spt5T1E_0min_term = cbind(Cdk9as_spt5_WT7_0min_DMSO_r1[,3], Cdk9as_spt5_WT7_0min_DMSO_r2[,3], Cdk9as_spt5_T1E_0min_DMSO_r1[,3], Cdk9as_spt5_T1E_0min_DMSO_r2[,3]) 
colnames(Spt5WT7_Spt5T1E_0min_term) <- c("D1", "D2", "M1", "M2")
rownames(Spt5WT7_Spt5T1E_0min_term) <- rownames(WT_D1)

#### 6th experiment: Comparisons between different Dis2 mutants #####
# Fisher WT vs dis2-11
FisherWT_dis2_11_gb = cbind(Fisher_WT_r1[,2], Fisher_WT_r2[,2], Fisher_dis2_11_r1[,2], Fisher_dis2_11_r2[,2]) 
colnames(FisherWT_dis2_11_gb) <- c("D1", "D2", "M1", "M2")
rownames(FisherWT_dis2_11_gb) <- rownames(WT_D1)
FisherWT_dis2_11_pr = cbind(Fisher_WT_r1[,1], Fisher_WT_r2[,1], Fisher_dis2_11_r1[,1], Fisher_dis2_11_r2[,1]) 
colnames(FisherWT_dis2_11_pr) <- c("D1", "D2", "M1", "M2")
rownames(FisherWT_dis2_11_pr) <- rownames(WT_D1)
FisherWT_dis2_11_term = cbind(Fisher_WT_r1[,3], Fisher_WT_r2[,3], Fisher_dis2_11_r1[,3], Fisher_dis2_11_r2[,3]) 
colnames(FisherWT_dis2_11_term) <- c("D1", "D2", "M1", "M2")
rownames(FisherWT_dis2_11_term) <- rownames(WT_D1)
# Fisher WT vs dis2-delete
FisherWT_dis2Del_gb = cbind(Fisher_WT_r1[,2], Fisher_WT_r2[,2], Fisher_dis2Del_r1[,2], Fisher_dis2Del_r2[,2]) 
colnames(FisherWT_dis2Del_gb) <- c("D1", "D2", "M1", "M2")
rownames(FisherWT_dis2Del_gb) <- rownames(WT_D1)
FisherWT_dis2Del_pr = cbind(Fisher_WT_r1[,1], Fisher_WT_r2[,1], Fisher_dis2Del_r1[,1], Fisher_dis2Del_r2[,1]) 
colnames(FisherWT_dis2Del_pr) <- c("D1", "D2", "M1", "M2")
rownames(FisherWT_dis2Del_pr) <- rownames(WT_D1)
FisherWT_dis2Del_term = cbind(Fisher_WT_r1[,3], Fisher_WT_r2[,3], Fisher_dis2Del_r1[,3], Fisher_dis2Del_r2[,3]) 
colnames(FisherWT_dis2Del_term) <- c("D1", "D2", "M1", "M2")
rownames(FisherWT_dis2Del_term) <- rownames(WT_D1)
# Fisher WT vs dis2-T316A
FisherWT_dis2_T316A_gb = cbind(Fisher_WT_r1[,2], Fisher_WT_r2[,2], Fisher_dis2_T316A_r1[,2], Fisher_dis2_T316A_r2[,2]) 
colnames(FisherWT_dis2_T316A_gb) <- c("D1", "D2", "M1", "M2")
rownames(FisherWT_dis2_T316A_gb) <- rownames(WT_D1)
FisherWT_dis2_T316A_pr = cbind(Fisher_WT_r1[,1], Fisher_WT_r2[,1], Fisher_dis2_T316A_r1[,1], Fisher_dis2_T316A_r2[,1]) 
colnames(FisherWT_dis2_T316A_pr) <- c("D1", "D2", "M1", "M2")
rownames(FisherWT_dis2_T316A_pr) <- rownames(WT_D1)
FisherWT_dis2_T316A_term = cbind(Fisher_WT_r1[,3], Fisher_WT_r2[,3], Fisher_dis2_T316A_r1[,3], Fisher_dis2_T316A_r2[,3]) 
colnames(FisherWT_dis2_T316A_term) <- c("D1", "D2", "M1", "M2")
rownames(FisherWT_dis2_T316A_term) <- rownames(WT_D1)
# Fisher WT vs dis2-T316D
FisherWT_dis2_T316D_gb = cbind(Fisher_WT_r1[,2], Fisher_WT_r2[,2], Fisher_dis2_T316D_r1[,2], Fisher_dis2_T316D_r2[,2]) 
colnames(FisherWT_dis2_T316D_gb) <- c("D1", "D2", "M1", "M2")
rownames(FisherWT_dis2_T316D_gb) <- rownames(WT_D1)
FisherWT_dis2_T316D_pr = cbind(Fisher_WT_r1[,1], Fisher_WT_r2[,1], Fisher_dis2_T316D_r1[,1], Fisher_dis2_T316D_r2[,1]) 
colnames(FisherWT_dis2_T316D_pr) <- c("D1", "D2", "M1", "M2")
rownames(FisherWT_dis2_T316D_pr) <- rownames(WT_D1)
FisherWT_dis2_T316D_term = cbind(Fisher_WT_r1[,3], Fisher_WT_r2[,3], Fisher_dis2_T316D_r1[,3], Fisher_dis2_T316D_r2[,3]) 
colnames(FisherWT_dis2_T316D_term) <- c("D1", "D2", "M1", "M2")
rownames(FisherWT_dis2_T316D_term) <- rownames(WT_D1)

# Make lists of the above data frames to apply function to.  
ExpList_gb = list(WT_gb, MCS6_gb, CDK9_gb, Isk1_gb, Mcs6_CDK9_gb, Cdk9as_18C_gb, dis2ts_18C_gb, dis2ts_Cdk9as_18C_gb, dis2ts_30C_gb, 
                  dis2ts_Cdk9as_30C_gb, dis2ts_30Cv18C_DMSO_gb, dis2ts_30Cv18C_3MBPP1_gb, dis2ts_Cdk9as_30Cv18C_DMSO_gb,dis2ts_Cdk9as_30Cv18C_3MBPP1_gb,
                  Cdk9as_30secVs0min_gb, Cdk9as_1minVs0min_gb, Cdk9as_2.5minVs0min_gb, Cdk9as_5minVs0min_gb, Cdk9as_7.5minVs0min_gb, Cdk9as_10minVs0min_gb, Cdk9as_20minVs0min_gb,
                  Lsk1_CDK9_gb, Spt5WT7_1minVs0min_gb, Spt5WT7_5minVs0min_gb, Spt5T1A_1minVs0min_gb, Spt5T1A_5minVs0min_gb, Spt5T1E_1minVs0min_gb, Spt5T1E_5minVs0min_gb,
                  Spt5WT7_Spt5T1A_0min_gb, Spt5WT7_Spt5T1E_0min_gb, FisherWT_dis2_11_gb, FisherWT_dis2Del_gb, FisherWT_dis2_T316A_gb, FisherWT_dis2_T316D_gb)
ExpList_pr = list(WT_pr, MCS6_pr, CDK9_pr, Isk1_pr, Mcs6_CDK9_pr, Cdk9as_18C_pr, dis2ts_18C_pr, dis2ts_Cdk9as_18C_pr, dis2ts_30C_pr, 
                  dis2ts_Cdk9as_30C_pr, dis2ts_30Cv18C_DMSO_pr, dis2ts_30Cv18C_3MBPP1_pr, dis2ts_Cdk9as_30Cv18C_DMSO_pr,dis2ts_Cdk9as_30Cv18C_3MBPP1_pr,
                  Cdk9as_30secVs0min_pr, Cdk9as_1minVs0min_pr, Cdk9as_2.5minVs0min_pr, Cdk9as_5minVs0min_pr, Cdk9as_7.5minVs0min_pr, Cdk9as_10minVs0min_pr, Cdk9as_20minVs0min_pr,
                  Lsk1_CDK9_pr, Spt5WT7_1minVs0min_pr, Spt5WT7_5minVs0min_pr, Spt5T1A_1minVs0min_pr, Spt5T1A_5minVs0min_pr, Spt5T1E_1minVs0min_pr, Spt5T1E_5minVs0min_pr,
                  Spt5WT7_Spt5T1A_0min_pr, Spt5WT7_Spt5T1E_0min_pr, FisherWT_dis2_11_pr, FisherWT_dis2Del_pr, FisherWT_dis2_T316A_pr, FisherWT_dis2_T316D_pr)
ExpList_term = list(WT_term, MCS6_term, CDK9_term, Isk1_term, Mcs6_CDK9_term, Cdk9as_18C_term, dis2ts_18C_term, dis2ts_Cdk9as_18C_term, dis2ts_30C_term, 
                    dis2ts_Cdk9as_30C_term, dis2ts_30Cv18C_DMSO_term, dis2ts_30Cv18C_3MBPP1_term, dis2ts_Cdk9as_30Cv18C_DMSO_term, dis2ts_Cdk9as_30Cv18C_3MBPP1_term,
                    Cdk9as_30secVs0min_term, Cdk9as_1minVs0min_term, Cdk9as_2.5minVs0min_term, Cdk9as_5minVs0min_term, Cdk9as_7.5minVs0min_term, Cdk9as_10minVs0min_term, Cdk9as_20minVs0min_term,
                    Lsk1_CDK9_term, Spt5WT7_1minVs0min_term, Spt5WT7_5minVs0min_term, Spt5T1A_1minVs0min_term, Spt5T1A_5minVs0min_term, Spt5T1E_1minVs0min_term, Spt5T1E_5minVs0min_term,
                    Spt5WT7_Spt5T1A_0min_term, Spt5WT7_Spt5T1E_0min_term, FisherWT_dis2_11_term, FisherWT_dis2Del_term, FisherWT_dis2_T316A_term, FisherWT_dis2_T316D_term)
## prepare a separate normFactor list for each experiment (i.e. WT, mcs6, CDK9)
NFlist = read.table(file = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/countData/replicate_spikeNormFactors.txt", head = T)
WT_NFs <- c(NFlist[c(1:4), 2])/100000
mcs6_NFs <- c(NFlist[c(5,6,8), 2])/100000  ## Note that we are skipping the first treated replicate here due to outlier effects
CDK9_NFs <- c(NFlist[c(9:12), 2])/100000
Isk1_NFs <- c(NFlist[c(13:16), 2])/100000
mcs6_CDK9_NFs <- c(NFlist[c(17:20), 2])/100000
CDK9_18C_NFs <- c(NFlist[c(21:24), 2])/100000
dis2ts_18C_NFs <- c(NFlist[c(25:28), 2])/100000
dis2ts_30C_NFs <- c(NFlist[c(29:32), 2])/100000
CDK9_dis2ts_18C_NFs <- c(NFlist[c(33:36), 2])/100000
CDK9_dis2ts_30C_NFs <- c(NFlist[c(37:40), 2])/100000
dis2ts_30Cv18C_DMSO_NFs <- c(NFlist[c(29,30,25,26), 2])/100000
dis2ts_30Cv18C_3MBPP1_NFs <- c(NFlist[c(31,32,27,28), 2])/100000
CDK9_dis2ts_30Cv18C_DMSO_NFs <- c(NFlist[c(37,38,33,34), 2])/100000
CDK9_dis2ts_30Cv18C_3MBPP1_NFs <- c(NFlist[c(39,40,35,36), 2])/100000
WTvCDK9_NFs <- c(NFlist[c(1,2,9,10), 2])/100000
WTvMcs6_NFs <- c(NFlist[c(1,2,5,6), 2])/100000
WTvMcs6_CDK9_NFs <- c(NFlist[c(1,2,17,18), 2])/100000
WTvIsk1_NFs <- c(NFlist[c(1,2,13,14), 2])/100000
CDK9_30secvs0min_NFs <- c(NFlist[c(51,52,53,54), 2])/100000
CDK9_1vs0min_NFs <- c(NFlist[c(51,52,55,56), 2])/100000
CDK9_2.5vs0min_NFs <- c(NFlist[c(51,52,57,58), 2])/100000
CDK9_5vs0min_NFs <- c(NFlist[c(51,52,59,60), 2])/100000
CDK9_7.5vs0min_NFs <- c(NFlist[c(51,52,61,62), 2])/100000
CDK9_10vs0min_NFs <- c(NFlist[c(51,52,63,64), 2])/100000
CDK9_20vs0min_NFs <- c(NFlist[c(51,52,65,66), 2])/100000
Lsk1as_Cdk9as_NFs <- c(NFlist[c(67,68,69,70), 2])/100000
Spt5WT7_1minVs0min_NFs <- c(NFlist[c(81,82,83,84), 2])/100000
Spt5WT7_5minVs0min_NFs <- c(NFlist[c(81,82,85,86), 2])/100000
Spt5T1A_1minVs0min_NFs <- c(NFlist[c(87,88,89,90), 2])/100000
Spt5T1A_5minVs0min_NFs <- c(NFlist[c(87,88,91,92), 2])/100000
Spt5T1E_1minVs0min_NFs <- c(NFlist[c(93,94,95,96), 2])/100000
Spt5T1E_5minVs0min_NFs <- c(NFlist[c(93,94,97,98), 2])/100000
Spt5WT7_Spt5T1A_0min_NFs <- c(NFlist[c(81,82,87,88), 2])/100000
Spt5WT7_Spt5T1E_0min_NFs <- c(NFlist[c(81,82,93,94), 2])/100000
FisherWT_dis2_11_NFs <- c(NFlist[c(71,72,73,74), 2])/100000
FisherWT_dis2Del_NFs <- c(NFlist[c(71,72,75,76), 2])/100000
FisherWT_dis2_T316A_NFs <- c(NFlist[c(71,72,77,78), 2])/100000
FisherWT_dis2_T316D_NFs <- c(NFlist[c(71,72,79,80), 2])/100000
#########################################################################################################
## Prepare Tables with columns of raw counts for samples you want to compare (format:  untreated, untreated, treated, treated)
# function to perform basic DESeq2 analysis on 4-column data frame

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
### apply function to proper experimental data and corresponding normalization factors 
#WT
WT_gb_res = DEseq2.analyze(df = ExpList_gb[[1]], sizeFactors = WT_NFs)
WT_pr_res = DEseq2.analyze(df = ExpList_pr[[1]], sizeFactors = WT_NFs)
WT_term_res = DEseq2.analyze(df = ExpList_term[[1]], sizeFactors = WT_NFs)
#mcs6-as. NOTE:  Not using mcs6_M1 because it is a huge outlier and affects the analysis. 
mcs6_gb_res = DEseq2.analyze(df = ExpList_gb[[2]], sizeFactors = mcs6_NFs, D_reps = 2, M_reps = 1)
mcs6_pr_res = DEseq2.analyze(df = ExpList_pr[[2]], sizeFactors = mcs6_NFs, D_reps = 2, M_reps = 1)
mcs6_term_res = DEseq2.analyze(df = ExpList_term[[2]], sizeFactors = mcs6_NFs, D_reps = 2, M_reps = 1)
#CDK9-as
CDK9_gb_res = DEseq2.analyze(df = ExpList_gb[[3]], sizeFactors = CDK9_NFs)
CDK9_pr_res = DEseq2.analyze(df = ExpList_pr[[3]], sizeFactors = CDK9_NFs)
CDK9_term_res = DEseq2.analyze(df = ExpList_term[[3]], sizeFactors = CDK9_NFs)
#Isk1-as
Isk1_gb_res = DEseq2.analyze(df = ExpList_gb[[4]], sizeFactors = Isk1_NFs)
Isk1_pr_res = DEseq2.analyze(df = ExpList_pr[[4]], sizeFactors = Isk1_NFs)
Isk1_term_res = DEseq2.analyze(df = ExpList_term[[4]], sizeFactors = Isk1_NFs)
#Mcs6_Cdk9as
mcs6_CDK9_gb_res = DEseq2.analyze(df = ExpList_gb[[5]], sizeFactors = mcs6_CDK9_NFs)
mcs6_CDK9_pr_res = DEseq2.analyze(df = ExpList_pr[[5]], sizeFactors = mcs6_CDK9_NFs)
mcs6_CDK9_term_res = DEseq2.analyze(df = ExpList_term[[5]], sizeFactors = mcs6_CDK9_NFs)
#Cdk9as_18C
Cdk9as_18C_gb_res = DEseq2.analyze(df = ExpList_gb[[6]], sizeFactors = CDK9_18C_NFs)
Cdk9as_18C_pr_res = DEseq2.analyze(df = ExpList_pr[[6]], sizeFactors = CDK9_18C_NFs)
Cdk9as_18C_term_res = DEseq2.analyze(df = ExpList_term[[6]], sizeFactors = CDK9_18C_NFs)
#dis2ts_18C
dis2ts_18C_gb_res = DEseq2.analyze(df = ExpList_gb[[7]], sizeFactors = dis2ts_18C_NFs)
dis2ts_18C_pr_res = DEseq2.analyze(df = ExpList_pr[[7]], sizeFactors = dis2ts_18C_NFs)
dis2ts_18C_term_res = DEseq2.analyze(df = ExpList_term[[7]], sizeFactors = dis2ts_18C_NFs)
#Cdk9as_dis2ts_18C
CDK9_dis2ts_18C_gb_res = DEseq2.analyze(df = ExpList_gb[[8]], sizeFactors = CDK9_dis2ts_18C_NFs)
CDK9_dis2ts_18C_pr_res = DEseq2.analyze(df = ExpList_pr[[8]], sizeFactors = CDK9_dis2ts_18C_NFs)
CDK9_dis2ts_18C_term_res = DEseq2.analyze(df = ExpList_term[[8]], sizeFactors = CDK9_dis2ts_18C_NFs)
#dis2ts_30C
dis2ts_30C_gb_res = DEseq2.analyze(df = ExpList_gb[[9]], sizeFactors = dis2ts_30C_NFs)
dis2ts_30C_pr_res = DEseq2.analyze(df = ExpList_pr[[9]], sizeFactors = dis2ts_30C_NFs)
dis2ts_30C_term_res = DEseq2.analyze(df = ExpList_term[[9]], sizeFactors = dis2ts_30C_NFs)
#Cdk9as_dis2ts_30C
CDK9_dis2ts_30C_gb_res = DEseq2.analyze(df = ExpList_gb[[10]], sizeFactors = CDK9_dis2ts_30C_NFs)
CDK9_dis2ts_30C_pr_res = DEseq2.analyze(df = ExpList_pr[[10]], sizeFactors = CDK9_dis2ts_30C_NFs)
CDK9_dis2ts_30C_term_res = DEseq2.analyze(df = ExpList_term[[10]], sizeFactors = CDK9_dis2ts_30C_NFs)
#dis2ts_30Cv18C_DMSO
dis2ts_30Cv18C_DMSO_gb_res = DEseq2.analyze(df = ExpList_gb[[11]], sizeFactors = dis2ts_30Cv18C_DMSO_NFs)
dis2ts_30Cv18C_DMSO_pr_res = DEseq2.analyze(df = ExpList_pr[[11]], sizeFactors = dis2ts_30Cv18C_DMSO_NFs)
dis2ts_30Cv18C_DMSO_term_res = DEseq2.analyze(df = ExpList_term[[11]], sizeFactors = dis2ts_30Cv18C_DMSO_NFs)
#dis2ts_30Cv18C_3MBPP1
dis2ts_30Cv18C_3MBPP1_gb_res = DEseq2.analyze(df = ExpList_gb[[12]], sizeFactors = dis2ts_30Cv18C_3MBPP1_NFs)
dis2ts_30Cv18C_3MBPP1_pr_res = DEseq2.analyze(df = ExpList_pr[[12]], sizeFactors = dis2ts_30Cv18C_3MBPP1_NFs)
dis2ts_30Cv18C_3MBPP1_term_res = DEseq2.analyze(df = ExpList_term[[12]], sizeFactors = dis2ts_30Cv18C_3MBPP1_NFs)
#CDK9_dis2ts_30Cv18C_DMSO
CDK9_dis2ts_30Cv18C_DMSO_gb_res = DEseq2.analyze(df = ExpList_gb[[13]], sizeFactors = CDK9_dis2ts_30Cv18C_DMSO_NFs)
CDK9_dis2ts_30Cv18C_DMSO_pr_res = DEseq2.analyze(df = ExpList_pr[[13]], sizeFactors = CDK9_dis2ts_30Cv18C_DMSO_NFs)
CDK9_dis2ts_30Cv18C_DMSO_term_res = DEseq2.analyze(df = ExpList_term[[13]], sizeFactors = CDK9_dis2ts_30Cv18C_DMSO_NFs)
#CDK9_dis2ts_30Cv18C_3MBPP1
CDK9_dis2ts_30Cv18C_3MBPP1_gb_res = DEseq2.analyze(df = ExpList_gb[[14]], sizeFactors = CDK9_dis2ts_30Cv18C_3MBPP1_NFs)
CDK9_dis2ts_30Cv18C_3MBPP1_pr_res = DEseq2.analyze(df = ExpList_pr[[14]], sizeFactors = CDK9_dis2ts_30Cv18C_3MBPP1_NFs)
CDK9_dis2ts_30Cv18C_3MBPP1_term_res = DEseq2.analyze(df = ExpList_term[[14]], sizeFactors = CDK9_dis2ts_30Cv18C_3MBPP1_NFs)

# WT vs mutants (untreated)
# WT vs Mcs6 DMSO
WTvMcs6_DMSO_gb_res = DEseq2.analyze(df = WTvMcs6_gb, sizeFactors = WTvMcs6_NFs)
WTvMcs6_DMSO_pr_res = DEseq2.analyze(df = WTvMcs6_pr, sizeFactors = WTvMcs6_NFs)
WTvMcs6_DMSO_term_res = DEseq2.analyze(df = WTvMcs6_term, sizeFactors = WTvMcs6_NFs)
# WT vs Cdk9 DMSO
WTvCDK9_DMSO_gb_res = DEseq2.analyze(df = WTvCDK9_gb, sizeFactors = WTvCDK9_NFs)
WTvCDK9_DMSO_pr_res = DEseq2.analyze(df = WTvCDK9_pr, sizeFactors = WTvCDK9_NFs)
WTvCDK9_DMSO_term_res = DEseq2.analyze(df = WTvCDK9_term, sizeFactors = WTvCDK9_NFs)
# WT vs Mcs6 Cdk9 DMSO
WTvMcs6_CDK9_DMSO_gb_res = DEseq2.analyze(df = WTvMcs6_CDK9_gb, sizeFactors = WTvMcs6_CDK9_NFs)
WTvMcs6_CDK9_DMSO_pr_res = DEseq2.analyze(df = WTvMcs6_CDK9_pr, sizeFactors = WTvMcs6_CDK9_NFs)
WTvMcs6_CDK9_DMSO_term_res = DEseq2.analyze(df = WTvMcs6_CDK9_term, sizeFactors = WTvMcs6_CDK9_NFs)
# WT vs Lsk1 DMSO
WTvIsk1_DMSO_gb_res = DEseq2.analyze(df = WTvIsk1_gb, sizeFactors = WTvIsk1_NFs)
WTvIsk1_DMSO_pr_res = DEseq2.analyze(df = WTvIsk1_pr, sizeFactors = WTvIsk1_NFs)
WTvIsk1_DMSO_term_res = DEseq2.analyze(df = WTvIsk1_term, sizeFactors = WTvIsk1_NFs)

# Cdk9as Time Course (Full) Experiment (New data)
# 0 min vs 30 seconds
CDK9_0.5vs0min_gb_res = DEseq2.analyze(df = Cdk9as_30secVs0min_gb, sizeFactors = CDK9_30secvs0min_NFs)
CDK9_0.5vs0min_pr_res = DEseq2.analyze(df = Cdk9as_30secVs0min_pr, sizeFactors = CDK9_30secvs0min_NFs)
CDK9_0.5vs0min_term_res = DEseq2.analyze(df = Cdk9as_30secVs0min_term, sizeFactors = CDK9_30secvs0min_NFs)
# 0 min vs 1 minute
CDK9_1vs0min_gb_res = DEseq2.analyze(df = Cdk9as_1minVs0min_gb, sizeFactors = CDK9_1vs0min_NFs)
CDK9_1vs0min_pr_res = DEseq2.analyze(df = Cdk9as_1minVs0min_pr, sizeFactors = CDK9_1vs0min_NFs)
CDK9_1vs0min_term_res = DEseq2.analyze(df = Cdk9as_1minVs0min_term, sizeFactors = CDK9_1vs0min_NFs)
# 0 min vs 2.5 minute
CDK9_2.5vs0min_gb_res = DEseq2.analyze(df = Cdk9as_2.5minVs0min_gb, sizeFactors = CDK9_2.5vs0min_NFs)
CDK9_2.5vs0min_pr_res = DEseq2.analyze(df = Cdk9as_2.5minVs0min_pr, sizeFactors = CDK9_2.5vs0min_NFs)
CDK9_2.5vs0min_term_res = DEseq2.analyze(df = Cdk9as_2.5minVs0min_term, sizeFactors = CDK9_2.5vs0min_NFs)
# 0 min vs 5 minute
CDK9_5vs0min_gb_res = DEseq2.analyze(df = Cdk9as_5minVs0min_gb, sizeFactors = CDK9_5vs0min_NFs)
CDK9_5vs0min_pr_res = DEseq2.analyze(df = Cdk9as_5minVs0min_pr, sizeFactors = CDK9_5vs0min_NFs)
CDK9_5vs0min_term_res = DEseq2.analyze(df = Cdk9as_5minVs0min_term, sizeFactors = CDK9_5vs0min_NFs)
# 0 min vs 7.5 minute
CDK9_7.5vs0min_gb_res = DEseq2.analyze(df = Cdk9as_7.5minVs0min_gb, sizeFactors = CDK9_7.5vs0min_NFs)
CDK9_7.5vs0min_pr_res = DEseq2.analyze(df = Cdk9as_7.5minVs0min_pr, sizeFactors = CDK9_7.5vs0min_NFs)
CDK9_7.5vs0min_term_res = DEseq2.analyze(df = Cdk9as_7.5minVs0min_term, sizeFactors = CDK9_7.5vs0min_NFs)
# 0 min vs 10 minute
CDK9_10vs0min_gb_res = DEseq2.analyze(df = Cdk9as_10minVs0min_gb, sizeFactors = CDK9_10vs0min_NFs)
CDK9_10vs0min_pr_res = DEseq2.analyze(df = Cdk9as_10minVs0min_pr, sizeFactors = CDK9_10vs0min_NFs)
CDK9_10vs0min_term_res = DEseq2.analyze(df = Cdk9as_10minVs0min_term, sizeFactors = CDK9_10vs0min_NFs)
# 0 min vs 20 minute
CDK9_20vs0min_gb_res = DEseq2.analyze(df = Cdk9as_20minVs0min_gb, sizeFactors = CDK9_20vs0min_NFs)
CDK9_20vs0min_pr_res = DEseq2.analyze(df = Cdk9as_20minVs0min_pr, sizeFactors = CDK9_20vs0min_NFs)
CDK9_20vs0min_term_res = DEseq2.analyze(df = Cdk9as_20minVs0min_term, sizeFactors = CDK9_20vs0min_NFs)
# Lsk1as Cdk9as double mutant
Lsk1_CDK9_gb_res = DEseq2.analyze(df = Lsk1_CDK9_gb, sizeFactors = Lsk1as_Cdk9as_NFs)
Lsk1_CDK9_pr_res = DEseq2.analyze(df = Lsk1_CDK9_pr, sizeFactors = Lsk1as_Cdk9as_NFs)
Lsk1_CDK9_term_res = DEseq2.analyze(df = Lsk1_CDK9_term, sizeFactors = Lsk1as_Cdk9as_NFs)

# Cdk9as Time Course in Spt5 CTR mutants 
# spt5WT7 1 minutes vs 0 minutes
Spt5WT7_1minVs0min_gb_res = DEseq2.analyze(df = Spt5WT7_1minVs0min_gb, sizeFactors = Spt5WT7_1minVs0min_NFs)
Spt5WT7_1minVs0min_pr_res = DEseq2.analyze(df = Spt5WT7_1minVs0min_pr, sizeFactors = Spt5WT7_1minVs0min_NFs)
Spt5WT7_1minVs0min_term_res = DEseq2.analyze(df = Spt5WT7_1minVs0min_term, sizeFactors = Spt5WT7_1minVs0min_NFs)
# spt5WT7 5 minutes vs 0 minutes
Spt5WT7_5minVs0min_gb_res = DEseq2.analyze(df = Spt5WT7_5minVs0min_gb, sizeFactors = Spt5WT7_5minVs0min_NFs)
Spt5WT7_5minVs0min_pr_res = DEseq2.analyze(df = Spt5WT7_5minVs0min_pr, sizeFactors = Spt5WT7_5minVs0min_NFs)
Spt5WT7_5minVs0min_term_res = DEseq2.analyze(df = Spt5WT7_5minVs0min_term, sizeFactors = Spt5WT7_5minVs0min_NFs)
# spt5T1A 1 minutes vs 0 minutes
Spt5T1A_1minVs0min_gb_res = DEseq2.analyze(df = Spt5T1A_1minVs0min_gb, sizeFactors = Spt5T1A_1minVs0min_NFs)
Spt5T1A_1minVs0min_pr_res = DEseq2.analyze(df = Spt5T1A_1minVs0min_pr, sizeFactors = Spt5T1A_1minVs0min_NFs)
Spt5T1A_1minVs0min_term_res = DEseq2.analyze(df = Spt5T1A_1minVs0min_term, sizeFactors = Spt5T1A_1minVs0min_NFs)
# spt5T1A 5 minutes vs 0 minutes
Spt5T1A_5minVs0min_gb_res = DEseq2.analyze(df = Spt5T1A_5minVs0min_gb, sizeFactors = Spt5T1A_5minVs0min_NFs)
Spt5T1A_5minVs0min_pr_res = DEseq2.analyze(df = Spt5T1A_5minVs0min_pr, sizeFactors = Spt5T1A_5minVs0min_NFs)
Spt5T1A_5minVs0min_term_res = DEseq2.analyze(df = Spt5T1A_5minVs0min_term, sizeFactors = Spt5T1A_5minVs0min_NFs)
# spt5T1E 1 minutes vs 0 minutes
Spt5T1E_1minVs0min_gb_res = DEseq2.analyze(df = Spt5T1A_1minVs0min_gb, sizeFactors = Spt5T1A_1minVs0min_NFs)
Spt5T1E_1minVs0min_pr_res = DEseq2.analyze(df = Spt5T1A_1minVs0min_pr, sizeFactors = Spt5T1A_1minVs0min_NFs)
Spt5T1E_1minVs0min_term_res = DEseq2.analyze(df = Spt5T1A_1minVs0min_term, sizeFactors = Spt5T1A_1minVs0min_NFs)
# spt5T1A 5 minutes vs 0 minutes
Spt5T1E_5minVs0min_gb_res = DEseq2.analyze(df = Spt5T1E_5minVs0min_gb, sizeFactors = Spt5T1E_5minVs0min_NFs)
Spt5T1E_5minVs0min_pr_res = DEseq2.analyze(df = Spt5T1E_5minVs0min_pr, sizeFactors = Spt5T1E_5minVs0min_NFs)
Spt5T1E_5minVs0min_term_res = DEseq2.analyze(df = Spt5T1E_5minVs0min_term, sizeFactors = Spt5T1E_5minVs0min_NFs)
# Comparison of Spt5WT7 with phosphomimetic and constitutively unphosphorylated Spt5 CTR
# spt5WT7 vs spt5T1A 0 min
Spt5WT7_Spt5T1A_0min_gb_res = DEseq2.analyze(df = Spt5WT7_Spt5T1A_0min_gb, sizeFactors = Spt5WT7_Spt5T1A_0min_NFs)
Spt5WT7_Spt5T1A_0min_pr_res = DEseq2.analyze(df = Spt5WT7_Spt5T1A_0min_pr, sizeFactors = Spt5WT7_Spt5T1A_0min_NFs)
Spt5WT7_Spt5T1A_0min_term_res = DEseq2.analyze(df = Spt5WT7_Spt5T1A_0min_term, sizeFactors = Spt5WT7_Spt5T1A_0min_NFs)
# spt5WT7 vs spt5T1A 0 min
Spt5WT7_Spt5T1E_0min_gb_res = DEseq2.analyze(df = Spt5WT7_Spt5T1E_0min_gb, sizeFactors = Spt5WT7_Spt5T1E_0min_NFs)
Spt5WT7_Spt5T1E_0min_pr_res = DEseq2.analyze(df = Spt5WT7_Spt5T1E_0min_pr, sizeFactors = Spt5WT7_Spt5T1E_0min_NFs)
Spt5WT7_Spt5T1E_0min_term_res = DEseq2.analyze(df = Spt5WT7_Spt5T1E_0min_term, sizeFactors = Spt5WT7_Spt5T1E_0min_NFs)

# Experiments comparing between different Dis2 mutants (Fisher)
# WT vs dis2-11
FisherWT_dis2_11_gb_res = DEseq2.analyze(df = FisherWT_dis2_11_gb, sizeFactors = FisherWT_dis2_11_NFs)
FisherWT_dis2_11_pr_res = DEseq2.analyze(df = FisherWT_dis2_11_pr, sizeFactors = FisherWT_dis2_11_NFs)
FisherWT_dis2_11_term_res = DEseq2.analyze(df = FisherWT_dis2_11_term, sizeFactors = FisherWT_dis2_11_NFs)
# WT vs dis2Delete
FisherWT_dis2Del_gb_res = DEseq2.analyze(df = FisherWT_dis2Del_gb, sizeFactors = FisherWT_dis2Del_NFs)
FisherWT_dis2Del_pr_res = DEseq2.analyze(df = FisherWT_dis2Del_pr, sizeFactors = FisherWT_dis2Del_NFs)
FisherWT_dis2Del_term_res = DEseq2.analyze(df = FisherWT_dis2Del_term, sizeFactors = FisherWT_dis2Del_NFs)
# WT vs dis2T316A
FisherWT_dis2_T316A_gb_res = DEseq2.analyze(df = FisherWT_dis2_T316A_gb, sizeFactors = FisherWT_dis2_T316A_NFs)
FisherWT_dis2_T316A_pr_res = DEseq2.analyze(df = FisherWT_dis2_T316A_pr, sizeFactors = FisherWT_dis2_T316A_NFs)
FisherWT_dis2_T316A_term_res = DEseq2.analyze(df = FisherWT_dis2_T316A_term, sizeFactors = FisherWT_dis2_T316A_NFs)
# WT vs dis2T316A
FisherWT_dis2_T316D_gb_res = DEseq2.analyze(df = FisherWT_dis2_T316D_gb, sizeFactors = FisherWT_dis2_T316D_NFs)
FisherWT_dis2_T316D_pr_res = DEseq2.analyze(df = FisherWT_dis2_T316D_pr, sizeFactors = FisherWT_dis2_T316D_NFs)
FisherWT_dis2_T316D_term_res = DEseq2.analyze(df = FisherWT_dis2_T316D_term, sizeFactors = FisherWT_dis2_T316D_NFs)


### write the DESeq2 results to files. 
write.table(WT_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_WT_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(WT_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_WT_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(WT_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_WT_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(mcs6_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_mcs6_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(mcs6_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_mcs6_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(mcs6_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_mcs6_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(CDK9_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(CDK9_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(CDK9_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(Isk1_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Isk1_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(Isk1_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Isk1_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(Isk1_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Isk1_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(mcs6_CDK9_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_mcs6_CDK9_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(mcs6_CDK9_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_mcs6_CDK9_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(mcs6_CDK9_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_mcs6_CDK9_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

write.table(Cdk9as_18C_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Cdk9as_18C_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(Cdk9as_18C_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Cdk9as_18C_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(Cdk9as_18C_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Cdk9as_18C_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(dis2ts_18C_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_dis2ts_18C_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(dis2ts_18C_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_dis2ts_18C_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(dis2ts_18C_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_dis2ts_18C_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(CDK9_dis2ts_18C_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_dis2ts_18C_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(CDK9_dis2ts_18C_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_dis2ts_18C_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(CDK9_dis2ts_18C_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_dis2ts_18C_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(dis2ts_30C_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_dis2ts_30C_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(dis2ts_30C_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_dis2ts_30C_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(dis2ts_30C_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_dis2ts_30C_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(CDK9_dis2ts_30C_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_dis2ts_30C_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(CDK9_dis2ts_30C_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_dis2ts_30C_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(CDK9_dis2ts_30C_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_dis2ts_30C_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(dis2ts_30Cv18C_DMSO_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_dis2ts_30Cv18C_DMSO_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(dis2ts_30Cv18C_DMSO_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_dis2ts_30Cv18C_DMSO_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(dis2ts_30Cv18C_DMSO_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_dis2ts_30Cv18C_DMSO_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(dis2ts_30Cv18C_3MBPP1_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_dis2ts_30Cv18C_3MBPP1_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(dis2ts_30Cv18C_3MBPP1_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_dis2ts_30Cv18C_3MBPP1_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(dis2ts_30Cv18C_3MBPP1_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_dis2ts_30Cv18C_3MBPP1_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(CDK9_dis2ts_30Cv18C_DMSO_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_dis2ts_30Cv18C_DMSO_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(CDK9_dis2ts_30Cv18C_DMSO_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_dis2ts_30Cv18C_DMSO_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(CDK9_dis2ts_30Cv18C_DMSO_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_dis2ts_30Cv18C_DMSO_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(CDK9_dis2ts_30Cv18C_3MBPP1_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_dis2ts_30Cv18C_3MBPP1_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(CDK9_dis2ts_30Cv18C_3MBPP1_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_dis2ts_30Cv18C_3MBPP1_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(CDK9_dis2ts_30Cv18C_3MBPP1_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_dis2ts_30Cv18C_3MBPP1_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
# 
write.table(WTvMcs6_DMSO_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_WTvsMcs6as_DMSO_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(WTvMcs6_DMSO_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_WTvsMcs6as_DMSO_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(WTvMcs6_DMSO_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_WTvsMcs6as_DMSO_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(WTvCDK9_DMSO_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_WTvsCDK9as_DMSO_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(WTvCDK9_DMSO_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_WTvsCDK9as_DMSO_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(WTvCDK9_DMSO_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_WTvsCDK9as_DMSO_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(WTvMcs6_CDK9_DMSO_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_WTvsMcs6_CDK9_DMSO_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(WTvMcs6_CDK9_DMSO_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_WTvsMcs6_CDK9_DMSO_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(WTvMcs6_CDK9_DMSO_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_WTvsMcs6_CDK9_DMSO_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(WTvIsk1_DMSO_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_WTvsIsk1_DMSO_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(WTvIsk1_DMSO_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_WTvsIsk1_DMSO_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(WTvIsk1_DMSO_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_WTvIsk1_DMSO_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
#30'' vs 0' 
write.table(CDK9_0.5vs0min_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_30secvs0min_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(CDK9_0.5vs0min_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_30secvs0min_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(CDK9_0.5vs0min_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_30secvs0min_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
#1' vs 0' 
write.table(CDK9_1vs0min_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_1vs0min_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(CDK9_1vs0min_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_1vs0min_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(CDK9_1vs0min_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_1vs0min_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
#2.5' vs 0' 
write.table(CDK9_2.5vs0min_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_2min30secVs0min_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(CDK9_2.5vs0min_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_2min30secVs0min_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(CDK9_2.5vs0min_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_2min30secVs0min_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
#5' vs 0' 
write.table(CDK9_5vs0min_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_5vs0min_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(CDK9_5vs0min_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_5vs0min_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(CDK9_5vs0min_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_5vs0min_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
#7.5' vs 0' 
write.table(CDK9_7.5vs0min_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_7min30secVs0min_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(CDK9_7.5vs0min_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_7min30secVs0min_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(CDK9_7.5vs0min_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_7min30secVs0min_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
#20 vs 0' 
write.table(CDK9_10vs0min_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_10vs0min_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(CDK9_10vs0min_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_10vs0min_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(CDK9_10vs0min_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_10vs0min_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
#20' vs 0' 
write.table(CDK9_20vs0min_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_20vs0min_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(CDK9_20vs0min_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_20vs0min_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(CDK9_20vs0min_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_CDK9_20vs0min_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
# Lsk1as Cdk9as double mutant
write.table(Lsk1_CDK9_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Lsk1_CDK9_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(Lsk1_CDK9_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Lsk1_CDK9_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(Lsk1_CDK9_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Lsk1_CDK9_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

# Cdk9as Time Course in Spt5 CTR mutants 
# spt5WT7 1 minutes vs 0 minutes
write.table(Spt5WT7_1minVs0min_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Spt5WT7_1minVs0min_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(Spt5WT7_1minVs0min_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Spt5WT7_1minVs0min_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(Spt5WT7_1minVs0min_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Spt5WT7_1minVs0min_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
# spt5WT7 5 minutes vs 0 minutes
write.table(Spt5WT7_5minVs0min_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Spt5WT7_5minVs0min_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(Spt5WT7_5minVs0min_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Spt5WT7_5minVs0min_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(Spt5WT7_5minVs0min_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Spt5WT7_5minVs0min_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
# Spt5T1A 1 minutes vs 0 minutes
write.table(Spt5T1A_1minVs0min_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Spt5T1A_1minVs0min_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(Spt5T1A_1minVs0min_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Spt5T1A_1minVs0min_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(Spt5T1A_1minVs0min_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Spt5T1A_1minVs0min_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
# Spt5T1A 5 minutes vs 0 minutes
write.table(Spt5T1A_5minVs0min_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Spt5T1A_5minVs0min_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(Spt5T1A_5minVs0min_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Spt5T1A_5minVs0min_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(Spt5T1A_5minVs0min_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Spt5T1A_5minVs0min_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
# Spt5T1E 1 minutes vs 0 minutes
write.table(Spt5T1E_1minVs0min_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Spt5T1E_1minVs0min_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(Spt5T1E_1minVs0min_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Spt5T1E_1minVs0min_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(Spt5T1E_1minVs0min_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Spt5T1E_1minVs0min_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
# Spt5T1E 5 minutes vs 0 minutes
write.table(Spt5T1E_5minVs0min_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Spt5T1E_5minVs0min_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(Spt5T1E_5minVs0min_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Spt5T1E_5minVs0min_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(Spt5T1E_5minVs0min_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Spt5T1E_5minVs0min_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
# Spt5WT7 vs spt5T1A 0 minutes
write.table(Spt5WT7_Spt5T1A_0min_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Spt5WT7_Spt5T1A_0min_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(Spt5WT7_Spt5T1A_0min_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Spt5WT7_Spt5T1A_0min_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(Spt5WT7_Spt5T1A_0min_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Spt5WT7_Spt5T1A_0min_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
# Spt5WT7 vs spt5T1E 0 minutes
write.table(Spt5WT7_Spt5T1E_0min_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Spt5WT7_Spt5T1E_0min_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(Spt5WT7_Spt5T1E_0min_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Spt5WT7_Spt5T1E_0min_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(Spt5WT7_Spt5T1E_0min_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_Spt5WT7_Spt5T1E_0min_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

# Comparisons between different Dis2 mutants:
# WT7 vs Dis2-11 
write.table(FisherWT_dis2_11_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_FisherWT_dis2_11_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(FisherWT_dis2_11_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_FisherWT_dis2_11_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(FisherWT_dis2_11_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_FisherWT_dis2_11_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
# WT7 vs Dis2Delete 
write.table(FisherWT_dis2Del_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_FisherWT_dis2Del_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(FisherWT_dis2Del_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_FisherWT_dis2Del_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(FisherWT_dis2Del_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_FisherWT_dis2Del_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
# WT7 vs Dis2-T316A
write.table(FisherWT_dis2_T316A_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_dis2_T316A_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(FisherWT_dis2_T316A_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_dis2_T316A_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(FisherWT_dis2_T316A_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_dis2_T316A_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
# WT7 vs Dis2-T316A
write.table(FisherWT_dis2_T316D_gb_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_dis2_T316D_gb_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(FisherWT_dis2_T316D_pr_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_dis2_T316D_pr_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(FisherWT_dis2_T316D_term_res, file = paste(DEseq_outpath ,"DEseq2_AllGenes_dis2_T316D_term_res.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")









