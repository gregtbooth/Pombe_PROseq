source("/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/scripts/GB_Functions.R")
library(gplots)
bwpath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/alignment/bw/"
infopath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/pombeinfo/"
fig_dir = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/figures/Dis2Analysis/MetaPlots/02-12-18/"
dir.create(fig_dir)
SPmappability = load.bigWig(paste(infopath, "SPombeMappable_SacCerNonOverLap_expanded_ChrFix.bw", sep = ""))
filteredGL = read.table(paste(infopath, "SP_CompleteFilteredGenes.bed", sep = ""))
SPGL_1kbsep = read.table(file = paste(infopath, "SP_PROcapObservedTSS_1000bpSep.bed", sep = ""), stringsAsFactors = F)
SPGL_1kbsep_filt = merge(filteredGL, SPGL_1kbsep, by.x = 4, by.y = 4)[, c(2,3,4,1,5,6)]
ObsTSS_NOoverlap = read.table(file = paste(infopath, "SP_PROcapObservedTSS_noOverlap.bed", sep = ""))
#pausedGenes = read.table(file = paste(infopath, "/SP_PausedGene_Filtered_bothReps.bed", sep = "")) #defined  in genome research paper
#NotpausedGenes = read.table(file = paste(infopath, "SP_nonPausedGene_Filtered_bothReps.bed", sep = "")) #defined  in genome research paper

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
RPM_NFs = read.table(file = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/countData/combined_RPM_NormFactors.txt", head = T)
RPM_repNFs = read.table(file = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/countData/replicate_RPM_NormFactors.txt", head = T)
wigset = rbind(c("5993_7157_26977_HJ72CBGXX_pombe_WT_COMBINED_DMSO_CAGATC_R1.fastq_pombe_plus.bw", "5993_7157_26977_HJ72CBGXX_pombe_WT_COMBINED_DMSO_CAGATC_R1.fastq_pombe_minus.bw", "WT_DMSO_combined", NFs[1,'Spikereads']/100000),
               c("5993_7157_26978_HJ72CBGXX_pombe_WT_COMBINED_3MB-PPI_ACTTGA_R1.fastq_pombe_plus.bw", "5993_7157_26978_HJ72CBGXX_pombe_WT_COMBINED_3MB-PPI_ACTTGA_R1.fastq_pombe_minus.bw", "WT_3MBPP1_combined", NFs[2,'Spikereads']/100000),
               c("7772_7157_43149_HVLYCBGXY_pombe_18C_Cdk9as_DMSO_COMBINED_GGCTAC_R1_pombe_plus.bw", "7772_7157_43149_HVLYCBGXY_pombe_18C_Cdk9as_DMSO_COMBINED_GGCTAC_R1_pombe_minus.bw", "Cdk9as_DMSO_18C_combined", NFs[11,'Spikereads']/100000),
               c("7772_7157_43150_HVLYCBGXY_pombe_18C_Cdk9as_3MBPP1_COMBINED_CTTGTA_R1_pombe_plus.bw", "7772_7157_43150_HVLYCBGXY_pombe_18C_Cdk9as_3MBPP1_COMBINED_CTTGTA_R1_pombe_minus.bw", "Cdk9as_3MBPP1_18C_combined", NFs[12,'Spikereads']/100000),
               c("7772_7157_43151_HVLYCBGXY_pombe_18C_dis2ts_DMSO_COMBINED_ATTCCT_R1_pombe_plus.bw", "7772_7157_43151_HVLYCBGXY_pombe_18C_dis2ts_DMSO_COMBINED_ATTCCT_R1_pombe_minus.bw", "dis2ts_DMSO_18C_combined", NFs[13,'Spikereads']/100000),
               c("7772_7157_43152_HVLYCBGXY_pombe_18C_dis2ts_3MBPP1_COMBINED_CAAAAG_R1_pombe_plus.bw", "7772_7157_43152_HVLYCBGXY_pombe_18C_dis2ts_3MBPP1_COMBINED_CAAAAG_R1_pombe_minus.bw", "dis2ts_3MBPP1_18C_combined", NFs[14,'Spikereads']/100000),
               c("7772_7157_43153_HVLYCBGXY_pombe_30C_dis2ts_DMSO_COMBINED_CAACTA_R1_pombe_plus.bw", "7772_7157_43153_HVLYCBGXY_pombe_30C_dis2ts_DMSO_COMBINED_CAACTA_R1_pombe_minus.bw", "dis2ts_DMSO_30C_combined", NFs[15,'Spikereads']/100000),
               c("7772_7157_43154_HVLYCBGXY_pombe_30C_dis2ts_3MBPP1_COMBINED_CACCGG_R1_pombe_plus.bw", "7772_7157_43154_HVLYCBGXY_pombe_30C_dis2ts_3MBPP1_COMBINED_CACCGG_R1_pombe_minus.bw", "dis2ts_3MBPP1_30C_combined", NFs[16,'Spikereads']/100000),
               c("7772_7157_43155_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_DMSO_COMBINED_CACGAT_R1_pombe_plus.bw", "7772_7157_43155_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_DMSO_COMBINED_CACGAT_R1_pombe_minus.bw", "CDK9as_dis2ts_DMSO_18C_combined", NFs[17,'Spikereads']/100000),
               c("7772_7157_43156_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_3MBPP1_COMBINED_CACTCA_R1_pombe_plus.bw", "7772_7157_43156_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_3MBPP1_COMBINED_CACTCA_R1_pombe_minus.bw", "CDK9as_dis2ts_3MBPP1_18C_combined", NFs[18,'Spikereads']/100000),
               c("7772_7157_43157_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_DMSO_COMBINED_CAGGCG_R1_pombe_plus.bw", "7772_7157_43157_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_DMSO_COMBINED_CAGGCG_R1_pombe_minus.bw", "CDK9as_dis2ts_DMSO_30C_combined", NFs[19,'Spikereads']/100000),
               c("7772_7157_43158_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_3MBPP1_COMBINED_CATGGC_R1_pombe_plus.bw", "7772_7157_43158_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_3MBPP1_COMBINED_CATGGC_R1_pombe_minus.bw", "CDK9as_dis2ts_3MBPP1_30C_combined", NFs[20,'Spikereads']/100000))

wigsetRPM = rbind(c("5993_7157_26977_HJ72CBGXX_pombe_WT_COMBINED_DMSO_CAGATC_R1.fastq_pombe_plus.bw", "5993_7157_26977_HJ72CBGXX_pombe_WT_COMBINED_DMSO_CAGATC_R1.fastq_pombe_minus.bw", "WT_DMSO_combined", RPM_NFs[1,'Mappedreads']/1000000),
               c("5993_7157_26978_HJ72CBGXX_pombe_WT_COMBINED_3MB-PPI_ACTTGA_R1.fastq_pombe_plus.bw", "5993_7157_26978_HJ72CBGXX_pombe_WT_COMBINED_3MB-PPI_ACTTGA_R1.fastq_pombe_minus.bw", "WT_3MBPP1_combined", RPM_NFs[2,'Mappedreads']/1000000),
               c("7772_7157_43149_HVLYCBGXY_pombe_18C_Cdk9as_DMSO_COMBINED_GGCTAC_R1_pombe_plus.bw", "7772_7157_43149_HVLYCBGXY_pombe_18C_Cdk9as_DMSO_COMBINED_GGCTAC_R1_pombe_minus.bw", "Cdk9as_DMSO_18C_combined", RPM_NFs[11,'Mappedreads']/1000000),
               c("7772_7157_43150_HVLYCBGXY_pombe_18C_Cdk9as_3MBPP1_COMBINED_CTTGTA_R1_pombe_plus.bw", "7772_7157_43150_HVLYCBGXY_pombe_18C_Cdk9as_3MBPP1_COMBINED_CTTGTA_R1_pombe_minus.bw", "Cdk9as_3MBPP1_18C_combined", RPM_NFs[12,'Mappedreads']/1000000),
               c("7772_7157_43151_HVLYCBGXY_pombe_18C_dis2ts_DMSO_COMBINED_ATTCCT_R1_pombe_plus.bw", "7772_7157_43151_HVLYCBGXY_pombe_18C_dis2ts_DMSO_COMBINED_ATTCCT_R1_pombe_minus.bw", "dis2ts_DMSO_18C_combined", RPM_NFs[13,'Mappedreads']/1000000),
               c("7772_7157_43152_HVLYCBGXY_pombe_18C_dis2ts_3MBPP1_COMBINED_CAAAAG_R1_pombe_plus.bw", "7772_7157_43152_HVLYCBGXY_pombe_18C_dis2ts_3MBPP1_COMBINED_CAAAAG_R1_pombe_minus.bw", "dis2ts_3MBPP1_18C_combined", RPM_NFs[14,'Mappedreads']/1000000),
               c("7772_7157_43153_HVLYCBGXY_pombe_30C_dis2ts_DMSO_COMBINED_CAACTA_R1_pombe_plus.bw", "7772_7157_43153_HVLYCBGXY_pombe_30C_dis2ts_DMSO_COMBINED_CAACTA_R1_pombe_minus.bw", "dis2ts_DMSO_30C_combined", RPM_NFs[15,'Mappedreads']/1000000),
               c("7772_7157_43154_HVLYCBGXY_pombe_30C_dis2ts_3MBPP1_COMBINED_CACCGG_R1_pombe_plus.bw", "7772_7157_43154_HVLYCBGXY_pombe_30C_dis2ts_3MBPP1_COMBINED_CACCGG_R1_pombe_minus.bw", "dis2ts_3MBPP1_30C_combined", RPM_NFs[16,'Mappedreads']/1000000),
               c("7772_7157_43155_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_DMSO_COMBINED_CACGAT_R1_pombe_plus.bw", "7772_7157_43155_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_DMSO_COMBINED_CACGAT_R1_pombe_minus.bw", "CDK9as_dis2ts_DMSO_18C_combined", RPM_NFs[17,'Mappedreads']/1000000),
               c("7772_7157_43156_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_3MBPP1_COMBINED_CACTCA_R1_pombe_plus.bw", "7772_7157_43156_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_3MBPP1_COMBINED_CACTCA_R1_pombe_minus.bw", "CDK9as_dis2ts_3MBPP1_18C_combined", RPM_NFs[18,'Mappedreads']/1000000),
               c("7772_7157_43157_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_DMSO_COMBINED_CAGGCG_R1_pombe_plus.bw", "7772_7157_43157_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_DMSO_COMBINED_CAGGCG_R1_pombe_minus.bw", "CDK9as_dis2ts_DMSO_30C_combined", RPM_NFs[19,'Mappedreads']/1000000),
               c("7772_7157_43158_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_3MBPP1_COMBINED_CATGGC_R1_pombe_plus.bw", "7772_7157_43158_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_3MBPP1_COMBINED_CATGGC_R1_pombe_minus.bw", "CDK9as_dis2ts_3MBPP1_30C_combined", RPM_NFs[20,'Mappedreads']/1000000))

wigsetReps = rbind(c("5993_7157_26971_HJ72CBGXX_pombe_WT_rep1_DMSO_ATCACG_R1.fastq_pombe_plus.bw", "5993_7157_26971_HJ72CBGXX_pombe_WT_rep1_DMSO_ATCACG_R1.fastq_pombe_minus.bw", "WT_DMSOr1", repNFs[1,'Spikereads']/100000),
               c("5993_7157_26977_HJ72CBGXX_pombe_WT_rep2_DMSO_CAGATC_R1.fastq_pombe_plus.bw", "5993_7157_26977_HJ72CBGXX_pombe_WT_rep2_DMSO_CAGATC_R1.fastq_pombe_minus.bw", "WT_DMSOr2", repNFs[2,'Spikereads']/100000),
               c("5993_7157_26972_HJ72CBGXX_pombe_WT_rep1_3MB-PPI_CGATGT_R1.fastq_pombe_plus.bw", "5993_7157_26972_HJ72CBGXX_pombe_WT_rep1_3MB-PPI_CGATGT_R1.fastq_pombe_minus.bw", "WT_3MBPP1r1", repNFs[3,'Spikereads']/100000),
               c("5993_7157_26978_HJ72CBGXX_pombe_WT_rep2_3MB-PPI_ACTTGA_R1.fastq_pombe_plus.bw", "5993_7157_26978_HJ72CBGXX_pombe_WT_rep2_3MB-PPI_ACTTGA_R1.fastq_pombe_minus.bw", "WT_3MBPP1r2", repNFs[4,'Spikereads']/100000),
               c("7772_7157_43139_HVLYCBGXY_pombe_18C_Cdk9as_DMSO_rep1_ATCACG_R1_pombe_plus.bw", "7772_7157_43139_HVLYCBGXY_pombe_18C_Cdk9as_DMSO_rep1_ATCACG_R1_pombe_minus.bw", "Cdk9as_DMSO_18Cr1", repNFs[21,'Spikereads']/100000),
               c("7772_7157_43149_HVLYCBGXY_pombe_18C_Cdk9as_DMSO_rep2_GGCTAC_R1_pombe_plus.bw", "7772_7157_43149_HVLYCBGXY_pombe_18C_Cdk9as_DMSO_rep2_GGCTAC_R1_pombe_minus.bw", "Cdk9as_DMSO_18Cr2", repNFs[22,'Spikereads']/100000),
               c("7772_7157_43140_HVLYCBGXY_pombe_18C_Cdk9as_3MBPP1_rep1_CGATGT_R1_pombe_plus.bw", "7772_7157_43140_HVLYCBGXY_pombe_18C_Cdk9as_3MBPP1_rep1_CGATGT_R1_pombe_minus.bw", "Cdk9as_3MBPP1_18Cr1", repNFs[23,'Spikereads']/100000),
               c("7772_7157_43150_HVLYCBGXY_pombe_18C_Cdk9as_3MBPP1_rep2_CTTGTA_R1_pombe_plus.bw", "7772_7157_43150_HVLYCBGXY_pombe_18C_Cdk9as_3MBPP1_rep2_CTTGTA_R1_pombe_minus.bw", "Cdk9as_3MBPP1_18Cr2", repNFs[24,'Spikereads']/100000),
               c("7772_7157_43141_HVLYCBGXY_pombe_18C_dis2ts_DMSO_rep1_TTAGGC_R1_pombe_plus.bw", "7772_7157_43141_HVLYCBGXY_pombe_18C_dis2ts_DMSO_rep1_TTAGGC_R1_pombe_minus.bw", "dis2ts_DMSO_18Cr1", repNFs[25,'Spikereads']/100000),
               c("7772_7157_43151_HVLYCBGXY_pombe_18C_dis2ts_DMSO_rep2_ATTCCT_R1_pombe_plus.bw", "7772_7157_43151_HVLYCBGXY_pombe_18C_dis2ts_DMSO_rep2_ATTCCT_R1_pombe_minus.bw", "dis2ts_DMSO_18Cr2", repNFs[26,'Spikereads']/100000),
               c("7772_7157_43142_HVLYCBGXY_pombe_18C_dis2ts_3MBPP1_rep1_TGACCA_R1_pombe_plus.bw", "7772_7157_43142_HVLYCBGXY_pombe_18C_dis2ts_3MBPP1_rep1_TGACCA_R1_pombe_minus.bw", "dis2ts_3MBPP1_18Cr1", repNFs[27,'Spikereads']/100000),
               c("7772_7157_43152_HVLYCBGXY_pombe_18C_dis2ts_3MBPP1_rep2_CAAAAG_R1_pombe_plus.bw", "7772_7157_43152_HVLYCBGXY_pombe_18C_dis2ts_3MBPP1_rep2_CAAAAG_R1_pombe_minus.bw", "dis2ts_3MBPP1_18Cr2", repNFs[28,'Spikereads']/100000),
               c("7772_7157_43143_HVLYCBGXY_pombe_30C_dis2ts_DMSO_rep1_ACAGTG_R1_pombe_plus.bw", "7772_7157_43143_HVLYCBGXY_pombe_30C_dis2ts_DMSO_rep1_ACAGTG_R1_pombe_minus.bw", "dis2ts_DMSO_30Cr1", repNFs[29,'Spikereads']/100000),
               c("7772_7157_43153_HVLYCBGXY_pombe_30C_dis2ts_DMSO_rep2_CAACTA_R1_pombe_plus.bw", "7772_7157_43153_HVLYCBGXY_pombe_30C_dis2ts_DMSO_rep2_CAACTA_R1_pombe_minus.bw", "dis2ts_DMSO_30Cr2", repNFs[30,'Spikereads']/100000),
               c("7772_7157_43144_HVLYCBGXY_pombe_30C_dis2ts_3MBPP1_rep1_GCCAAT_R1_pombe_plus.bw", "7772_7157_43144_HVLYCBGXY_pombe_30C_dis2ts_3MBPP1_rep1_GCCAAT_R1_pombe_minus.bw", "dis2ts_3MBPP1_30Cr1", repNFs[31,'Spikereads']/100000),
               c("7772_7157_43154_HVLYCBGXY_pombe_30C_dis2ts_3MBPP1_rep2_CACCGG_R1_pombe_plus.bw", "7772_7157_43154_HVLYCBGXY_pombe_30C_dis2ts_3MBPP1_rep2_CACCGG_R1_pombe_minus.bw", "dis2ts_3MBPP1_30Cr2", repNFs[32,'Spikereads']/100000),
               c("7772_7157_43145_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_DMSO_rep1_CAGATC_R1_pombe_plus.bw", "7772_7157_43145_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_DMSO_rep1_CAGATC_R1_pombe_minus.bw", "CDK9as_dis2ts_DMSO_18Cr1", repNFs[33,'Spikereads']/100000),
               c("7772_7157_43155_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_DMSO_rep2_CACGAT_R1_pombe_plus.bw", "7772_7157_43155_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_DMSO_rep2_CACGAT_R1_pombe_minus.bw", "CDK9as_dis2ts_DMSO_18Cr2", repNFs[34,'Spikereads']/100000),
               c("7772_7157_43146_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_3MBPP1_rep1_ACTTGA_R1_pombe_plus.bw", "7772_7157_43146_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_3MBPP1_rep1_ACTTGA_R1_pombe_minus.bw", "CDK9as_dis2ts_3MBPP1_18Cr1", repNFs[35,'Spikereads']/100000),
               c("7772_7157_43156_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_3MBPP1_rep2_CACTCA_R1_pombe_plus.bw", "7772_7157_43156_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_3MBPP1_rep2_CACTCA_R1_pombe_minus.bw", "CDK9as_dis2ts_3MBPP1_18Cr2", repNFs[36,'Spikereads']/100000),
               c("7772_7157_43147_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_DMSO_rep1_GATCAG_R1_pombe_plus.bw", "7772_7157_43147_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_DMSO_rep1_GATCAG_R1_pombe_minus.bw", "CDK9as_dis2ts_DMSO_30Cr1", repNFs[37,'Spikereads']/100000),
               c("7772_7157_43157_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_DMSO_rep2_CAGGCG_R1_pombe_plus.bw", "7772_7157_43157_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_DMSO_rep2_CAGGCG_R1_pombe_minus.bw", "CDK9as_dis2ts_DMSO_30Cr2", repNFs[38,'Spikereads']/100000),
               c("7772_7157_43148_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_3MBPP1_rep1_TAGCTT_R1_pombe_plus.bw", "7772_7157_43148_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_3MBPP1_rep1_TAGCTT_R1_pombe_minus.bw", "CDK9as_dis2ts_3MBPP1_30Cr1", repNFs[39,'Spikereads']/100000),
               c("7772_7157_43158_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_3MBPP1_rep2_CATGGC_R1_pombe_plus.bw", "7772_7157_43158_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_3MBPP1_rep2_CATGGC_R1_pombe_minus.bw", "CDK9as_dis2ts_3MBPP1_30Cr2", repNFs[40,'Spikereads']/100000))

wigsetReps_RPM = rbind(c("5993_7157_26971_HJ72CBGXX_pombe_WT_rep1_DMSO_ATCACG_R1.fastq_pombe_plus.bw", "5993_7157_26971_HJ72CBGXX_pombe_WT_rep1_DMSO_ATCACG_R1.fastq_pombe_minus.bw", "WT_DMSOr1", RPM_repNFs[1,'Mappedreads']/1000000),
                   c("5993_7157_26977_HJ72CBGXX_pombe_WT_rep2_DMSO_CAGATC_R1.fastq_pombe_plus.bw", "5993_7157_26977_HJ72CBGXX_pombe_WT_rep2_DMSO_CAGATC_R1.fastq_pombe_minus.bw", "WT_DMSOr2", RPM_repNFs[2,'Mappedreads']/1000000),
                   c("5993_7157_26972_HJ72CBGXX_pombe_WT_rep1_3MB-PPI_CGATGT_R1.fastq_pombe_plus.bw", "5993_7157_26972_HJ72CBGXX_pombe_WT_rep1_3MB-PPI_CGATGT_R1.fastq_pombe_minus.bw", "WT_3MBPP1r1", RPM_repNFs[3,'Mappedreads']/1000000),
                   c("5993_7157_26978_HJ72CBGXX_pombe_WT_rep2_3MB-PPI_ACTTGA_R1.fastq_pombe_plus.bw", "5993_7157_26978_HJ72CBGXX_pombe_WT_rep2_3MB-PPI_ACTTGA_R1.fastq_pombe_minus.bw", "WT_3MBPP1r2", RPM_repNFs[4,'Mappedreads']/1000000),
                   c("7772_7157_43139_HVLYCBGXY_pombe_18C_Cdk9as_DMSO_rep1_ATCACG_R1_pombe_plus.bw", "7772_7157_43139_HVLYCBGXY_pombe_18C_Cdk9as_DMSO_rep1_ATCACG_R1_pombe_minus.bw", "Cdk9as_DMSO_18Cr1", RPM_repNFs[21,'Mappedreads']/1000000),
                   c("7772_7157_43149_HVLYCBGXY_pombe_18C_Cdk9as_DMSO_rep2_GGCTAC_R1_pombe_plus.bw", "7772_7157_43149_HVLYCBGXY_pombe_18C_Cdk9as_DMSO_rep2_GGCTAC_R1_pombe_minus.bw", "Cdk9as_DMSO_18Cr2", RPM_repNFs[22,'Mappedreads']/1000000),
                   c("7772_7157_43140_HVLYCBGXY_pombe_18C_Cdk9as_3MBPP1_rep1_CGATGT_R1_pombe_plus.bw", "7772_7157_43140_HVLYCBGXY_pombe_18C_Cdk9as_3MBPP1_rep1_CGATGT_R1_pombe_minus.bw", "Cdk9as_3MBPP1_18Cr1", RPM_repNFs[23,'Mappedreads']/1000000),
                   c("7772_7157_43150_HVLYCBGXY_pombe_18C_Cdk9as_3MBPP1_rep2_CTTGTA_R1_pombe_plus.bw", "7772_7157_43150_HVLYCBGXY_pombe_18C_Cdk9as_3MBPP1_rep2_CTTGTA_R1_pombe_minus.bw", "Cdk9as_3MBPP1_18Cr2", RPM_repNFs[24,'Mappedreads']/1000000),
                   c("7772_7157_43141_HVLYCBGXY_pombe_18C_dis2ts_DMSO_rep1_TTAGGC_R1_pombe_plus.bw", "7772_7157_43141_HVLYCBGXY_pombe_18C_dis2ts_DMSO_rep1_TTAGGC_R1_pombe_minus.bw", "dis2ts_DMSO_18Cr1", RPM_repNFs[25,'Mappedreads']/1000000),
                   c("7772_7157_43151_HVLYCBGXY_pombe_18C_dis2ts_DMSO_rep2_ATTCCT_R1_pombe_plus.bw", "7772_7157_43151_HVLYCBGXY_pombe_18C_dis2ts_DMSO_rep2_ATTCCT_R1_pombe_minus.bw", "dis2ts_DMSO_18Cr2", RPM_repNFs[26,'Mappedreads']/1000000),
                   c("7772_7157_43142_HVLYCBGXY_pombe_18C_dis2ts_3MBPP1_rep1_TGACCA_R1_pombe_plus.bw", "7772_7157_43142_HVLYCBGXY_pombe_18C_dis2ts_3MBPP1_rep1_TGACCA_R1_pombe_minus.bw", "dis2ts_3MBPP1_18Cr1", RPM_repNFs[27,'Mappedreads']/1000000),
                   c("7772_7157_43152_HVLYCBGXY_pombe_18C_dis2ts_3MBPP1_rep2_CAAAAG_R1_pombe_plus.bw", "7772_7157_43152_HVLYCBGXY_pombe_18C_dis2ts_3MBPP1_rep2_CAAAAG_R1_pombe_minus.bw", "dis2ts_3MBPP1_18Cr2", RPM_repNFs[28,'Mappedreads']/1000000),
                   c("7772_7157_43143_HVLYCBGXY_pombe_30C_dis2ts_DMSO_rep1_ACAGTG_R1_pombe_plus.bw", "7772_7157_43143_HVLYCBGXY_pombe_30C_dis2ts_DMSO_rep1_ACAGTG_R1_pombe_minus.bw", "dis2ts_DMSO_30Cr1", RPM_repNFs[29,'Mappedreads']/1000000),
                   c("7772_7157_43153_HVLYCBGXY_pombe_30C_dis2ts_DMSO_rep2_CAACTA_R1_pombe_plus.bw", "7772_7157_43153_HVLYCBGXY_pombe_30C_dis2ts_DMSO_rep2_CAACTA_R1_pombe_minus.bw", "dis2ts_DMSO_30Cr2", RPM_repNFs[30,'Mappedreads']/1000000),
                   c("7772_7157_43144_HVLYCBGXY_pombe_30C_dis2ts_3MBPP1_rep1_GCCAAT_R1_pombe_plus.bw", "7772_7157_43144_HVLYCBGXY_pombe_30C_dis2ts_3MBPP1_rep1_GCCAAT_R1_pombe_minus.bw", "dis2ts_3MBPP1_30Cr1", RPM_repNFs[31,'Mappedreads']/1000000),
                   c("7772_7157_43154_HVLYCBGXY_pombe_30C_dis2ts_3MBPP1_rep2_CACCGG_R1_pombe_plus.bw", "7772_7157_43154_HVLYCBGXY_pombe_30C_dis2ts_3MBPP1_rep2_CACCGG_R1_pombe_minus.bw", "dis2ts_3MBPP1_30Cr2", RPM_repNFs[32,'Mappedreads']/1000000),
                   c("7772_7157_43145_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_DMSO_rep1_CAGATC_R1_pombe_plus.bw", "7772_7157_43145_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_DMSO_rep1_CAGATC_R1_pombe_minus.bw", "CDK9as_dis2ts_DMSO_18Cr1", RPM_repNFs[33,'Mappedreads']/1000000),
                   c("7772_7157_43155_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_DMSO_rep2_CACGAT_R1_pombe_plus.bw", "7772_7157_43155_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_DMSO_rep2_CACGAT_R1_pombe_minus.bw", "CDK9as_dis2ts_DMSO_18Cr2", RPM_repNFs[34,'Mappedreads']/1000000),
                   c("7772_7157_43146_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_3MBPP1_rep1_ACTTGA_R1_pombe_plus.bw", "7772_7157_43146_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_3MBPP1_rep1_ACTTGA_R1_pombe_minus.bw", "CDK9as_dis2ts_3MBPP1_18Cr1", RPM_repNFs[35,'Mappedreads']/1000000),
                   c("7772_7157_43156_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_3MBPP1_rep2_CACTCA_R1_pombe_plus.bw", "7772_7157_43156_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_3MBPP1_rep2_CACTCA_R1_pombe_minus.bw", "CDK9as_dis2ts_3MBPP1_18Cr2", RPM_repNFs[36,'Mappedreads']/1000000),
                   c("7772_7157_43147_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_DMSO_rep1_GATCAG_R1_pombe_plus.bw", "7772_7157_43147_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_DMSO_rep1_GATCAG_R1_pombe_minus.bw", "CDK9as_dis2ts_DMSO_30Cr1", RPM_repNFs[37,'Mappedreads']/1000000),
                   c("7772_7157_43157_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_DMSO_rep2_CAGGCG_R1_pombe_plus.bw", "7772_7157_43157_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_DMSO_rep2_CAGGCG_R1_pombe_minus.bw", "CDK9as_dis2ts_DMSO_30Cr2", RPM_repNFs[38,'Mappedreads']/1000000),
                   c("7772_7157_43148_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_3MBPP1_rep1_TAGCTT_R1_pombe_plus.bw", "7772_7157_43148_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_3MBPP1_rep1_TAGCTT_R1_pombe_minus.bw", "CDK9as_dis2ts_3MBPP1_30Cr1", RPM_repNFs[39,'Mappedreads']/1000000),
                   c("7772_7157_43158_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_3MBPP1_rep2_CATGGC_R1_pombe_plus.bw", "7772_7157_43158_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_3MBPP1_rep2_CATGGC_R1_pombe_minus.bw", "CDK9as_dis2ts_3MBPP1_30Cr2", RPM_repNFs[40,'Mappedreads']/1000000))

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
  #result <- xyplot(mean ~ x | factor(background), data = df,
  result <- xyplot(mean ~ x | factor(strain) + factor(treatment), data = df, 
                   #group = factor(treatment, labels = c("DMSO", "3MB-PP1")), 
                   group = factor(temp, labels = c("30 deg C", "18 deg C")), 
                   scales = list(tck=c(1,0),alternating = c(1,1),x=list(relation='free',axs='i'),
                                 y=list(relation='free',axs='i')),
                   key=list(corner=c(0.98,0.85),padding.text=3,
                            #text=list(c("DMSO", "3MB-PP1"),col=c("black", "black"), cex=0.8, font=2),
                            text=list(c("30 deg C", "18 deg C"),col=c("black", "black"), cex=0.8, font=2),
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
## Same as above function, but used when you want to compare two samples directly. 
lattice_meta.proSeq.directCompare = function(filename= "test.pdf", df, ylim = c(0, 10), main = NULL, 
                               xlab = "Distance to TSS (bp)",  ylab = "Median PRO-seq intensity", 
                               labs = c("sample1", "sample2")){
  pdf(paste(fig_dir, filename, sep=""), width = 5, height = 5)
  result <- xyplot(mean ~ x | factor(treatment), data = df,
                   group = factor(sample, labels = labs),
                   scales = list(tck=c(1,0),alternating = c(1,1),x=list(relation='free',axs='i'),
                                 y=list(relation='free',axs='i')),
                   key=list(corner=c(0.98,0.85),padding.text=3,
                              text=list(labs ,col=c("black", "black"), cex=0.8, font=2),
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

lattice_Scaledmeta.proSeq = function(filename= "test.pdf", df, ylim = c(0, 10), main = NULL, groupFact = "strain",
                               xlab = "Distance to gene boundaries (bp)",  ylab = "Median PRO-seq intensity", labs = c("sample1", "sample2")){
  pdf(paste(fig_dir, filename, sep=""), width = 8, height = 5)
  result <- xyplot(mean ~ x | factor(treatment), data = df,
  #result <- xyplot(mean ~ x | factor(strain) + factor(treatment), data = df, 
                   group = factor(sample, labels = labs),
                   #group = factor(background, labels = labs), 
                   #group = factor(treatment, labels = c("DMSO", "3MB-PP1")), 
                   #group = factor(temp, labels = c("30 deg C", "18 deg C")),
                   scales = list(tck=c(1,0),alternating = c(1,1),
                                 x=list(relation='free',axs='i', labels = c('-1000', 'TSS', "+300", "-300", "CPS", "1000"), at = c(1, 101, 131, 191, 221, 320)),
                                 y=list(relation='free',axs='i')),
                   key=list(corner=c(0.98,0.85),padding.text=3,
                            #text=list(c("DMSO", "3MB-PP1"),col=c("black", "black"), cex=0.8, font=2),
                            text=list(labs ,col=c("black", "black"), cex=0.8, font=2),
                            rectangles=list(col=c(rgb(0,0,0.5), rgb(0.5,0,0)), size=2.7, height=0.8, border='white')),
                   type = 'l', 
                   ylim = ylim,
                   col = c( rgb(0,0,0.5,0.9), rgb(0.5,0,0,0.9)),# rgb(0,0.5,0,0.9)),  ## just add or take away colors to match number of samples in each panel
                   fill = c(rgb(0,0,0.5,0.3), rgb(0.5,0,0,0.3)),# rgb(0,0.5,0,0.3)),
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
  #BG_Indx = 0 # counter for tracking genotype, which changes every other sample.
  for (i in 1:N){
    if (i == ceiling(N/2)){cat("* 50% complete ... \n")}
    cat("* loading", i, "\n")
    wigs = load.wigset(wigset, i) #this should give all the bigwig files
    cat("* generating and combining all sample scaled meta-plot data ...\n")
    meta= meta.subsample.scaled(bed, bigWig.plus = wigs[[1]], bigWig.minus = wigs[[2]], 1000, 300, 10, do.sum = T)
    metaNorm = meta.normalize(result = meta, scaleFactor = 1/as.numeric(wigs[[4]]))
    sampleNames = c(sampleNames, wigs[[3]])
    sampleVal = i ### index for later retrieval of sample name from sampleNames list.
    #if ((-1)**i == -1){BG_Indx = BG_Indx + 1}
    #else{BG_Indx = BG_Indx}
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
    if(grepl("r1", wigs[[3]])){Rep = 1}
    else if (grepl("r2", wigs[[3]])){Rep = 2}
    else {Rep = 3}
    print(strain)
    print(temp)
    xAxis = seq(from = 1, to = 320)
    sample_df = cbind(xAxis, metaNorm[[4]], metaNorm[[3]], metaNorm[[2]], sampleVal, treatment, strain, temp, Rep)
    df = rbind.data.frame(df, sample_df)
  }
  colnames(df) = c("x", "mean", "lower", "upper", "sample", 'treatment', "strain", "temp", "replicate")
  # lastly I will replace values in the 'sample', 'treatment', and 'background' columns with actual names
  for (ii in 1:length(sampleNames)){
    df$sample[df$sample == ii] <- sampleNames[ii]
  }
  df$treatment[df$treatment == 0] <- "DMSO"
  df$treatment[df$treatment == 1] <- "3-MB_PP1"
  df$strain[df$strain == 1] <- "WT"
  df$strain[df$strain == 2] <- "CDK9_as"
  df$strain[df$strain == 6] <- "dis2_ts"
  df$strain[df$strain == 7] <- "CDK9as_dis2ts"
  return(df)
}

## formatting meta plots for lattice (not scaled Metaplots)
FormatMetaData_lattice <- function(wigset, bed, halfWindow, step, FractMax = FALSE){
  N = dim(wigset)[1]
  df = data.frame(matrix(ncol = 9, nrow = 0), stringsAsFactors = F)
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
    if(grepl("r1", wigs[[3]])){Rep = 1}
    else if (grepl("r2", wigs[[3]])){Rep = 2}
    else {Rep = 3}
    print(strain)
    print(temp)
    xAxis = seq(from = -halfWindow, to = halfWindow, by = step)
    if (FractMax){
      FracMax_mean = metaNorm[[4]]/(max(metaNorm[[4]]))
      FracMax_uQ = metaNorm[[3]]/(max(metaNorm[[4]]))
      FracMax_lQ = metaNorm[[2]]/(max(metaNorm[[4]]))
      sample_df = cbind(xAxis, FracMax_mean, FracMax_uQ, FracMax_lQ, sampleVal, treatment, strain, temp, Rep)
    }
    else{
      sample_df = cbind(xAxis, metaNorm[[4]], metaNorm[[3]], metaNorm[[2]], sampleVal, treatment, strain, temp, Rep)
    }
    df = rbind.data.frame(df, sample_df)
  }
  colnames(df) = c("x", "mean", "lower", "upper", "sample", 'treatment', "strain", "temp", "replicate")
  # lastly I will replace values in the 'sample', 'treatment', and 'background' columns with actual names
  for (ii in 1:length(sampleNames)){
    df$sample[df$sample == ii] <- sampleNames[ii]
  }
  df$treatment[df$treatment == 0] <- "DMSO"
  df$treatment[df$treatment == 1] <- "3-MB_PP1"
  df$strain[df$strain == 1] <- "WT"
  df$strain[df$strain == 2] <- "CDK9_as"
  df$strain[df$strain == 6] <- "dis2_ts"
  df$strain[df$strain == 7] <- "CDK9as_dis2ts"
  return(df)
}
#####################################################################################################################
#plot scaled meta-plots with lattice. 
scaledMetaData = FormatScaledMetaData_lattice(wigset, bed)

## from the full table, pull out only the samples I want to compare
# compare Cdk9as_DMSO @18
dir.create(paste(fig_dir, "spikeNorm/", sep = ""))
# Cdk9as (18C) vs Cdk9as_Dis2-11 (18C)
cdk9as_Vs_cdk9asDis2_18C_DMSO = (scaledMetaData$sample == "Cdk9as_DMSO_18C_combined" |  scaledMetaData$sample == "CDK9as_dis2ts_DMSO_18C_combined")
lattice_Scaledmeta.proSeq(filename = "spikeNorm/Cdk9as_vs_Cdk9asDis2ts_18C_DMSO_scaledMeta_kbsep_SpikeNorm.pdf", 
                          df = scaledMetaData[cdk9as_Vs_cdk9asDis2_18C_DMSO,], ylim = c(0, 60), main = "scaled Gene Length", 
                          xlab = "Distance to Gene boundaries (bp)", ylab = "Median PRO-seq intensity", labs = c("cdk9as_dis2-11_18C", "ckd9as_18C"))
# Cdk9as (18C) vs Cdk9as_Dis2-11 (30C)
cdk9as_Vs_cdk9asDis2_30C_DMSO = (scaledMetaData$sample == "Cdk9as_DMSO_18C_combined" |  scaledMetaData$sample == "CDK9as_dis2ts_DMSO_30C_combined")
lattice_Scaledmeta.proSeq(filename = "spikeNorm/Cdk9as_18C_vs_Cdk9asDis2ts_30C_DMSO_scaledMeta_kbsep_SpikeNorm.pdf", 
                          df = scaledMetaData[cdk9as_Vs_cdk9asDis2_30C_DMSO,], ylim = c(0, 60), main = "scaled Gene Length", 
                          xlab = "Distance to Gene boundaries (bp)", ylab = "Median PRO-seq intensity", labs = c("cdk9as_dis2-11_30C", "ckd9as_18C"))

# WT (30C) vs Dis2-11 (18C)
WT_30CvsDis2_18C_DMSO = (scaledMetaData$sample == "WT_DMSO_combined" |  scaledMetaData$sample == "dis2ts_DMSO_18C_combined")
lattice_Scaledmeta.proSeq(filename = "spikeNorm/WT_30CvsDis2_18C_DMSO_scaledMeta_kbsep_SpikeNorm.pdf", 
                          df = scaledMetaData[WT_30CvsDis2_18C_DMSO,], ylim = c(0, 45), main = "scaled Gene Length", 
                          xlab = "Distance to Gene boundaries (bp)", ylab = "Median PRO-seq intensity", labs = c("Dis2-11_18C", "WT_30C"))
#  Cdk9as (18C) vs Cdk9as_Dis2-11 (18C) both treated with 3MBPP1
cdk9as_Vs_cdk9asDis2_18C_3MBPP1 = (scaledMetaData$sample == "Cdk9as_3MBPP1_18C_combined" |  scaledMetaData$sample == "CDK9as_dis2ts_3MBPP1_18C_combined")
lattice_Scaledmeta.proSeq(filename = "spikeNorm/Cdk9as_vs_Cdk9asDis2ts_18C_3MBPP1_scaledMeta_kbsep_SpikeNorm.pdf", 
                          df = scaledMetaData[cdk9as_Vs_cdk9asDis2_18C_3MBPP1,], ylim = c(0, 45), main = "scaled Gene Length", 
                          xlab = "Distance to Gene boundaries (bp)", ylab = "Median PRO-seq intensity", labs = c("cdk9as_dis2-11_18C_3MBPP1", "ckd9as_18C_3MBPP1"))
#  Cdk9as (18C) vs Dis2-11 (18C) both DMSO
cdk9as_Vs_Dis2_18C_DMSO = (scaledMetaData$sample == "Cdk9as_DMSO_18C_combined" |  scaledMetaData$sample == "dis2ts_DMSO_18C_combined")
lattice_Scaledmeta.proSeq(filename = "spikeNorm/Cdk9as_vs_Dis2ts_18C_DMSO_scaledMeta_kbsep_SpikeNorm.pdf", 
                          df = scaledMetaData[cdk9as_Vs_Dis2_18C_DMSO,], ylim = c(0, 45), main = "scaled Gene Length", 
                          xlab = "Distance to Gene boundaries (bp)", ylab = "Median PRO-seq intensity", labs = c("dis2ts_DMSO_18C_combined", "Cdk9as_DMSO_18C_combined"))


#####################################################################################################################
### Plots around TSS and CPS (for KBsep genes) 
# TSS
TSSmetaData = FormatMetaData_lattice(wigset, bed, halfWindow = 1000, step = 10)
lattice_meta.proSeq(filename = "SP_TSS_experiments_byTemp_Meta.pdf", df = TSSmetaData, ylim = c(0, 75), main = 'Distance to TSS', 
                    xlab = "Distance to TSS (bp)",  ylab = "Median PRO-seq intensity")
# CPS
CPSmetaData = FormatMetaData_lattice(wigset, bed[, c(1,3,2,4,5,6)], halfWindow = 1000, step = 10)
lattice_meta.proSeq(filename = "SP_CPS_experiments_byTemp_Meta.pdf", df = CPSmetaData, ylim = c(0, 75), main = 'Distance to CPS', 
                    xlab = "Distance to CPS (bp)",  ylab = "Median PRO-seq intensity")
########################
#section added 06-22-17 
## same as above, but selected comparisons using fraction of maximum as a normalization.
CPSmetaData_FM = FormatMetaData_lattice(wigset, bed[, c(1,3,2,4,5,6)], halfWindow = 1000, step = 10, FractMax = T)
WT_30CvsDis2_18C_FM = (CPSmetaData_FM$sample == "WT_DMSO_combined" |  CPSmetaData_FM$sample == "dis2ts_DMSO_18C_combined")
lattice_meta.proSeq.directCompare(filename = "SP_CPS_WTvsDis2_11_18C_CombinedReps_Meta.pdf", df = CPSmetaData_FM[WT_30CvsDis2_18C_FM,], ylim = c(0, 1.2), main = 'Distance to CPS', 
                                  xlab = "Distance to CPS (bp)",  ylab = "Fraction of Max Signal", labs = c("WT", "dis2-11 18C"))
WT_30CvsDis2_30C_FM = (CPSmetaData_FM$sample == "WT_DMSO_combined" |  CPSmetaData_FM$sample == "dis2ts_DMSO_30C_combined")
lattice_meta.proSeq.directCompare(filename = "SP_CPS_WTvsDis2_11_30C_CombinedReps_Meta.pdf", df = CPSmetaData_FM[WT_30CvsDis2_30C_FM,], ylim = c(0, 1.2), main = 'Distance to CPS', 
                                  xlab = "Distance to CPS (bp)",  ylab = "Fraction of Max Signal", labs = c("WT", "dis2-11 30C"))

#####################################################################################################################
## prepare plots for individual replicates separately 
#####################################################################################################################
## make metaplot comparisons of desired strains with separate plots for each replicated
lattice_ScaledmetaReps.proSeq = function(filename= "test.pdf", df, ylim = c(0, 10), main = NULL, groupFact = "strain",
                                     xlab = "Distance to gene boundaries (bp)",  ylab = "Median PRO-seq intensity", labs = c("sample1", "sample2")){
  pdf(paste(fig_dir, filename, sep=""), width = 5, height = 8)
  result <- xyplot(mean ~ x | factor(replicate, labels = c("rep 1", "rep 2")), data = df,
                   #result <- xyplot(mean ~ x | factor(strain) + factor(treatment), data = df, 
                   group = factor(df[[groupFact]]), 
                   #group = factor(treatment, labels = c("DMSO", "3MB-PP1")), 
                   #group = factor(temp, labels = c("30 deg C", "18 deg C")),
                   scales = list(tck=c(1,0),alternating = c(1,1),
                                 x=list(relation='free',axs='i', labels = c('-1000', 'TSS', "+300", "-300", "CPS", "1000"), at = c(1, 101, 131, 191, 221, 320)),
                                 y=list(relation='free',axs='i')),
                   key=list(corner=c(0.98,0.85),padding.text=3,
                            #text=list(c("DMSO", "3MB-PP1"),col=c("black", "black"), cex=0.8, font=2),
                            text=list(labs ,col=c("black", "black"), cex=0.8, font=2),
                            rectangles=list(col=c(rgb(0,0,0.5), rgb(0.5,0,0)), size=2.7, height=0.8, border='white')),
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

lattice_metaReps.proSeq = function(filename= "test.pdf", df, ylim = c(0, 10), main = NULL,  groupFact = "strain",
                               xlab = "Distance to TSS (bp)",  ylab = "Median PRO-seq intensity", labs = c("sample1", "sample2")){
  pdf(paste(fig_dir, filename, sep = ""), width = 3, height = 6)
  result <- xyplot(mean ~ x | factor(replicate, labels = c("rep 1", "rep 2")), data = df,
                   #result <- xyplot(mean ~ x | factor(strain) + factor(treatment), data = df, 
                   group = factor(df[[groupFact]]), 
                   #group = factor(temp, labels = c("30 deg C", "18 deg C")), 
                   scales = list(tck=c(1,0),alternating = c(1,1),x=list(relation='free',axs='i'),
                                 y=list(relation='free',axs='i')),
                   key=list(corner=c(0.98,0.85),padding.text=3,
                            text=list(labs, col=c("black", "black"), cex=0.8, font=2),
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
######################################################################################################
## functions below are for preparing Metaplots that show data as a fraction of peak height within the window. 
## This serves as a nice qualitative way of comparing two samples, without the need for normalization. 
## want to use to determine if, relative to amount of pol II signal, dis2-11 strains have longer elongation beyond CPS (i.e. termination defects)

WTseqplusObs = meta.subsample(bed, SCPbw, SCMbw, 300, 5, do.sum = T, at.TSS=T)
WTcapplusObs = meta.subsample(bed, SCcapPbw, SCcapMbw, 300, 5, do.sum = T, at.TSS=T)
WTseqplusAnn = meta.subsample(MatchedFullgList, SCPbw, SCMbw, 300, 5, do.sum = T, at.TSS=T)
WTcapplusAnn = meta.subsample(MatchedFullgList, SCcapPbw, SCcapMbw, 300, 5, do.sum = T, at.TSS=T)
WT_PROseq_FracMax_Obs_pl = WTseqplusObs[[4]]/(max(WTseqplusObs[[4]]))
WT_PROcap_FracMax_Obs_pl = WTcapplusObs[[4]]/(max(WTcapplusObs[[4]]))
WT_PROseq_FracMax_Ann_pl = WTseqplusAnn[[4]]/(max(WTseqplusAnn[[4]]))
WT_PROcap_FracMax_Ann_pl = WTcapplusAnn[[4]]/(max(WTcapplusAnn[[4]]))
base = seq(-295, 300,5) 
FracMaxTable = cbind(base, WT_PROseq_FracMax_Obs_pl, WT_PROcap_FracMax_Obs_pl, WT_PROseq_FracMax_Ann_pl, WT_PROcap_FracMax_Ann_pl)
# FracMax around Observed TSS
pdf(file = paste(directory, "PROcapVsSeq_300bp_ObservedTSS_FractionOfMaximum_meta.pdf", sep = ""), width = 5, height = 5)
plot(FracMaxTable[,1], FracMaxTable[,3], type = 'l', lwd = 4, col = 'red4',
     main="SC Observed TSS",  ylab = "Fraction of Maximum", xlab = "Distance to TSS",
     cex.main=2, cex.lab=1.5, cex.axis=1.5)
points(FracMaxTable[,1], FracMaxTable[,2], type = 'l', lwd = 4, col = 'blue4')
abline(v=c(0), lwd=2, lty=2)
legend('topleft', c("PRO-Cap", "PRO-seq"), text.col= c("sienna4", "blue4"), cex = 1.5)
dev.off()





#####################################################################################################################
scaledMetaData_IndivReplicates = FormatScaledMetaData_lattice(wigsetReps, bed)
dir.create(paste(fig_dir, "spikeNorm/individualReps/", sep = ""))
# Cdk9as (18C) vs Cdk9as_Dis2-11 (18C)
cdk9as_Vs_cdk9asDis2_18C_DMSO_reps = (scaledMetaData_IndivReplicates$sample == "Cdk9as_DMSO_18Cr1" | scaledMetaData_IndivReplicates$sample == "Cdk9as_DMSO_18Cr2" |  scaledMetaData_IndivReplicates$sample == "CDK9as_dis2ts_DMSO_18Cr1" | scaledMetaData_IndivReplicates$sample == "CDK9as_dis2ts_DMSO_18Cr2")
lattice_ScaledmetaReps.proSeq(filename = "spikeNorm/individualReps/Cdk9as_vs_Cdk9asDis2ts_18C_DMSO_scaledMeta_kbsep_SpikeNorm.pdf", 
                          df = scaledMetaData_IndivReplicates[cdk9as_Vs_cdk9asDis2_18C_DMSO_reps,], ylim = c(0, 45), main = "scaled Gene Length", 
                          xlab = "Distance to Gene boundaries (bp)", ylab = "Median PRO-seq intensity", labs = c("cdk9as_dis2-11_18C","ckd9as_18C"))
# Cdk9as (18C) vs Cdk9as_Dis2-11 (30C) (Both DMSO treated)
cdk9as18C_Vs_cdk9asDis2_30C_DMSO_reps = (scaledMetaData_IndivReplicates$sample == "Cdk9as_DMSO_18Cr1" | scaledMetaData_IndivReplicates$sample == "Cdk9as_DMSO_18Cr2" |  scaledMetaData_IndivReplicates$sample == "CDK9as_dis2ts_DMSO_30Cr1" | scaledMetaData_IndivReplicates$sample == "CDK9as_dis2ts_DMSO_30Cr2")
lattice_ScaledmetaReps.proSeq(filename = "spikeNorm/individualReps/Cdk9as18C_vs_Cdk9asDis2ts_30C_DMSO_scaledMeta_kbsep_SpikeNorm.pdf", 
                              df = scaledMetaData_IndivReplicates[cdk9as18C_Vs_cdk9asDis2_30C_DMSO_reps,], ylim = c(0, 45), main = "scaled Gene Length", 
                              xlab = "Distance to Gene boundaries (bp)", ylab = "Median PRO-seq intensity", labs = c("cdk9as_dis2-11_30C","ckd9as_18C"))
# Cdk9as (18C) vs Cdk9as_Dis2-11 (18C) (Both 3MBPP1 treated)
cdk9as18C_Vs_cdk9asDis2_18C_3MBPP1_reps = (scaledMetaData_IndivReplicates$sample == "Cdk9as_3MBPP1_18Cr1" | scaledMetaData_IndivReplicates$sample == "Cdk9as_3MBPP1_18Cr2" |  scaledMetaData_IndivReplicates$sample == "CDK9as_dis2ts_3MBPP1_18Cr1" | scaledMetaData_IndivReplicates$sample == "CDK9as_dis2ts_3MBPP1_18Cr2")
lattice_ScaledmetaReps.proSeq(filename = "spikeNorm/individualReps/Cdk9as18C_vs_Cdk9asDis2ts_18C_3MBPP1_scaledMeta_kbsep_SpikeNorm.pdf", 
                              df = scaledMetaData_IndivReplicates[cdk9as18C_Vs_cdk9asDis2_18C_3MBPP1_reps,], ylim = c(0, 45), main = "scaled Gene Length", 
                              xlab = "Distance to Gene boundaries (bp)", ylab = "Median PRO-seq intensity", labs = c("cdk9as_dis2-11_18C_3MBPP1","ckd9as_18C_3MBPP1"))
# Cdk9as_Dis2-11 (18C) vs Cdk9as_Dis2-11 (30C)
cdk9asDis2_18C_Vs_cdk9asDis2_30C_DMSO_reps = (scaledMetaData_IndivReplicates$sample == "CDK9as_dis2ts_DMSO_18Cr1" | scaledMetaData_IndivReplicates$sample == "CDK9as_dis2ts_DMSO_18Cr2" |  scaledMetaData_IndivReplicates$sample == "CDK9as_dis2ts_DMSO_30Cr1" | scaledMetaData_IndivReplicates$sample == "CDK9as_dis2ts_DMSO_30Cr2")
lattice_ScaledmetaReps.proSeq(filename = "spikeNorm/individualReps/Cdk9asDis2ts_18C_vs_Cdk9asDis2ts_30C_DMSO_scaledMeta_kbsep_SpikeNorm.pdf", 
                              df = scaledMetaData_IndivReplicates[cdk9asDis2_18C_Vs_cdk9asDis2_30C_DMSO_reps,], ylim = c(0, 45), main = "scaled Gene Length", groupFact = "temp",
                              xlab = "Distance to Gene boundaries (bp)", ylab = "Median PRO-seq intensity", labs = c("cdk9as_dis2-11_30C","cdk9as_dis2-11_18C"))
# Cdk9as_Dis2-11 (18C) vs Cdk9as_Dis2-11 (30C)
cdk9asDis2_18C_Vs_cdk9asDis2_30C_3MBPP1_reps = (scaledMetaData_IndivReplicates$sample == "CDK9as_dis2ts_3MBPP1_18Cr1" | scaledMetaData_IndivReplicates$sample == "CDK9as_dis2ts_3MBPP1_18Cr2" |  scaledMetaData_IndivReplicates$sample == "CDK9as_dis2ts_3MBPP1_30Cr1" | scaledMetaData_IndivReplicates$sample == "CDK9as_dis2ts_3MBPP1_30Cr2")
lattice_ScaledmetaReps.proSeq(filename = "spikeNorm/individualReps/Cdk9asDis2ts_18C_vs_Cdk9asDis2ts_30C_3MBPP1_scaledMeta_kbsep_SpikeNorm.pdf", 
                              df = scaledMetaData_IndivReplicates[cdk9asDis2_18C_Vs_cdk9asDis2_30C_3MBPP1_reps,], ylim = c(0, 45), main = "scaled Gene Length", groupFact = "temp",
                              xlab = "Distance to Gene boundaries (bp)", ylab = "Median PRO-seq intensity", labs = c("cdk9as_dis2-11_30C_3MBPP1","cdk9as_dis2-11_18C_3MBPP1"))
# Dis2-11 (18C) vs Dis2-11 (30C) (Both DMSO treated)
Dis2_18C_Vs_Dis2_30C_DMSO_reps = (scaledMetaData_IndivReplicates$sample == "dis2ts_DMSO_18Cr1" | scaledMetaData_IndivReplicates$sample == "dis2ts_DMSO_18Cr2" |  scaledMetaData_IndivReplicates$sample == "dis2ts_DMSO_30Cr1" | scaledMetaData_IndivReplicates$sample == "dis2ts_DMSO_30Cr2")
lattice_ScaledmetaReps.proSeq(filename = "spikeNorm/individualReps/Dis2ts_18C_vs_Dis2ts_30C_DMSO_scaledMeta_kbsep_SpikeNorm.pdf", 
                              df = scaledMetaData_IndivReplicates[Dis2_18C_Vs_Dis2_30C_DMSO_reps,], ylim = c(0, 45), main = "scaled Gene Length", groupFact = "temp",
                              xlab = "Distance to Gene boundaries (bp)", ylab = "Median PRO-seq intensity", labs = c("dis2-11_30C","dis2-11_18C"))
# Dis2-11 (18C) vs Dis2-11 (30C) (Both 3MBPP1 treated)
Dis2_18C_Vs_Dis2_30C_3MBPP1_reps = (scaledMetaData_IndivReplicates$sample == "dis2ts_3MBPP1_18Cr1" | scaledMetaData_IndivReplicates$sample == "dis2ts_3MBPP1_18Cr2" |  scaledMetaData_IndivReplicates$sample == "dis2ts_3MBPP1_30Cr1" | scaledMetaData_IndivReplicates$sample == "dis2ts_3MBPP1_30Cr2")
lattice_ScaledmetaReps.proSeq(filename = "spikeNorm/individualReps/Dis2ts_18C_vs_Dis2ts_30C_3MBPP1_scaledMeta_kbsep_SpikeNorm.pdf", 
                              df = scaledMetaData_IndivReplicates[Dis2_18C_Vs_Dis2_30C_3MBPP1_reps,], ylim = c(0, 45), main = "scaled Gene Length", groupFact = "temp",
                              xlab = "Distance to Gene boundaries (bp)", ylab = "Median PRO-seq intensity", labs = c("dis2-11_30C_3MBPP1","dis2-11_18C_3MBPP1"))
#####################################################################################################################
## plot individual rep data (unscaled) around TSS and CPS. 
TSSMetaData_IndivReplicates_FracMax = FormatMetaData_lattice(wigsetReps, bed, halfWindow = 1000, step = 10, FractMax = TRUE)
CPSMetaData_IndivReplicates_FracMax = FormatMetaData_lattice(wigsetReps, bed[, c(1,3,2,4,5,6)], halfWindow = 1000, step = 10, FractMax = TRUE)
TSSMetaData_IndivReplicates = FormatMetaData_lattice(wigsetReps, bed, halfWindow = 1000, step = 10)
CPSMetaData_IndivReplicates = FormatMetaData_lattice(wigsetReps, bed[, c(1,3,2,4,5,6)], halfWindow = 1000, step = 10)
dir.create(paste(fig_dir, "spikeNorm/individualReps/centeredMetaplots/", sep = ""))#
dir.create(paste(fig_dir, "spikeNorm/individualReps/centeredMetaplots/FractionOfMax/", sep = ""))#", sep = ""))
# Cdk9as (18C) vs Cdk9as_Dis2-11 (18C)
cdk9as_Vs_cdk9asDis2_18C_DMSO_reps_TSS = (TSSMetaData_IndivReplicates$sample == "Cdk9as_DMSO_18Cr1" | TSSMetaData_IndivReplicates$sample == "Cdk9as_DMSO_18Cr2" |  TSSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_DMSO_18Cr1" | TSSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_DMSO_18Cr2")
lattice_metaReps.proSeq(filename = "spikeNorm/individualReps/centeredMetaplots/Cdk9as_vs_Cdk9asDis2ts_18C_DMSO_TSSMeta_kbsep_SpikeNorm.pdf", 
                              df = TSSMetaData_IndivReplicates[cdk9as_Vs_cdk9asDis2_18C_DMSO_reps_TSS,], ylim = c(0, 45), main = "Distance to TSS", 
                              xlab = "Distance to TSS (bp)", ylab = "Median PRO-seq intensity", labs = c("cdk9as_dis2-11_18C","ckd9as_18C"))
cdk9as_Vs_cdk9asDis2_18C_DMSO_reps_CPS = (CPSMetaData_IndivReplicates$sample == "Cdk9as_DMSO_18Cr1" | CPSMetaData_IndivReplicates$sample == "Cdk9as_DMSO_18Cr2" |  CPSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_DMSO_18Cr1" | CPSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_DMSO_18Cr2")
lattice_metaReps.proSeq(filename = "spikeNorm/individualReps/centeredMetaplots/Cdk9as_vs_Cdk9asDis2ts_18C_DMSO_CPSMeta_kbsep_SpikeNorm.pdf", 
                        df = CPSMetaData_IndivReplicates[cdk9as_Vs_cdk9asDis2_18C_DMSO_reps_CPS,], ylim = c(0, 45), main = "Distance to CPS", 
                        xlab = "Distance to CPS (bp)", ylab = "Median PRO-seq intensity", labs = c("cdk9as_dis2-11_18C","ckd9as_18C"))
# Cdk9as (18C) vs Cdk9as_Dis2-11 (30C) (Both DMSO treated)
cdk9as18C_Vs_cdk9asDis2_30C_DMSO_reps_TSS = (TSSMetaData_IndivReplicates$sample == "Cdk9as_DMSO_18Cr1" | TSSMetaData_IndivReplicates$sample == "Cdk9as_DMSO_18Cr2" |  TSSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_DMSO_30Cr1" | TSSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_DMSO_18Cr2")
lattice_metaReps.proSeq(filename = "spikeNorm/individualReps/centeredMetaplots/Cdk9as18C_vs_Cdk9asDis2ts_30C_DMSO_TSSMeta_kbsep_SpikeNorm.pdf", 
                        df = TSSMetaData_IndivReplicates[cdk9as18C_Vs_cdk9asDis2_30C_DMSO_reps_TSS,], ylim = c(0, 45), main = "Distance to TSS", 
                        xlab = "Distance to TSS (bp)", ylab = "Median PRO-seq intensity", labs = c("cdk9as_dis2-11_30C","ckd9as_18C"))
cdk9as18C_Vs_cdk9asDis2_30C_DMSO_reps_CPS = (CPSMetaData_IndivReplicates$sample == "Cdk9as_DMSO_18Cr1" | CPSMetaData_IndivReplicates$sample == "Cdk9as_DMSO_18Cr2" |  CPSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_DMSO_30Cr1" | CPSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_DMSO_30Cr2")
lattice_metaReps.proSeq(filename = "spikeNorm/individualReps/centeredMetaplots/Cdk9as18C_vs_Cdk9asDis2ts_30C_DMSO_CPSMeta_kbsep_SpikeNorm.pdf", 
                        df = CPSMetaData_IndivReplicates[cdk9as18C_Vs_cdk9asDis2_30C_DMSO_reps_CPS,], ylim = c(0, 45), main = "Distance to CPS", 
                        xlab = "Distance to CPS (bp)", ylab = "Median PRO-seq intensity", labs = c("cdk9as_dis2-11_30C","ckd9as_18C"))
# Cdk9as (18C) vs Cdk9as_Dis2-11 (18C) (Both 3MBPP1 treated)
cdk9as18C_Vs_cdk9asDis2_18C_3MBPP1_reps_TSS = (TSSMetaData_IndivReplicates$sample == "Cdk9as_3MBPP1_18Cr1" | TSSMetaData_IndivReplicates$sample == "Cdk9as_3MBPP1_18Cr2" |  TSSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_3MBPP1_18Cr1" | TSSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_3MBPP1_18Cr2")
lattice_metaReps.proSeq(filename = "spikeNorm/individualReps/centeredMetaplots/Cdk9as18C_vs_Cdk9asDis2ts_18C_3MBPP1_TSSMeta_kbsep_SpikeNorm.pdf", 
                        df = TSSMetaData_IndivReplicates[cdk9as18C_Vs_cdk9asDis2_18C_3MBPP1_reps_TSS,], ylim = c(0, 45), main = "Distance to TSS", 
                        xlab = "Distance to TSS (bp)", ylab = "Median PRO-seq intensity", labs = c("cdk9as_dis2-11_18C_3MBPP1","ckd9as_18C_3MBPP1"))
cdk9as18C_Vs_cdk9asDis2_18C_3MBPP1_reps_CPS = (CPSMetaData_IndivReplicates$sample == "Cdk9as_3MBPP1_18Cr1" | CPSMetaData_IndivReplicates$sample == "Cdk9as_3MBPP1_18Cr2" |  CPSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_3MBPP1_18Cr1" | CPSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_3MBPP1_18Cr2")
lattice_metaReps.proSeq(filename = "spikeNorm/individualReps/centeredMetaplots/Cdk9as18C_vs_Cdk9asDis2ts_18C_3MBPP1_CPSMeta_kbsep_SpikeNorm.pdf", 
                        df = CPSMetaData_IndivReplicates[cdk9as18C_Vs_cdk9asDis2_18C_3MBPP1_reps_CPS,], ylim = c(0, 45), main = "Distance to CPS", 
                        xlab = "Distance to CPS (bp)", ylab = "Median PRO-seq intensity", labs = c("cdk9as_dis2-11_18C_3MBPP1","ckd9as_18C_3MBPP1"))
# Cdk9as_Dis2-11 (18C) vs Cdk9as_Dis2-11 (30C) (Both DMSO treated)
cdk9asDis2_18C_Vs_cdk9asDis2_30C_DMSO_reps_TSS = (TSSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_DMSO_18Cr1" | TSSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_DMSO_18Cr2" |  TSSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_DMSO_30Cr1" | TSSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_DMSO_30Cr2")
lattice_metaReps.proSeq(filename = "spikeNorm/individualReps/centeredMetaplots/Cdk9asDis2ts_18C_vs_Cdk9asDis2ts_30C_DMSO_TSSMeta_kbsep_SpikeNorm.pdf", 
                        df = TSSMetaData_IndivReplicates[cdk9asDis2_18C_Vs_cdk9asDis2_30C_DMSO_reps_TSS,], ylim = c(0, 45), main = "Distance to TSS",  groupFact = "temp", 
                        xlab = "Distance to TSS (bp)", ylab = "Median PRO-seq intensity", labs = c("cdk9as_dis2-11_30C_DMSO","cdk9as_dis2-11_18C_DMSO"))
cdk9asDis2_18C_Vs_cdk9asDis2_30C_DMSO_reps_CPS = (CPSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_DMSO_18Cr1" | CPSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_DMSO_18Cr2" |  CPSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_DMSO_30Cr1" | CPSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_DMSO_30Cr2")
lattice_metaReps.proSeq(filename = "spikeNorm/individualReps/centeredMetaplots/Cdk9asDis2ts_18C_vs_Cdk9asDis2ts_30C_DMSO_CPSMeta_kbsep_SpikeNorm.pdf", 
                        df = CPSMetaData_IndivReplicates[cdk9asDis2_18C_Vs_cdk9asDis2_30C_DMSO_reps_CPS,], ylim = c(0, 45), main = "Distance to CPS",  groupFact = "temp", 
                        xlab = "Distance to CPS (bp)", ylab = "Median PRO-seq intensity", labs = c("cdk9as_dis2-11_30C_DMSO","cdk9as_dis2-11_18C_DMSO"))
# Cdk9as_Dis2-11 (18C) vs Cdk9as_Dis2-11 (30C) (Both 3MBPP1 treated)
cdk9asDis2_18C_Vs_cdk9asDis2_30C_3MBPP1_reps_TSS = (TSSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_3MBPP1_18Cr1" | TSSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_3MBPP1_18Cr2" |  TSSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_3MBPP1_30Cr1" | TSSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_3MBPP1_30Cr2")
lattice_metaReps.proSeq(filename = "spikeNorm/individualReps/centeredMetaplots/Cdk9asDis2ts_18C_vs_Cdk9asDis2ts_30C_3MBPP1_TSSMeta_kbsep_SpikeNorm.pdf", 
                        df = TSSMetaData_IndivReplicates[cdk9asDis2_18C_Vs_cdk9asDis2_30C_3MBPP1_reps_TSS,], ylim = c(0, 45), main = "Distance to TSS",  groupFact = "temp", 
                        xlab = "Distance to TSS (bp)", ylab = "Median PRO-seq intensity", labs = c("cdk9as_dis2-11_30C_3MBPP1","cdk9as_dis2-11_18C_3MBPP1"))
cdk9asDis2_18C_Vs_cdk9asDis2_30C_3MBPP1_reps_CPS = (CPSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_3MBPP1_18Cr1" | CPSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_3MBPP1_18Cr2" |  CPSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_3MBPP1_30Cr1" | CPSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_3MBPP1_30Cr2")
lattice_metaReps.proSeq(filename = "spikeNorm/individualReps/centeredMetaplots/Cdk9asDis2ts_18C_vs_Cdk9asDis2ts_30C_3MBPP1_CPSMeta_kbsep_SpikeNorm.pdf", 
                        df = CPSMetaData_IndivReplicates[cdk9asDis2_18C_Vs_cdk9asDis2_30C_3MBPP1_reps_CPS,], ylim = c(0, 45), main = "Distance to CPS",  groupFact = "temp", 
                        xlab = "Distance to CPS (bp)", ylab = "Median PRO-seq intensity", labs = c("cdk9as_dis2-11_30C_3MBPP1","cdk9as_dis2-11_18C_3MBPP1"))
# Cdk9as (18C) vs Cdk9as_Dis2-11 (18C) (Both 3MBPP1 treated)
cdk9asDis2_18C_Vs_cdk9asDis2_30C_3MBPP1_reps_TSS = (TSSMetaData_IndivReplicates$sample == "Cdk9as_3MBPP1_18Cr1" | TSSMetaData_IndivReplicates$sample == "Cdk9as_3MBPP1_18Cr2" |  TSSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_3MBPP1_18Cr1" | TSSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_3MBPP1_18Cr2")
lattice_metaReps.proSeq(filename = "spikeNorm/individualReps/centeredMetaplots/Cdk9as_18C_vs_Cdk9asDis2ts_18C_3MBPP1_TSSMeta_kbsep_SpikeNorm.pdf", 
                        df = TSSMetaData_IndivReplicates[cdk9asDis2_18C_Vs_cdk9asDis2_30C_3MBPP1_reps_TSS,], ylim = c(0, 45), main = "Distance to TSS", 
                        xlab = "Distance to TSS (bp)", ylab = "Median PRO-seq intensity", labs = c("cdk9as_dis2-11_18C_3MBPP1", "cdk9as_18C_3MBPP1"))
cdk9asDis2_18C_Vs_cdk9asDis2_30C_3MBPP1_reps_CPS = (CPSMetaData_IndivReplicates$sample == "Cdk9as_3MBPP1_18Cr1" | CPSMetaData_IndivReplicates$sample == "Cdk9as_3MBPP1_18Cr2" |  CPSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_3MBPP1_18Cr1" | CPSMetaData_IndivReplicates$sample == "CDK9as_dis2ts_3MBPP1_18Cr2")
lattice_metaReps.proSeq(filename = "spikeNorm/individualReps/centeredMetaplots/Cdk9as_18C_vs_Cdk9asDis2ts_18C_3MBPP1_CPSMeta_kbsep_SpikeNorm.pdf", 
                        df = CPSMetaData_IndivReplicates[cdk9asDis2_18C_Vs_cdk9asDis2_30C_3MBPP1_reps_CPS,], ylim = c(0, 45), main = "Distance to CPS", 
                        xlab = "Distance to CPS (bp)", ylab = "Median PRO-seq intensity", labs = c("cdk9as_dis2-11_18C_3MBPP1", "cdk9as_18C_3MBPP1"))
#####################################################################################################################
# Plots of Fraction of Maximum:
#####################################################################################################################
# Cdk9as (18C) vs Cdk9as_Dis2-11 (18C)
cdk9as_Vs_cdk9asDis2_18C_DMSO_reps_CPS_FractMax = (CPSMetaData_IndivReplicates_FracMax$sample == "Cdk9as_DMSO_18Cr1" | CPSMetaData_IndivReplicates_FracMax$sample == "Cdk9as_DMSO_18Cr2" |  CPSMetaData_IndivReplicates_FracMax$sample == "CDK9as_dis2ts_DMSO_18Cr1" | CPSMetaData_IndivReplicates_FracMax$sample == "CDK9as_dis2ts_DMSO_18Cr2")
lattice_metaReps.proSeq(filename = "spikeNorm/individualReps/centeredMetaplots/FractionOfMax/Cdk9as_vs_Cdk9asDis2ts_18C_DMSO_CPSMeta_kbsep_FractionOfMaximum.pdf", 
                        df = CPSMetaData_IndivReplicates_FracMax[cdk9as_Vs_cdk9asDis2_18C_DMSO_reps_CPS_FractMax,], ylim = c(0, 1.2), main = "Distance to CPS", 
                        xlab = "Distance to CPS (bp)", ylab = "Fraction of Maximum", labs = c("cdk9as_dis2-11_18C","ckd9as_18C"))
# WT (30C) vs Dis2-11 (18C)
WT_30CvsDis2_18C_DMSO_RPM = (CPSMetaData_IndivReplicates_FracMax$sample == "WT_DMSOr1" |  CPSMetaData_IndivReplicates_FracMax$sample == "WT_DMSOr2" | CPSMetaData_IndivReplicates_FracMax$sample == "dis2ts_DMSO_18Cr1" | CPSMetaData_IndivReplicates_FracMax$sample == "dis2ts_DMSO_18Cr2")
lattice_metaReps.proSeq(filename = "spikeNorm/individualReps/centeredMetaplots/FractionOfMax/WT_30CvsDis2_18C_DMSO_CPSMeta_kbsep_FractionMax.pdf", 
                          df = CPSMetaData_IndivReplicates_FracMax[WT_30CvsDis2_18C_DMSO_RPM,], ylim = c(0, 1.2), main = "Distance to CPS", 
                          xlab = "Distance to CPS (bp)", ylab = "Fraction of Maximum", labs = c("WT_30C", "Dis2-11_18C"))
# WT (30C) vs Dis2-11 (30C)
WT_30CvsDis2_30C_DMSO_RPM = (CPSMetaData_IndivReplicates_FracMax$sample == "WT_DMSOr1" |  CPSMetaData_IndivReplicates_FracMax$sample == "WT_DMSOr2" | CPSMetaData_IndivReplicates_FracMax$sample == "dis2ts_DMSO_30Cr1" | CPSMetaData_IndivReplicates_FracMax$sample == "dis2ts_DMSO_30Cr2")
lattice_metaReps.proSeq(filename = "spikeNorm/individualReps/centeredMetaplots/FractionOfMax/WT_30CvsDis2_30C_DMSO_CPSMeta_kbsep_FractionMax.pdf", 
                        df = CPSMetaData_IndivReplicates_FracMax[WT_30CvsDis2_30C_DMSO_RPM,], ylim = c(0, 1.2), main = "Distance to CPS", 
                        xlab = "Distance to CPS (bp)", ylab = "Fraction of Maximum", labs = c("WT_30C", "Dis2-11_30C"))


#####################################################################################################################
#####################################################################################################################
# RPM Normalized Plots
#####################################################################################################################
#####################################################################################################################
scaledMetaData.RPM = FormatScaledMetaData_lattice(wigsetRPM, bed)
# generate same plots, but normalizing using RPM instead of Spike-Ins
#  Cdk9as (18C) vs Cdk9as_Dis2-11 (18C) both treated with DMSO
dir.create(paste(fig_dir, "RPM_Norm/", sep = ""))
cdk9as_Vs_cdk9asDis2_18C_DMSO_RPM = (scaledMetaData.RPM$sample == "Cdk9as_DMSO_18C_combined" |  scaledMetaData.RPM$sample == "CDK9as_dis2ts_DMSO_18C_combined")
lattice_Scaledmeta.proSeq(filename = "RPM_Norm/Cdk9as_vs_Cdk9asDis2ts_18C_DMSO_scaledMeta_kbsep_RPMnorm.pdf", 
                          df = scaledMetaData.RPM[cdk9as_Vs_cdk9asDis2_18C_DMSO_RPM,], ylim = c(0, 2), main = "scaled Gene Length", 
                          xlab = "Distance to Gene boundaries (bp)", ylab = "Median PRO-seq intensity", labs = c("cdk9as_dis2-11_18C", "ckd9as_18C"))
# WT (30C) vs Dis2-11 (18C)
WT_30CvsDis2_18C_DMSO_RPM = (scaledMetaData.RPM$sample == "WT_DMSO_combined" |  scaledMetaData.RPM$sample == "dis2ts_DMSO_18C_combined")
lattice_Scaledmeta.proSeq(filename = "RPM_Norm/WT_30CvsDis2_18C_DMSO_scaledMeta_kbsep_RPMNorm.pdf", 
                          df = scaledMetaData.RPM[WT_30CvsDis2_18C_DMSO_RPM,], ylim = c(0, 2), main = "scaled Gene Length", 
                          xlab = "Distance to Gene boundaries (bp)", ylab = "Median PRO-seq intensity", labs = c("WT_30C", "Dis2-11_18C"))

#  Cdk9as (18C) vs Cdk9as_Dis2-11 (18C) both treated with 3MBPP1
cdk9as_Vs_cdk9asDis2_18C_3MBPP1_RPM = (scaledMetaData.RPM$sample == "Cdk9as_3MBPP1_18C_combined" |  scaledMetaData.RPM$sample == "CDK9as_dis2ts_3MBPP1_18C_combined")
lattice_Scaledmeta.proSeq(filename = "RPM_Norm/Cdk9as_vs_Cdk9asDis2ts_18C_3MBPP1_scaledMeta_kbsep_RPMnorm.pdf", 
                          df = scaledMetaData.RPM[cdk9as_Vs_cdk9asDis2_18C_3MBPP1_RPM,], ylim = c(0, 2), main = "scaled Gene Length", 
                          xlab = "Distance to Gene boundaries (bp)", ylab = "Median PRO-seq intensity", labs = c("cdk9as_dis2-11_18C_3MBPP1", "ckd9as_18C_3MBPP1"))
#  Cdk9as (18C) vs Dis2-11 (18C) both DMSO
cdk9as_Vs_Dis2_18C_DMSO_RPM = (scaledMetaData.RPM$sample == "Cdk9as_DMSO_18C_combined" |  scaledMetaData.RPM$sample == "dis2ts_DMSO_18C_combined")
lattice_Scaledmeta.proSeq(filename = "RPM_Norm/Cdk9as_vs_Dis2ts_18C_DMSO_scaledMeta_kbsep_RPMnorm.pdf", 
                          df = scaledMetaData.RPM[cdk9as_Vs_Dis2_18C_DMSO_RPM,], ylim = c(0, 2), main = "scaled Gene Length", 
                          xlab = "Distance to Gene boundaries (bp)", ylab = "Median PRO-seq intensity", labs = c("dis2ts_DMSO_18C_combined", "Cdk9as_DMSO_18C_combined"))
#####################################################################################################################
scaledMetaData_IndivReplicates_RPM = FormatScaledMetaData_lattice(wigsetReps_RPM, bed)
dir.create(paste(fig_dir, "RPM_Norm/individualReps/", sep = ""))
# WT (DMSO, 30 C) vs Cdk9as (DMSO 18 C) 
WT_30C_Vs_cdk9as_18C_DMSO_reps_RPM = (scaledMetaData_IndivReplicates_RPM$sample == "WT_DMSOr1" | scaledMetaData_IndivReplicates_RPM$sample == "WT_DMSOr2" |  scaledMetaData_IndivReplicates_RPM$sample == "Cdk9as_DMSO_18Cr1" | scaledMetaData_IndivReplicates_RPM$sample == "Cdk9as_DMSO_18Cr2")
lattice_ScaledmetaReps.proSeq(filename = "RPM_Norm/individualReps/WT_30C_vs_Cdk9as_18C_DMSO_scaledMeta_kbsep_RPM.pdf", 
                              df = scaledMetaData_IndivReplicates_RPM[WT_30C_Vs_cdk9as_18C_DMSO_reps_RPM,], ylim = c(0, 2), main = "scaled Gene Length", 
                              xlab = "Distance to Gene boundaries (bp)", ylab = "Median PRO-seq intensity", labs = c("WT_30C","ckd9as_18C"))
# WT (DMSO, 30 C) vs Dis2ts (DMSO 18 C) 
WT_30C_Vs_dis218C_DMSO_reps_RPM = (scaledMetaData_IndivReplicates_RPM$sample == "WT_DMSOr1" | scaledMetaData_IndivReplicates_RPM$sample == "WT_DMSOr2" |  scaledMetaData_IndivReplicates_RPM$sample == "dis2ts_DMSO_18Cr1" | scaledMetaData_IndivReplicates_RPM$sample == "dis2ts_DMSO_18Cr2")
lattice_ScaledmetaReps.proSeq(filename = "RPM_Norm/individualReps/WT_30C_vs_dis2ts_18C_DMSO_scaledMeta_kbsep_RPM.pdf", 
                              df = scaledMetaData_IndivReplicates_RPM[WT_30C_Vs_dis218C_DMSO_reps_RPM,], ylim = c(0, 2), main = "scaled Gene Length", 
                              xlab = "Distance to Gene boundaries (bp)", ylab = "Median PRO-seq intensity", labs = c("WT_30C","dis2ts_18C"))
# WT (DMSO, 30 C) vs Dis2ts (DMSO 30 C) 
WT_30C_Vs_dis2_30C_DMSO_reps_RPM = (scaledMetaData_IndivReplicates_RPM$sample == "WT_DMSOr1" | scaledMetaData_IndivReplicates_RPM$sample == "WT_DMSOr2" |  scaledMetaData_IndivReplicates_RPM$sample == "dis2ts_DMSO_30Cr1" | scaledMetaData_IndivReplicates_RPM$sample == "dis2ts_DMSO_30Cr2")
lattice_ScaledmetaReps.proSeq(filename = "RPM_Norm/individualReps/WT_30C_vs_dis2ts_30C_DMSO_scaledMeta_kbsep_RPM.pdf", 
                              df = scaledMetaData_IndivReplicates_RPM[WT_30C_Vs_dis2_30C_DMSO_reps_RPM,], ylim = c(0, 2), main = "scaled Gene Length", 
                              xlab = "Distance to Gene boundaries (bp)", ylab = "Median PRO-seq intensity", labs = c("WT_30C","dis2ts_30C"))

# Cdk9as (18C) vs Cdk9as_Dis2-11 (18C) (Both DMSO treated)
cdk9as_Vs_cdk9asDis2_18C_DMSO_reps_RPM = (scaledMetaData_IndivReplicates_RPM$sample == "Cdk9as_DMSO_18Cr1" | scaledMetaData_IndivReplicates_RPM$sample == "Cdk9as_DMSO_18Cr2" |  scaledMetaData_IndivReplicates_RPM$sample == "CDK9as_dis2ts_DMSO_18Cr1" | scaledMetaData_IndivReplicates_RPM$sample == "CDK9as_dis2ts_DMSO_18Cr2")
lattice_ScaledmetaReps.proSeq(filename = "RPM_Norm/individualReps/Cdk9as_vs_Cdk9asDis2ts_18C_DMSO_scaledMeta_kbsep_RPM.pdf", 
                              df = scaledMetaData_IndivReplicates_RPM[cdk9as_Vs_cdk9asDis2_18C_DMSO_reps_RPM,], ylim = c(0, 2), main = "scaled Gene Length", 
                              xlab = "Distance to Gene boundaries (bp)", ylab = "Median PRO-seq intensity", labs = c("cdk9as_dis2-11_18C","ckd9as_18C"))
# Cdk9as (18C) vs Cdk9as_Dis2-11 (30C) (Both DMSO treated)
cdk9as18C_Vs_cdk9asDis2_30C_DMSO_reps_RPM = (scaledMetaData_IndivReplicates_RPM$sample == "Cdk9as_DMSO_18Cr1" | scaledMetaData_IndivReplicates_RPM$sample == "Cdk9as_DMSO_18Cr2" |  scaledMetaData_IndivReplicates_RPM$sample == "CDK9as_dis2ts_DMSO_30Cr1" | scaledMetaData_IndivReplicates_RPM$sample == "CDK9as_dis2ts_DMSO_30Cr2")
lattice_ScaledmetaReps.proSeq(filename = "RPM_Norm/individualReps/Cdk9as18C_vs_Cdk9asDis2ts_30C_DMSO_scaledMeta_kbsep_RPM.pdf", 
                              df = scaledMetaData_IndivReplicates_RPM[cdk9as18C_Vs_cdk9asDis2_30C_DMSO_reps_RPM,], ylim = c(0, 2), main = "scaled Gene Length", 
                              xlab = "Distance to Gene boundaries (bp)", ylab = "Median PRO-seq intensity", labs = c("cdk9as_dis2-11_30C","ckd9as_18C"))
# Cdk9as (18C) vs Cdk9as_Dis2-11 (18C) (Both 3MBPP1 treated)
cdk9as18C_Vs_cdk9asDis2_18C_3MBPP1_reps_RPM = (scaledMetaData_IndivReplicates_RPM$sample == "Cdk9as_3MBPP1_18Cr1" | scaledMetaData_IndivReplicates_RPM$sample == "Cdk9as_3MBPP1_18Cr2" |  scaledMetaData_IndivReplicates_RPM$sample == "CDK9as_dis2ts_3MBPP1_18Cr1" | scaledMetaData_IndivReplicates_RPM$sample == "CDK9as_dis2ts_3MBPP1_18Cr2")
lattice_ScaledmetaReps.proSeq(filename = "RPM_Norm/individualReps/Cdk9as18C_vs_Cdk9asDis2ts_18C_3MBPP1_scaledMeta_kbsep_RPM.pdf", 
                              df = scaledMetaData_IndivReplicates_RPM[cdk9as18C_Vs_cdk9asDis2_18C_3MBPP1_reps_RPM,], ylim = c(0, 2), main = "scaled Gene Length", 
                              xlab = "Distance to Gene boundaries (bp)", ylab = "Median PRO-seq intensity", labs = c("cdk9as_dis2-11_18C_3MBPP1","ckd9as_18C_3MBPP1"))
# Cdk9as (18C) vs Cdk9as_Dis2-11 (30C) (Both 3MBPP1 treated)
cdk9as18C_Vs_cdk9asDis2_30C_3MBPP1_reps_RPM = (scaledMetaData_IndivReplicates_RPM$sample == "Cdk9as_3MBPP1_18Cr1" | scaledMetaData_IndivReplicates_RPM$sample == "Cdk9as_3MBPP1_18Cr2" |  scaledMetaData_IndivReplicates_RPM$sample == "CDK9as_dis2ts_3MBPP1_30Cr1" | scaledMetaData_IndivReplicates_RPM$sample == "CDK9as_dis2ts_3MBPP1_30Cr2")
lattice_ScaledmetaReps.proSeq(filename = "RPM_Norm/individualReps/Cdk9as18C_vs_Cdk9asDis2ts_30C_3MBPP1_scaledMeta_kbsep_RPM.pdf", 
                              df = scaledMetaData_IndivReplicates_RPM[cdk9as18C_Vs_cdk9asDis2_30C_3MBPP1_reps_RPM,], ylim = c(0, 2), main = "scaled Gene Length", 
                              xlab = "Distance to Gene boundaries (bp)", ylab = "Median PRO-seq intensity", labs = c("cdk9as_dis2-11_30C_3MBPP1","ckd9as_18C_3MBPP1"))

# Cdk9as_Dis2-11 (18C) vs Cdk9as_Dis2-11 (30C)
cdk9asDis2_18C_Vs_cdk9asDis2_30C_DMSO_reps_RPM = (scaledMetaData_IndivReplicates_RPM$sample == "CDK9as_dis2ts_DMSO_18Cr1" | scaledMetaData_IndivReplicates_RPM$sample == "CDK9as_dis2ts_DMSO_18Cr2" |  scaledMetaData_IndivReplicates_RPM$sample == "CDK9as_dis2ts_DMSO_30Cr1" | scaledMetaData_IndivReplicates_RPM$sample == "CDK9as_dis2ts_DMSO_30Cr2")
lattice_ScaledmetaReps.proSeq(filename = "RPM_Norm/individualReps/Cdk9asDis2ts_18C_vs_Cdk9asDis2ts_30C_DMSO_scaledMeta_kbsep_RPM.pdf", 
                              df = scaledMetaData_IndivReplicates_RPM[cdk9asDis2_18C_Vs_cdk9asDis2_30C_DMSO_reps_RPM,], ylim = c(0, 2), main = "scaled Gene Length", groupFact = "temp",
                              xlab = "Distance to Gene boundaries (bp)", ylab = "Median PRO-seq intensity", labs = c("cdk9as_dis2-11_30C","cdk9as_dis2-11_18C"))
# Cdk9as_Dis2-11 (18C) vs Cdk9as_Dis2-11 (30C)
cdk9asDis2_18C_Vs_cdk9asDis2_30C_3MBPP1_reps_RPM = (scaledMetaData_IndivReplicates_RPM$sample == "CDK9as_dis2ts_3MBPP1_18Cr1" | scaledMetaData_IndivReplicates_RPM$sample == "CDK9as_dis2ts_3MBPP1_18Cr2" |  scaledMetaData_IndivReplicates_RPM$sample == "CDK9as_dis2ts_3MBPP1_30Cr1" | scaledMetaData_IndivReplicates_RPM$sample == "CDK9as_dis2ts_3MBPP1_30Cr2")
lattice_ScaledmetaReps.proSeq(filename = "RPM_Norm/individualReps/Cdk9asDis2ts_18C_vs_Cdk9asDis2ts_30C_3MBPP1_scaledMeta_kbsep_RPM.pdf", 
                              df = scaledMetaData_IndivReplicates_RPM[cdk9asDis2_18C_Vs_cdk9asDis2_30C_3MBPP1_reps_RPM,], ylim = c(0, 2), main = "scaled Gene Length", groupFact = "temp",
                              xlab = "Distance to Gene boundaries (bp)", ylab = "Median PRO-seq intensity", labs = c("cdk9as_dis2-11_30C_3MBPP1","cdk9as_dis2-11_18C_3MBPP1"))
# Dis2-11 (18C) vs Dis2-11 (30C) (Both DMSO treated)
Dis2_18C_Vs_Dis2_30C_DMSO_reps_RPM = (scaledMetaData_IndivReplicates_RPM$sample == "dis2ts_DMSO_18Cr1" | scaledMetaData_IndivReplicates_RPM$sample == "dis2ts_DMSO_18Cr2" |  scaledMetaData_IndivReplicates_RPM$sample == "dis2ts_DMSO_30Cr1" | scaledMetaData_IndivReplicates_RPM$sample == "dis2ts_DMSO_30Cr2")
lattice_ScaledmetaReps.proSeq(filename = "RPM_Norm/individualReps/Dis2ts_18C_vs_Dis2ts_30C_DMSO_scaledMeta_kbsep_RPM.pdf", 
                              df = scaledMetaData_IndivReplicates_RPM[Dis2_18C_Vs_Dis2_30C_DMSO_reps_RPM,], ylim = c(0, 2), main = "scaled Gene Length", groupFact = "temp",
                              xlab = "Distance to Gene boundaries (bp)", ylab = "Median PRO-seq intensity", labs = c("dis2-11_30C","dis2-11_18C"))
# Dis2-11 (18C) vs Dis2-11 (30C) (Both 3MBPP1 treated)
Dis2_18C_Vs_Dis2_30C_3MBPP1_reps_RPM = (scaledMetaData_IndivReplicates_RPM$sample == "dis2ts_3MBPP1_18Cr1" | scaledMetaData_IndivReplicates_RPM$sample == "dis2ts_3MBPP1_18Cr2" |  scaledMetaData_IndivReplicates_RPM$sample == "dis2ts_3MBPP1_30Cr1" | scaledMetaData_IndivReplicates_RPM$sample == "dis2ts_3MBPP1_30Cr2")
lattice_ScaledmetaReps.proSeq(filename = "RPM_Norm/individualReps/Dis2ts_18C_vs_Dis2ts_30C_3MBPP1_scaledMeta_kbsep_RPM.pdf", 
                              df = scaledMetaData_IndivReplicates_RPM[Dis2_18C_Vs_Dis2_30C_3MBPP1_reps_RPM,], ylim = c(0, 2), main = "scaled Gene Length", groupFact = "temp",
                              xlab = "Distance to Gene boundaries (bp)", ylab = "Median PRO-seq intensity", labs = c("dis2-11_30C_3MBPP1","dis2-11_18C_3MBPP1"))


