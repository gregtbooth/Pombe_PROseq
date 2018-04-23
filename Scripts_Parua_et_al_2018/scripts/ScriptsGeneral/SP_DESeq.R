
library(DESeq2)
bwpath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/alignment/bw/"
fig_dir = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/figures/DESeq2/12-14-17/"
dir.create(fig_dir)
infopath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/pombeinfo/"
countpath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/countData/"
DEseq_Table_path = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/countData/AllGenesObsTSS/DESeqOutput/"
filteredGL = read.table(paste(infopath, "SP_CompleteFilteredGenes.bed", sep = ""))
#########################################################################################################
## Note: Rather than running DESeq here, I'm now just loading the DESeq2 tables that were generated 
## from the various comparisons for all genes (with spike in normalization) and then using this script 
## to make MA plots. 
## the following script: /Volumes/SEAGATE_EXP/Fisher_collaboration/scripts/SP_write_PItables_DEseqTables_Allgenes.R
## was used to make the tables used here.
WT_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_WT_gb_res.txt", sep = ""), head = T)
WT_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_WT_pr_res.txt", sep = ""), head = T)
WT_term_res= read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_WT_term_res.txt", sep = ""), head = T)
mcs6_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_mcs6_gb_res.txt", sep = ""), head = T)
mcs6_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_mcs6_pr_res.txt", sep = ""), head = T)
mcs6_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_mcs6_term_res.txt", sep = ""), head = T)
CDK9_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_gb_res.txt", sep = ""), head = T)
CDK9_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_pr_res.txt", sep = ""), head = T)
CDK9_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_term_res.txt", sep = ""), head = T)
Isk1_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Isk1_gb_res.txt", sep = ""), head = T)
Isk1_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Isk1_pr_res.txt", sep = ""), head = T)
Isk1_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Isk1_term_res.txt", sep = ""), head = T)
mcs6_CDK9_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_mcs6_CDK9_gb_res.txt", sep = ""), head = T)
mcs6_CDK9_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_mcs6_CDK9_pr_res.txt", sep = ""), head = T)
mcs6_CDK9_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_mcs6_CDK9_term_res.txt", sep = ""), head = T)
Cdk9as_18C_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Cdk9as_18C_gb_res.txt", sep = ""), head = T)
Cdk9as_18C_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Cdk9as_18C_pr_res.txt", sep = ""), head = T)
Cdk9as_18C_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Cdk9as_18C_term_res.txt", sep = ""), head = T)
dis2ts_18C_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_dis2ts_18C_gb_res.txt", sep = ""), head = T)
dis2ts_18C_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_dis2ts_18C_pr_res.txt", sep = ""), head = T)
dis2ts_18C_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_dis2ts_18C_term_res.txt", sep = ""), head = T)
CDK9_dis2ts_18C_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_dis2ts_18C_gb_res.txt", sep = ""), head = T)
CDK9_dis2ts_18C_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_dis2ts_18C_pr_res.txt", sep = ""), head = T)
CDK9_dis2ts_18C_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_dis2ts_18C_term_res.txt", sep = ""), head = T)
dis2ts_30C_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_dis2ts_30C_gb_res.txt", sep = ""), head = T)
dis2ts_30C_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_dis2ts_30C_pr_res.txt", sep = ""), head = T)
dis2ts_30C_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_dis2ts_30C_term_res.txt", sep = ""), head = T)
CDK9_dis2ts_30C_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_dis2ts_30C_gb_res.txt", sep = ""), head = T)
CDK9_dis2ts_30C_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_dis2ts_30C_pr_res.txt", sep = ""), head = T)
CDK9_dis2ts_30C_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_dis2ts_30C_term_res.txt", sep = ""), head = T)
dis2ts_30Cv18C_DMSO_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_dis2ts_30Cv18C_DMSO_gb_res.txt", sep = ""), head = T)
dis2ts_30Cv18C_DMSO_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_dis2ts_30Cv18C_DMSO_pr_res.txt", sep = ""), head = T)
dis2ts_30Cv18C_DMSO_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_dis2ts_30Cv18C_DMSO_term_res.txt", sep = ""), head = T)
dis2ts_30Cv18C_3MBPP1_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_dis2ts_30Cv18C_3MBPP1_gb_res.txt", sep = ""), head = T)
dis2ts_30Cv18C_3MBPP1_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_dis2ts_30Cv18C_3MBPP1_pr_res.txt", sep = ""), head = T)
dis2ts_30Cv18C_3MBPP1_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_dis2ts_30Cv18C_3MBPP1_term_res.txt", sep = ""), head = T)
CDK9_dis2ts_30Cv18C_DMSO_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_dis2ts_30Cv18C_DMSO_gb_res.txt", sep = ""), head = T)
CDK9_dis2ts_30Cv18C_DMSO_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_dis2ts_30Cv18C_DMSO_pr_res.txt", sep = ""), head = T)
CDK9_dis2ts_30Cv18C_DMSO_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_dis2ts_30Cv18C_DMSO_term_res.txt", sep = ""), head = T)
CDK9_dis2ts_30Cv18C_3MBPP1_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_dis2ts_30Cv18C_3MBPP1_gb_res.txt", sep = ""), head = T)
CDK9_dis2ts_30Cv18C_3MBPP1_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_dis2ts_30Cv18C_3MBPP1_pr_res.txt", sep = ""), head = T)
CDK9_dis2ts_30Cv18C_3MBPP1_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_dis2ts_30Cv18C_3MBPP1_term_res.txt", sep = ""), head = T)

WTvCDK9_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_WTvsCDK9as_DMSO_gb_res.txt", sep = ""), head = T)
WTvCDK9_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_WTvsCDK9as_DMSO_pr_res.txt", sep = ""), head = T)
WTvCDK9_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_WTvsCDK9as_DMSO_term_res.txt", sep = ""), head = T)
WTvMcs6_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_WTvsMcs6as_DMSO_gb_res.txt", sep = ""), head = T)
WTvMcs6_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_WTvsMcs6as_DMSO_pr_res.txt", sep = ""), head = T)
WTvMcs6_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_WTvsMcs6as_DMSO_term_res.txt", sep = ""), head = T)
WTvIsk1_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_WTvsIsk1_DMSO_gb_res.txt", sep = ""), head = T)
WTvIsk1_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_WTvsIsk1_DMSO_pr_res.txt", sep = ""), head = T)
WTvIsk1_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_WTvIsk1_DMSO_term_res.txt", sep = ""), head = T)
WTvsMcs6_CDK9_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_WTvsMcs6_CDK9_DMSO_gb_res.txt", sep = ""), head = T)
WTvsMcs6_CDK9_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_WTvsMcs6_CDK9_DMSO_pr_res.txt", sep = ""), head = T)
WTvsMcs6_CDK9_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_WTvsMcs6_CDK9_DMSO_term_res.txt", sep = ""), head = T)

#30'' vs 0' 
CDK9_0.5vs0min_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_30secvs0min_gb_res.txt", sep = ""), head = T)
CDK9_0.5vs0min_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_30secvs0min_pr_res.txt", sep = ""), head = T)
CDK9_0.5vs0min_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_30secvs0min_term_res.txt", sep = ""), head = T)
#1' vs 0' 
CDK9_1vs0min_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_1vs0min_gb_res.txt", sep = ""), head = T)
CDK9_1vs0min_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_1vs0min_pr_res.txt", sep = ""), head = T)
CDK9_1vs0min_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_1vs0min_term_res.txt", sep = ""), head = T)
#2.5' vs 0' 
CDK9_2.5vs0min_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_2min30secVs0min_gb_res.txt", sep = ""), head = T)
CDK9_2.5vs0min_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_2min30secVs0min_pr_res.txt", sep = ""), head = T)
CDK9_2.5vs0min_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_2min30secVs0min_term_res.txt", sep = ""), head = T)
#5' vs 0' 
CDK9_5vs0min_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_5vs0min_gb_res.txt", sep = ""), head = T)
CDK9_5vs0min_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_5vs0min_pr_res.txt", sep = ""), head = T)
CDK9_5vs0min_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_5vs0min_term_res.txt", sep = ""), head = T)
#7.5' vs 0' 
CDK9_7.5vs0min_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_7min30secVs0min_gb_res.txt", sep = ""), head = T)
CDK9_7.5vs0min_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_7min30secVs0min_pr_res.txt", sep = ""), head = T)
CDK9_7.5vs0min_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_7min30secVs0min_term_res.txt", sep = ""), head = T)
#20 vs 0' 
CDK9_10vs0min_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_10vs0min_gb_res.txt", sep = ""), head = T)
CDK9_10vs0min_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_10vs0min_pr_res.txt", sep = ""), head = T)
CDK9_10vs0min_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_10vs0min_term_res.txt", sep = ""), head = T)
#20' vs 0' 
CDK9_20vs0min_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_20vs0min_gb_res.txt", sep = ""), head = T)
CDK9_20vs0min_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_20vs0min_pr_res.txt", sep = ""), head = T)
CDK9_20vs0min_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_CDK9_20vs0min_term_res.txt", sep = ""), head = T)
# Lsk1as Cdk9as double mutant
Lsk1_CDK9_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Lsk1_CDK9_gb_res.txt", sep = ""), head = T)
Lsk1_CDK9_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Lsk1_CDK9_pr_res.txt", sep = ""), head = T)
Lsk1_CDK9_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Lsk1_CDK9_term_res.txt", sep = ""), head = T)
# Cdk9as Time Course in Spt5 CTR mutants 
# spt5WT7 1 minutes vs 0 minutes
Spt5WT7_1minVs0min_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Spt5WT7_1minVs0min_gb_res.txt", sep = ""), head = T)
Spt5WT7_1minVs0min_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Spt5WT7_1minVs0min_pr_res.txt", sep = ""), head = T)
Spt5WT7_1minVs0min_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Spt5WT7_1minVs0min_term_res.txt", sep = ""), head = T)
# spt5WT7 5 minutes vs 0 minutes
Spt5WT7_5minVs0min_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Spt5WT7_5minVs0min_gb_res.txt", sep = ""), head = T)
Spt5WT7_5minVs0min_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Spt5WT7_5minVs0min_pr_res.txt", sep = ""), head = T)
Spt5WT7_5minVs0min_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Spt5WT7_5minVs0min_term_res.txt", sep = ""), head = T)
# Spt5T1A 1 minutes vs 0 minutes
Spt5T1A_1minVs0min_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Spt5T1A_1minVs0min_gb_res.txt", sep = ""), head = T)
Spt5T1A_1minVs0min_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Spt5T1A_1minVs0min_pr_res.txt", sep = ""), head = T)
Spt5T1A_1minVs0min_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Spt5T1A_1minVs0min_term_res.txt", sep = ""), head = T)
# Spt5T1A 5 minutes vs 0 minutes
Spt5T1A_5minVs0min_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Spt5T1A_5minVs0min_gb_res.txt", sep = ""), head = T)
Spt5T1A_5minVs0min_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Spt5T1A_5minVs0min_pr_res.txt", sep = ""), head = T)
Spt5T1A_5minVs0min_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Spt5T1A_5minVs0min_term_res.txt", sep = ""), head = T)
# Spt5T1E 1 minutes vs 0 minutes
Spt5T1E_1minVs0min_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Spt5T1E_1minVs0min_gb_res.txt", sep = ""), head = T)
Spt5T1E_1minVs0min_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Spt5T1E_1minVs0min_pr_res.txt", sep = ""), head = T)
Spt5T1E_1minVs0min_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Spt5T1E_1minVs0min_term_res.txt", sep = ""), head = T)
# Spt5T1E 5 minutes vs 0 minutes
Spt5T1E_5minVs0min_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Spt5T1E_5minVs0min_gb_res.txt", sep = ""), head = T)
Spt5T1E_5minVs0min_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Spt5T1E_5minVs0min_pr_res.txt", sep = ""), head = T)
Spt5T1E_5minVs0min_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Spt5T1E_5minVs0min_term_res.txt", sep = ""), head = T)
# Spt5WT7 vs spt5T1A 0 minutes
Spt5WT7_Spt5T1A_0min_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Spt5WT7_Spt5T1A_0min_gb_res.txt", sep = ""), head = T)
Spt5WT7_Spt5T1A_0min_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Spt5WT7_Spt5T1A_0min_pr_res.txt", sep = ""), head = T)
Spt5WT7_Spt5T1A_0min_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Spt5WT7_Spt5T1A_0min_term_res.txt", sep = ""), head = T)
# Spt5WT7 vs spt5T1E 0 minutes
Spt5WT7_Spt5T1E_0min_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Spt5WT7_Spt5T1E_0min_gb_res.txt", sep = ""), head = T)
Spt5WT7_Spt5T1E_0min_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Spt5WT7_Spt5T1E_0min_pr_res.txt", sep = ""), head = T)
Spt5WT7_Spt5T1E_0min_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_Spt5WT7_Spt5T1E_0min_term_res.txt", sep = ""), head = T)
# Comparisons between different Dis2 mutants:
# WT7 vs Dis2-11 
FisherWT_dis2_11_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_FisherWT_dis2_11_gb_res.txt", sep = ""), head = T)
FisherWT_dis2_11_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_FisherWT_dis2_11_pr_res.txt", sep = ""), head = T)
FisherWT_dis2_11_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_FisherWT_dis2_11_term_res.txt", sep = ""), head = T)
# WT7 vs Dis2Delete 
FisherWT_dis2Del_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_FisherWT_dis2Del_gb_res.txt", sep = ""), head = T)
FisherWT_dis2Del_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_FisherWT_dis2Del_pr_res.txt", sep = ""), head = T)
FisherWT_dis2Del_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_FisherWT_dis2Del_term_res.txt", sep = ""), head = T)
# WT7 vs Dis2-T316A
FisherWT_dis2_T316A_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_dis2_T316A_gb_res.txt", sep = ""), head = T)
FisherWT_dis2_T316A_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_dis2_T316A_pr_res.txt", sep = ""), head = T)
FisherWT_dis2_T316A_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_dis2_T316A_term_res.txt", sep = ""), head = T)
# WT7 vs Dis2-T316A
FisherWT_dis2_T316D_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_dis2_T316D_gb_res.txt", sep = ""), head = T)
FisherWT_dis2_T316D_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_dis2_T316D_pr_res.txt", sep = ""), head = T)
FisherWT_dis2_T316D_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_dis2_T316D_term_res.txt", sep = ""), head = T)



#################################################################################################
# function for plotting MA-plots by parsing out the sig. changed genes 
makeMA = function(DESeq_res, pval = 0.01, filename = 'Experiment.pdf', main = NULL, xlab = "Mean PRO-seq Counts (log10)", ylab = "fold change MB3-PP1/DMSO (log2)", geneList = NULL){
  if (is.null(geneList)){
    DESeq_res = DESeq_res
  }
  else{
    DESeq_res_filt = merge(DESeq_res, geneList, by.x = 0, by.y = 4)[,c(1:7)]
    rownames(DESeq_res_filt) = DESeq_res_filt[,1]
    DESeq_res = DESeq_res_filt[,c(2:7)]
  }
  unchangedGenes = DESeq_res[complete.cases(DESeq_res) & DESeq_res$padj>= pval, ]
  SigUpgenes = DESeq_res[complete.cases(DESeq_res) & DESeq_res$padj<= pval & DESeq_res$log2FoldChange > 0 , ]
  SigDownGenes = DESeq_res[complete.cases(DESeq_res) & DESeq_res$padj<= pval & DESeq_res$log2FoldChange < 0 , ]
  cat("Unchanged Genes = ", length(unchangedGenes[,1]), "\n", "Sig Up genes = ", length(SigUpgenes[,1]), "\n", "SigDownGenes = ", length(SigDownGenes[,1]), "\n" )
  pdf(file = paste(fig_dir, filename, sep = ""), width = 7, height = 7)
  plot(log10(unchangedGenes$baseMean), unchangedGenes$log2FoldChange, ylim = c(-3, 3), pch = 16, 
       main = main, xlab = xlab, ylab = ylab, panel.first = grid( nx = NULL, lty = 1, lwd = 2),
       cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, col = adjustcolor("black", alpha = 0.5))
  points(log10(SigUpgenes$baseMean), SigUpgenes$log2FoldChange, col = adjustcolor("steelblue", alpha = 0.5), pch = 16)
  points(log10(SigDownGenes$baseMean), SigDownGenes$log2FoldChange, col = adjustcolor("forestgreen",alpha = 0.5), pch = 16)
  legend("topleft", c(sprintf("Sig Up: N = %i", length(SigUpgenes[,1])), sprintf("Unchanged: N = %i", length(unchangedGenes[,1])), sprintf("Sig Down: N = %i", length(SigDownGenes[,1]))), text.col = c("steelblue", "black", "forestgreen"), cex = 1.5)
  abline(h=0, lty = 2, lwd = 2)
  dev.off()
}
## Apply the makeMA function to results dataframes generated above
# WT
makeMA(WT_gb_res, filename = "SP_WTgb_DESeqResults.pdf", main = "WT Change in Gene-body Counts", geneList = filteredGL)
makeMA(WT_pr_res, filename = "SP_WTpr_DESeqResults.pdf", main = "WT Change in Promoter-prox Counts", geneList = filteredGL)
makeMA(WT_term_res, filename = "SP_WTterm_DESeqResults.pdf", main = "WT Change in Post-CPS Counts", geneList = filteredGL)
# Mcs6-as
makeMA(mcs6_gb_res, filename = "SP_mcs6gb_DESeqResults.pdf", main = "mcs6-as Change in Gene-body Counts", geneList = filteredGL)
makeMA(mcs6_pr_res, filename = "SP_mcs6pr_DESeqResults.pdf", main = "mcs6-as Change in Promoter-prox Counts", geneList = filteredGL)
makeMA(mcs6_term_res, filename = "SP_mcs6term_DESeqResults.pdf", main = "mcs6-as Change in Post-CPS Counts", geneList = filteredGL)
# CDK9-as
makeMA(CDK9_gb_res, filename = "SP_CDK9gb_DESeqResults.pdf", main = "CDK9-as Change in Gene-body Counts", geneList = filteredGL)
makeMA(CDK9_pr_res, filename = "SP_CDK9pr_DESeqResults.pdf", main = "CDK9-as Change in Promoter-prox Counts", geneList = filteredGL)
makeMA(CDK9_term_res, filename = "SP_CDK9term_DESeqResults.pdf", main = "CDK9-as Change in Post-CPS Counts", geneList = filteredGL)
# Isk1-as
makeMA(Isk1_gb_res, filename = "SP_Isk1gb_DESeqResults.pdf", main = "Isk1-as Change in Gene-body Counts", geneList = filteredGL)
makeMA(Isk1_pr_res, filename = "SP_Isk1pr_DESeqResults.pdf", main = "Isk1-as Change in Promoter-prox Counts", geneList = filteredGL)
makeMA(Isk1_term_res, filename = "SP_Isk1term_DESeqResults.pdf", main = "Isk1-as Change in Post-CPS Counts", geneList = filteredGL)
# mcs6-as CDK9-as (double mutant)
makeMA(mcs6_CDK9_gb_res, filename = "SP_mcs6_CDK9gb_DESeqResults.pdf", main = "mcs6-as CDK9-as Change in Gene-body Counts", geneList = filteredGL)
makeMA(mcs6_CDK9_pr_res, filename = "SP_mcs6_CDK9pr_DESeqResults.pdf", main = "mcs6-as CDK9-as Change in Promoter-prox Counts", geneList = filteredGL)
makeMA(mcs6_CDK9_term_res, filename = "SP_mcs6_CDK9term_DESeqResults.pdf", main = "mcs6-as CDK9-as Change in Post-CPS Counts", geneList = filteredGL)

# WT vs CDK9as Untreated
makeMA(WTvCDK9_gb_res, filename = "SP_WTvsCDK9_Untreated_gb_DESeqResults.pdf", main = "WT vs CDK9as Untreated gb Counts", geneList = filteredGL)
makeMA(WTvCDK9_pr_res, filename = "SP_WTvsCDK9_Untreated_pr_DESeqResults.pdf", main = "WT vs CDK9as Untreated pr Counts", geneList = filteredGL)
makeMA(WTvCDK9_term_res, filename = "SP_WTvsCDK9_Untreated_term_DESeqResults.pdf", main = "WT vs CDK9as Untreated Post-CPS Counts", geneList = filteredGL)
# WT vs mcs6as Untreated
makeMA(WTvMcs6_gb_res, filename = "SP_WTvsMcs6_Untreated_gb_DESeqResults.pdf", main = "WT vs Mcs6as Untreated gb Counts", geneList = filteredGL)
makeMA(WTvMcs6_pr_res, filename = "SP_WTvsMcs6_Untreated_pr_DESeqResults.pdf", main = "WT vs Mcs6as Untreated pr Counts", geneList = filteredGL)
makeMA(WTvMcs6_term_res, filename = "SP_WTvsMcs6_Untreated_term_DESeqResults.pdf", main = "WT vs Mcs6as Untreated Post-CPS Counts", geneList = filteredGL)
# WT vs Isk1as Untreated
makeMA(WTvIsk1_gb_res, filename = "SP_WTvsIsk1_Untreated_gb_DESeqResults.pdf", main = "WT vs Isk1as Untreated gb Counts", geneList = filteredGL)
makeMA(WTvIsk1_pr_res, filename = "SP_WTvsIsk1_Untreated_pr_DESeqResults.pdf", main = "WT vs Isk1as Untreated pr Counts", geneList = filteredGL)
makeMA(WTvIsk1_term_res, filename = "SP_WTvsIsk1_Untreated_term_DESeqResults.pdf", main = "WT vs Isk1as Untreated Post-CPS Counts", geneList = filteredGL)
# WT vs mcs6as_CDK9as (double mutant) Untreated
makeMA(WTvsMcs6_CDK9_gb_res, filename = "SP_WTvsMcs6_CDK9_Untreated_gb_DESeqResults.pdf", main = "WT vs Mcs6as CDK9as Untreated gb Counts", geneList = filteredGL)
makeMA(WTvsMcs6_CDK9_pr_res, filename = "SP_WTvsMcs6_CDK9Untreated_pr_DESeqResults.pdf", main = "WT vs Mcs6as CDK9as Untreated pr Counts", geneList = filteredGL)
makeMA(WTvsMcs6_CDK9_term_res, filename = "SP_WTvsMcs6_CDK9_Untreated_term_DESeqResults.pdf", main = "WT vs Mcs6as CDK9as Untreated Post-CPS Counts", geneList = filteredGL)

# Cdk9as_18C
makeMA(Cdk9as_18C_gb_res, filename = "Cdk9as_18C_gb_DESeqResults.pdf", main = "Cdk9as 18C gb Counts", geneList = filteredGL)
makeMA(Cdk9as_18C_pr_res, filename = "Cdk9as_18C_pr_DESeqResults.pdf", main = "Cdk9as 18C pr Counts", geneList = filteredGL)
makeMA(Cdk9as_18C_term_res, filename = "Cdk9as_18C_term_DESeqResults.pdf", main = "Cdk9as 18C Post-CPS Counts", geneList = filteredGL)
# dis2ts_18C
makeMA(dis2ts_18C_gb_res, filename = "dis2ts_18C_gb_DESeqResults.pdf", main = "dis2ts 18C gb Counts", geneList = filteredGL)
makeMA(dis2ts_18C_pr_res, filename = "dis2ts_18C_pr_DESeqResults.pdf", main = "dis2ts 18C pr Counts", geneList = filteredGL)
makeMA(dis2ts_18C_term_res, filename = "dis2ts_18C_term_DESeqResults.pdf", main = "dis2ts 18C Post-CPS Counts", geneList = filteredGL)
# dis2ts_30C
makeMA(dis2ts_30C_gb_res, filename = "dis2ts_30C_gb_DESeqResults.pdf", main = "dis2ts_30C gb Counts", geneList = filteredGL)
makeMA(dis2ts_30C_pr_res, filename = "dis2ts_30C_pr_DESeqResults.pdf", main = "dis2ts_30C pr Counts", geneList = filteredGL)
makeMA(dis2ts_30C_term_res, filename = "dis2ts_30C_term_DESeqResults.pdf", main = "dis2ts_30C Post-CPS Counts", geneList = filteredGL)
# Cdk9as_dis2ts_18C
makeMA(CDK9_dis2ts_18C_gb_res, filename = "Cdk9as_dis2ts_18C_gb_DESeqResults.pdf", main = "Cdk9as_dis2ts 18C gb Counts", geneList = filteredGL)
makeMA(CDK9_dis2ts_18C_pr_res, filename = "Cdk9as_dis2ts_18C_pr_DESeqResults.pdf", main = "Cdk9as_dis2ts 18C pr Counts", geneList = filteredGL)
makeMA(CDK9_dis2ts_18C_term_res, filename = "Cdk9as_dis2ts_18C_term_DESeqResults.pdf", main = "Cdk9as_dis2ts 18C Post-CPS Counts", geneList = filteredGL)
# Cdk9as_dis2ts_30C
makeMA(CDK9_dis2ts_30C_gb_res, filename = "Cdk9as_dis2ts_30C_gb_DESeqResults.pdf", main = "Cdk9as_dis2ts 30C gb Counts", geneList = filteredGL)
makeMA(CDK9_dis2ts_30C_pr_res, filename = "Cdk9as_dis2ts_30C_pr_DESeqResults.pdf", main = "Cdk9as_dis2ts 30C pr Counts", geneList = filteredGL)
makeMA(CDK9_dis2ts_30C_term_res, filename = "Cdk9as_dis2ts_30C_term_DESeqResults.pdf", main = "Cdk9as_dis2ts 30C Post-CPS Counts", geneList = filteredGL)
# dis2ts_30Cv18C_DMSO
makeMA(dis2ts_30Cv18C_DMSO_gb_res, filename = "dis2ts_30Cv18C_DMSO_gb_DESeqResults.pdf", main = "dis2ts 30C vs 18C DMSO gb Counts", geneList = filteredGL)
makeMA(dis2ts_30Cv18C_DMSO_pr_res, filename = "dis2ts_30Cv18C_DMSO_pr_DESeqResults.pdf", main = "dis2ts 30C vs 18C DMSO  pr Counts", geneList = filteredGL)
makeMA(dis2ts_30Cv18C_DMSO_term_res, filename = "dis2ts_30Cv18C_DMSO_term_DESeqResults.pdf", main = "dis2ts 30C vs 18C DMSO  Post-CPS Counts", geneList = filteredGL)
# dis2ts_30Cv18C_3MBPP1
makeMA(dis2ts_30Cv18C_3MBPP1_gb_res, filename = "dis2ts_30Cv18C_3MBPP1_gb_DESeqResults.pdf", main = "dis2ts 30C vs 18C 3MBPP1 gb Counts", geneList = filteredGL)
makeMA(dis2ts_30Cv18C_3MBPP1_pr_res, filename = "dis2ts_30Cv18C_3MBPP1_pr_DESeqResults.pdf", main = "dis2ts 30C vs 18C 3MBPP1 pr Counts", geneList = filteredGL)
makeMA(dis2ts_30Cv18C_3MBPP1_term_res, filename = "dis2ts_30Cv18C_3MBPP1_term_DESeqResults.pdf", main = "dis2ts 30C vs 18C 3MBPP1  Post-CPS Counts", geneList = filteredGL)
# Cdk9as_dis2ts_30Cv18C_DMSO
makeMA(CDK9_dis2ts_30Cv18C_DMSO_gb_res, filename = "Cdk9as_dis2ts_30Cv18C_DMSO_gb_DESeqResults.pdf", main = "Cdk9as_dis2ts 30C vs 18C DMSO gb Counts", geneList = filteredGL)
makeMA(CDK9_dis2ts_30Cv18C_DMSO_pr_res, filename = "Cdk9as_dis2ts_30Cv18C_DMSO_pr_DESeqResults.pdf", main = "Cdk9as_dis2ts 30C vs 18C DMSO pr Counts", geneList = filteredGL)
makeMA(CDK9_dis2ts_30Cv18C_DMSO_term_res, filename = "Cdk9as_dis2ts_30Cv18C_DMSO_term_DESeqResults.pdf", main = "Cdk9as_dis2ts 30C vs 18C DMSO Post-CPS Counts", geneList = filteredGL)
# Cdk9as_dis2ts_30Cv18C_3MBPP1
makeMA(CDK9_dis2ts_30Cv18C_3MBPP1_gb_res, filename = "Cdk9as_dis2ts_30Cv18C_3MBPP1_gb_DESeqResults.pdf", main = "Cdk9as_dis2ts 30C vs 18C 3MBPP1 gb Counts", geneList = filteredGL)
makeMA(CDK9_dis2ts_30Cv18C_3MBPP1_pr_res, filename = "Cdk9as_dis2ts_30Cv18C_3MBPP1_pr_DESeqResults.pdf", main = "Cdk9as_dis2ts 30C vs 18C 3MBPP1 pr Counts", geneList = filteredGL)
makeMA(CDK9_dis2ts_30Cv18C_3MBPP1_term_res, filename = "Cdk9as_dis2ts_30Cv18C_3MBPP1_term_DESeqResults.pdf", main = "Cdk9as_dis2ts 30C vs 18C 3MBPP1 Post-CPS Counts", geneList = filteredGL)


## WT vs Mutants (untreated)
#WT vs Mcs6as 
makeMA(WTvMcs6_gb_res, filename = "WTvMcs6_DMSO_gb_DESeqResults.pdf", main = "WT vs Mcs6as DMSO gb Counts", geneList = filteredGL)
makeMA(WTvMcs6_pr_res, filename = "WTvMcs6_DMSO_pr_DESeqResults.pdf", main = "WT vs Mcs6as DMSO pr Counts", geneList = filteredGL)
makeMA(WTvMcs6_term_res, filename = "WTvMcs6_DMSO_term_DESeqResults.pdf", main = "WT vs Mcs6as DMSO Post-CPS Counts", geneList = filteredGL)
#WT vs Cdk9as
makeMA(WTvCDK9_gb_res, filename = "WTvCdk9_DMSO_gb_DESeqResults.pdf", main = "WT vs Cdk9as DMSO gb Counts", geneList = filteredGL)
makeMA(WTvCDK9_pr_res, filename = "WTvCdk9_DMSO_pr_DESeqResults.pdf", main = "WT vs Cdk9as DMSO pr Counts", geneList = filteredGL)
makeMA(WTvCDK9_term_res, filename = "WTvCdk9_DMSO_term_DESeqResults.pdf", main = "WT vs Cdk9as DMSO Post-CPS Counts", geneList = filteredGL)
#WT vs Mcs6as Cdk9as
makeMA(WTvsMcs6_CDK9_gb_res, filename = "WTvMcs6_Cdk9_DMSO_gb_DESeqResults.pdf", main = "WT vs Mcs6_Cdk9as DMSO gb Counts", geneList = filteredGL)
makeMA(WTvsMcs6_CDK9_pr_res, filename = "WTvMcs6_Cdk9_DMSO_pr_DESeqResults.pdf", main = "WT vs Mcs6_Cdk9as DMSO pr Counts", geneList = filteredGL)
makeMA(WTvsMcs6_CDK9_term_res, filename = "WTvMcs6_Cdk9_DMSO_term_DESeqResults.pdf", main = "WT vs Mcs6_Cdk9as DMSO Post-CPS Counts", geneList = filteredGL)
#WT vs Isk1as
makeMA(WTvIsk1_gb_res, filename = "WTvIsk1_DMSO_gb_DESeqResults.pdf", main = "WT vs Isk1as DMSO gb Counts", geneList = filteredGL)
makeMA(WTvIsk1_pr_res, filename = "WTvIsk1_DMSO_pr_DESeqResults.pdf", main = "WT vs Isk1as DMSO pr Counts", geneList = filteredGL)
makeMA(WTvIsk1_term_res, filename = "WTvIsk1_DMSO_term_DESeqResults.pdf", main = "WT vs Isk1as DMSO Post-CPS Counts", geneList = filteredGL)

### Time Course of Cdk9as inhibition.  (treatments vs. 0 mins)
# 30 sec vs 0 min gb res
makeMA(CDK9_0.5vs0min_gb_res, filename = "CDK9_30secvs0min_gb_DESeqResults.pdf", main = "CDK9as 30sec vs 0min gb Counts", geneList = filteredGL)
makeMA(CDK9_0.5vs0min_pr_res, filename = "CDK9_30secvs0min_pr_DESeqResults.pdf", main = "CDK9as 30sec vs 0min pr Counts", geneList = filteredGL)
makeMA(CDK9_0.5vs0min_term_res, filename = "CDK9_30secvs0min_term_DESeqResults.pdf", main = "CDK9as 30sec vs 0min Post-CPS Counts", geneList = filteredGL)
# 1 min vs 0 min gb res
makeMA(CDK9_1vs0min_gb_res, filename = "CDK9_1vs0min_gb_DESeqResults.pdf", main = "CDK9as 1min vs 0min gb Counts", geneList = filteredGL)
makeMA(CDK9_1vs0min_pr_res, filename = "CDK9_1vs0min_pr_DESeqResults.pdf", main = "CDK9as 1min vs 0min pr Counts", geneList = filteredGL)
makeMA(CDK9_1vs0min_term_res, filename = "CDK9_1vs0min_term_DESeqResults.pdf", main = "CDK9as 1min vs 0min Post-CPS Counts", geneList = filteredGL)
# 2.5 min vs 0 min gb res
makeMA(CDK9_2.5vs0min_gb_res, filename = "CDK9_2min30secvs0min_gb_DESeqResults.pdf", main = "CDK9as 2.5min vs 0min gb Counts", geneList = filteredGL)
makeMA(CDK9_2.5vs0min_pr_res, filename = "CDK9_2min30secvs0min_pr_DESeqResults.pdf", main = "CDK9as 2.5min vs 0min pr Counts", geneList = filteredGL)
makeMA(CDK9_2.5vs0min_term_res, filename = "CDK9_2min30secvs0min_term_DESeqResults.pdf", main = "CDK9as 2.5min vs 0min Post-CPS Counts", geneList = filteredGL)
# 5 min vs 0 min gb res
makeMA(CDK9_5vs0min_gb_res, filename = "CDK9_5vs0min_gb_DESeqResults.pdf", main = "CDK9as 5min vs 0min gb Counts", geneList = filteredGL)
makeMA(CDK9_5vs0min_pr_res, filename = "CDK9_5vs0min_pr_DESeqResults.pdf", main = "CDK9as 5min vs 0min pr Counts", geneList = filteredGL)
makeMA(CDK9_5vs0min_term_res, filename = "CDK9_5vs0min_term_DESeqResults.pdf", main = "CDK9as 5min vs 0min Post-CPS Counts", geneList = filteredGL)
# 1 min vs 0 min gb res
makeMA(CDK9_7.5vs0min_gb_res, filename = "CDK9_7min30secvs0min_gb_DESeqResults.pdf", main = "CDK9as 7.5min vs 0min gb Counts", geneList = filteredGL)
makeMA(CDK9_7.5vs0min_pr_res, filename = "CDK9_7min30secvs0min_pr_DESeqResults.pdf", main = "CDK9as 7.5min vs 0min pr Counts", geneList = filteredGL)
makeMA(CDK9_7.5vs0min_term_res, filename = "CDK9_7min30secvs0min_term_DESeqResults.pdf", main = "CDK9as 7.5min vs 0min Post-CPS Counts", geneList = filteredGL)
# 10 min vs 0 min gb res
makeMA(CDK9_10vs0min_gb_res, filename = "CDK9_10vs0min_gb_DESeqResults.pdf", main = "CDK9as 10min vs 0min gb Counts", geneList = filteredGL)
makeMA(CDK9_10vs0min_pr_res, filename = "CDK9_10vs0min_pr_DESeqResults.pdf", main = "CDK9as 10min vs 0min pr Counts", geneList = filteredGL)
makeMA(CDK9_10vs0min_term_res, filename = "CDK9_10vs0min_term_DESeqResults.pdf", main = "CDK9as 10min vs 0min Post-CPS Counts", geneList = filteredGL)
# 20 min vs 0 min gb res
makeMA(CDK9_20vs0min_gb_res, filename = "CDK9_20vs0min_gb_DESeqResults.pdf", main = "CDK9as 20min vs 0min gb Counts", geneList = filteredGL)
makeMA(CDK9_20vs0min_pr_res, filename = "CDK9_20vs0min_pr_DESeqResults.pdf", main = "CDK9as 20min vs 0min pr Counts", geneList = filteredGL)
makeMA(CDK9_20vs0min_term_res, filename = "CDK9_20vs0min_term_DESeqResults.pdf", main = "CDK9as 20min vs 0min Post-CPS Counts", geneList = filteredGL)
# Lsk1-as CDK9-as (double mutant)
makeMA(Lsk1_CDK9_gb_res, filename = "SP_Lsk1_CDK9gb_DESeqResults.pdf", main = "Lsk1-as CDK9-as Change in Gene-body Counts", geneList = filteredGL)
makeMA(Lsk1_CDK9_pr_res, filename = "SP_Lsk1_CDK9pr_DESeqResults.pdf", main = "Lsk1-as CDK9-as Change in Promoter-prox Counts", geneList = filteredGL)
makeMA(Lsk1_CDK9_term_res, filename = "SP_Lsk1_CDK9term_DESeqResults.pdf", main = "Lsk1-as CDK9-as Change in Post-CPS Counts", geneList = filteredGL)

# Cdk9as Time Course in Spt5 CTR mutants 
# spt5WT7 1 minutes vs 0 minutes
makeMA(Spt5WT7_1minVs0min_gb_res, filename = "SP_Spt5WT7_1minVs0min_gb_DESeqResults.pdf", main = "Spt5WT7 1min Vs 0min Change in Gene-body Counts", geneList = filteredGL)
makeMA(Spt5WT7_1minVs0min_pr_res, filename = "SP_Spt5WT7_1minVs0min_pr_DESeqResults.pdf", main = "Spt5WT7 1min Vs 0min Change in Promoter-prox Counts", geneList = filteredGL)
makeMA(Spt5WT7_1minVs0min_term_res, filename = "SP_Spt5WT7_1minVs0min_term_DESeqResults.pdf", main = "Spt5WT7 1min Vs 0min Change in Post-CPS Counts", geneList = filteredGL)
# spt5WT7 5 minutes vs 0 minutes
makeMA(Spt5WT7_5minVs0min_gb_res, filename = "SP_Spt5WT7_5minVs0min_gb_DESeqResults.pdf", main = "Spt5WT7 5min Vs 0min Change in Gene-body Counts", geneList = filteredGL)
makeMA(Spt5WT7_5minVs0min_pr_res, filename = "SP_Spt5WT7_5minVs0min_pr_DESeqResults.pdf", main = "Spt5WT7 5min Vs 0min Change in Promoter-prox Counts", geneList = filteredGL)
makeMA(Spt5WT7_5minVs0min_term_res, filename = "SP_Spt5WT7_5minVs0min_term_DESeqResults.pdf", main = "Spt5WT7 5min Vs 0min Change in Post-CPS Counts", geneList = filteredGL)
# Spt5T1A 1 minutes vs 0 minutes
makeMA(Spt5T1A_1minVs0min_gb_res, filename = "SP_Spt5T1A_1minVs0min_gb_DESeqResults.pdf", main = "Spt5T1A 1min Vs 0min Change in Gene-body Counts", geneList = filteredGL)
makeMA(Spt5T1A_1minVs0min_pr_res, filename = "SP_Spt5T1A_1minVs0min_pr_DESeqResults.pdf", main = "Spt5T1A 1min Vs 0min Change in Promoter-prox Counts", geneList = filteredGL)
makeMA(Spt5T1A_1minVs0min_term_res, filename = "SP_Spt5T1A_1minVs0min_term_DESeqResults.pdf", main = "Spt5T1A 1min Vs 0min Change in Post-CPS Counts", geneList = filteredGL)
# Spt5T1A 5 minutes vs 0 minutes
makeMA(Spt5T1A_5minVs0min_gb_res, filename = "SP_Spt5T1A_5minVs0min_gb_DESeqResults.pdf", main = "Spt5T1A 5min Vs 0min Change in Gene-body Counts", geneList = filteredGL)
makeMA(Spt5T1A_5minVs0min_pr_res, filename = "SP_Spt5T1A_5minVs0min_pr_DESeqResults.pdf", main = "Spt5T1A 5min Vs 0min Change in Promoter-prox Counts", geneList = filteredGL)
makeMA(Spt5T1A_5minVs0min_term_res, filename = "SP_Spt5T1A_5minVs0min_term_DESeqResults.pdf", main = "Spt5T1A 5min Vs 0min Change in Post-CPS Counts", geneList = filteredGL)
# Spt5T1E 1 minutes vs 0 minutes
makeMA(Spt5T1E_1minVs0min_gb_res, filename = "SP_Spt5T1E_1minVs0min_gb_DESeqResults.pdf", main = "Spt5T1E 1min Vs 0min Change in Gene-body Counts", geneList = filteredGL)
makeMA(Spt5T1E_1minVs0min_pr_res, filename = "SP_Spt5T1E_1minVs0min_pr_DESeqResults.pdf", main = "Spt5T1E 1min Vs 0min Change in Promoter-prox Counts", geneList = filteredGL)
makeMA(Spt5T1E_1minVs0min_term_res, filename = "SP_Spt5T1E_1minVs0min_term_DESeqResults.pdf", main = "Spt5T1E 1min Vs 0min Change in Post-CPS Counts", geneList = filteredGL)
# Spt5T1E 5 minutes vs 0 minutes
makeMA(Spt5T1E_5minVs0min_gb_res, filename = "SP_Spt5T1E_5minVs0min_gb_DESeqResults.pdf", main = "Spt5T1E 5min Vs 0min Change in Gene-body Counts", geneList = filteredGL)
makeMA(Spt5T1E_5minVs0min_pr_res, filename = "SP_Spt5T1E_5minVs0min_pr_DESeqResults.pdf", main = "Spt5T1E 5min Vs 0min Change in Promoter-prox Counts", geneList = filteredGL)
makeMA(Spt5T1E_5minVs0min_term_res, filename = "SP_Spt5T1E_5minVs0min_term_DESeqResults.pdf", main = "Spt5T1E 5min Vs 0min Change in Post-CPS Counts", geneList = filteredGL)
# Spt5WT7 vs spt5T1A 0 minutes
makeMA(Spt5WT7_Spt5T1A_0min_gb_res, filename = "SP_Spt5WT7_Spt5T1A_0min_gb_DESeqResults.pdf", main = "Spt5T1A vs Spt5WT7 Change in Gene-body Counts", geneList = filteredGL)
makeMA(Spt5WT7_Spt5T1A_0min_pr_res, filename = "SP_Spt5WT7_Spt5T1A_0min_pr_DESeqResults.pdf", main = "Spt5T1A vs Spt5WT7 Change in Promoter-prox Counts", geneList = filteredGL)
makeMA(Spt5WT7_Spt5T1A_0min_term_res, filename = "SP_Spt5WT7_Spt5T1A_0min_term_DESeqResults.pdf", main = "Spt5T1A vs Spt5WT7 Change in Post-CPS Counts", geneList = filteredGL)
# Spt5WT7 vs spt5T1E 0 minutes
makeMA(Spt5WT7_Spt5T1E_0min_gb_res, filename = "SP_Spt5WT7_Spt5T1E_0min_gb_DESeqResults.pdf", main = "Spt5T1E vs Spt5WT7 Change in Gene-body Counts", geneList = filteredGL)
makeMA(Spt5WT7_Spt5T1E_0min_pr_res, filename = "SP_Spt5WT7_Spt5T1E_0min_pr_DESeqResults.pdf", main = "Spt5T1E vs Spt5WT7 Change in Promoter-prox Counts", geneList = filteredGL)
makeMA(Spt5WT7_Spt5T1E_0min_term_res, filename = "SP_Spt5WT7_Spt5T1E_0min_term_DESeqResults.pdf", main = "Spt5T1E vs Spt5WT7 Change in Post-CPS Counts", geneList = filteredGL)
# Comparisons between different Dis2 mutants:
# WT7 vs Dis2-11 
makeMA(FisherWT_dis2_11_gb_res, filename = "SP_FisherWT_dis2_11_gb_DESeqResults.pdf", main = "WT vs dis2-11 Change in Gene-body Counts", geneList = filteredGL)
makeMA(FisherWT_dis2_11_pr_res, filename = "SP_FisherWT_dis2_11_pr_DESeqResults.pdf", main = "WT vs dis2-11 Change in Promoter-prox Counts", geneList = filteredGL)
makeMA(FisherWT_dis2_11_term_res, filename = "SP_FisherWT_dis2_11_term_DESeqResults.pdf", main = "WT vs dis2-11 Change in Post-CPS Counts", geneList = filteredGL)
# WT7 vs Dis2Delete 
makeMA(FisherWT_dis2Del_gb_res, filename = "SP_FisherWT_dis2Del_gb_DESeqResults.pdf", main = "WT vs dis2Del Change in Gene-body Counts", geneList = filteredGL)
makeMA(FisherWT_dis2Del_pr_res, filename = "SP_FisherWT_dis2Del_pr_DESeqResults.pdf", main = "WT vs dis2Del Change in Promoter-prox Counts", geneList = filteredGL)
makeMA(FisherWT_dis2Del_term_res, filename = "SP_FisherWT_dis2Del_term_DESeqResults.pdf", main = "WT vs dis2Del Change in Post-CPS Counts", geneList = filteredGL)
# WT7 vs Dis2-T316A
makeMA(FisherWT_dis2_T316A_gb_res, filename = "SP_FisherWT_dis2_T316A_gb_DESeqResults.pdf", main = "WT vs dis2_T316A Change in Gene-body Counts", geneList = filteredGL)
makeMA(FisherWT_dis2_T316A_pr_res, filename = "SP_FisherWT_dis2_T316A_pr_DESeqResults.pdf", main = "WT vs dis2_T316A Change in Promoter-prox Counts", geneList = filteredGL)
makeMA(FisherWT_dis2_T316A_term_res, filename = "SP_FisherWT_dis2_T316A_term_DESeqResults.pdf", main = "WT vs dis2_T316A Change in Post-CPS Counts", geneList = filteredGL)
# WT7 vs Dis2-T316D
makeMA(FisherWT_dis2_T316D_gb_res, filename = "SP_FisherWT_dis2_T316D_gb_DESeqResults.pdf", main = "WT vs dis2_T316D Change in Gene-body Counts", geneList = filteredGL)
makeMA(FisherWT_dis2_T316D_pr_res, filename = "SP_FisherWT_dis2_T316D_pr_DESeqResults.pdf", main = "WT vs dis2_T316D Change in Promoter-prox Counts", geneList = filteredGL)
makeMA(FisherWT_dis2_T316D_term_res, filename = "SP_FisherWT_dis2_T316D_term_DESeqResults.pdf", main = "WT vs dis2_T316D Change in Post-CPS Counts", geneList = filteredGL)






