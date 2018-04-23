
#source('/Volumes/SEAGATE_EXP/Greg_YeastData/YeastAnalysisScripts/Pombe/ObservedTSS/Calc_SP_PI_data_07_2015.r')
source("/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Greg_YeastData/YeastAnalysisScripts/GB_Functions.R")
library(hexbin)
library(bigWig)
library(lattice)
library(ggplot2)
library(gplots)
library(MASS)
fig_dir = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/figures/ReplicateCorrelation/12-14-17/"
dir.create(fig_dir)
countpath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/countData/AllGenesObsTSS/"

#########################################################################################################
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
# 4th experiment: Cdk9as Time course (Full TC from new data)
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
NFs = read.table(file = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/countData/replicate_spikeNormFactors.txt", head = T)
WT_DMSOr1NF = NFs[1, "Spikereads"]/100000
WT_DMSOr2NF = NFs[2, "Spikereads"]/100000
WT_3MBPP1r1NF = NFs[3, "Spikereads"]/100000
WT_3MBPP1r2NF = NFs[4, "Spikereads"]/100000
mcs6as_DMSOr1NF = NFs[5, "Spikereads"]/100000
mcs6as_DMSOr2NF = NFs[6, "Spikereads"]/100000
mcs6as_3MBPP1r1NF = NFs[7, "Spikereads"]/100000
mcs6as_3MBPP1r2NF = NFs[8, "Spikereads"]/100000
Cdk9as_DMSOr1NF = NFs[9, "Spikereads"]/100000
Cdk9as_DMSOr2NF = NFs[10, "Spikereads"]/100000
Cdk9as_3MBPP1r1NF = NFs[11, "Spikereads"]/100000
Cdk9as_3MBPP1r2NF = NFs[12, "Spikereads"]/100000
Isk1as_DMSOr1NF = NFs[13, "Spikereads"]/100000
Isk1as_DMSOr2NF = NFs[14, "Spikereads"]/100000
Isk1as_3MBPP1r1NF = NFs[15, "Spikereads"]/100000
Isk1as_3MBPP1r2NF = NFs[16, "Spikereads"]/100000
mcs6as_Cdk9as_DMSOr1NF = NFs[17, "Spikereads"]/100000
mcs6as_Cdk9as_DMSOr2NF = NFs[18, "Spikereads"]/100000
mcs6as_Cdk9as_3MBPP1r1NF = NFs[19, "Spikereads"]/100000
mcs6as_Cdk9as_3MBPP1r2NF = NFs[20, "Spikereads"]/100000
CDK9_18C_DMSOr1NF = NFs[21, "Spikereads"]/100000
CDK9_18C_DMSOr2NF = NFs[22, "Spikereads"]/100000
CDK9_18C_3MBPP1r1NF = NFs[23, "Spikereads"]/100000
CDK9_18C_3MBPP1r2NF = NFs[24, "Spikereads"]/100000
dis2_18C_DMSOr1NF = NFs[25, "Spikereads"]/100000
dis2_18C_DMSOr2NF = NFs[26, "Spikereads"]/100000
dis2_18C_3MBPP1r1NF = NFs[27, "Spikereads"]/100000
dis2_18C_3MBPP1r2NF = NFs[28, "Spikereads"]/100000
dis2_30C_DMSOr1NF = NFs[29, "Spikereads"]/100000
dis2_30C_DMSOr2NF = NFs[30, "Spikereads"]/100000
dis2_30C_3MBPP1r1NF = NFs[31, "Spikereads"]/100000
dis2_30C_3MBPP1r2NF = NFs[32, "Spikereads"]/100000
CDK9_dis2_18C_DMSOr1NF = NFs[33, "Spikereads"]/100000
CDK9_dis2_18C_DMSOr2NF = NFs[34, "Spikereads"]/100000
CDK9_dis2_18C_3MBPP1r1NF = NFs[35, "Spikereads"]/100000
CDK9_dis2_18C_3MBPP1r2NF = NFs[36, "Spikereads"]/100000
CDK9_dis2_30C_DMSOr1NF = NFs[37, "Spikereads"]/100000
CDK9_dis2_30C_DMSOr2NF = NFs[38, "Spikereads"]/100000
CDK9_dis2_30C_3MBPP1r1NF = NFs[39, "Spikereads"]/100000
CDK9_dis2_30C_3MBPP1r2NF = NFs[40, "Spikereads"]/100000
CDK9_0minr1NF = NFs[51, "Spikereads"]/100000
CDK9_0minr2NF = NFs[52, "Spikereads"]/100000
CDK9_0.5minr1NF = NFs[53, "Spikereads"]/100000
CDK9_0.5minr2NF = NFs[54, "Spikereads"]/100000
CDK9_1minr1NF = NFs[55, "Spikereads"]/100000
CDK9_1minr2NF = NFs[56, "Spikereads"]/100000
CDK9_2.5minr1NF = NFs[57, "Spikereads"]/100000
CDK9_2.5minr2NF = NFs[58, "Spikereads"]/100000
CDK9_5minr1NF = NFs[59, "Spikereads"]/100000
CDK9_5minr2NF = NFs[60, "Spikereads"]/100000
CDK9_7.5minr1NF = NFs[61, "Spikereads"]/100000
CDK9_7.5minr2NF = NFs[62, "Spikereads"]/100000
CDK9_10minr1NF = NFs[63, "Spikereads"]/100000
CDK9_10minr2NF = NFs[64, "Spikereads"]/100000
CDK9_20minr1NF = NFs[65, "Spikereads"]/100000
CDK9_20minr2NF = NFs[66, "Spikereads"]/100000
Lsk1as_Cdk9as_DMSOr1NF = NFs[67, "Spikereads"]/100000
Lsk1as_Cdk9as_DMSOr2NF = NFs[68, "Spikereads"]/100000
Lsk1as_Cdk9as_3MBPP1r1NF = NFs[69, "Spikereads"]/100000
Lsk1as_Cdk9as_3MBPP1r2NF = NFs[70, "Spikereads"]/100000
Cdk9as_spt5_WT7_0min_DMSO_r1NF = NFs[81, "Spikereads"]/100000
Cdk9as_spt5_WT7_0min_DMSO_r2NF = NFs[82, "Spikereads"]/100000
Cdk9as_spt5_WT7_1min_3MBPP1_r1NF = NFs[83, "Spikereads"]/100000
Cdk9as_spt5_WT7_1min_3MBPP1_r2NF = NFs[84, "Spikereads"]/100000
Cdk9as_spt5_WT7_5min_3MBPP1_r1NF = NFs[85, "Spikereads"]/100000
Cdk9as_spt5_WT7_5min_3MBPP1_r2NF = NFs[86, "Spikereads"]/100000
Cdk9as_spt5_T1A_0min_DMSO_r1NF = NFs[87, "Spikereads"]/100000
Cdk9as_spt5_T1A_0min_DMSO_r2NF = NFs[88, "Spikereads"]/100000
Cdk9as_spt5_T1A_1min_3MBPP1_r1NF = NFs[89, "Spikereads"]/100000
Cdk9as_spt5_T1A_1min_3MBPP1_r2NF = NFs[90, "Spikereads"]/100000
Cdk9as_spt5_T1A_5min_3MBPP1_r1NF = NFs[91, "Spikereads"]/100000
Cdk9as_spt5_T1A_5min_3MBPP1_r2NF = NFs[92, "Spikereads"]/100000
Cdk9as_spt5_T1E_0min_DMSO_r1NF = NFs[93, "Spikereads"]/100000
Cdk9as_spt5_T1E_0min_DMSO_r2NF = NFs[94, "Spikereads"]/100000
Cdk9as_spt5_T1E_1min_3MBPP1_r1NF = NFs[95, "Spikereads"]/100000
Cdk9as_spt5_T1E_1min_3MBPP1_r2NF = NFs[96, "Spikereads"]/100000
Cdk9as_spt5_T1E_5min_3MBPP1_r1NF = NFs[97, "Spikereads"]/100000
Cdk9as_spt5_T1E_5min_3MBPP1_r2NF = NFs[98, "Spikereads"]/100000
Fisher_WT_r1NF = NFs[71, "Spikereads"]/100000
Fisher_WT_r2NF = NFs[72, "Spikereads"]/100000
Fisher_dis2_11_r1NF = NFs[73, "Spikereads"]/100000
Fisher_dis2_11_r2NF = NFs[74, "Spikereads"]/100000
Fisher_dis2Del_r1NF = NFs[75, "Spikereads"]/100000
Fisher_dis2Del_r2NF = NFs[76, "Spikereads"]/100000
Fisher_dis2_T316A_r1NF = NFs[77, "Spikereads"]/100000
Fisher_dis2_T316A_r2NF = NFs[78, "Spikereads"]/100000
Fisher_dis2_T316D_r1NF = NFs[79, "Spikereads"]/100000
Fisher_dis2_T316D_r2NF = NFs[80, "Spikereads"]/100000


df_List_GB = list(WT_D1, WT_D2, WT_M1, WT_M2, MCS6_D1, MCS6_D2,
               MCS6_M1, MCS6_M2, CDK9_D1, CDK9_D2,  CDK9_M1, 
               CDK9_M2, Isk1as_D1, Isk1as_D2, Isk1as_M1, Isk1as_M2,
               mcs6as_Cdk9as_D1, mcs6as_Cdk9as_D2, mcs6as_Cdk9as_M1, mcs6as_Cdk9as_M2,
               ePRO_Cdk9as_0min_r1, ePRO_Cdk9as_0min_r2, ePRO_Cdk9as_30sec_r1, ePRO_Cdk9as_30sec_r2, ePRO_Cdk9as_1min_r1, ePRO_Cdk9as_1min_r2, 
               ePRO_Cdk9as_2min30sec_r1, ePRO_Cdk9as_2min30sec_r2, ePRO_Cdk9as_5min_r1, ePRO_Cdk9as_5min_r2, ePRO_Cdk9as_7min30sec_r1, ePRO_Cdk9as_7min30sec_r2, 
               ePRO_Cdk9as_10min_r1, ePRO_Cdk9as_10min_r2, ePRO_Cdk9as_20min_r1, ePRO_Cdk9as_20min_r2, Lsk1as_Cdk9as_D1, Lsk1as_Cdk9as_D2, Lsk1as_Cdk9as_M1, Lsk1as_Cdk9as_M2)

df_List_Spt5Muts = list(Cdk9as_spt5_WT7_0min_DMSO_r1, Cdk9as_spt5_WT7_0min_DMSO_r2, Cdk9as_spt5_WT7_1min_3MBPP1_r1, Cdk9as_spt5_WT7_1min_3MBPP1_r2, 
                        Cdk9as_spt5_WT7_5min_3MBPP1_r1, Cdk9as_spt5_WT7_5min_3MBPP1_r2, Cdk9as_spt5_T1A_0min_DMSO_r1, Cdk9as_spt5_T1A_0min_DMSO_r2, 
                        Cdk9as_spt5_T1A_1min_3MBPP1_r1, Cdk9as_spt5_T1A_1min_3MBPP1_r2, Cdk9as_spt5_T1A_5min_3MBPP1_r1, Cdk9as_spt5_T1A_5min_3MBPP1_r2, 
                        Cdk9as_spt5_T1E_0min_DMSO_r1, Cdk9as_spt5_T1E_0min_DMSO_r2, Cdk9as_spt5_T1E_1min_3MBPP1_r1, Cdk9as_spt5_T1E_1min_3MBPP1_r2, 
                        Cdk9as_spt5_T1E_5min_3MBPP1_r1, Cdk9as_spt5_T1E_5min_3MBPP1_r2)

df_List_Pabitra = list(CDK9_18C_D1, CDK9_18C_D2, CDK9_18C_M1, CDK9_18C_M2, Dis2ts_18_D1, Dis2ts_18_D2, Dis2ts_18_M1, Dis2ts_18_M2,
                      Dis2ts_30_D1, Dis2ts_30_D2, Dis2ts_30_M1, Dis2ts_30_M2, Cdk9as_dis2ts_18_D1, Cdk9as_dis2ts_18_D2, Cdk9as_dis2ts_18_M1, Cdk9as_dis2ts_18_M2, 
                      Cdk9as_dis2ts_30_D1, Cdk9as_dis2ts_30_D2, Cdk9as_dis2ts_30_M1, Cdk9as_dis2ts_30_M2)

df_List_Fisher = list(Fisher_WT_r1, Fisher_WT_r2, Fisher_dis2_11_r1, Fisher_dis2_11_r2,Fisher_dis2Del_r1,
                      Fisher_dis2Del_r2, Fisher_dis2_T316A_r1, Fisher_dis2_T316A_r2, Fisher_dis2_T316D_r1, 
                      Fisher_dis2_T316D_r2)

NF_List_GB = list(WT_DMSOr1NF, WT_DMSOr2NF, WT_3MBPP1r1NF, WT_3MBPP1r2NF, 
               mcs6as_DMSOr1NF, mcs6as_DMSOr2NF, mcs6as_3MBPP1r1NF, mcs6as_3MBPP1r2NF,
               Cdk9as_DMSOr1NF, Cdk9as_DMSOr2NF, Cdk9as_3MBPP1r1NF, Cdk9as_3MBPP1r2NF,
               Isk1as_DMSOr1NF, Isk1as_DMSOr2NF, Isk1as_3MBPP1r1NF, Isk1as_3MBPP1r2NF, 
               mcs6as_Cdk9as_DMSOr1NF, mcs6as_Cdk9as_DMSOr2NF, mcs6as_Cdk9as_3MBPP1r1NF, mcs6as_Cdk9as_3MBPP1r2NF,
               CDK9_0minr1NF, CDK9_0minr2NF, CDK9_0.5minr1NF, CDK9_0.5minr2NF, CDK9_1minr1NF, CDK9_1minr2NF,
               CDK9_2.5minr1NF, CDK9_2.5minr2NF, CDK9_5minr1NF, CDK9_5minr2NF, CDK9_7.5minr1NF, CDK9_7.5minr2NF,
               CDK9_10minr1NF, CDK9_10minr2NF, CDK9_20minr1NF, CDK9_20minr2NF, Lsk1as_Cdk9as_DMSOr1NF, Lsk1as_Cdk9as_DMSOr2NF, Lsk1as_Cdk9as_3MBPP1r1NF, Lsk1as_Cdk9as_3MBPP1r2NF)

NF_List_Spt5Muts = list(Cdk9as_spt5_WT7_0min_DMSO_r1NF, Cdk9as_spt5_WT7_0min_DMSO_r2NF, Cdk9as_spt5_WT7_1min_3MBPP1_r1NF, Cdk9as_spt5_WT7_1min_3MBPP1_r2NF, 
                        Cdk9as_spt5_WT7_5min_3MBPP1_r1NF, Cdk9as_spt5_WT7_5min_3MBPP1_r2NF, Cdk9as_spt5_T1A_0min_DMSO_r1NF, Cdk9as_spt5_T1A_0min_DMSO_r2NF, 
                        Cdk9as_spt5_T1A_1min_3MBPP1_r1NF, Cdk9as_spt5_T1A_1min_3MBPP1_r2NF, Cdk9as_spt5_T1A_5min_3MBPP1_r1NF, Cdk9as_spt5_T1A_5min_3MBPP1_r2NF, 
                        Cdk9as_spt5_T1E_0min_DMSO_r1NF, Cdk9as_spt5_T1E_0min_DMSO_r2NF, Cdk9as_spt5_T1E_1min_3MBPP1_r1NF,  Cdk9as_spt5_T1E_1min_3MBPP1_r2NF, 
                        Cdk9as_spt5_T1E_5min_3MBPP1_r1NF, Cdk9as_spt5_T1E_5min_3MBPP1_r2NF)

NF_List_Pabitra = list(CDK9_18C_DMSOr1NF, CDK9_18C_DMSOr2NF, CDK9_18C_3MBPP1r1NF, CDK9_18C_3MBPP1r2NF,
                  dis2_18C_DMSOr1NF, dis2_18C_DMSOr2NF, dis2_18C_3MBPP1r1NF, dis2_18C_3MBPP1r2NF,
                  dis2_30C_DMSOr1NF, dis2_30C_DMSOr2NF, dis2_30C_3MBPP1r1NF, dis2_30C_3MBPP1r2NF,
                  CDK9_dis2_18C_DMSOr1NF, CDK9_dis2_18C_DMSOr2NF, CDK9_dis2_18C_3MBPP1r1NF, CDK9_dis2_18C_3MBPP1r2NF,
                  CDK9_dis2_30C_DMSOr1NF, CDK9_dis2_30C_DMSOr2NF, CDK9_dis2_30C_3MBPP1r1NF, CDK9_dis2_30C_3MBPP1r2NF)

NF_List_Fisher = list(Fisher_WT_r1NF, Fisher_WT_r2NF, Fisher_dis2_11_r1NF, Fisher_dis2_11_r2NF,Fisher_dis2Del_r1NF,
                      Fisher_dis2Del_r2NF, Fisher_dis2_T316A_r1NF, Fisher_dis2_T316A_r2NF, Fisher_dis2_T316D_r1NF, 
                      Fisher_dis2_T316D_r2NF)

SampleList_GB = c("WT_DMSO", "WT_DMSO", "WT_3MBPP1", "WT_3MBPP1","mcs6as_DMSO", "mcs6as_DMSO", "mcs6as_3MBPP1", "mcs6as_3MBPP1", 
               "Cdk9as_DMSO", "Cdk9as_DMSO", "Cdk9_3MBPP1", "Cdk9_3MBPP1", "Isk1as_DMSO", "Isk1as_DMSO", "Isk1as_3MBPP1", 
               "Isk1as_3MBPP1", "mcs6as_Cdk9as_DMSO",  "mcs6as_Cdk9as_DMSO",  "mcs6as_Cdk9as_3MBPP1", "mcs6as_Cdk9as_3MBPP1",
               "ePRO_Cdk9as_0min_r1", "ePRO_Cdk9as_0min_r2", "ePRO_Cdk9as_30sec_r1", "ePRO_Cdk9as_30sec_r2", "ePRO_Cdk9as_1min_r1", "ePRO_Cdk9as_1min_r2", 
               "ePRO_Cdk9as_2min30sec_r1", "ePRO_Cdk9as_2min30sec_r2", "ePRO_Cdk9as_5min_r1", "ePRO_Cdk9as_5min_r2", "ePRO_Cdk9as_7min30sec_r1", "ePRO_Cdk9as_7min30sec_r2", 
               "ePRO_Cdk9as_10min_r1", "ePRO_Cdk9as_10min_r2", "ePRO_Cdk9as_20min_r1", "ePRO_Cdk9as_20min_r2", "Lsk1as_Cdk9as_D1", "Lsk1as_Cdk9as_D2", "Lsk1as_Cdk9as_M1", "Lsk1as_Cdk9as_M2")

SampleList_Spt5Muts = c("Cdk9as_spt5_WT7_0min_DMSO", "Cdk9as_spt5_WT7_0min_DMSO", "Cdk9as_spt5_WT7_1min_3MBPP1", "Cdk9as_spt5_WT7_1min_3MBPP1", 
                        "Cdk9as_spt5_WT7_5min_3MBPP1", "Cdk9as_spt5_WT7_5min_3MBPP1", "Cdk9as_spt5_T1A_0min_DMSO", "Cdk9as_spt5_T1A_0min_DMSO", 
                        "Cdk9as_spt5_T1A_1min_3MBPP1", "Cdk9as_spt5_T1A_1min_3MBPP1", "Cdk9as_spt5_T1A_5min_3MBPP1", "Cdk9as_spt5_T1A_5min_3MBPP1", 
                        'Cdk9as_spt5_T1E_0min_DMSO', "Cdk9as_spt5_T1E_0min_DMSO", "Cdk9as_spt5_T1E_1min_3MBPP1", "Cdk9as_spt5_T1E_1min_3MBPP1", 
                        "Cdk9as_spt5_T1E_5min_3MBPP1", "Cdk9as_spt5_T1E_5min_3MBPP1")

SampleList_Pabitra = c("CDK9_18C_DMSO", "CDK9_18C_DMSO", "CDK9_18C_3MBPP1", "CDK9_18C_3MBPP1",
                    "dis2_18C_DMSO", "dis2_18C_DMSO", "dis2_18C_3MBPP1", "dis2_18C_3MBPP1",
                    "dis2_30C_DMSO", "dis2_30C_DMSO", "dis2_30C_3MBPP1", "dis2_30C_3MBPP1",
                    "CDK9_dis2_18C_DMSO", "CDK9_dis2_18C_DMSO", "CDK9_dis2_18C_3MBPP1", "CDK9_dis2_18C_3MBPP1",
                    "CDK9_dis2_30C_DMSO", "CDK9_dis2_30C_DMSO", "CDK9_dis2_30C_3MBPP1", "CDK9_dis2_30C_3MBPP1")

SampleList_Fisher = c("Fisher_WT", "Fisher_WT", "Fisher_dis2_11", "Fisher_dis2_11", "Fisher_dis2Del",
                      "Fisher_dis2Del", "Fisher_dis2_T316A", "Fisher_dis2_T316A", "Fisher_dis2_T316D", 
                      "Fisher_dis2_T316D")

ConstructCountTables = function(df_List = df_List,  NF_List = NF_List, Samplenames = SampleList){
  N = length(df_List)
  UnNormed = cbind(df_List[[1]][,1], df_List[[2]][,1], df_List[[1]][,2], df_List[[2]][,2], rep(Samplenames[1], dim(df_List[[1]])[1]))
  RPKMnormed = cbind(df_List[[1]][,8], df_List[[2]][,8], df_List[[1]][,9], df_List[[2]][,9], rep(Samplenames[1], dim(df_List[[1]])[1]))
  SpikeNormed = cbind(df_List[[1]][,1]/as.numeric(NF_List[1]), df_List[[2]][,1]/as.numeric(NF_List[2]), df_List[[1]][,2]/as.numeric(NF_List[1]), df_List[[2]][,2]/as.numeric(NF_List[2]), rep(Samplenames[1], dim(df_List[[1]])[1]))
  for (i in seq(3, N, 2)){
    Un = cbind(df_List[[i]][,1], df_List[[i+1]][,1], df_List[[i]][,2], df_List[[i+1]][,2], rep(Samplenames[i], dim(df_List[[i]])[1]))
    RPKM = cbind(df_List[[i]][,8], df_List[[i+1]][,8], df_List[[i]][,9], df_List[[i+1]][,9], rep(Samplenames[i], dim(df_List[[i]])[1]))
    Spike = cbind(df_List[[i]][,1]/as.numeric(NF_List[i]), df_List[[i+1]][,1]/as.numeric(NF_List[i+1]), df_List[[i]][,2]/as.numeric(NF_List[i]), df_List[[i+1]][,2]/as.numeric(NF_List[i+1]), rep(Samplenames[i], dim(df_List[[i]])[1]))
    UnNormed = rbind(UnNormed, Un)
    RPKMnormed = rbind(RPKMnormed, RPKM)
    SpikeNormed = rbind(SpikeNormed, Spike)
  }
  UNrow_sub = apply(UnNormed, 1, function(row) all(row !=0 ))  ### function for removing rows with 0's
  RPKMrow_sub = apply(RPKMnormed, 1, function(row) all(row !=0 ))
  Spikerow_sub = apply(SpikeNormed, 1, function(row) all(row !=0 ))
  result = list(UnNormed[UNrow_sub,], RPKMnormed[RPKMrow_sub,], SpikeNormed[Spikerow_sub,])
  return(result)
}

replicateSample_data_GB = ConstructCountTables(df_List = df_List_GB,  NF_List = NF_List_GB, Samplenames = SampleList_GB)
replicateSample_data_Pabitra = ConstructCountTables(df_List = df_List_Pabitra,  NF_List = NF_List_Pabitra, Samplenames = SampleList_Pabitra)
replicateSample_data_Spt5Muts = ConstructCountTables(df_List = df_List_Spt5Muts,  NF_List = NF_List_Spt5Muts, Samplenames = SampleList_Spt5Muts)
replicateSample_data_Fisher = ConstructCountTables(df_List = df_List_Fisher,  NF_List = NF_List_Fisher, Samplenames = SampleList_Fisher)

pdf(file = paste(fig_dir, "geneBody_spikeNorm_BoothSamples.pdf", sep = ""), width = 20, height = 20)
hexbinplot( log10(as.numeric(replicateSample_data_GB[[3]][,3])) ~ log10(as.numeric(replicateSample_data_GB[[3]][,4])) | replicateSample_data_GB[[3]][,5], 
           xbins = 100, colramp = rich.colors,
           as.table = T, # plots right to left, bottom to top
           #index.cond=list(c(10,8,6,4,2,9,7,5,3,1)), # provides the order of the panels
           aspect = 1,
           ylab = list(label = "GB reads(log10) Rep1", fontsize = 20), 
           xlab = list(label = "GB reads(log10) Rep2", fontsize =20), type = 'g',
           scales = list(cex = 1.5)) + latticeExtra::layer(panel.abline(a = 0, b = 1, lty = 2, lwd = 3, col = 'violet'))# + latticeExtra::layer(panel.key(c(sprintf("Rho = %f", cor$estimate), sprintf("cor_Pval = %f", cor$p.value), sprintf("N = %i", Ngenes)), corner = c(0,.98), lines = FALSE, points = FALSE, cex = 1.5), packets = 1)
dev.off()

pdf(file = paste(fig_dir, "geneBody_spikeNorm_ParuaSamples.pdf", sep = ""), width = 20, height = 20)
hexbinplot( log10(as.numeric(replicateSample_data_Pabitra[[3]][,3])) ~ log10(as.numeric(replicateSample_data_Pabitra[[3]][,4])) | replicateSample_data_Pabitra[[3]][,5], 
            xbins = 100, colramp = rich.colors,
            as.table = T, # plots right to left, bottom to top
            #index.cond=list(c(10,8,6,4,2,9,7,5,3,1)), # provides the order of the panels
            aspect = 1,
            ylab = list(label = "GB reads(log10) Rep1", fontsize = 20), 
            xlab = list(label = "GB reads(log10) Rep2", fontsize =20), type = 'g',
            scales = list(cex = 1.5)) + latticeExtra::layer(panel.abline(a = 0, b = 1, lty = 2, lwd = 3, col = 'violet'))# + latticeExtra::layer(panel.key(c(sprintf("Rho = %f", cor$estimate), sprintf("cor_Pval = %f", cor$p.value), sprintf("N = %i", Ngenes)), corner = c(0,.98), lines = FALSE, points = FALSE, cex = 1.5), packets = 1)
dev.off()

pdf(file = paste(fig_dir, "geneBody_spikeNorm_Cdk9as_Spt5Mut_TC_Samples.pdf", sep = ""), width = 10, height = 10)
hexbinplot( log10(as.numeric(replicateSample_data_Spt5Muts[[3]][,3])) ~ log10(as.numeric(replicateSample_data_Spt5Muts[[3]][,4])) | replicateSample_data_Spt5Muts[[3]][,5], 
            xbins = 100, colramp = rich.colors,
            as.table = T, # plots right to left, bottom to top
            #index.cond=list(c(10,8,6,4,2,9,7,5,3,1)), # provides the order of the panels
            aspect = 1,
            ylab = list(label = "GB reads(log10) Rep1", fontsize = 20), 
            xlab = list(label = "GB reads(log10) Rep2", fontsize =20), type = 'g',
            scales = list(cex = 1.5)) + latticeExtra::layer(panel.abline(a = 0, b = 1, lty = 2, lwd = 3, col = 'violet'))# + latticeExtra::layer(panel.key(c(sprintf("Rho = %f", cor$estimate), sprintf("cor_Pval = %f", cor$p.value), sprintf("N = %i", Ngenes)), corner = c(0,.98), lines = FALSE, points = FALSE, cex = 1.5), packets = 1)
dev.off()

pdf(file = paste(fig_dir, "geneBody_spikeNorm_FisherSamples.pdf", sep = ""), width = 10, height = 10)
hexbinplot( log10(as.numeric(replicateSample_data_Fisher[[3]][,3])) ~ log10(as.numeric(replicateSample_data_Fisher[[3]][,4])) | replicateSample_data_Fisher[[3]][,5], 
            xbins = 100, colramp = rich.colors,
            as.table = T, # plots right to left, bottom to top
            #index.cond=list(c(10,8,6,4,2,9,7,5,3,1)), # provides the order of the panels
            aspect = 1,
            ylab = list(label = "GB reads(log10) Rep1", fontsize = 20), 
            xlab = list(label = "GB reads(log10) Rep2", fontsize =20), type = 'g',
            scales = list(cex = 1.5)) + latticeExtra::layer(panel.abline(a = 0, b = 1, lty = 2, lwd = 3, col = 'violet'))# + latticeExtra::layer(panel.key(c(sprintf("Rho = %f", cor$estimate), sprintf("cor_Pval = %f", cor$p.value), sprintf("N = %i", Ngenes)), corner = c(0,.98), lines = FALSE, points = FALSE, cex = 1.5), packets = 1)
dev.off()
###############################################################################
## want to get all the correlation coeffs. 
## Get Correlation values for samples used in Booth paper 
sink(file = paste(fig_dir, "correlationCoefficientsBooth.txt", sep = ""))
cat("Spearmans correlation coefficients for each sample:\n")
for (sample in levels(as.factor(replicateSample_data_GB[[3]][,5]))){
  cat(paste("Correlation values for sample: ", sample, "\n"))
  SampleDat = replicateSample_data_GB[[3]][replicateSample_data_GB[[3]][,5]==sample, ]
  x = cor.test(log10(as.numeric(SampleDat[,3])), log10(as.numeric(SampleDat[,4])), method = "spearman")
  print(x)}
sink()
## Get Correlation values for samples used in Parua paper 
sink(file = paste(fig_dir, "correlationCoefficientsParua.txt", sep = ""))
cat("Spearmans correlation coefficients for each sample:\n")
for (sample in levels(as.factor(replicateSample_data_Pabitra[[3]][,5]))){
  cat(paste("Correlation values for sample: ", sample, "\n"))
  SampleDat = replicateSample_data_Pabitra[[3]][replicateSample_data_Pabitra[[3]][,5]==sample, ]
  x = cor.test(log10(as.numeric(SampleDat[,3])), log10(as.numeric(SampleDat[,4])), method = "spearman")
  print(x)}
sink()

