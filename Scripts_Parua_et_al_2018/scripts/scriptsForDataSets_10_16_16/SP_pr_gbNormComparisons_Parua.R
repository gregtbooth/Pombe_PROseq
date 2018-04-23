
#source('/Volumes/SEAGATE_EXP/Greg_YeastData/YeastAnalysisScripts/Pombe/ObservedTSS/Calc_SP_PI_data_07_2015.r')
source("/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/scripts/GB_Functions.R")
library(hexbin)
library(bigWig)
library(lattice)
library(ggplot2)
library(gplots)
library(MASS)
fig_dir = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/figures/Dis2Analysis/ReplicateCorrelation/03-01-18/"
dir.create(fig_dir)
countpath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/countData/AllGenesObsTSS/"

#########################################################################################################
# 02-13-18 Reduced the input files to only those relevant to the Dis2 project. 
#########################################################################################################

# 1st experiment
WT_D1 = read.table(file = paste(countpath, "WT_DMSOr1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
WT_D2 = read.table(file = paste(countpath, "WT_DMSOr2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
# 3rd experiment (not using any 3-MB-PP1 treated samples)
CDK9_18C_D1 = read.table(file = paste(countpath, "CDK9_18C_DMSOr1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
CDK9_18C_D2 = read.table(file = paste(countpath, "CDK9_18C_DMSOr2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Dis2ts_18_D1 = read.table(file = paste(countpath, "Dis2ts_18_DMSOr1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Dis2ts_18_D2 = read.table(file = paste(countpath, "Dis2ts_18_DMSOr2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Dis2ts_30_D1 = read.table(file = paste(countpath, "Dis2ts_30_DMSOr1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Dis2ts_30_D2 = read.table(file = paste(countpath, "Dis2ts_30_DMSOr2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Cdk9as_dis2ts_18_D1 = read.table(file = paste(countpath, "Cdk9as_dis2ts_18_DMSOr1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Cdk9as_dis2ts_18_D2 = read.table(file = paste(countpath, "Cdk9as_dis2ts_18_DMSOr2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Cdk9as_dis2ts_30_D1 = read.table(file = paste(countpath, "Cdk9as_dis2ts_30_DMSOr1PIcountData_Fix.txt", sep = ""), head = T, row.names=1)
Cdk9as_dis2ts_30_D2 = read.table(file = paste(countpath, "Cdk9as_dis2ts_30_DMSOr2PIcountData_Fix.txt", sep = ""), head = T, row.names=1)

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
CDK9_18C_DMSOr1NF = NFs[21, "Spikereads"]/100000
CDK9_18C_DMSOr2NF = NFs[22, "Spikereads"]/100000
dis2_18C_DMSOr1NF = NFs[25, "Spikereads"]/100000
dis2_18C_DMSOr2NF = NFs[26, "Spikereads"]/100000
dis2_30C_DMSOr1NF = NFs[29, "Spikereads"]/100000
dis2_30C_DMSOr2NF = NFs[30, "Spikereads"]/100000
CDK9_dis2_18C_DMSOr1NF = NFs[33, "Spikereads"]/100000
CDK9_dis2_18C_DMSOr2NF = NFs[34, "Spikereads"]/100000
CDK9_dis2_30C_DMSOr1NF = NFs[37, "Spikereads"]/100000
CDK9_dis2_30C_DMSOr2NF = NFs[38, "Spikereads"]/100000
#
#Fisher_WT_r1NF = NFs[71, "Spikereads"]/100000
#Fisher_WT_r2NF = NFs[72, "Spikereads"]/100000
#Fisher_dis2_11_r1NF = NFs[73, "Spikereads"]/100000
#Fisher_dis2_11_r2NF = NFs[74, "Spikereads"]/100000
#Fisher_dis2Del_r1NF = NFs[75, "Spikereads"]/100000
#Fisher_dis2Del_r2NF = NFs[76, "Spikereads"]/100000
#Fisher_dis2_T316A_r1NF = NFs[77, "Spikereads"]/100000
#Fisher_dis2_T316A_r2NF = NFs[78, "Spikereads"]/100000
#Fisher_dis2_T316D_r1NF = NFs[79, "Spikereads"]/100000
#Fisher_dis2_T316D_r2NF = NFs[80, "Spikereads"]/100000

# NFs now based on the ratio of ratios ((spikeRNA/sampleRNA)/(spikeDNA/sampleDNA))
WGS_NFs = read.table(file = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/scripts/scriptsForPabitra/countData/WGSrepNormFactorTable.txt", head = T)
# try dividing spike counts by the ratio of spike/sample DNA
WGS_NFs$spikeRNA_DNAcorr = WGS_NFs[,3]/WGS_NFs[,7]

Fisher_WT_r1NF = 1/WGS_NFs[1, 8]
Fisher_WT_r2NF = 1/WGS_NFs[2, 8]
Fisher_dis2_11_r1NF = 1/WGS_NFs[3, 8]
Fisher_dis2_11_r2NF = 1/WGS_NFs[4, 8]
Fisher_dis2Del_r1NF = 1/WGS_NFs[5, 8]
Fisher_dis2Del_r2NF = 1/WGS_NFs[6, 8]
Fisher_dis2_T316A_r1NF = 1/WGS_NFs[7, 8]
Fisher_dis2_T316A_r2NF = 1/WGS_NFs[8, 8]
Fisher_dis2_T316D_r1NF = 1/WGS_NFs[9, 8]
Fisher_dis2_T316D_r2NF = 1/WGS_NFs[10, 8]


df_List_Pabitra = list(WT_D1, WT_D2, CDK9_18C_D1, CDK9_18C_D2, Dis2ts_18_D1, Dis2ts_18_D2,
                      Dis2ts_30_D1, Dis2ts_30_D2, Cdk9as_dis2ts_18_D1, Cdk9as_dis2ts_18_D2,
                      Cdk9as_dis2ts_30_D1, Cdk9as_dis2ts_30_D2)

df_List_Fisher = list(Fisher_WT_r1, Fisher_WT_r2, Fisher_dis2_11_r1, Fisher_dis2_11_r2,Fisher_dis2Del_r1,
                      Fisher_dis2Del_r2, Fisher_dis2_T316A_r1, Fisher_dis2_T316A_r2, Fisher_dis2_T316D_r1, 
                      Fisher_dis2_T316D_r2)

NF_List_Pabitra = list(WT_DMSOr1NF, WT_DMSOr2NF,CDK9_18C_DMSOr1NF, CDK9_18C_DMSOr2NF, 
                  dis2_18C_DMSOr1NF, dis2_18C_DMSOr2NF,
                  dis2_30C_DMSOr1NF, dis2_30C_DMSOr2NF,
                  CDK9_dis2_18C_DMSOr1NF, CDK9_dis2_18C_DMSOr2NF,
                  CDK9_dis2_30C_DMSOr1NF, CDK9_dis2_30C_DMSOr2NF)

NF_List_Fisher = list(Fisher_WT_r1NF, Fisher_WT_r2NF, Fisher_dis2_11_r1NF, Fisher_dis2_11_r2NF,Fisher_dis2Del_r1NF,
                      Fisher_dis2Del_r2NF, Fisher_dis2_T316A_r1NF, Fisher_dis2_T316A_r2NF, Fisher_dis2_T316D_r1NF, 
                      Fisher_dis2_T316D_r2NF)

SampleList_Pabitra = c("WT_DMSO", "WT_DMSO", "CDK9_18C_DMSO", "CDK9_18C_DMSO",
                    "dis2_18C_DMSO", "dis2_18C_DMSO", 
                    "dis2_30C_DMSO", "dis2_30C_DMSO", 
                    "CDK9_dis2_18C_DMSO", "CDK9_dis2_18C_DMSO", 
                    "CDK9_dis2_30C_DMSO", "CDK9_dis2_30C_DMSO")

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

replicateSample_data_Pabitra = ConstructCountTables(df_List = df_List_Pabitra,  NF_List = NF_List_Pabitra, Samplenames = SampleList_Pabitra)
replicateSample_data_Fisher = ConstructCountTables(df_List = df_List_Fisher,  NF_List = NF_List_Fisher, Samplenames = SampleList_Fisher)


pdf(file = paste(fig_dir, "geneBody_spikeNorm_ParuaSamples.pdf", sep = ""), width = 10, height = 20)
hexbinplot( log10(as.numeric(replicateSample_data_Pabitra[[3]][,3])) ~ log10(as.numeric(replicateSample_data_Pabitra[[3]][,4])) | replicateSample_data_Pabitra[[3]][,5], 
            xbins = 100, colramp = rich.colors,
            as.table = T, # plots right to left, bottom to top
            #index.cond=list(c(10,8,6,4,2,9,7,5,3,1)), # provides the order of the panels
            aspect = 1,
            ylab = list(label = "GB reads(log10) Rep1", fontsize = 20), 
            xlab = list(label = "GB reads(log10) Rep2", fontsize =20), type = 'g',
            scales = list(cex = 1.5)) + latticeExtra::layer(panel.abline(a = 0, b = 1, lty = 2, lwd = 3, col = 'violet'))# + latticeExtra::layer(panel.key(c(sprintf("Rho = %f", cor$estimate), sprintf("cor_Pval = %f", cor$p.value), sprintf("N = %i", Ngenes)), corner = c(0,.98), lines = FALSE, points = FALSE, cex = 1.5), packets = 1)
dev.off()


pdf(file = paste(fig_dir, "geneBody_spikeNorm_FisherSamples_WGSnormed.pdf", sep = ""), width = 10, height = 10)
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
## Get Correlation values for samples used in Parua paper 
sink(file = paste(fig_dir, "correlationCoefficientsParua.txt", sep = ""))
cat("Spearmans correlation coefficients for each sample:\n")
for (sample in levels(as.factor(replicateSample_data_Pabitra[[3]][,5]))){
  cat(paste("Correlation values for sample: ", sample, "\n"))
  SampleDat = replicateSample_data_Pabitra[[3]][replicateSample_data_Pabitra[[3]][,5]==sample, ]
  x = cor.test(log10(as.numeric(SampleDat[,3])), log10(as.numeric(SampleDat[,4])), method = "spearman")
  print(x)}
sink()

## Get Correlation values for samples from new (Dis2 Mutant) Fisher data
sink(file = paste(fig_dir, "correlationCoefficientsFisher.txt", sep = ""))
cat("Spearmans correlation coefficients for each sample:\n")
for (sample in levels(as.factor(replicateSample_data_Fisher[[3]][,5]))){
  cat(paste("Correlation values for sample: ", sample, "\n"))
  SampleDat = replicateSample_data_Fisher[[3]][replicateSample_data_Fisher[[3]][,5]==sample, ]
  x = cor.test(log10(as.numeric(SampleDat[,3])), log10(as.numeric(SampleDat[,4])), method = "spearman")
  print(x)}
sink()

#########################################################################################################
#########################################################################################################
# 03-01-18
#########################################################################################################
## XY plots between different samples for all replicate comparisons. 
## Want to assess reproducibility of sample effects. 
df_List_WTvs_dis2_11 = list(Fisher_WT_r1, Fisher_dis2_11_r1, Fisher_WT_r2,  Fisher_dis2_11_r2, Fisher_WT_r1, Fisher_dis2_11_r2, Fisher_WT_r2, Fisher_dis2_11_r1)
NF_List_WTvs_dis2_11 = list(Fisher_WT_r1NF, Fisher_dis2_11_r1NF, Fisher_WT_r2NF,  Fisher_dis2_11_r2NF, Fisher_WT_r1NF, Fisher_dis2_11_r2NF, Fisher_WT_r2NF, Fisher_dis2_11_r1NF)
SampleList = c("Fisher_WT", "Fisher_WT", "Fisher_dis2_11", "Fisher_dis2_11", "Fisher_dis2Del",
                            "Fisher_dis2Del", "Fisher_dis2_T316A", "Fisher_dis2_T316A")
replicateSample_WTvs_dis2_11 = ConstructCountTables(df_List = df_List_WTvs_dis2_11,  NF_List = NF_List_WTvs_dis2_11, Samplenames = SampleList)
pdf(file = paste(fig_dir, "geneBody_spikeNorm_WT_vs_dis2_11reps_WGSnormed.pdf", sep = ""), width = 10, height = 10)
hexbinplot( log10(as.numeric(replicateSample_WTvs_dis2_11[[3]][,3])) ~ log10(as.numeric(replicateSample_WTvs_dis2_11[[3]][,4])) | replicateSample_WTvs_dis2_11[[3]][,5], 
            xbins = 100, colramp = rich.colors,
            as.table = T, # plots right to left, bottom to top
            #index.cond=list(c(10,8,6,4,2,9,7,5,3,1)), # provides the order of the panels
            aspect = 1,
            ylab = list(label = "GB reads(log10) Rep1", fontsize = 20), 
            xlab = list(label = "GB reads(log10) Rep2", fontsize =20), type = 'g',
            scales = list(cex = 1.5)) + latticeExtra::layer(panel.abline(a = 0, b = 1, lty = 2, lwd = 3, col = 'violet'))# + latticeExtra::layer(panel.key(c(sprintf("Rho = %f", cor$estimate), sprintf("cor_Pval = %f", cor$p.value), sprintf("N = %i", Ngenes)), corner = c(0,.98), lines = FALSE, points = FALSE, cex = 1.5), packets = 1)
dev.off()

df_List_WTvs_dis2_del = list(Fisher_WT_r1, Fisher_dis2Del_r1, Fisher_WT_r2, Fisher_dis2Del_r2, Fisher_WT_r1, Fisher_dis2Del_r2, Fisher_WT_r2, Fisher_dis2Del_r1)
NF_List_WTvs_dis2_del = list(Fisher_WT_r1NF, Fisher_dis2Del_r1NF, Fisher_WT_r2NF, Fisher_dis2Del_r2NF, Fisher_WT_r1NF, Fisher_dis2Del_r2NF, Fisher_WT_r2NF, Fisher_dis2Del_r1NF)
replicateSampleWTvs_dis2_del = ConstructCountTables(df_List = df_List_WTvs_dis2_del,  NF_List = NF_List_WTvs_dis2_del, Samplenames = SampleList_Fisher)
pdf(file = paste(fig_dir, "geneBody_spikeNorm_WT_vs_dis2Delreps_WGSnormed.pdf", sep = ""), width = 10, height = 10)
hexbinplot( log10(as.numeric(replicateSampleWTvs_dis2_del[[3]][,3])) ~ log10(as.numeric(replicateSampleWTvs_dis2_del[[3]][,4])) | replicateSampleWTvs_dis2_del[[3]][,5], 
            xbins = 100, colramp = rich.colors,
            as.table = T, # plots right to left, bottom to top
            #index.cond=list(c(10,8,6,4,2,9,7,5,3,1)), # provides the order of the panels
            aspect = 1,
            ylab = list(label = "GB reads(log10) Rep1", fontsize = 20), 
            xlab = list(label = "GB reads(log10) Rep2", fontsize =20), type = 'g',
            scales = list(cex = 1.5)) + latticeExtra::layer(panel.abline(a = 0, b = 1, lty = 2, lwd = 3, col = 'violet'))# + latticeExtra::layer(panel.key(c(sprintf("Rho = %f", cor$estimate), sprintf("cor_Pval = %f", cor$p.value), sprintf("N = %i", Ngenes)), corner = c(0,.98), lines = FALSE, points = FALSE, cex = 1.5), packets = 1)
dev.off()

df_List_WTvs_dis2_T316D = list(Fisher_WT_r1, Fisher_dis2_T316D_r1, Fisher_WT_r2, Fisher_dis2_T316D_r2, Fisher_WT_r1, Fisher_dis2_T316D_r2, Fisher_WT_r2, Fisher_dis2_T316D_r1)
NF_List_WTvs_dis2_T316D = list(Fisher_WT_r1NF, Fisher_dis2_T316D_r1NF, Fisher_WT_r2NF, Fisher_dis2_T316D_r2NF, Fisher_WT_r1NF, Fisher_dis2_T316D_r2NF, Fisher_WT_r2NF, Fisher_dis2_T316D_r1NF)
replicateSample_WTvs_dis2_T316D = ConstructCountTables(df_List = df_List_WTvs_dis2_T316D,  NF_List = NF_List_WTvs_dis2_T316D, Samplenames = SampleList_Fisher)
pdf(file = paste(fig_dir, "geneBody_spikeNorm_WT_vs_dis2T316Dreps_WGSnormed.pdf", sep = ""), width = 10, height = 10)
hexbinplot( log10(as.numeric(replicateSample_WTvs_dis2_T316D[[3]][,3])) ~ log10(as.numeric(replicateSample_WTvs_dis2_T316D[[3]][,4])) | replicateSample_WTvs_dis2_T316D[[3]][,5], 
            xbins = 100, colramp = rich.colors,
            as.table = T, # plots right to left, bottom to top
            #index.cond=list(c(10,8,6,4,2,9,7,5,3,1)), # provides the order of the panels
            aspect = 1,
            ylab = list(label = "GB reads(log10) Rep1", fontsize = 20), 
            xlab = list(label = "GB reads(log10) Rep2", fontsize =20), type = 'g',
            scales = list(cex = 1.5)) + latticeExtra::layer(panel.abline(a = 0, b = 1, lty = 2, lwd = 3, col = 'violet'))# + latticeExtra::layer(panel.key(c(sprintf("Rho = %f", cor$estimate), sprintf("cor_Pval = %f", cor$p.value), sprintf("N = %i", Ngenes)), corner = c(0,.98), lines = FALSE, points = FALSE, cex = 1.5), packets = 1)
dev.off()

