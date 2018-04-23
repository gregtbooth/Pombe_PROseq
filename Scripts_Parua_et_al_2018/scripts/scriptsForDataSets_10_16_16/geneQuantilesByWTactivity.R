### This Script is just for generating a list of genes with thier activity quantiles. 
source("/Volumes/SEAGATE_EXP/Fisher_collaboration/scripts/GB_Functions.R")
library(bigWig)
infopath = "/Volumes/SEAGATE_EXP/Fisher_collaboration/pombeinfo/"
countpath = "/Volumes/SEAGATE_EXP/Fisher_collaboration/analysis/countData/"
outpath = "/Volumes/SEAGATE_EXP/Fisher_collaboration/scripts/scriptsForPabitra/OutPut/GenesByQuantile/05-17-19/"
dir.create(outpath)
filteredGL = read.table(paste(infopath, "SP_CompleteFilteredGenes.bed", sep = ""))
SPGL_1kbsep = read.table(file = paste(infopath, "SP_PROcapObservedTSS_1000bpSep.bed", sep = ""), stringsAsFactors = F)
SPGL_1kbsep_filt = merge(filteredGL, SPGL_1kbsep, by.x = 4, by.y = 4)[, c(2,3,4,1,5,6)]
WT_DMSO = read.table(file = paste(countpath, "AllGenesObsTSS/WT_DMSO_combinedPIcountData.txt", sep = ""), head = T)

# function returns a list of gene lists based on their genbody readDensity Quantile (based on WT_DMSO combined Reps)
activity_quantile.genes = function(bed, WT_df, quantiles = 4){
  bed = merge(x = bed, y = WT_df[,c(8,9)], by.x = 'V4', by.y = 'row.names')
  bed = bed[is.finite(bed[,8]) & bed[,8] > 0,][,c(2,3,4,1,5,6,7,8)]
  quants = quantile(bed[,8], prob = seq(0,1, length = quantiles+1))
  bed$quantile = 1
  for (i in 1:quantiles){
    bed[bed[,8] < quants[[i+1]] & bed[,8] >= quants[[i]],]$quantile = i
  }
  bed = bed[,c(1,2,3,4,5,6,8,9)]
  colnames(bed) = c("chr", "start", "end", "gene", "CapSig", "strand", "gbPROseqDensity", "quantile")
  return(bed)
}

## run separately on all filtered genes, and Kb sep genes (also filtered)
kbsepGenes_quants = activity_quantile.genes(SPGL_1kbsep_filt, WT_df = WT_DMSO, quantiles = 4)
filtGenes_quants = activity_quantile.genes(filteredGL, WT_df = WT_DMSO, quantiles = 4)
## write gene lists to file
write.table(kbsepGenes_quants, file = paste(outpath, "Genes_Filtered_KbSep_activityQuantile.bed", sep = ""), quote = F, sep ="\t", row.names = F, col.names = T)
write.table(filtGenes_quants, file = paste(outpath, "Genes_Filtered_activityQuantile.bed", sep = ""), quote = F, sep ="\t", row.names = F, col.names = T)

