source("/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/scripts/GB_Functions.R")
library(lattice)
library(grid)

bwpath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/alignment/bw/"
count_outpath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/countData/"

SPmappability = load.bigWig("/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/pombeInfo/SPombeMappable_SacCerNonOverLap_expanded_ChrFix.bw")
ObservedTSS = read.table(file = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/pombeInfo/SP_PROcapObservedTSS_noOverlap.bed", sep = "\t")
onlyObservedTSS = ObservedTSS[ObservedTSS[,5] > 0.457, ] # 0.457 is equal to 5 (RPM) normalized reads for WT PRO-cap data
bed = onlyObservedTSS

# get background counts and rate estimate 
pbw = load.bigWig(paste(bwpath, "5993_7157_26977_HJ72CBGXX_pombe_WT_COMBINED_DMSO_CAGATC_R1.fastq_pombe_plus.bw", sep = ""))
mbw = load.bigWig(paste(bwpath, "5993_7157_26977_HJ72CBGXX_pombe_WT_COMBINED_DMSO_CAGATC_R1.fastq_pombe_minus.bw", sep = ""))
BGregions = read.table("/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/pombeInfo/SP_backGround_regions.bed")
BGreg = BGregions[,c(1,6,2,3,4,5)]
BGcountTable = regionCountTable(pbw, mbw, map = SPmappability, genes = BGreg)
BGrate = mean(BGcountTable$region_density) ## will be used to define a poisson distribution for comparing GB counts against to determine activity

# Filter out genes with upstream Run-through and genes that are inactive (not sig different from background)
WTPItable = pause.index(pbw, mbw, SPmappability, bed[,c(1,6,2,3,4)])  ## calculated for identifying active genes
activity = apply(WTPItable, MARGIN = 1, FUN = ActiveGeneProb, lambda = BGrate) # get an probability of activity for each gene based on GB counts and Background rate
bed$P_activity = activity
cat("number of genes = ", length(bed[,1]), "number active = ", length(bed[bed$P_activity < 0.01, 1]))
bed$upstreamRatio = apply(bed, 1, Upstream.Ratio, pbw = pbw, mbw = mbw, map = SPmappability)
bedFilt1 = bed[bed$upstreamRatio < 1 & bed$P_activity < 0.01 , ] # if upstream to downstream ratio > 1 gene is omitted.
bedFilt = bedFilt1[complete.cases(bedFilt1),]
#SPGL = bedFilt[,c(1,6,2,3,4)]
SPGL = bedFilt[,c(1,6,2,3,4)]
# save active and un-influence gene list for downstream analysis. Write in standard Bed format.
write.table(bedFilt[, c(1,2,3,4,5,6)], "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/pombeInfo/SP_CompleteFilteredGenes.bed", quote = F, sep = "\t", row.names =F, col.names =F)
cat("number of genes being analyzed (filtered List) = ", length(SPGL[,1]), "\n")
####################################################################################################
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
load.CombinedspikeIn_wigset<- function(wigset_row) {
  file = CombinedspikeIn_wigset[wigset_row, 1]
  wig.p = NULL
  if (file != "")
    wig.p = load.bigWig(paste(bwpath, file, sep=''))
  file = CombinedspikeIn_wigset[wigset_row, 2]
  wig.m = NULL
  if (file != "")
    wig.m = load.bigWig(paste(bwpath, file, sep=''))
  return(list(wig.p, wig.m, CombinedspikeIn_wigset[wigset_row, 3]))
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

spikeIn_wigset = rbind(
  c("SpikeIn_bw/5993_7157_26971_HJ72CBGXX_pombe_WT_rep1_DMSO_ATCACG_R1.fastq_cerevisiae_plus.bw", "SpikeIn_bw/5993_7157_26971_HJ72CBGXX_pombe_WT_rep1_DMSO_ATCACG_R1.fastq_cerevisiae_minus.bw", "WT_DMSOr1"),
  c("SpikeIn_bw/5993_7157_26977_HJ72CBGXX_pombe_WT_rep2_DMSO_CAGATC_R1.fastq_cerevisiae_plus.bw", "SpikeIn_bw/5993_7157_26977_HJ72CBGXX_pombe_WT_rep2_DMSO_CAGATC_R1.fastq_cerevisiae_minus.bw", "WT_DMSOr2"),
  c("SpikeIn_bw/5993_7157_26972_HJ72CBGXX_pombe_WT_rep1_3MB-PPI_CGATGT_R1.fastq_cerevisiae_plus.bw", "SpikeIn_bw/5993_7157_26972_HJ72CBGXX_pombe_WT_rep1_3MB-PPI_CGATGT_R1.fastq_cerevisiae_minus.bw", "WT_3MBPP1r1"),
  c("SpikeIn_bw/5993_7157_26978_HJ72CBGXX_pombe_WT_rep2_3MB-PPI_ACTTGA_R1.fastq_cerevisiae_plus.bw", "SpikeIn_bw/5993_7157_26978_HJ72CBGXX_pombe_WT_rep2_3MB-PPI_ACTTGA_R1.fastq_cerevisiae_minus.bw", "WT_3MBPP1r2"),
  c("SpikeIn_bw/5993_7157_26973_HJ72CBGXX_pombe_mcs6as_rep1_DMSO_TTAGGC_R1.fastq_cerevisiae_plus.bw", "SpikeIn_bw/5993_7157_26973_HJ72CBGXX_pombe_mcs6as_rep1_DMSO_TTAGGC_R1.fastq_cerevisiae_minus.bw", "mcs6as_DMSOr1"),
  c("SpikeIn_bw/5993_7157_26979_HJ72CBGXX_pombe_mcs6as_rep2_DMSO_GATCAG_R1.fastq_cerevisiae_plus.bw", "SpikeIn_bw/5993_7157_26979_HJ72CBGXX_pombe_mcs6as_rep2_DMSO_GATCAG_R1.fastq_cerevisiae_minus.bw", "mcs6as_DMSOr2"),
  c("SpikeIn_bw/5993_7157_26974_HJ72CBGXX_pombe_mcs6as_rep1_3MB-PPI_TGACCA_R1.fastq_cerevisiae_plus.bw", "SpikeIn_bw/5993_7157_26974_HJ72CBGXX_pombe_mcs6as_rep1_3MB-PPI_TGACCA_R1.fastq_cerevisiae_minus.bw", "mcs6as_3MBPP1r1"),
  c("SpikeIn_bw/5993_7157_26980_HJ72CBGXX_pombe_mcs6as_rep2_3MB-PPI_TAGCTT_R1.fastq_cerevisiae_plus.bw", "SpikeIn_bw/5993_7157_26980_HJ72CBGXX_pombe_mcs6as_rep2_3MB-PPI_TAGCTT_R1.fastq_cerevisiae_minus.bw", "mcs6as_3MBPP1r2"),
  c("SpikeIn_bw/5993_7157_26975_HJ72CBGXX_pombe_CDK9as_rep1_DMSO_ACAGTG_R1.fastq_cerevisiae_plus.bw", "SpikeIn_bw/5993_7157_26975_HJ72CBGXX_pombe_CDK9as_rep1_DMSO_ACAGTG_R1.fastq_cerevisiae_minus.bw", "Cdk9as_DMSOr1"),
  c("SpikeIn_bw/5993_7157_26981_HJ72CBGXX_pombe_CDK9as_rep2_DMSO_GGCTAC_R1.fastq_cerevisiae_plus.bw", "SpikeIn_bw/5993_7157_26981_HJ72CBGXX_pombe_CDK9as_rep2_DMSO_GGCTAC_R1.fastq_cerevisiae_minus.bw", "Cdk9as_DMSOr2"),
  c("SpikeIn_bw/5993_7157_26976_HJ72CBGXX_pombe_CDK9as_rep1_3MB-PPI_GCCAAT_R1.fastq_cerevisiae_plus.bw", "SpikeIn_bw/5993_7157_26976_HJ72CBGXX_pombe_CDK9as_rep1_3MB-PPI_GCCAAT_R1.fastq_cerevisiae_minus.bw", "Cdk9as_3MBPP1r1"),
  c("SpikeIn_bw/5993_7157_26982_HJ72CBGXX_pombe_CDK9as_rep2_3MB-PPI_CTTGTA_R1.fastq_cerevisiae_plus.bw", "SpikeIn_bw/5993_7157_26982_HJ72CBGXX_pombe_CDK9as_rep2_3MB-PPI_CTTGTA_R1.fastq_cerevisiae_minus.bw", "Cdk9as_3MBPP1r2"), 
  c("SpikeIn_bw/6713_7157_31627_HWV7YBGXX_pombe_Isk1as_DMSO_rep1_ACAGTG_R1.fastq_cerevisiae_plus.bw", "SpikeIn_bw/6713_7157_31627_HWV7YBGXX_pombe_Isk1as_DMSO_rep1_ACAGTG_R1.fastq_cerevisiae_minus.bw", "Isk1as_DMSOr1"),
  c("SpikeIn_bw/6713_7157_31629_HWV7YBGXX_pombe_Isk1as_DMSO_rep2_CAGATC_R1.fastq_cerevisiae_plus.bw", "SpikeIn_bw/6713_7157_31629_HWV7YBGXX_pombe_Isk1as_DMSO_rep2_CAGATC_R1.fastq_cerevisiae_minus.bw", "Isk1as_DMSOr2"),
  c("SpikeIn_bw/6713_7157_31628_HWV7YBGXX_pombe_Isk1as_3MB-PP1_rep1_GCCAAT_R1.fastq_cerevisiae_plus.bw", "SpikeIn_bw/6713_7157_31628_HWV7YBGXX_pombe_Isk1as_3MB-PP1_rep1_GCCAAT_R1.fastq_cerevisiae_minus.bw", "Isk1as_3MBPP1r1"),
  c("SpikeIn_bw/6713_7157_31630_HWV7YBGXX_pombe_Isk1as_3MB-PP1_rep2_ACTTGA_R1.fastq_cerevisiae_plus.bw", "SpikeIn_bw/6713_7157_31630_HWV7YBGXX_pombe_Isk1as_3MB-PP1_rep2_ACTTGA_R1.fastq_cerevisiae_minus.bw", "Isk1as_3MBPP1r2"),
  c("SpikeIn_bw/6713_7157_31623_HWV7YBGXX_pombe_Mcs6as_CDK9as_DMSO_rep1_ATCACG_R1.fastq_cerevisiae_plus.bw", "SpikeIn_bw/6713_7157_31623_HWV7YBGXX_pombe_Mcs6as_CDK9as_DMSO_rep1_ATCACG_R1.fastq_cerevisiae_minus.bw", "Mcs6as_Cdk9as_DMSOr1"),
  c("SpikeIn_bw/6713_7157_31625_HWV7YBGXX_pombe_Mcs6as_CDK9as_DMSO_rep2_TTAGGC_R1.fastq_cerevisiae_plus.bw", "SpikeIn_bw/6713_7157_31625_HWV7YBGXX_pombe_Mcs6as_CDK9as_DMSO_rep2_TTAGGC_R1.fastq_cerevisiae_minus.bw", "Mcs6as_Cdk9as_DMSOr2"),
  c("SpikeIn_bw/6713_7157_31624_HWV7YBGXX_pombe_Mcs6as_CDK9as_3MB-PP1_rep1_CGATGT_R1.fastq_cerevisiae_plus.bw", "SpikeIn_bw/6713_7157_31624_HWV7YBGXX_pombe_Mcs6as_CDK9as_3MB-PP1_rep1_CGATGT_R1.fastq_cerevisiae_minus.bw", "Mcs6as_Cdk9as_3MBPP1r1"),
  c("SpikeIn_bw/6713_7157_31626_HWV7YBGXX_pombe_Mcs6as_CDK9as_3MB-PP1_rep2_TGACCA_R1.fastq_cerevisiae_plus.bw", "SpikeIn_bw/6713_7157_31626_HWV7YBGXX_pombe_Mcs6as_CDK9as_3MB-PP1_rep2_TGACCA_R1.fastq_cerevisiae_minus.bw", "Mcs6as_Cdk9as_3MBPP1r2"))
CombinedspikeIn_wigset = rbind(
  c("SpikeIn_bw/5993_7157_26977_HJ72CBGXX_pombe_WT_COMBINED_DMSO_CAGATC_R1.fastq_cerevisiae_plus.bw", "SpikeIn_bw/5993_7157_26977_HJ72CBGXX_pombe_WT_COMBINED_DMSO_CAGATC_R1.fastq_cerevisiae_minus.bw", "WT_DMSO_combined"),
  c("SpikeIn_bw/5993_7157_26978_HJ72CBGXX_pombe_WT_COMBINED_3MB-PPI_ACTTGA_R1.fastq_cerevisiae_plus.bw", "SpikeIn_bw/5993_7157_26978_HJ72CBGXX_pombe_WT_COMBINED_3MB-PPI_ACTTGA_R1.fastq_cerevisiae_minus.bw", "WT_3MBPP1_combined"),
  c("SpikeIn_bw/5993_7157_26979_HJ72CBGXX_pombe_mcs6as_COMBINED_DMSO_GATCAG_R1.fastq_cerevisiae_plus.bw", "SpikeIn_bw/5993_7157_26979_HJ72CBGXX_pombe_mcs6as_COMBINED_DMSO_GATCAG_R1.fastq_cerevisiae_minus.bw", "mcs6as_DMSO_combined"),
  c("SpikeIn_bw/5993_7157_26980_HJ72CBGXX_pombe_mcs6as_COMBINED_3MB-PPI_TAGCTT_R1.fastq_cerevisiae_plus.bw", "SpikeIn_bw/5993_7157_26980_HJ72CBGXX_pombe_mcs6as_COMBINED_3MB-PPI_TAGCTT_R1.fastq_cerevisiae_minus.bw", "mcs6as_3MBPP1_combined"),
  c("SpikeIn_bw/5993_7157_26981_HJ72CBGXX_pombe_CDK9as_COMBINED_DMSO_GGCTAC_R1.fastq_cerevisiae_plus.bw", "SpikeIn_bw/5993_7157_26981_HJ72CBGXX_pombe_CDK9as_COMBINED_DMSO_GGCTAC_R1.fastq_cerevisiae_minus.bw", "Cdk9as_DMSO_combined"),
  c("SpikeIn_bw/5993_7157_26982_HJ72CBGXX_pombe_CDK9as_COMBINED_3MB-PPI_CTTGTA_R1.fastq_cerevisiae_plus.bw", "SpikeIn_bw/5993_7157_26982_HJ72CBGXX_pombe_CDK9as_COMBINED_3MB-PPI_CTTGTA_R1.fastq_cerevisiae_minus.bw", "Cdk9as_3MBPP1_combined"),
  c("SpikeIn_bw/6713_7157_31629_HWV7YBGXX_pombe_Isk1as_DMSO_COMBINED_CAGATC_R1.fastq_cerevisiae_plus.bw", "SpikeIn_bw/6713_7157_31629_HWV7YBGXX_pombe_Isk1as_DMSO_COMBINED_CAGATC_R1.fastq_cerevisiae_minus.bw", "Isk1as_DMSO_combined"),
  c("SpikeIn_bw/6713_7157_31630_HWV7YBGXX_pombe_Isk1as_3MB-PP1_COMBINED_ACTTGA_R1.fastq_cerevisiae_plus.bw", "SpikeIn_bw/6713_7157_31630_HWV7YBGXX_pombe_Isk1as_3MB-PP1_COMBINED_ACTTGA_R1.fastq_cerevisiae_minus.bw", "Isk1as_3MBPP1_combined"),
  c("SpikeIn_bw/6713_7157_31625_HWV7YBGXX_pombe_Mcs6as_CDK9as_DMSO_COMBINED_TTAGGC_R1.fastq_cerevisiae_plus.bw", "SpikeIn_bw/6713_7157_31625_HWV7YBGXX_pombe_Mcs6as_CDK9as_DMSO_COMBINED_TTAGGC_R1.fastq_cerevisiae_minus.bw", "Mcs6as_Cdk9as_DMSO_combined"),
  c("SpikeIn_bw/6713_7157_31626_HWV7YBGXX_pombe_Mcs6as_CDK9as_3MB-PP1_COMBINED_TGACCA_R1.fastq_cerevisiae_plus.bw", "SpikeIn_bw/6713_7157_31626_HWV7YBGXX_pombe_Mcs6as_CDK9as_3MB-PP1_COMBINED_TGACCA_R1.fastq_cerevisiae_minus.bw", "Mcs6as_Cdk9as_3MBPP1_combined"))



get_spikeIn_NFs = function(wigset){
  result = matrix(nrow = 0, ncol = 3)
  countList = c()
  N = dim(wigset)[1]
  for (i in 1:N){
    wigs = load.wigset(wigset, i)
    reads = total.reads(wigs[[1]], wigs[[2]])
    sample = wigset[i, 3]
    countList = c(countList, reads)
    result = rbind(result, c(sample, reads, reads))
  }
  maxReads = max(countList)
  result[,3] = as.numeric(result[,3])/maxReads
  colnames(result) = c("sample", "Spikereads", "NormFactor")
  return(result)
}
NF_list = get_spikeIn_NFs(spikeIn_wigset)
write.table(NF_list, file = paste(count_outpath, "replicate_spikeNormFactors.txt", sep = ""), quote = F, sep = "\t", row.names = F)
combinedNF_list = get_spikeIn_NFs(CombinedspikeIn_wigset)
write.table(combinedNF_list, file = paste(count_outpath, "combined_spikeNormFactors.txt", sep = ""), quote = F, sep = "\t", row.names = F)

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
  c("5993_7157_26977_HJ72CBGXX_pombe_WT_COMBINED_DMSO_CAGATC_R1.fastq_pombe_plus.bw", "5993_7157_26977_HJ72CBGXX_pombe_WT_COMBINED_DMSO_CAGATC_R1.fastq_pombe_minus.bw", "WT_DMSO_combined"),
  c("5993_7157_26978_HJ72CBGXX_pombe_WT_COMBINED_3MB-PPI_ACTTGA_R1.fastq_pombe_plus.bw", "5993_7157_26978_HJ72CBGXX_pombe_WT_COMBINED_3MB-PPI_ACTTGA_R1.fastq_pombe_minus.bw", "WT_3MBPP1_combined"),
  c("5993_7157_26979_HJ72CBGXX_pombe_mcs6as_COMBINED_DMSO_GATCAG_R1.fastq_pombe_plus.bw", "5993_7157_26979_HJ72CBGXX_pombe_mcs6as_COMBINED_DMSO_GATCAG_R1.fastq_pombe_minus.bw", "mcs6as_DMSO_combined"),
  c("5993_7157_26980_HJ72CBGXX_pombe_mcs6as_COMBINED_3MB-PPI_TAGCTT_R1.fastq_pombe_plus.bw", "5993_7157_26980_HJ72CBGXX_pombe_mcs6as_COMBINED_3MB-PPI_TAGCTT_R1.fastq_pombe_minus.bw", "mcs6as_3MBPP1_combined"),
  c("5993_7157_26981_HJ72CBGXX_pombe_CDK9as_COMBINED_DMSO_GGCTAC_R1.fastq_pombe_plus.bw", "5993_7157_26981_HJ72CBGXX_pombe_CDK9as_COMBINED_DMSO_GGCTAC_R1.fastq_pombe_minus.bw", "Cdk9as_DMSO_combined"),
  c("5993_7157_26982_HJ72CBGXX_pombe_CDK9as_COMBINED_3MB-PPI_CTTGTA_R1.fastq_pombe_plus.bw", "5993_7157_26982_HJ72CBGXX_pombe_CDK9as_COMBINED_3MB-PPI_CTTGTA_R1.fastq_pombe_minus.bw", "Cdk9as_3MBPP1_combined"),
  c("6713_7157_31629_HWV7YBGXX_pombe_Isk1as_DMSO_COMBINED_CAGATC_R1.fastq_pombe_plus.bw", "6713_7157_31629_HWV7YBGXX_pombe_Isk1as_DMSO_COMBINED_CAGATC_R1.fastq_pombe_minus.bw", "Isk1as_DMSO_combined"),
  c("6713_7157_31630_HWV7YBGXX_pombe_Isk1as_3MB-PP1_COMBINED_ACTTGA_R1.fastq_pombe_plus.bw", "6713_7157_31630_HWV7YBGXX_pombe_Isk1as_3MB-PP1_COMBINED_ACTTGA_R1.fastq_pombe_minus.bw", "Isk1as_3MBPP1_combined"),
  c("6713_7157_31625_HWV7YBGXX_pombe_Mcs6as_CDK9as_DMSO_COMBINED_TTAGGC_R1.fastq_pombe_plus.bw", "6713_7157_31625_HWV7YBGXX_pombe_Mcs6as_CDK9as_DMSO_COMBINED_TTAGGC_R1.fastq_pombe_minus.bw", "Mcs6as_Cdk9as_DMSO_combined"),
  c("6713_7157_31626_HWV7YBGXX_pombe_Mcs6as_CDK9as_3MB-PP1_COMBINED_TGACCA_R1.fastq_pombe_plus.bw", "6713_7157_31626_HWV7YBGXX_pombe_Mcs6as_CDK9as_3MB-PP1_COMBINED_TGACCA_R1.fastq_pombe_minus.bw", "Mcs6as_Cdk9as_3MBPP1_combined"))

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
# for testing functions
#test = calc.pause.index(pbw, mbw, mappability = SPmappability, genes = SPGL, name = 'why')


