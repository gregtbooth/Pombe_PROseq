## This script is intended for analysing changes in terminating Pol II between two samples.  
## One Measure will be "termination index", which is the ratio of reads upstream of CPS to Downstream. 

source("/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/scripts/GB_Functions.R")
library(gplots)
bwpath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/alignment/bw/"
infopath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/pombeinfo/"
fig_dir = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/figures/Dis2Analysis/Termination/03-29-18/"
dir.create(fig_dir)
SPmappability = load.bigWig(paste(infopath, "SPombeMappable_SacCerNonOverLap_expanded_ChrFix.bw", sep = ""))
filteredGL = read.table(paste(infopath, "SP_CompleteFilteredGenes.bed", sep = ""))
SPGL_1kbsep = read.table(file = paste(infopath, "SP_PROcapObservedTSS_1000bpSep.bed", sep = ""), stringsAsFactors = F)
SPGL_1kbsep_filt = merge(filteredGL, SPGL_1kbsep, by.x = 4, by.y = 4)[, c(2,3,4,1,5,6)]
ObsTSS_NOoverlap = read.table(file = paste(infopath, "SP_PROcapObservedTSS_noOverlap.bed", sep = ""))
pausedGenes = read.table(file = paste(infopath, "/SP_PausedGene_Filtered_bothReps.bed", sep = "")) #defined  in genome research paper
NotpausedGenes = read.table(file = paste(infopath, "SP_nonPausedGene_Filtered_bothReps.bed", sep = "")) #defined  in genome research paper

bed = SPGL_1kbsep_filt
cat("number of genes being used =", length(bed[,1]), "\n")
bed3p = bed[, c(1,3,2,4,5,6)]
Ngenes = length(bed[,1])

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
wigsetReps = rbind(c("5993_7157_26971_HJ72CBGXX_pombe_WT_rep1_DMSO_ATCACG_R1.fastq_pombe_plus.bw", "5993_7157_26971_HJ72CBGXX_pombe_WT_rep1_DMSO_ATCACG_R1.fastq_pombe_minus.bw", "WT_DMSOr1", RPM_repNFs[1,'Spikereads']/100000),
                   c("5993_7157_26977_HJ72CBGXX_pombe_WT_rep2_DMSO_CAGATC_R1.fastq_pombe_plus.bw", "5993_7157_26977_HJ72CBGXX_pombe_WT_rep2_DMSO_CAGATC_R1.fastq_pombe_minus.bw", "WT_DMSOr2", RPM_repNFs[2,'Spikereads']/100000),
                   c("5993_7157_26972_HJ72CBGXX_pombe_WT_rep1_3MB-PPI_CGATGT_R1.fastq_pombe_plus.bw", "5993_7157_26972_HJ72CBGXX_pombe_WT_rep1_3MB-PPI_CGATGT_R1.fastq_pombe_minus.bw", "WT_3MBPP1r1", RPM_repNFs[3,'Spikereads']/100000),
                   c("5993_7157_26978_HJ72CBGXX_pombe_WT_rep2_3MB-PPI_ACTTGA_R1.fastq_pombe_plus.bw", "5993_7157_26978_HJ72CBGXX_pombe_WT_rep2_3MB-PPI_ACTTGA_R1.fastq_pombe_minus.bw", "WT_3MBPP1r2", RPM_repNFs[4,'Spikereads']/100000),
                   c("7772_7157_43139_HVLYCBGXY_pombe_18C_Cdk9as_DMSO_rep1_ATCACG_R1_pombe_plus.bw", "7772_7157_43139_HVLYCBGXY_pombe_18C_Cdk9as_DMSO_rep1_ATCACG_R1_pombe_minus.bw", "Cdk9as_DMSO_18Cr1", RPM_repNFs[21,'Spikereads']/100000),
                   c("7772_7157_43149_HVLYCBGXY_pombe_18C_Cdk9as_DMSO_rep2_GGCTAC_R1_pombe_plus.bw", "7772_7157_43149_HVLYCBGXY_pombe_18C_Cdk9as_DMSO_rep2_GGCTAC_R1_pombe_minus.bw", "Cdk9as_DMSO_18Cr2", RPM_repNFs[22,'Spikereads']/100000),
                   c("7772_7157_43140_HVLYCBGXY_pombe_18C_Cdk9as_3MBPP1_rep1_CGATGT_R1_pombe_plus.bw", "7772_7157_43140_HVLYCBGXY_pombe_18C_Cdk9as_3MBPP1_rep1_CGATGT_R1_pombe_minus.bw", "Cdk9as_3MBPP1_18Cr1", RPM_repNFs[23,'Spikereads']/100000),
                   c("7772_7157_43150_HVLYCBGXY_pombe_18C_Cdk9as_3MBPP1_rep2_CTTGTA_R1_pombe_plus.bw", "7772_7157_43150_HVLYCBGXY_pombe_18C_Cdk9as_3MBPP1_rep2_CTTGTA_R1_pombe_minus.bw", "Cdk9as_3MBPP1_18Cr2", RPM_repNFs[24,'Spikereads']/100000),
                   c("7772_7157_43141_HVLYCBGXY_pombe_18C_dis2ts_DMSO_rep1_TTAGGC_R1_pombe_plus.bw", "7772_7157_43141_HVLYCBGXY_pombe_18C_dis2ts_DMSO_rep1_TTAGGC_R1_pombe_minus.bw", "dis2ts_DMSO_18Cr1", RPM_repNFs[25,'Spikereads']/100000),
                   c("7772_7157_43151_HVLYCBGXY_pombe_18C_dis2ts_DMSO_rep2_ATTCCT_R1_pombe_plus.bw", "7772_7157_43151_HVLYCBGXY_pombe_18C_dis2ts_DMSO_rep2_ATTCCT_R1_pombe_minus.bw", "dis2ts_DMSO_18Cr2", RPM_repNFs[26,'Spikereads']/100000),
                   c("7772_7157_43142_HVLYCBGXY_pombe_18C_dis2ts_3MBPP1_rep1_TGACCA_R1_pombe_plus.bw", "7772_7157_43142_HVLYCBGXY_pombe_18C_dis2ts_3MBPP1_rep1_TGACCA_R1_pombe_minus.bw", "dis2ts_3MBPP1_18Cr1", RPM_repNFs[27,'Spikereads']/100000),
                   c("7772_7157_43152_HVLYCBGXY_pombe_18C_dis2ts_3MBPP1_rep2_CAAAAG_R1_pombe_plus.bw", "7772_7157_43152_HVLYCBGXY_pombe_18C_dis2ts_3MBPP1_rep2_CAAAAG_R1_pombe_minus.bw", "dis2ts_3MBPP1_18Cr2", RPM_repNFs[28,'Spikereads']/100000),
                   c("7772_7157_43143_HVLYCBGXY_pombe_30C_dis2ts_DMSO_rep1_ACAGTG_R1_pombe_plus.bw", "7772_7157_43143_HVLYCBGXY_pombe_30C_dis2ts_DMSO_rep1_ACAGTG_R1_pombe_minus.bw", "dis2ts_DMSO_30Cr1", RPM_repNFs[29,'Spikereads']/100000),
                   c("7772_7157_43153_HVLYCBGXY_pombe_30C_dis2ts_DMSO_rep2_CAACTA_R1_pombe_plus.bw", "7772_7157_43153_HVLYCBGXY_pombe_30C_dis2ts_DMSO_rep2_CAACTA_R1_pombe_minus.bw", "dis2ts_DMSO_30Cr2", RPM_repNFs[30,'Spikereads']/100000),
                   c("7772_7157_43144_HVLYCBGXY_pombe_30C_dis2ts_3MBPP1_rep1_GCCAAT_R1_pombe_plus.bw", "7772_7157_43144_HVLYCBGXY_pombe_30C_dis2ts_3MBPP1_rep1_GCCAAT_R1_pombe_minus.bw", "dis2ts_3MBPP1_30Cr1", RPM_repNFs[31,'Spikereads']/100000),
                   c("7772_7157_43154_HVLYCBGXY_pombe_30C_dis2ts_3MBPP1_rep2_CACCGG_R1_pombe_plus.bw", "7772_7157_43154_HVLYCBGXY_pombe_30C_dis2ts_3MBPP1_rep2_CACCGG_R1_pombe_minus.bw", "dis2ts_3MBPP1_30Cr2", RPM_repNFs[32,'Spikereads']/100000),
                   c("7772_7157_43145_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_DMSO_rep1_CAGATC_R1_pombe_plus.bw", "7772_7157_43145_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_DMSO_rep1_CAGATC_R1_pombe_minus.bw", "CDK9as_dis2ts_DMSO_18Cr1", RPM_repNFs[33,'Spikereads']/100000),
                   c("7772_7157_43155_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_DMSO_rep2_CACGAT_R1_pombe_plus.bw", "7772_7157_43155_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_DMSO_rep2_CACGAT_R1_pombe_minus.bw", "CDK9as_dis2ts_DMSO_18Cr2", RPM_repNFs[34,'Spikereads']/100000),
                   c("7772_7157_43146_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_3MBPP1_rep1_ACTTGA_R1_pombe_plus.bw", "7772_7157_43146_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_3MBPP1_rep1_ACTTGA_R1_pombe_minus.bw", "CDK9as_dis2ts_3MBPP1_18Cr1", RPM_repNFs[35,'Spikereads']/100000),
                   c("7772_7157_43156_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_3MBPP1_rep2_CACTCA_R1_pombe_plus.bw", "7772_7157_43156_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_3MBPP1_rep2_CACTCA_R1_pombe_minus.bw", "CDK9as_dis2ts_3MBPP1_18Cr2", RPM_repNFs[36,'Spikereads']/100000),
                   c("7772_7157_43147_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_DMSO_rep1_GATCAG_R1_pombe_plus.bw", "7772_7157_43147_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_DMSO_rep1_GATCAG_R1_pombe_minus.bw", "CDK9as_dis2ts_DMSO_30Cr1", RPM_repNFs[37,'Spikereads']/100000),
                   c("7772_7157_43157_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_DMSO_rep2_CAGGCG_R1_pombe_plus.bw", "7772_7157_43157_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_DMSO_rep2_CAGGCG_R1_pombe_minus.bw", "CDK9as_dis2ts_DMSO_30Cr2", RPM_repNFs[38,'Spikereads']/100000),
                   c("7772_7157_43148_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_3MBPP1_rep1_TAGCTT_R1_pombe_plus.bw", "7772_7157_43148_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_3MBPP1_rep1_TAGCTT_R1_pombe_minus.bw", "CDK9as_dis2ts_3MBPP1_30Cr1", RPM_repNFs[39,'Spikereads']/100000),
                   c("7772_7157_43158_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_3MBPP1_rep2_CATGGC_R1_pombe_plus.bw", "7772_7157_43158_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_3MBPP1_rep2_CATGGC_R1_pombe_minus.bw", "CDK9as_dis2ts_3MBPP1_30Cr2", RPM_repNFs[40,'Spikereads']/100000))



#####################################################################################
# Function designed to identify genes whose Pol II density profiles are influenced by upstream genes
# Purpose:  Calculate ratios of upstream region (boundaries, all positive values, chosen relative to center) and downstream region (relative to TSS or CPS)
# Function is to be applied to rows of a bed file containg gene information
# requires bigWig package
Downstream.Ratio = function(x, pbw, mbw, map, center = "CPS", UpBound1 = 500, UpBound2 = 0, 
                          DownBound1 = 0, DownBound2 = 500){ # note X is the row (i.e. gene) to apply this function to
  pbw = load.bigWig(paste(bwpath, pbw, sep = ""))
  mbw = load.bigWig(paste(bwpath, mbw, sep = ""))
  if (!(center %in% c("CPS" , "TSS"))){
    print("center must be set to either 'CPS' or 'TSS'.\n")
  }
  if (center == "CPS"){
    cp = 3
    cm = 2}
  else{
    cp = 2
    cm = 3}
  chro = x[1]
  if (x[6] == "+"){
    bw = pbw
    upstart = as.numeric(x[cp]) - UpBound1
    upend = as.numeric(x[cp]) - UpBound2
    downstart = as.numeric(x[cp]) + DownBound1
    downend = as.numeric(x[cp]) + DownBound2
    upmap = sum(query.bigWig(map, chro, upstart, upend)[,3])
    upQuery = query.bigWig(bw, chro, upstart, upend)
    upDens = sum(upQuery[,3])/upmap
    downmap = sum(query.bigWig(map, chro, downstart, downend)[,3])
    downQuery = query.bigWig(bw, chro, downstart, downend)
    downDens = sum(downQuery[,3])/downmap
    downupRatio = (downDens/upDens)
  }
  else{
    bw = mbw
    upstart = as.numeric(x[cm]) + UpBound2
    upend = as.numeric(x[cm]) + UpBound1
    downstart = as.numeric(x[cm]) - DownBound2
    downend = as.numeric(x[cm]) - DownBound1
    upmap = sum(query.bigWig(map, chro, upstart, upend)[,3])
    upQuery = query.bigWig(bw, chro, upstart, upend)
    upDens = sum(upQuery[,3])/upmap
    downmap = sum(query.bigWig(map, chro, downstart, downend)[,3])
    downQuery = query.bigWig(bw, chro, downstart, downend)
    downDens = sum(downQuery[,3])/downmap
    downupRatio = (downDens/upDens)
  }
  return(downupRatio)
}
### example of use: 
#bed$upDownRatio = apply(bed, 1, Upstream.Ratio, pbw = pbw, mbw = mbw, map = mappability)
#####################################################################################
# function below applies the Upstream.Ratio function to all samples in a "wigset" 
# It generates a meta table with two columns for each sample.  
## Termination Index = (-500 to CPS)/ (CPS to +500)
## Termination Elongation Index = (-500 to CPS)/ (CPS+250 to +750)
CPS_TI_table = function(bed, wigset){
  res = data.frame(bed[,"V4"])
  for (i in 1:dim(wigset)[1]){
    cat("calculating Termination index for sample", wigset[i,3], "...\n")
    TI = apply(bed, 1, Downstream.Ratio, pbw = wigset[i,1], mbw = wigset[i,2], map = SPmappability, 
               center = "CPS", UpBound1 = 500, UpBound2 = 0, DownBound1 = 0, DownBound2 = 500)
    cat("calculating Termination elongation index for sample", wigset[i,3], "...\n")
    TEI = apply(bed, 1, Downstream.Ratio, pbw = wigset[i,1], mbw = wigset[i,2], map = SPmappability,
                center = "CPS", UpBound1 = 500, UpBound2 = 0, DownBound1 = 250, DownBound2 = 750)
    df = cbind(TI, TEI)
    colnames(df) <- c(paste(wigset[i,3], "_TI", sep = ""), paste(wigset[i,3], "_TEI", sep = ""))
    res = cbind(res, df)
  }
return(res)
  }
  
# Function for comparing boxplots of TI and TEI between samples. 
# it will draw data from the table produced as output from CPS_TI_table.
## must provide a list of sample names (as in wigset) for samples to compare.
library(ggplot2)
makeboxplots.TI_TEI <- function(postCPS_data, sampleList, sampleLabs = c("sample1", "sample2"), filename = "PI_quantile_PRdensity_boxplot.pdf", TEIonly=FALSE){
  reformatRes = matrix(ncol = 3, nrow =0)
  for(i in 1:length(sampleList)){
    TIindx = paste(sampleList[i], "_TI", sep = "")
    TEIindx =paste(sampleList[i], "_TEI", sep = "")
    TIdat = cbind(postCPS_data[TIindx], 1, i)
    colnames(TIdat) = c("ratio", "TIorTEI", 'sample')
    TEIdat = cbind(postCPS_data[TEIindx], 2, i)
    colnames(TEIdat) = c("ratio", "TIorTEI", 'sample')
    if (TEIonly){reformatRes = rbind(reformatRes, TEIdat)}
    else {reformatRes = rbind(reformatRes, TIdat, TEIdat)}
  }
  if(TEIonly){p <- ggplot(data = reformatRes, aes(x = factor(sample, labels = sampleLabs), y = log10(ratio), 
                                                  fill = factor(TIorTEI, labels = c("Term. Elong. Idx.")))) + 
    geom_boxplot() + 
    xlab("sample") + 
    ylab("log10 (postCPS / PreCPS)") +
    labs(fill = "") +# gets rid of legend name
    scale_fill_brewer(palette = "Accent")}
  else{p <- ggplot(data = reformatRes, aes(x = factor(sample, labels = sampleLabs), y = log10(ratio), 
                                           fill = factor(TIorTEI, labels = c("Term. Idx.", "Term. Elong. Idx.")))) + 
    geom_boxplot() + 
    xlab("sample") + 
    ylab("log10 (postCPS / PreCPS)") +
    labs(fill = "") +# gets rid of legend name
    scale_fill_brewer(palette = "Accent")}
  ggsave(filename = paste(fig_dir, filename, sep = ""), plot = p, width = 4, height = 4)
  # print box plot parameters to outputfile
  x = ggplot_build(p)$data
  box_info = cbind(x[[1]]$group, x[[1]]$ymin_final, x[[1]]$ymin, x[[1]]$lower, 
                   x[[1]]$middle, x[[1]]$upper, x[[1]]$ymax, x[[1]]$ymax_final)
  box_info <- box_info[order(box_info[,1]),]
  colnames(box_info) = c("boxNum", "ymin", "lowerWhisker", "lowBoxEdge_25percentile",
                         "Median", "upBoxEdge_75percentile", "upperWhisker", "ymax")
  sink(file = paste(fig_dir, filename, "_info", ".txt", sep = ""))
  cat(paste("box info for ", filename, "\n", sep=""))
  print(box_info)
  sink()
  return(reformatRes) 
}

is.finite.data.frame <- function(obj){
  apply(obj, 1, FUN = function(x) all(is.finite(x)))
}
# function below returns a list of t-test results for comparison TI and TEI between two samples
diffTest_TI_TEI <- function(postCPS_data, sampleList, method = "ttest"){
  dat = postCPS_data[complete.cases(postCPS_data),]
  TI_s1_indx = paste(sampleList[1], "_TI", sep = "")
  TEI_s1indx =paste(sampleList[1], "_TEI", sep = "")
  TI_s2_indx = paste(sampleList[2], "_TI", sep = "")
  TEI_s2indx =paste(sampleList[2], "_TEI", sep = "")
  comp1 = cbind(dat[TI_s1_indx],dat[TI_s2_indx])
  comp1_fin = comp1[is.finite.data.frame(comp1),]
  print(dim(comp1_fin))
  comp2 = cbind(dat[TEI_s1indx],dat[TEI_s2indx])
  comp2_fin = comp2[is.finite.data.frame(comp2),]
  TItest = t.test(comp1_fin[,1], comp1_fin[,2])
  TEItest = t.test(comp2_fin[,1], comp2_fin[,2])
  result = list()
  result[["TI_comparison"]] = TItest
  result[["TEI_comparison"]] = TEItest
  return(result)
}

postCPS_ratios = CPS_TI_table(bed, wigset)
### prepare boxplots of different comparisons simply by using different lists of samples. 
WT_Dis2ts_SampleList = c("WT_DMSO_combined", "dis2ts_DMSO_30C_combined")
Cdk9_Cdk9Dis2ts_SampleList = c("Cdk9as_DMSO_18C_combined", "CDK9as_dis2ts_DMSO_18C_combined")

WT_Dis2ts_res = makeboxplots.TI_TEI(postCPS_data = postCPS_ratios, sampleList = WT_Dis2ts_SampleList, sampleLabs = c("WT", "Dis2-11"), filename = "WT_Dis2ts_30C_DMSO_TerminationRatios_boxplot.pdf")
WT_Dis2ts_res_TEI= makeboxplots.TI_TEI(postCPS_data = postCPS_ratios, sampleList = WT_Dis2ts_SampleList, sampleLabs = c("WT", "Dis2-11"), filename = "WT_Dis2ts_30C_DMSO_TerminationRatios_boxplot_TEIonly.pdf", TEIonly = T)
WT_Dis2ts_diffTest = diffTest_TI_TEI(postCPS_data = postCPS_ratios, sampleList = WT_Dis2ts_SampleList)
print(WT_Dis2ts_diffTest)
Cdk9_Cdk9Dis2ts_res = makeboxplots.TI_TEI(postCPS_data = postCPS_ratios, sampleList = Cdk9_Cdk9Dis2ts_SampleList, sampleLabs = c("CDK9as", "Cdk9as_Dis2-11"), filename = "Cdk9_Cdk9Dis2ts_18C_DMSO_TerminationRatios_boxplot.pdf")
Cdk9_Cdk9Dis2ts_res_TEI = makeboxplots.TI_TEI(postCPS_data = postCPS_ratios, sampleList = Cdk9_Cdk9Dis2ts_SampleList, sampleLabs = c("CDK9as", "Cdk9as_Dis2-11"), filename = "Cdk9_Cdk9Dis2ts_18C_DMSO_TerminationRatios_boxplot_TEIonly.pdf", TEIonly = T)
Cdk9_Cdk9Dis2ts_diffTest = diffTest_TI_TEI(postCPS_data = postCPS_ratios, sampleList = Cdk9_Cdk9Dis2ts_SampleList)
print(Cdk9_Cdk9Dis2ts_diffTest)
#####################################################################################
## For comparisons that can be normalized, prepare FC heatmops from -250 to +1000 around CPS
# function for sorting genes by Termination Index or Term Elongation Index
orderGenesbyTI = function(genes, CPS_TI_Table, samplename = "CDK9as_dis2ts_DMSO_18C_combined", TIorTEI = "TEI"){
  if (TIorTEI == "TEI"){indx = paste(samplename, "_TEI", sep = "")}
  else if (TIorTEI == "TI"){indx = paste(samplename, "_TI", sep = "")}
  else {print ("TIorTEI must be set as 'TI' or 'TEI'\n")}
  GL = cbind(genes, CPS_TI_Table[indx])
  res = GL[order(GL[indx]),]
  return(res)
}
# sort genes based on TEI of Cdk9Dis2_18C_DMSO
GL_TEIsort_Cdk9Dis2_18C_DMSO = orderGenesbyTI(genes = bed, CPS_TI_Table = postCPS_ratios, samplename = "CDK9as_dis2ts_DMSO_18C_combined", TIorTEI = "TEI")
# prepare fold change heat data 
Cdk9_Cdk9Dis2_18C_DMSO_heat_CPS_TEIsort = compareHeats(GL_TEIsort_Cdk9Dis2_18C_DMSO[,c(1,3,2,4,5,6)], BWset1 = wigset[3,], BWset2 = wigset[9,], maxLength = 1000, step = 10, BWpath = bwpath)
# plot the fold change heat maps
Plot3Heat(Cdk9_Cdk9Dis2_18C_DMSO_heat_CPS_TEIsort, filename = "Cdk9_Cdk9Dis2_18C_DMSO_FCheat_aroundCPS_TEISort.pdf", distToTSS = 1000)# RAWbreaks = seq(0, 3, length.out = 100), FCbreaks = seq(-2,2, length.out = 100))

#######################################################################################
# get box plot parameters, 

WT_Dis2ts_box = ggplot(data = WT_Dis2ts_res, aes(x = factor(sample), y = log10(ratio), 
                                                fill = factor(TIorTEI))) + 
                                                geom_boxplot()
x = ggplot_build(WT_Dis2ts_box)$data
y = cbind(x[[1]]$group, x[[1]]$ymin_final, x[[1]]$ymin, x[[1]]$lower, x[[1]]$middle, x[[1]]$upper, x[[1]]$ymax, x[[1]]$ymax_final)
colnames(y) = c("boxNum", "ymin", "lowerWhisker", "lowBoxEdge_25percentile", "Median", "upBoxEdge_75percentile", "upperWhisker", "ymax")




