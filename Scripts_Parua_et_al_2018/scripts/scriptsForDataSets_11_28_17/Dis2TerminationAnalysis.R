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

bed = SPGL_1kbsep_filt
cat("number of genes being used =", length(bed[,1]), "\n")
bed3p = bed[, c(1,3,2,4,5,6)]
Ngenes = length(bed[,1])

NFs = read.table(file = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/countData/combined_spikeNormFactors.txt", head = T)
repNFs = read.table(file = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/countData/replicate_spikeNormFactors.txt", head = T)
wigset = rbind(c("ePROseq_Pombe_Fisher_WT_COMBINED_pombe_plus.bw", "ePROseq_Pombe_Fisher_WT_COMBINED_pombe_minus.bw", "Fisher_WT_combined", NFs[36,'Spikereads']/100000),
               c("ePROseq_Pombe_Fisher_dis2_11_COMBINED_pombe_plus.bw", "ePROseq_Pombe_Fisher_dis2_11_COMBINED_pombe_minus.bw", "Fisher_dis2_11_combined", NFs[37,'Spikereads']/100000),
               c("ePROseq_Pombe_Fisher_dis2Del_COMBINED_pombe_plus.bw", "ePROseq_Pombe_Fisher_dis2Del_COMBINED_pombe_minus.bw", "Fisher_dis2Del_combined", NFs[38,'Spikereads']/100000),
               c("ePROseq_Pombe_Fisher_dis2_T316A_COMBINED_pombe_plus.bw", "ePROseq_Pombe_Fisher_dis2_T316A_COMBINED_pombe_minus.bw", "Fisher_dis2_T316A_combined", NFs[39,'Spikereads']/100000),
               c("ePROseq_Pombe_Fisher_dis2_T316D_COMBINED_pombe_plus.bw", "ePROseq_Pombe_Fisher_dis2_T316D_COMBINED_pombe_minus.bw", "Fisher_dis2_T316D_combined", NFs[40,'Spikereads']/100000))

wigsetReps = rbind(c("ePROseq_Pombe_Fisher_WT_rep1_pombe_plus.bw", "ePROseq_Pombe_Fisher_WT_rep1_pombe_minus.bw", "Fisher_WT_r1", repNFs[71,'Spikereads']/100000),
                   c("ePROseq_Pombe_Fisher_WT_rep2_pombe_plus.bw", "ePROseq_Pombe_Fisher_WT_rep2_pombe_minus.bw", "Fisher_WT_r2", repNFs[72,'Spikereads']/100000),
                   c("ePROseq_Pombe_Fisher_dis2_11_rep1_pombe_plus.bw", "ePROseq_Pombe_Fisher_dis2_11_rep1_pombe_minus.bw", "Fisher_dis2_11_r1", repNFs[73,'Spikereads']/100000),
                   c("ePROseq_Pombe_Fisher_dis2_11_rep2_pombe_plus.bw", "ePROseq_Pombe_Fisher_dis2_11_rep2_pombe_minus.bw", "Fisher_dis2_11_r2", repNFs[74,'Spikereads']/100000),
                   c("ePROseq_Pombe_Fisher_dis2Del_rep1_pombe_plus.bw", "ePROseq_Pombe_Fisher_dis2Del_rep1_pombe_minus.bw", "Fisher_dis2Del_r1", repNFs[75,'Spikereads']/100000),
                   c("ePROseq_Pombe_Fisher_dis2Del_rep2_pombe_plus.bw", "ePROseq_Pombe_Fisher_dis2Del_rep2_pombe_minus.bw", "Fisher_dis2Del_r2", repNFs[76,'Spikereads']/100000),
                   c("ePROseq_Pombe_Fisher_dis2_T316A_rep1_pombe_plus.bw", "ePROseq_Pombe_Fisher_dis2_T316A_rep1_pombe_minus.bw", "Fisher_dis2_T316A_r1", repNFs[77,'Spikereads']/100000),
                   c("ePROseq_Pombe_Fisher_dis2_T316A_rep2_pombe_plus.bw", "ePROseq_Pombe_Fisher_dis2_T316A_rep2_pombe_minus.bw", "Fisher_dis2_T316A_r2", repNFs[78,'Spikereads']/100000),
                   c("ePROseq_Pombe_Fisher_dis2_T316D_rep1_pombe_plus.bw", "ePROseq_Pombe_Fisher_dis2_T316D_rep1_pombe_minus.bw", "Fisher_dis2_T316D_r1", repNFs[79,'Spikereads']/100000),
                   c("ePROseq_Pombe_Fisher_dis2_T316D_rep2_pombe_plus.bw", "ePROseq_Pombe_Fisher_dis2_T316D_rep2_pombe_minus.bw", "Fisher_dis2_T316D_r2", repNFs[80,'Spikereads']/100000))

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
WT_dis2_11_SampleList = c("Fisher_WT_combined", "Fisher_dis2_11_combined")
WT_dis2Del_SampleList = c("Fisher_WT_combined", "Fisher_dis2Del_combined")
WT_dis2_T316A_SampleList = c("Fisher_WT_combined", "Fisher_dis2_T316A_combined")
WT_dis2_T316D_SampleList = c("Fisher_WT_combined", "Fisher_dis2_T316D_combined")
dis2_T316A_dis2_T316D_SampleList = c("Fisher_dis2_T316A_combined", "Fisher_dis2_T316D_combined")
Allsamples = c("Fisher_WT_combined", "Fisher_dis2_11_combined", "Fisher_dis2Del_combined", "Fisher_dis2_T316A_combined", "Fisher_dis2_T316D_combined")

#Analyze comparisons between different strains
# WT vs dis2-11
WT_dis2_11_res = makeboxplots.TI_TEI(postCPS_data = postCPS_ratios, sampleList = WT_dis2_11_SampleList, sampleLabs = c("WT", "dis2-11"), filename = "WT_di2-11_TerminationRatios_boxplot.pdf")
WT_dis2_11_diffTest = diffTest_TI_TEI(postCPS_data = postCPS_ratios, sampleList = WT_dis2_11_SampleList)
print(WT_dis2_11_diffTest)
# WT vs dis2-delete
WT_dis2Del_res = makeboxplots.TI_TEI(postCPS_data = postCPS_ratios, sampleList = WT_dis2Del_SampleList, sampleLabs = c("WT", "dis2Del"), filename = "WT_di2Del_TerminationRatios_boxplot.pdf")
WT_dis2Del_diffTest = diffTest_TI_TEI(postCPS_data = postCPS_ratios, sampleList = WT_dis2Del_SampleList)
print(WT_dis2Del_diffTest)
# WT vs dis2-T316A
WT_dis2_T316A_res = makeboxplots.TI_TEI(postCPS_data = postCPS_ratios, sampleList = WT_dis2_T316A_SampleList, sampleLabs = c("WT", "dis2-T316A"), filename = "WT_dis2_T316A_TerminationRatios_boxplot.pdf")
WT_dis2_T316A_diffTest = diffTest_TI_TEI(postCPS_data = postCPS_ratios, sampleList = WT_dis2_T316A_SampleList)
print(WT_dis2_T316A_diffTest)
# WT vs dis2-T316D
WT_dis2_T316D_res = makeboxplots.TI_TEI(postCPS_data = postCPS_ratios, sampleList = WT_dis2_T316D_SampleList, sampleLabs = c("WT", "dis2-T316D"), filename = "WT_dis2_T316D_TerminationRatios_boxplot.pdf")
WT_dis2_T316D_diffTest = diffTest_TI_TEI(postCPS_data = postCPS_ratios, sampleList = WT_dis2_T316D_SampleList)
print(WT_dis2_T316D_diffTest)
# dis2-T316A vs dis2-T316D
dis2_T316A_dis2_T316D_res = makeboxplots.TI_TEI(postCPS_data = postCPS_ratios, sampleList = dis2_T316A_dis2_T316D_SampleList, sampleLabs = c("dis2-T316A", "dis2-T316D"), filename = "dis2_T316A_dis2_T316D_TerminationRatios_boxplot.pdf")
dis2_T316A_dis2_T316D_diffTest = diffTest_TI_TEI(postCPS_data = postCPS_ratios, sampleList = dis2_T316A_dis2_T316D_SampleList)
print(dis2_T316A_dis2_T316D_diffTest)
# Plot only TEI boxplots for all samples on same plot
All_TEI_res = makeboxplots.TI_TEI(postCPS_data = postCPS_ratios, sampleList = Allsamples, sampleLabs = c("WT", "dis2-11", "dis2Del", "dis2-T316A", "dis2-T316D"), filename = "All_sample_TEIonly_boxplot.pdf", TEIonly = T)

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
###########
# for each comparison:
## resort genes based on TEI of the WT strain to keep sorting for all plots
## prepare heatmap comparisons
## plot 3 heatmaps
##########
# WT vs dis2-11
GL_TEIsort_WT = orderGenesbyTI(genes = bed, CPS_TI_Table = postCPS_ratios, samplename = "Fisher_WT_combined", TIorTEI = "TEI")
#GL_TEIsort_dis2_11 = orderGenesbyTI(genes = bed, CPS_TI_Table = postCPS_ratios, samplename = "Fisher_dis2_11_combined", TIorTEI = "TEI")
WT_dis2_11_heat_CPS_TEIsort = compareHeats(GL_TEIsort_WT[,c(1,3,2,4,5,6)], BWset1 = wigset[1,], BWset2 = wigset[2,], maxLength = 1000, step = 10, BWpath = bwpath)
Plot3Heat(WT_dis2_11_heat_CPS_TEIsort, filename = "WT_vs_dis2-11_FCheat_aroundCPS_WT_TEISort.pdf", distToTSS = 1000)# RAWbreaks = seq(0, 3, length.out = 100), FCbreaks = seq(-2,2, length.out = 100))
# WT vs dis2-delete
#GL_TEIsort_dis2Del = orderGenesbyTI(genes = bed, CPS_TI_Table = postCPS_ratios, samplename = "Fisher_dis2Del_combined", TIorTEI = "TEI")
WT_dis2Del_heat_CPS_TEIsort = compareHeats(GL_TEIsort_WT[,c(1,3,2,4,5,6)], BWset1 = wigset[1,], BWset2 = wigset[3,], maxLength = 1000, step = 10, BWpath = bwpath)
Plot3Heat(WT_dis2Del_heat_CPS_TEIsort, filename = "WT_vs_dis2Del_FCheat_aroundCPS_WT_TEISort.pdf", distToTSS = 1000)# RAWbreaks = seq(0, 3, length.out = 100), FCbreaks = seq(-2,2, length.out = 100))
# WT vs dis2-T316A
#GL_TEIsort_dis2_T316A = orderGenesbyTI(genes = bed, CPS_TI_Table = postCPS_ratios, samplename = "Fisher_dis2_T316A_combined", TIorTEI = "TEI")
WT_dis2_T316A_heat_CPS_TEIsort = compareHeats(GL_TEIsort_WT[,c(1,3,2,4,5,6)], BWset1 = wigset[1,], BWset2 = wigset[4,], maxLength = 1000, step = 10, BWpath = bwpath)
Plot3Heat(WT_dis2_T316A_heat_CPS_TEIsort, filename = "WT_vs_dis2_T316A_FCheat_aroundCPS_WT_TEISort.pdf", distToTSS = 1000)# RAWbreaks = seq(0, 3, length.out = 100), FCbreaks = seq(-2,2, length.out = 100))
# WT vs dis2-T316D
#GL_TEIsort_dis2_T316D = orderGenesbyTI(genes = bed, CPS_TI_Table = postCPS_ratios, samplename = "Fisher_dis2_T316D_combined", TIorTEI = "TEI")
WT_dis2_T316D_heat_CPS_TEIsort = compareHeats(GL_TEIsort_WT[,c(1,3,2,4,5,6)], BWset1 = wigset[1,], BWset2 = wigset[5,], maxLength = 1000, step = 10, BWpath = bwpath)
Plot3Heat(WT_dis2_T316D_heat_CPS_TEIsort, filename = "WT_vs_dis2_T316D_FCheat_aroundCPS_WT_TEISort.pdf", distToTSS = 1000)# RAWbreaks = seq(0, 3, length.out = 100), FCbreaks = seq(-2,2, length.out = 100))
# dis2-T316A vs dis2-T316D (Sorted by dis2_T316A)
GL_TEIsort_dis2_T316A = orderGenesbyTI(genes = bed, CPS_TI_Table = postCPS_ratios, samplename = "Fisher_dis2_T316A_combined", TIorTEI = "TEI")
dis2_T316A_dis2_T316D_heat_CPS_TEIsort = compareHeats(GL_TEIsort_dis2_T316A[,c(1,3,2,4,5,6)], BWset1 = wigset[4,], BWset2 = wigset[5,], maxLength = 1000, step = 10, BWpath = bwpath)
Plot3Heat(dis2_T316A_dis2_T316D_heat_CPS_TEIsort, filename = "dis2_T316A_vs_dis2_T316D_FCheat_aroundCPS_TEISort.pdf", distToTSS = 1000)# RAWbreaks = seq(0, 3, length.out = 100), FCbreaks = seq(-2,2, length.out = 100))





