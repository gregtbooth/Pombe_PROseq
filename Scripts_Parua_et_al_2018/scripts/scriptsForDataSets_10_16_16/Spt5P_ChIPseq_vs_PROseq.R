source("/Volumes/SEAGATE_EXP/Fisher_collaboration/scripts/GB_Functions.R")
library(gplots)
library(ggplot2)
library(lattice)
bwpath = "/Volumes/SEAGATE_EXP/Fisher_collaboration/alignment/bw/"
chip_bwpath = "/Volumes/SEAGATE_EXP/GEO_datasets/Fisher_data/ChIPseq/chrNumeric/"
infopath = "/Volumes/SEAGATE_EXP/Fisher_collaboration/pombeinfo/"
fig_dir = "/Volumes/SEAGATE_EXP/Fisher_collaboration/analysis/figures/Pabitra/PROseqvsSpt5chipSeq/06-21-17/"
dir.create(fig_dir)
SPmappability = load.bigWig(paste(infopath, "SPombeMappable_SacCerNonOverLap_expanded_ChrFix.bw", sep = ""))
ObsTSS_NOoverlap = read.table(file = paste(infopath, "SP_PROcapObservedTSS_noOverlap.bed", sep = ""))
filteredGL = read.table(paste(infopath, "SP_CompleteFilteredGenes.bed", sep = ""))
SPGL_1kbsep = read.table(file = paste(infopath, "SP_PROcapObservedTSS_1000bpSep.bed", sep = ""), stringsAsFactors = F)
SPGL_1kbsep_filt = merge(filteredGL, SPGL_1kbsep, by.x = 4, by.y = 4)[, c(2,3,4,1,5,6)]

###################################################################################################################
wigset = rbind(c("5993_7157_26977_HJ72CBGXX_pombe_WT_COMBINED_DMSO_CAGATC_R1.fastq_pombe_plus.bw", "5993_7157_26977_HJ72CBGXX_pombe_WT_COMBINED_DMSO_CAGATC_R1.fastq_pombe_minus.bw", "WT_DMSO_combined", "PROseq"),
               c("5993_7157_26978_HJ72CBGXX_pombe_WT_COMBINED_3MB-PPI_ACTTGA_R1.fastq_pombe_plus.bw", "5993_7157_26978_HJ72CBGXX_pombe_WT_COMBINED_3MB-PPI_ACTTGA_R1.fastq_pombe_minus.bw", "WT_3MBPP1_combined", "PROseq"),
               c("7772_7157_43149_HVLYCBGXY_pombe_18C_Cdk9as_DMSO_COMBINED_GGCTAC_R1_pombe_plus.bw", "7772_7157_43149_HVLYCBGXY_pombe_18C_Cdk9as_DMSO_COMBINED_GGCTAC_R1_pombe_minus.bw", "Cdk9as_DMSO_18C_combined", "PROseq"),
               c("7772_7157_43150_HVLYCBGXY_pombe_18C_Cdk9as_3MBPP1_COMBINED_CTTGTA_R1_pombe_plus.bw", "7772_7157_43150_HVLYCBGXY_pombe_18C_Cdk9as_3MBPP1_COMBINED_CTTGTA_R1_pombe_minus.bw", "Cdk9as_3MBPP1_18C_combined", "PROseq"),
               c("7772_7157_43151_HVLYCBGXY_pombe_18C_dis2ts_DMSO_COMBINED_ATTCCT_R1_pombe_plus.bw", "7772_7157_43151_HVLYCBGXY_pombe_18C_dis2ts_DMSO_COMBINED_ATTCCT_R1_pombe_minus.bw", "dis2ts_DMSO_18C_combined", "PROseq"),
               c("7772_7157_43152_HVLYCBGXY_pombe_18C_dis2ts_3MBPP1_COMBINED_CAAAAG_R1_pombe_plus.bw", "7772_7157_43152_HVLYCBGXY_pombe_18C_dis2ts_3MBPP1_COMBINED_CAAAAG_R1_pombe_minus.bw", "dis2ts_3MBPP1_18C_combined", "PROseq"),
               c("7772_7157_43153_HVLYCBGXY_pombe_30C_dis2ts_DMSO_COMBINED_CAACTA_R1_pombe_plus.bw", "7772_7157_43153_HVLYCBGXY_pombe_30C_dis2ts_DMSO_COMBINED_CAACTA_R1_pombe_minus.bw", "dis2ts_DMSO_30C_combined", "PROseq"),
               c("7772_7157_43154_HVLYCBGXY_pombe_30C_dis2ts_3MBPP1_COMBINED_CACCGG_R1_pombe_plus.bw", "7772_7157_43154_HVLYCBGXY_pombe_30C_dis2ts_3MBPP1_COMBINED_CACCGG_R1_pombe_minus.bw", "dis2ts_3MBPP1_30C_combined", "PROseq"),
               c("7772_7157_43155_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_DMSO_COMBINED_CACGAT_R1_pombe_plus.bw", "7772_7157_43155_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_DMSO_COMBINED_CACGAT_R1_pombe_minus.bw", "CDK9as_dis2ts_DMSO_18C_combined", "PROseq"),
               c("7772_7157_43156_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_3MBPP1_COMBINED_CACTCA_R1_pombe_plus.bw", "7772_7157_43156_HVLYCBGXY_pombe_18C_Cdk9as_dis2ts_3MBPP1_COMBINED_CACTCA_R1_pombe_minus.bw", "CDK9as_dis2ts_3MBPP1_18C_combined", "PROseq"),
               c("7772_7157_43157_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_DMSO_COMBINED_CAGGCG_R1_pombe_plus.bw", "7772_7157_43157_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_DMSO_COMBINED_CAGGCG_R1_pombe_minus.bw", "CDK9as_dis2ts_DMSO_30C_combined", "PROseq"),
               c("7772_7157_43158_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_3MBPP1_COMBINED_CATGGC_R1_pombe_plus.bw", "7772_7157_43158_HVLYCBGXY_pombe_30C_Cdk9as_dis2ts_3MBPP1_COMBINED_CATGGC_R1_pombe_minus.bw", "CDK9as_dis2ts_3MBPP1_30C_combined", "PROseq"),
               c("pSpt5-AFixChr.bw", "pSpt5-AFixChr.bw", "pSpt5_r1", "ChIPseq"),
               c("pSpt5-BFixChr.bw", "pSpt5-BFixChr.bw", "pSpt5_r2", "ChIPseq"),
               c("Spt5-V5-AFixChr.bw", "Spt5-V5-AFixChr.bw", "Spt5-V5_r1", "ChIPseq"),
               c("Spt5-V5-BFixChr.bw", "Spt5-V5-BFixChr.bw", "Spt5-V5_r2", "ChIPseq"),
               c("Spt5-Myc-AFixChr.bw", "Spt5-Myc-AFixChr.bw", "Spt5-Myc_r1", "ChIPseq"),
               c("Spt5-Myc-BFixChr.bw", "Spt5-Myc-BFixChr.bw", "Spt5-Myc_r2", "ChIPseq"))

###################################################################################################################
load.wigset <- function(wigset, wigset_row) {
  expmt = wigset[wigset_row, 4]
  if(expmt == "ChIPseq"){basepath = chip_bwpath}
  else if(expmt == "PROseq"){basepath = bwpath}
  file = wigset[wigset_row, 1]
  wig.p = NULL
  if (file != "")
    wig.p = load.bigWig(paste(basepath, file, sep=''))
  file = wigset[wigset_row, 2]
  wig.m = NULL
  if (file != "")
    wig.m = load.bigWig(paste(basepath, file, sep=''))
  return(list(wig.p, wig.m, wigset[wigset_row, 3], wigset[wigset_row, 4]))
}
## formatting meta plots for lattice (not scaled Metaplots)
## modified to work with PROseq or ChIPseq samples (needs to be specified in the wigset row)
## only gives fraction of max data & samples must be deduced from the sample number (corresponds to wigset row #)
FormatMetaData_lattice_FracMax <- function(wigset, bed, halfWindow, step, RPM = F){
  N = dim(wigset)[1]
  df = data.frame(matrix(ncol = 5, nrow = 0), stringsAsFactors = F)
  for (i in 1:N){
    sample = i
    if (i == ceiling(N/2)){cat("* 50% complete ... \n")}
    cat("* loading", i, "\n")
    wigs = load.wigset(wigset, i) #this should give all the bigwig files
    cat("* generating and combining all sample meta-plot data ...\n")
    if (wigs[[4]] == "PROseq"){
      meta= meta.subsample(bed, bigWig.plus = wigs[[1]], bigWig.minus = wigs[[2]], 
                          step = step, at.TSS = T, halfWindow = halfWindow, do.sum = T)
      tc = total.reads(wigs[[1]], wigs[[2]])
      nf = 1/(tc/1000000)
    }
    else if (wigs[[4]] == "ChIPseq"){
      meta= meta.subsample.ChIPchip(bed, bigWig.plus = wigs[[1]], bigWig.minus = wigs[[2]], 
                           step = step, at.TSS = T, halfWindow = halfWindow, do.sum = T)
      tc = (total.reads(wigs[[1]], wigs[[2]]) / 2) # need to devide by two because only have one bw for chIP seq samples 
      nf = 1/(tc/1000000)
      }
    xAxis = seq(from = -halfWindow, to = halfWindow, by = step)
    #if(RPM){
    #  MetaNorm = meta.normalize(meta, nf)
    #  sample_df = cbind(xAxis, MetaNorm[[4]], MetaNorm[[3]], MetaNorm[[2]], sample)
    #}
    #else{
      FracMax_mean = meta[[4]]/(max(meta[[4]]))
      FracMax_uQ = meta[[3]]/(max(meta[[4]]))
      FracMax_lQ = meta[[2]]/(max(meta[[4]]))
      sample_df = cbind(xAxis, FracMax_mean, FracMax_uQ, FracMax_lQ, sample)
     # }
    df = rbind.data.frame(df, sample_df)
  }
  colnames(df) = c("x", "mean", "lower", "upper", "sample")
  return(df)
}

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
                               xlab = "Distance to TSS (bp)",  ylab = "Median PRO-seq intensity",
                               labs = c("WT PROseq", "Spt5 ChIPseq")){
  pdf(paste(fig_dir, filename), width = 5, height = 5)
  result <- xyplot(mean ~ x, data = df,
                   group = factor(sample, labels = labs), 
                   scales = list(tck=c(1,0),alternating = c(1,1),x=list(relation='free',axs='i'),
                                 y=list(relation='free',axs='i')),
                   key=list(corner=c(0.98,0.85),padding.text=3,
                            text=list(labs,col=c("black", "black"), cex=0.8, font=2),
                            #text=list(c("30 deg C", "18 deg C"),col=c("black", "black"), cex=0.8, font=2),
                            rectangles=list(col=c(rgb(0,0,0.5), rgb(0.5,0,0), rgb(0,0.5,0)), size=2.7, height=0.8, border='white')),
                   type = 'l', 
                   ylim = ylim,
                   col = c(rgb(0.5,0,0,0.9), rgb(0,0,0.5,0.9), rgb(0,0.5,0,0.9)),  ## just add or take away colors to match number of samples in each panel
                   fill = c(rgb(0.5,0,0,0.3), rgb(0,0,0.5,0.3), rgb(0,0.5,0,0.3)),
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


######################
# plots around TSS
#####################
metaDataTSS = FormatMetaData_lattice_FracMax(wigset, bed = SPGL_1kbsep_filt, halfWindow = 1000, step = 10)

# WT (DMSO) PROseq vs Spt5-V5A ChipSeq
WTDMSOPROseq_Vs_Spt5V5A = (metaDataTSS$sample == 1 |  metaDataTSS$sample == 15)
lattice_meta.proSeq(filename = "WT_DMSO_PROseq_vs_Spt5V5A_ChIPseq_TSS_Meta.pdf", 
                          df = metaDataTSS[WTDMSOPROseq_Vs_Spt5V5A,], ylim = c(0, 1.2), main = NULL, 
                          xlab = "Distance to TSS (bp)", ylab = "Fraction of Max Signal", labs = c("WT PROseq", "Spt5 ChIPseq"))
# WT (DMSO) PROseq vs Spt5-V5B ChipSeq
WTDMSOPROseq_Vs_Spt5V5B = (metaDataTSS$sample == 1 |  metaDataTSS$sample == 16)
lattice_meta.proSeq(filename = "WT_DMSO_PROseq_vs_Spt5V5B_ChIPseq_TSS_Meta.pdf", 
                    df = metaDataTSS[WTDMSOPROseq_Vs_Spt5V5B,], ylim = c(0, 1.2), main = NULL, 
                    xlab = "Distance to TSS (bp)", ylab = "Fraction of Max Signal", labs = c("WT PROseq", "Spt5 ChIPseq"))
# WT (DMSO) PROseq vs pSpt5 ChipSeq
WTDMSOPROseq_Vs_pSpt5A = (metaDataTSS$sample == 1 |  metaDataTSS$sample == 13)
lattice_meta.proSeq(filename = "WT_DMSO_PROseq_vs_pSpt5A_ChIPseq_TSS_Meta.pdf", 
                    df = metaDataTSS[WTDMSOPROseq_Vs_pSpt5A,], ylim = c(0, 1.2), main = NULL, 
                    xlab = "Distance to TSS (bp)", ylab = "Fraction of Max Signal", labs = c("WT PROseq", "pSpt5 ChIPseq"))
# WT (DMSO) PROseq vs pSpt5B ChipSeq
WTDMSOPROseq_Vs_pSpt5B = (metaDataTSS$sample == 1 |  metaDataTSS$sample == 14)
lattice_meta.proSeq(filename = "WT_DMSO_PROseq_vs_pSpt5B_ChIPseq_TSS_Meta.pdf", 
                    df = metaDataTSS[WTDMSOPROseq_Vs_pSpt5B,], ylim = c(0, 1.2), main = NULL, 
                    xlab = "Distance to TSS (bp)", ylab = "Fraction of Max Signal", labs = c("WT PROseq", "pSpt5 ChIPseq"))

######################
# plots around CPS
#####################
metaDataCPS = FormatMetaData_lattice_FracMax(wigset, bed = filteredGL[,c(1,3,2,4,5,6)], halfWindow = 1000, step = 10)
#metaDataCPS_RPM = FormatMetaData_lattice_FracMax(wigset, bed = SPGL_1kbsep_filt[,c(1,3,2,4,5,6)], halfWindow = 1000, step = 10, RPM = T)

# WT (DMSO) PROseq vs Spt5-V5A ChipSeq
WTDMSOPROseq_Vs_Spt5V5A_cps = (metaDataCPS$sample == 1 |  metaDataCPS$sample == 15)
lattice_meta.proSeq(filename = "WT_DMSO_PROseq_vs_Spt5V5A_ChIPseq_CPS_Meta_AllFiltgenes.pdf", 
                    df = metaDataCPS[WTDMSOPROseq_Vs_Spt5V5A_cps,], ylim = c(0, 1.2), main = NULL, 
                    xlab = "Distance to CPS (bp)", ylab = "Normalized Signal", labs = c("WT PROseq", "Spt5 ChIPseq"))
# WT (DMSO) PROseq vs Spt5-V5B ChipSeq
WTDMSOPROseq_Vs_Spt5V5B_cps = (metaDataCPS$sample == 1 |  metaDataCPS$sample == 16)
lattice_meta.proSeq(filename = "WT_DMSO_PROseq_vs_Spt5V5B_ChIPseq_CPS_Meta_AllFiltgenes.pdf", 
                    df = metaDataCPS[WTDMSOPROseq_Vs_Spt5V5B_cps,], ylim = c(0, 1.2), main = NULL, 
                    xlab = "Distance to CPS (bp)", ylab = "Fraction of Max Signal", labs = c("WT PROseq", "Spt5 ChIPseq"))
# WT (DMSO) PROseq vs pSpt5 ChipSeq
WTDMSOPROseq_Vs_pSpt5A_cps = (metaDataCPS$sample == 1 |  metaDataCPS$sample == 13)
lattice_meta.proSeq(filename = "WT_DMSO_PROseq_vs_pSpt5A_ChIPseq_CPS_Meta_AllFiltgenes.pdf", 
                    df = metaDataCPS[WTDMSOPROseq_Vs_pSpt5A_cps,], ylim = c(0, 1.2), main = NULL, 
                    xlab = "Distance to CPS (bp)", ylab = "Fraction of Max Signal", labs = c("WT PROseq", "pSpt5 ChIPseq"))
# WT (DMSO) PROseq vs pSpt5B ChipSeq
WTDMSOPROseq_Vs_pSpt5B_cps = (metaDataCPS$sample == 1 |  metaDataCPS$sample == 14)
lattice_meta.proSeq(filename = "WT_DMSO_PROseq_vs_pSpt5B_ChIPseq_CPS_Meta_AllFiltgenes.pdf", 
                    df = metaDataCPS[WTDMSOPROseq_Vs_pSpt5B_cps,], ylim = c(0, 1.2), main = NULL, 
                    xlab = "Distance to CPS (bp)", ylab = "Fraction of Max Signal", labs = c("WT PROseq", "pSpt5 ChIPseq"))

#######################
# plot PRO-seq and ChIP-seq for both pSpt5 and Spt5 all on same plot
WTDMSOPROseq_Vs_pSpt5A_Spt5V5A_cps = (metaDataCPS$sample == 1 |  metaDataCPS$sample == 13 |  metaDataCPS$sample == 15)
lattice_meta.proSeq(filename = "WT_DMSO_PROseq_vs_pSpt5A_Spt5V5A_ChIPseq_CPS_Meta_AllFiltgenes.pdf", 
                    df = metaDataCPS[WTDMSOPROseq_Vs_pSpt5A_Spt5V5A_cps,], ylim = c(0, 1.2), main = NULL, 
                    xlab = "Distance to CPS (bp)", ylab = "Normalized Signal", labs = c( "Spt5p ChIPseq", "WT PROseq", "Spt5 ChIPseq"))

metaDataCPS_kbsep = FormatMetaData_lattice_FracMax(wigset, bed = SPGL_1kbsep_filt[,c(1,3,2,4,5,6)], halfWindow = 1000, step = 10)
WTDMSOPROseq_Vs_pSpt5A_Spt5V5A_cps_kbsep = (metaDataCPS_kbsep$sample == 1 |  metaDataCPS_kbsep$sample == 13 |  metaDataCPS_kbsep$sample == 15)
lattice_meta.proSeq(filename = "WT_DMSO_PROseq_vs_pSpt5A_Spt5V5A_ChIPseq_CPS_Meta.pdf", 
                    df = metaDataCPS_kbsep[WTDMSOPROseq_Vs_pSpt5A_Spt5V5A_cps,], ylim = c(0, 1.2), main = NULL, 
                    xlab = "Distance to CPS (bp)", ylab = "Normalized Signal", labs = c( "Spt5p ChIPseq", "WT PROseq", "Spt5 ChIPseq"))

