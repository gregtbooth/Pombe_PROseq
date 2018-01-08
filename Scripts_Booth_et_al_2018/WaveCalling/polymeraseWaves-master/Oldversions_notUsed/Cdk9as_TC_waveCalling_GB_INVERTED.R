## 05-11-17
## Using Charles Danko's wave caller to try to identify Clearing Wave backends after inhibiting Cdk9as in S. pombe. 
#
## NOTE:  This is a HACKED attempt to just flip genes strands as well as the data to try to start wave calls from 3' to 5' end
#
# the groHMM package was largely developed by Charles Danko and downloaded from bioconductor. 
# This lager package contains the scripts for wave calling (polymeraseWave) that was first demonstrated by Hah et al (2011) and  Danko et al. 2013,
# but used in Mahat et al. 2016
## Note:  Also need to source Charles Danko's script (downloaded from GitHub) that has modified the polymeraseWave function to be compatible with bigwigs.
## Note:  wave calling function depends on functions from a more recent version of the bigWig package. Had to use lunchroom Mac to run scripts. 
library(groHMM)
library(bigWig)
source("/Volumes/SEAGATE_EXP/Fisher_collaboration/scripts/GB_Functions.R")
source("/Volumes/SEAGATE_EXP/Fisher_collaboration/scripts/WaveCalling/polymeraseWaves-master/polymeraseWave_gbEdit.bw.R")

fig_dir = "/Volumes/SEAGATE_EXP/Fisher_collaboration/analysis/figures/cdk9asTC_waves/05-12-17/ClearingWaves_Hacked/"
dir.create(fig_dir)
infopath = "/Volumes/SEAGATE_EXP/Fisher_collaboration/pombeinfo/"
normedBWpath = "/Volumes/SEAGATE_EXP/Fisher_collaboration/alignment_04-12-17/SpikeNormed_bw/CombinedReps/"
normedBWpath_indivReps = "/Volumes/SEAGATE_EXP/Fisher_collaboration/alignment_04-12-17/SpikeNormed_bw/"
SPmappability = load.bigWig(paste(infopath, "SPombeMappable_SacCerNonOverLap_expanded_ChrFix.bw", sep = ""))
SPGL_1kbsep = read.table(file = paste(infopath, "SP_PROcapObservedTSS_1000bpSep.bed", sep = ""), stringsAsFactors = F)
ObsTSS_NOoverlap = read.table(file = paste(infopath, "SP_PROcapObservedTSS_noOverlap.bed", sep = ""))
filteredGL = read.table(paste(infopath, "SP_CompleteFilteredGenes.bed", sep = ""))
# fix parameters here for all runs. 
size = 500
upstreamDist = 500
## function below will just place the CPS farther downstream, since we want the script to know where to expect signal to begin.
## it will also swap the strands of the genes, because we are going to be scanning from 3' to 5' 
AdjCPS <- function(genes, downStreamDist = 500){
  plusStrand = genes[,6] == "+"
  minusStrand = genes[,6]== "-"
  res = genes
  res[plusStrand,]$V3 = res[plusStrand,]$V3 + downStreamDist
  res[minusStrand,]$V2 = res[minusStrand,]$V2 - downStreamDist
  res[plusStrand,]$V6 = "-"
  res[minusStrand,]$V6 = "+"
  return(res)
}

# load bigwigs for each time point
#zeroMinBW = list(load.bigWig(paste(normedBWpath, "ePROseq_Pombe_CDK9as_0min_DMSO_COMBINED_pombe_SpikeNormed_plus.bw", sep = "" )), load.bigWig(paste(normedBWpath, "ePROseq_Pombe_CDK9as_0min_DMSO_COMBINED_pombe_SpikeNormed_minus.bw", sep = "" )), "Cdk9as_0min_combined")

##  Combined reps
zeroMinBW.pl = paste(normedBWpath, "ePROseq_Pombe_CDK9as_0min_DMSO_COMBINED_pombe_SpikeNormed_plus.bw", sep = "" )
zeroMinBW.mn = paste(normedBWpath, "ePROseq_Pombe_CDK9as_0min_DMSO_COMBINED_pombe_SpikeNormed_minus.bw", sep = "" )
halfMinBW.pl = paste(normedBWpath, "ePROseq_Pombe_CDK9as_30sec_3MBPP1_10uM_COMBINED_pombe_SpikeNormed_plus.bw", sep = "" )
halfMinBW.mn = paste(normedBWpath, "ePROseq_Pombe_CDK9as_30sec_3MBPP1_10uM_COMBINED_pombe_SpikeNormed_minus.bw", sep = "" )
oneMinBW.pl = paste(normedBWpath, "ePROseq_Pombe_CDK9as_1min_3MBPP1_10uM_COMBINED_pombe_SpikeNormed_plus.bw", sep = "" )
oneMinBW.mn = paste(normedBWpath, "ePROseq_Pombe_CDK9as_1min_3MBPP1_10uM_COMBINED_pombe_SpikeNormed_minus.bw", sep = "" )
two.5MinBW.pl = paste(normedBWpath, "ePROseq_Pombe_CDK9as_2min30sec_3MBPP1_10uM_COMBINED_pombe_SpikeNormed_plus.bw", sep = "" )
two.5MinBW.mn = paste(normedBWpath, "ePROseq_Pombe_CDK9as_2min30sec_3MBPP1_10uM_COMBINED_pombe_SpikeNormed_minus.bw", sep = "" )
fiveMinBW.pl = paste(normedBWpath, "ePROseq_Pombe_CDK9as_5min_3MBPP1_10uM_COMBINED_pombe_SpikeNormed_plus.bw", sep = "" )
fiveMinBW.mn = paste(normedBWpath, "ePROseq_Pombe_CDK9as_5min_3MBPP1_10uM_COMBINED_pombe_SpikeNormed_minus.bw", sep = "" )


## Test with single gene (Swap the strand)
## Also, because of transcription beyond 3' CPS, place gene 3' ends 500 bp farther downstream
#SPAPB1E7.07 <- data.frame(chr = "chr01", start = 3302825, end = 3310869, 
#                          str = "-", SYMBOL = "SPAPB1E7.07", ID = "testGene",
#                          stringsAsFactors = F)
# Reverse strands of bigWig files 
#test_1min_pw = polymeraseWaveBW(reads1_plus = fiveMinBW.mn, reads1_minus = fiveMinBW.pl, 
#                                reads2_plus = zeroMinBW.mn, reads2_minus = zeroMinBW.pl, 
#                                finterWindowSize = 1000, TSmooth = NA, genes = SPAPB1E7.07, 
#                                approxDist = 1500, size = 100, upstreamDist = 500, 
#                                prefix = paste(fig_dir, "5minClearWave", sep = ""), 
#                                emissionDistAssumption = "gamma", returnVal = "alldata")


## run on all filtered genes of at least a given length (i.e. 3 KB)
## reorder genelist columns to match required 
filteredGL_adjCPS = AdjCPS(filteredGL, downStreamDist = 500)
GL = filteredGL_adjCPS[,c(1,2,3,6,4,5)]
idx = (GL[,3] - GL[,2] >= 4500) # increased to 4500 because of the extra length added to 3' ends
dir.create(paste(fig_dir, "/CombinedReps/", sep = ""))
dir.create(paste(fig_dir, "/CombinedReps/", size, "bpWindowing/", sep = ""))
# Note that I've set all of the ApproxDist = 1000
## this is because the clearing wave transition point, relative to CPS, will be influenced by gene length in addition to time point. 
pw_30sec= polymeraseWaveBW(reads1_plus = halfMinBW.mn, reads1_minus = halfMinBW.pl, reads2_plus = zeroMinBW.mn, reads2_minus = zeroMinBW.pl, finterWindowSize = 1000,
                                 TSmooth = NA, genes = GL[idx,], approxDist = 1000, size = size, upstreamDist = upstreamDist, prefix = paste(fig_dir, "/CombinedReps/", size, "bpWindowing/", "30secClearWave", sep = ""), emissionDistAssumption = "gamma",
                                 returnVal = "alldata")
pw_1min = polymeraseWaveBW(reads1_plus = oneMinBW.mn, reads1_minus = oneMinBW.pl, reads2_plus = zeroMinBW.mn, reads2_minus = zeroMinBW.pl, finterWindowSize = 1000,
                                TSmooth = NA, genes = GL[idx,], approxDist = 1000, size = size, upstreamDist = upstreamDist, prefix = paste(fig_dir, "/CombinedReps/", size, "bpWindowing/", "1minClearWave", sep = ""), emissionDistAssumption = "gamma",
                                returnVal = "alldata")
pw_2.5min = polymeraseWaveBW(reads1_plus = two.5MinBW.mn, reads1_minus = two.5MinBW.pl, reads2_plus = zeroMinBW.mn, reads2_minus = zeroMinBW.pl, finterWindowSize = 1000,
                                  TSmooth = NA, genes = GL[idx,], approxDist = 1000, size = size, upstreamDist = upstreamDist, prefix = paste(fig_dir, "/CombinedReps/", size, "bpWindowing/", "2.5minClearWave", sep = ""), emissionDistAssumption = "gamma",
                                  returnVal = "alldata")
pw_5min = polymeraseWaveBW(reads1_plus = fiveMinBW.mn, reads1_minus = fiveMinBW.pl, reads2_plus = zeroMinBW.mn, reads2_minus = zeroMinBW.pl, finterWindowSize = 1000,
                                TSmooth = NA, genes = GL[idx,], approxDist = 1000, size = size, upstreamDist = upstreamDist, prefix = paste(fig_dir, "/CombinedReps/", size, "bpWindowing/", "5minClearWave", sep = ""), emissionDistAssumption = "gamma",
                                returnVal = "alldata")
# save the wave calling output above to an ".RData" file in the fig_dir directory
save(pw_30sec, pw_1min, pw_2.5min, pw_5min, file = paste(fig_dir, "/CombinedReps/", size, "bpWindowing/", "cdk9_TC_ClearingWaveCalldat.RData", sep=""))
# if you want to specifically load the output data (pw_30sec, pw_1min, etc.) generated while calling waves
#load(file = paste(fig_dir, "cdk9_TC_WaveCalldat.RData", sep = ""))
# get just the simplified wave coordinate data
#simple_0.5min <- data.frame(pw_30sec[225])
#simple_1min <- data.frame(pw_1min[225])
#simple_2.5min <- data.frame(pw_2.5min[225])
#simple_5min <- data.frame(pw_5min[225])

##############################################################################################################################################
# run same analysis for only bio replicate 1
##  Combined reps
zeroMinBWr1.pl = paste(normedBWpath_indivReps, "ePROseq_Pombe_CDK9as_0min_DMSO_rep1_pombe_SpikeNormed_plus.bw", sep = "" )
zeroMinBWr1.mn = paste(normedBWpath_indivReps, "ePROseq_Pombe_CDK9as_0min_DMSO_rep1_pombe_SpikeNormed_minus.bw", sep = "" )
halfMinBWr1.pl = paste(normedBWpath_indivReps, "ePROseq_Pombe_CDK9as_30sec_3MBPP1_10uM_rep1_pombe_SpikeNormed_plus.bw", sep = "" )
halfMinBWr1.mn = paste(normedBWpath_indivReps, "ePROseq_Pombe_CDK9as_30sec_3MBPP1_10uM_rep1_pombe_SpikeNormed_minus.bw", sep = "" )
oneMinBWr1.pl = paste(normedBWpath_indivReps, "ePROseq_Pombe_CDK9as_1min_3MBPP1_10uM_rep1_pombe_SpikeNormed_plus.bw", sep = "" )
oneMinBWr1.mn = paste(normedBWpath_indivReps, "ePROseq_Pombe_CDK9as_1min_3MBPP1_10uM_rep1_pombe_SpikeNormed_minus.bw", sep = "" )
two.5MinBWr1.pl = paste(normedBWpath_indivReps, "ePROseq_Pombe_CDK9as_2min30sec_3MBPP1_10uM_rep1_pombe_SpikeNormed_plus.bw", sep = "" )
two.5MinBWr1.mn = paste(normedBWpath_indivReps, "ePROseq_Pombe_CDK9as_2min30sec_3MBPP1_10uM_rep1_pombe_SpikeNormed_minus.bw", sep = "" )
fiveMinBWr1.pl = paste(normedBWpath_indivReps, "ePROseq_Pombe_CDK9as_5min_3MBPP1_10uM_rep1_pombe_SpikeNormed_plus.bw", sep = "" )
fiveMinBWr1.mn = paste(normedBWpath_indivReps, "ePROseq_Pombe_CDK9as_5min_3MBPP1_10uM_rep1_pombe_SpikeNormed_minus.bw", sep = "" )

dir.create(paste(fig_dir, "/Rep1/", sep = ""))
dir.create(paste(fig_dir, "/Rep1/", size, "bpWindowing/", sep = ""))
pw_30sec_r1 = polymeraseWaveBW(reads1_plus = halfMinBWr1.mn, reads1_minus = halfMinBWr1.pl, reads2_plus = zeroMinBWr1.mn, reads2_minus = zeroMinBWr1.pl, finterWindowSize = 1000,
                           TSmooth = NA, genes = GL[idx,], approxDist = 1000, size = size, upstreamDist = upstreamDist, emissionDistAssumption = "gamma",
                           returnVal = "alldata")
pw_1min_r1 = polymeraseWaveBW(reads1_plus = oneMinBWr1.mn, reads1_minus = oneMinBWr1.pl, reads2_plus = zeroMinBWr1.mn, reads2_minus = zeroMinBWr1.pl, finterWindowSize = 1000,
                           TSmooth = NA, genes = GL[idx,], approxDist = 1000, size = size, upstreamDist = upstreamDist, emissionDistAssumption = "gamma",
                           returnVal = "alldata")
pw_2.5min_r1 = polymeraseWaveBW(reads1_plus = two.5MinBWr1.mn, reads1_minus = two.5MinBWr1.pl, reads2_plus = zeroMinBWr1.mn, reads2_minus = zeroMinBWr1.pl, finterWindowSize = 1000,
                             TSmooth = NA, genes = GL[idx,], approxDist = 1000, size = size, upstreamDist = upstreamDist, emissionDistAssumption = "gamma",
                             returnVal = "alldata")
pw_5min_r1 = polymeraseWaveBW(reads1_plus = fiveMinBWr1.mn, reads1_minus = fiveMinBWr1.pl, reads2_plus = zeroMinBWr1.mn, reads2_minus = zeroMinBWr1.pl, finterWindowSize = 1000,
                           TSmooth = NA, genes = GL[idx,], approxDist = 1000, size = size, upstreamDist = upstreamDist, emissionDistAssumption = "gamma",
                           returnVal = "alldata")
# save the wave calling output above to an ".RData" file in the fig_dir directory
save(pw_30sec_r1, pw_1min_r1, pw_2.5min_r1, pw_5min_r1, file = paste(fig_dir, "/Rep1/", size, "bpWindowing/", "cdk9_TCr1_ClearingWaveCalldat.RData", sep=""))

##############################################################################################################################################
# run same analysis for only bio replicate 1
##  Combined reps
zeroMinBWr2.pl = paste(normedBWpath_indivReps, "ePROseq_Pombe_CDK9as_0min_DMSO_rep2_pombe_SpikeNormed_plus.bw", sep = "" )
zeroMinBWr2.mn = paste(normedBWpath_indivReps, "ePROseq_Pombe_CDK9as_0min_DMSO_rep2_pombe_SpikeNormed_minus.bw", sep = "" )
halfMinBWr2.pl = paste(normedBWpath_indivReps, "ePROseq_Pombe_CDK9as_30sec_3MBPP1_10uM_rep2_pombe_SpikeNormed_plus.bw", sep = "" )
halfMinBWr2.mn = paste(normedBWpath_indivReps, "ePROseq_Pombe_CDK9as_30sec_3MBPP1_10uM_rep2_pombe_SpikeNormed_minus.bw", sep = "" )
oneMinBWr2.pl = paste(normedBWpath_indivReps, "ePROseq_Pombe_CDK9as_1min_3MBPP1_10uM_rep2_pombe_SpikeNormed_plus.bw", sep = "" )
oneMinBWr2.mn = paste(normedBWpath_indivReps, "ePROseq_Pombe_CDK9as_1min_3MBPP1_10uM_rep2_pombe_SpikeNormed_minus.bw", sep = "" )
two.5MinBWr2.pl = paste(normedBWpath_indivReps, "ePROseq_Pombe_CDK9as_2min30sec_3MBPP1_10uM_rep2_pombe_SpikeNormed_plus.bw", sep = "" )
two.5MinBWr2.mn = paste(normedBWpath_indivReps, "ePROseq_Pombe_CDK9as_2min30sec_3MBPP1_10uM_rep2_pombe_SpikeNormed_minus.bw", sep = "" )
fiveMinBWr2.pl = paste(normedBWpath_indivReps, "ePROseq_Pombe_CDK9as_5min_3MBPP1_10uM_rep2_pombe_SpikeNormed_plus.bw", sep = "" )
fiveMinBWr2.mn = paste(normedBWpath_indivReps, "ePROseq_Pombe_CDK9as_5min_3MBPP1_10uM_rep2_pombe_SpikeNormed_minus.bw", sep = "" )

dir.create(paste(fig_dir, "/Rep2/", sep = ""))
dir.create(paste(fig_dir, "/Rep2/", size, "bpWindowing/", sep = ""))
pw_30sec_r2 = polymeraseWaveBW(reads1_plus = halfMinBWr2.mn, reads1_minus = halfMinBWr2.pl, reads2_plus = zeroMinBWr2.mn, reads2_minus = zeroMinBWr2.pl, finterWindowSize = 1000,
                               TSmooth = NA, genes = GL[idx,], approxDist = 1000, size = size, upstreamDist = upstreamDist, emissionDistAssumption = "gamma",
                               returnVal = "alldata")
pw_1min_r2 = polymeraseWaveBW(reads1_plus = oneMinBWr2.mn, reads1_minus = oneMinBWr2.pl, reads2_plus = zeroMinBWr2.mn, reads2_minus = zeroMinBWr2.pl, finterWindowSize = 1000,
                              TSmooth = NA, genes = GL[idx,], approxDist = 1000, size = size, upstreamDist = upstreamDist, emissionDistAssumption = "gamma",
                              returnVal = "alldata")
pw_2.5min_r2 = polymeraseWaveBW(reads1_plus = two.5MinBWr2.mn, reads1_minus = two.5MinBWr2.pl, reads2_plus = zeroMinBWr2.mn, reads2_minus = zeroMinBWr2.pl, finterWindowSize = 1000,
                                TSmooth = NA, genes = GL[idx,], approxDist = 1000, size = size, upstreamDist = upstreamDist, emissionDistAssumption = "gamma",
                                returnVal = "alldata")
pw_5min_r2 = polymeraseWaveBW(reads1_plus = fiveMinBWr2.mn, reads1_minus = fiveMinBWr2.pl, reads2_plus = zeroMinBWr2.mn, reads2_minus = zeroMinBWr2.pl, finterWindowSize = 1000,
                              TSmooth = NA, genes = GL[idx,], approxDist = 1000, size = size, upstreamDist = upstreamDist, emissionDistAssumption = "gamma",
                              returnVal = "alldata")
# save the wave calling output above to an ".RData" file in the fig_dir directory
save(pw_30sec_r2, pw_1min_r2, pw_2.5min_r2, pw_5min_r2, file = paste(fig_dir, "/Rep2/", size, "bpWindowing/", "cdk9_TCr2_ClearingWaveCalldat.RData", sep=""))









