## 05-16-17
## Using Charles Danko's wave caller to try to identify the cleared genebody region between the advancing and clearing waves after inhibiting Cdk9as in S. pombe. 
## Such boundaries, should in effect, simultaneously give us the front of the advancing wave and back of the clearing wave. 
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
source("/Volumes/SEAGATE_EXP/Fisher_collaboration/scripts/WaveCalling/polymeraseWaves-master/polymeraseWave_gbEdit_callBOTH.bw.v2.R")

fig_dir = "/Volumes/SEAGATE_EXP/Fisher_collaboration/analysis/figures/cdk9asTC_waves/06-01-17/BothWaves/TryTSmooth20_6Kbgenes/"#
dir.create(fig_dir)
infopath = "/Volumes/SEAGATE_EXP/Fisher_collaboration/pombeinfo/"
normedBWpath = "/Volumes/SEAGATE_EXP/Fisher_collaboration/alignment_04-12-17/SpikeNormed_bw/CombinedReps/"
normedBWpath_indivReps = "/Volumes/SEAGATE_EXP/Fisher_collaboration/alignment_04-12-17/SpikeNormed_bw/"
SPmappability = load.bigWig(paste(infopath, "SPombeMappable_SacCerNonOverLap_expanded_ChrFix.bw", sep = ""))
SPGL_1kbsep = read.table(file = paste(infopath, "SP_PROcapObservedTSS_1000bpSep.bed", sep = ""), stringsAsFactors = F)
ObsTSS_NOoverlap = read.table(file = paste(infopath, "SP_PROcapObservedTSS_noOverlap.bed", sep = ""))
filteredGL = read.table(paste(infopath, "SP_CompleteFilteredGenes.bed", sep = ""))
## function below will just place the CPS farther downstream, since we want the script to know where to expect signal to begin.
## it will also swap the strands of the genes, because we are going to be scanning from 3' to 5' 
AdjTSS <- function(genes, downStreamDist = 500){
  plusStrand = genes[,6] == "+"
  minusStrand = genes[,6]== "-"
  res = genes
  res[plusStrand,]$V2 = res[plusStrand,]$V2 + downStreamDist
  res[minusStrand,]$V3 = res[minusStrand,]$V3 - downStreamDist
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

## Test with single gene 
# set the TSS downstream of the actual TSS.  We want the first state to be initiated within the first wave (Trying to call the second state:  between two waves)
#SPAPB1E7.07 <- data.frame(chr = "chr01", start = 3302825, end = 3310369, 
#                               str = "+", SYMBOL = "SPAPB1E7.07", ID = "testGene",
#                               stringsAsFactors = F)
#SPAPB1E7.07_offset1000 <- data.frame(chr = "chr01", start = 3303825, end = 3310369, 
#                          str = "+", SYMBOL = "SPAPB1E7.07", ID = "testGene",
#                          stringsAsFactors = F)
#test_1min_pw = polymeraseWaveBW(reads1_plus = halfMinBW.pl, reads1_minus = halfMinBW.mn, 
#                                reads2_plus = zeroMinBW.pl, reads2_minus = zeroMinBW.mn, 
#                                finterWindowSize = 1000, TSmooth = 5, genes = SPAPB1E7.07_offset1000, 
#                                approxDist = 3000, size = 25, upstreamDist = 1000, 
#                                prefix = paste(fig_dir, "1minBothWave_25BPwin_apx3000_upstream1000_TSmooth5", sep = ""), 
#                                emissionDistAssumption = "2waves", returnVal = "alldata")


## run on all filtered genes of at least a given length (i.e. 4 KB)
## reorder genelist columns to match required 
## NOTE:  For calling the center of advancing and clearing waves, we have to hack the parameters quite a bit. 
### We want to initialize the HMM within the advancing wave (This will be slightly different for each time point): 
#######  To do this, we will adjust the TSS close to the estimated end of the advancing wave (i.e. downstream by 1kb for 30sec, 1.5 kb for 1 min., etc)
#######  The upstream distance should not go all the way up to the TSS, but about half of this distance (thus it will also differ between timepoints)
#######  approxDist now represents where we think the start of the clearing wave is, and relative to the newly set TSS.  
#######  Therefore, we will need to back calculate where the "wave calls" are relative to the actual TSS (needs to take into account the upstreamDist and the TSS offset)
idx = (filteredGL[,3] - filteredGL[,2] >= 6000) 
#size = 500
##########
## Loop entire wave calling script for varying size and TSmooth values
for (size in c(25,50,100,250,500)){
  for (Tsmooth in c(NA, 2, 5, 10, 20)){
    cat("determining space between advancing waves for all time points with combined replicates using window size: ", size, ", and TSmooth = ", Tsmooth, "\n")
dir.create(paste(fig_dir, "CombinedReps/", sep = ""))
dir.create(paste(fig_dir, "CombinedReps/", size, "bpWindowing/", sep = ""))
dir.create(paste(fig_dir, "CombinedReps/", size, "bpWindowing/", Tsmooth, "_TSmooth/",sep = ""))
# Note that I've set all of the ApproxDist = 1000
## this is because the clearing wave transition point, relative to CPS, will be influenced by gene length in addition to time point. 
# Time-point specific parameters: 
GL_adjTSS_30sec = AdjTSS(filteredGL[idx,], downStreamDist = 1000) # downstream approximately the distance of the advancing wave (avg.)
GL_adjTSS_1min = AdjTSS(filteredGL[idx,], downStreamDist = 1000) # downstream approximately the distance of the advancing wave (avg.)
GL_adjTSS_2.5min = AdjTSS(filteredGL[idx,], downStreamDist = 1000) # downstream approximately the distance of the advancing wave (avg.)
#GL_adjTSS_5min = AdjTSS(filteredGL[idx,], downStreamDist = 2500) # downstream approximately the distance of the advancing wave (avg.)
GL_30sec = GL_adjTSS_30sec[,c(1,2,3,6,4,5)]
GL_1min = GL_adjTSS_1min[,c(1,2,3,6,4,5)]
GL_2.5min = GL_adjTSS_2.5min[,c(1,2,3,6,4,5)]
#GL_5min = GL_adjTSS_5min[,c(1,2,3,6,4,5)]
upstreamDist_30sec = 1000
upstreamDist_1min = 1000
upstreamDist_2.5min = 1000
#upstreamDist_5min = 1000
approxClearStart_30sec = 2000 # note:  relative to the shifted TSS (actually 3 kb from TSS)
approxClearStart_1min = 3000 # note:  relative to the shifted TSS (actually 3 kb from TSS)
approxClearStart2.5min = 4500 # note:  relative to the shifted TSS (actually 3 kb from TSS).  Actually farther, but some genes are only 4kb long (i.e. max).
#approxClearStart5min = 1000 # note:  relative to the shifted TSS (actually 3.5 kb from TSS).  Actually farther, but some genes are only 4kb long (i.e. max).
## run 30 second caller using above parameters
pw_30sec= polymeraseWaveBW(reads1_plus = halfMinBW.pl, reads1_minus = halfMinBW.mn, reads2_plus = zeroMinBW.pl, reads2_minus = zeroMinBW.mn, finterWindowSize = 1000,
                                 TSmooth = Tsmooth, genes = GL_30sec, approxDist = approxClearStart_30sec, size = size, upstreamDist = upstreamDist_30sec, prefix = paste(fig_dir, "/CombinedReps/", size, "bpWindowing/", Tsmooth, "_TSmooth/","30secBothWaves", sep = ""), emissionDistAssumption = "2waves",
                                 returnVal = "alldata")
## run 1 minute caller using above parameters
pw_1min = polymeraseWaveBW(reads1_plus = oneMinBW.pl, reads1_minus = oneMinBW.mn, reads2_plus = zeroMinBW.pl, reads2_minus = zeroMinBW.mn, finterWindowSize = 1000,
                                TSmooth = Tsmooth, genes = GL_1min, approxDist = approxClearStart_1min, size = size, upstreamDist = upstreamDist_1min, prefix = paste(fig_dir, "/CombinedReps/", size, "bpWindowing/", Tsmooth, "_TSmooth/","1minBothWaves", sep = ""), emissionDistAssumption = "2waves",
                                returnVal = "alldata")
# run 2.5 minute specific parameters: 
pw_2_5min = polymeraseWaveBW(reads1_plus = two.5MinBW.pl, reads1_minus = two.5MinBW.mn, reads2_plus = zeroMinBW.pl, reads2_minus = zeroMinBW.mn, finterWindowSize = 1000,
                                  TSmooth = Tsmooth, genes = GL_2.5min, approxDist = approxClearStart2.5min, size = size, upstreamDist = upstreamDist_2.5min, prefix = paste(fig_dir, "/CombinedReps/", size, "bpWindowing/", Tsmooth, "_TSmooth/","2.5minBothWaves", sep = ""), emissionDistAssumption = "2waves",
                                  returnVal = "alldata")
# run 5 min specific parameters (use the same as 2.5 min (probably not going to use data at all)): 
#pw_5min = polymeraseWaveBW(reads1_plus = fiveMinBW.pl, reads1_minus = fiveMinBW.mn, reads2_plus = zeroMinBW.pl, reads2_minus = zeroMinBW.mn, finterWindowSize = 1000,
#                                TSmooth = NA, genes = GL_5min, approxDist = approxClearStart5min, size = size, upstreamDist = upstreamDist_5min, prefix = paste(fig_dir, "/CombinedReps/", size, "bpWindowing/", "5minBothWaves", sep = ""), emissionDistAssumption = "2waves",
#                                returnVal = "alldata")
# save the wave calling output above to an ".RData" file in the fig_dir directory
save(pw_30sec, pw_1min, pw_2_5min, file = paste(fig_dir, "CombinedReps/", size, "bpWindowing/", Tsmooth, "_TSmooth/", "cdk9_TC_BothWavesCalldat.RData", sep=""))
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
cat("determining space between advancing waves for all time points with replicate 1 using window size: ", size, ", and TSmooth = ", Tsmooth, "\n")
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

dir.create(paste(fig_dir, "Rep1/", sep = ""))
dir.create(paste(fig_dir, "Rep1/", size, "bpWindowing/", sep = ""))
dir.create(paste(fig_dir, "Rep1/", size, "bpWindowing/", Tsmooth, "_TSmooth/", sep = ""))
pw_30sec_r1 = polymeraseWaveBW(reads1_plus = halfMinBWr1.pl, reads1_minus = halfMinBWr1.mn, reads2_plus = zeroMinBWr1.pl, reads2_minus = zeroMinBWr1.mn, finterWindowSize = 1000,
                           TSmooth = Tsmooth, genes = GL_30sec, approxDist = approxClearStart_30sec, size = size, upstreamDist = upstreamDist_30sec, emissionDistAssumption = "2waves",
                           returnVal = "alldata")
pw_1min_r1 = polymeraseWaveBW(reads1_plus = oneMinBWr1.pl, reads1_minus = oneMinBWr1.mn, reads2_plus = zeroMinBWr1.pl, reads2_minus = zeroMinBWr1.mn, finterWindowSize = 1000,
                           TSmooth = Tsmooth, genes = GL_1min, approxDist = approxClearStart_1min, size = size, upstreamDist = upstreamDist_1min, emissionDistAssumption = "2waves",
                           returnVal = "alldata")
pw_2_5min_r1 = polymeraseWaveBW(reads1_plus = two.5MinBWr1.pl, reads1_minus = two.5MinBWr1.mn, reads2_plus = zeroMinBWr1.pl, reads2_minus = zeroMinBWr1.mn, finterWindowSize = 1000,
                             TSmooth = Tsmooth, genes = GL_2.5min, approxDist = approxClearStart2.5min, size = size, upstreamDist = upstreamDist_2.5min, emissionDistAssumption = "2waves",
                             returnVal = "alldata")
#pw_5min_r1 = polymeraseWaveBW(reads1_plus = fiveMinBWr1.pl, reads1_minus = fiveMinBWr1.mn, reads2_plus = zeroMinBWr1.pl, reads2_minus = zeroMinBWr1.mn, finterWindowSize = 1000,
#                           TSmooth = NA, genes = GL_5min, approxDist = approxClearStart5min, size = size, upstreamDist = upstreamDist_5min, emissionDistAssumption = "2waves",
#                           returnVal = "alldata")
# save the wave calling output above to an ".RData" file in the fig_dir directory
save(pw_30sec_r1, pw_1min_r1, pw_2_5min_r1, file = paste(fig_dir, "Rep1/", size, "bpWindowing/", Tsmooth, "_TSmooth/", "cdk9_TCr1_BothWavesCalldat.RData", sep=""))

##############################################################################################################################################
# run same analysis for only bio replicate 2
##  Combined reps
cat("determining space between advancing waves for all time points with replicate two using window size: ", size, ", and TSmooth = ", Tsmooth, "\n")
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

dir.create(paste(fig_dir, "Rep2/", sep = ""))
dir.create(paste(fig_dir, "Rep2/", size, "bpWindowing/", sep = ""))
dir.create(paste(fig_dir, "Rep2/", size, "bpWindowing/", Tsmooth, "_TSmooth/", sep = ""))
pw_30sec_r2 = polymeraseWaveBW(reads1_plus = halfMinBWr2.pl, reads1_minus = halfMinBWr2.mn, reads2_plus = zeroMinBWr2.pl, reads2_minus = zeroMinBWr1.mn, finterWindowSize = 1000,
                               TSmooth = Tsmooth, genes = GL_30sec, approxDist = approxClearStart_30sec, size = size, upstreamDist = upstreamDist_30sec, emissionDistAssumption = "2waves",
                               returnVal = "alldata")
pw_1min_r2 = polymeraseWaveBW(reads1_plus = oneMinBWr2.pl, reads1_minus = oneMinBWr2.mn, reads2_plus = zeroMinBWr2.pl, reads2_minus = zeroMinBWr1.mn, finterWindowSize = 1000,
                              TSmooth = Tsmooth, genes = GL_1min, approxDist = approxClearStart_1min, size = size, upstreamDist = upstreamDist_1min, emissionDistAssumption = "2waves",
                              returnVal = "alldata")
pw_2_5min_r2 = polymeraseWaveBW(reads1_plus = two.5MinBWr2.pl, reads1_minus = two.5MinBWr2.mn, reads2_plus = zeroMinBWr2.pl, reads2_minus = zeroMinBWr1.mn, finterWindowSize = 1000,
                                TSmooth = Tsmooth, genes = GL_2.5min, approxDist = approxClearStart2.5min, size = size, upstreamDist = upstreamDist_2.5min, emissionDistAssumption = "2waves",
                                returnVal = "alldata")
#pw_5min_r2 = polymeraseWaveBW(reads1_plus = fiveMinBWr2.pl, reads1_minus = fiveMinBWr2.mn, reads2_plus = zeroMinBWr2.pl, reads2_minus = zeroMinBWr1.mn, finterWindowSize = 1000,
#                              TSmooth = NA, genes = GL_5min, approxDist = approxClearStart5min, size = size, upstreamDist = upstreamDist_5min, emissionDistAssumption = "2waves",
#                              returnVal = "alldata")
# save the wave calling output above to an ".RData" file in the fig_dir directory
save(pw_30sec_r2, pw_1min_r2, pw_2_5min_r2, file = paste(fig_dir, "Rep2/", size, "bpWindowing/", Tsmooth, "_TSmooth/", "cdk9_TCr2_BothWavesCalldat.RData", sep=""))
}
}







