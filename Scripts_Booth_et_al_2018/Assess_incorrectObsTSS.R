## Noticed an issure with how I was calling Observed TSSs on the Minus strand (found mistake on 09-27-17). Mistake was a result of how I was subtracting background on the minus strand. 
## This bit of code is just trying to determine how many genes were called incorrectly, and generally how incorrect the calls were (how far from correct site).
library(ggplot2)
fig_dir = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/figures/AssesssIncorrectTSS/10-02-17/"
dir.create(fig_dir)

AllGenesObs = read.table(file = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/pombeInfo/pombe.ASM294v1.16.cleaned_sorted_PROcapObservedTSS.bed", sep = "\t")
AllGenesObs_old = read.table(file = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/pombeInfo/incorrectObsTSS/pombe.ASM294v1.16.cleaned_sorted_PROcapObservedTSS.bed", sep = "\t")

GLcompare = merge(AllGenesObs[, c(1:7)], AllGenesObs_old[, c(1:7)], by.x = "V4", by.y = "V4")
startDiff = abs(GLcompare[,"V2.x"] - GLcompare[,"V2.y"])
endDiff = abs(GLcompare[,"V3.x"] - GLcompare[,"V3.y"])

GLcompare1 = cbind(GLcompare, startDiff, endDiff)
badIndx = GLcompare1[,"endDiff"] > 0

nBadCalls = dim(GLcompare1[badIndx,])[1]
percentBadCalls = (nBadCalls / dim(GLcompare)[1])*100
sink(file = paste(fig_dir, "incorrectTSSinformation.txt", sep = ""))
cat("total number of genes with a previously miscalled TSS = ", nBadCalls, "\n")
cat("percent of all genes with a previously miscalled TSS = ", percentBadCalls, "\n")
cat("summary of all distances = \n")
print(summary(GLcompare1[badIndx,]$endDiff))

p = qplot(GLcompare1[badIndx,]$endDiff,
          geom="histogram",
          binwidth = 10,  
          main = "Incorrect obsTSS offset", 
          xlab = "Dist. to Actual TSS")
ggsave(p, filename = paste(fig_dir, "AllIncorrectTSS_distance_Hist.pdf", sep = ""), device = "pdf")

##########################################################################################################################################
## now do the same, but for only the filtered genes

FiltGenesObs = read.table(file = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/pombeInfo/SP_CompleteFilteredGenes.bed", sep = "\t")
FiltGenesObs_old = read.table(file = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/pombeInfo/incorrectObsTSS/SP_CompleteFilteredGenes.bed", sep = "\t")

FiltGLcompare = merge(FiltGenesObs[, c(1:6)], FiltGenesObs_old[, c(1:6)], by.x = "V4", by.y = "V4")
FiltstartDiff = abs(FiltGLcompare[,"V2.x"] - FiltGLcompare[,"V2.y"])
FiltendDiff = abs(FiltGLcompare[,"V3.x"] - FiltGLcompare[,"V3.y"])

FiltGLcompare1 = cbind(FiltGLcompare, FiltstartDiff, FiltendDiff)
FiltbadIndx = FiltGLcompare1[,"FiltendDiff"] > 0

nBadCalls_filt = dim(FiltGLcompare1[FiltbadIndx,])[1]
percentBadCalls_filt = (nBadCalls_filt / dim(FiltGLcompare)[1])*100
cat("number of filtered genes with a previously miscalled TSS = ", nBadCalls_filt, "\n")
cat("percent of filtered genes with a previously miscalled TSS = ", percentBadCalls_filt, "\n")
cat("summary of distances for filtered Genes= \n")
print(summary(FiltGLcompare1[FiltbadIndx,]$FiltendDiff))
sink()

p1 = qplot(FiltGLcompare1[FiltbadIndx,]$FiltendDiff,
          geom="histogram",
          binwidth = 10,  
          main = "Incorrect obsTSS offset", 
          xlab = "Dist. to Actual TSS")
ggsave(p1, filename = paste(fig_dir, "FiltIncorrectTSS_distance_Hist.pdf", sep = ""), device = "pdf")



