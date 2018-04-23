library(DESeq2)
library(ggplot2)
bwpath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/alignment/bw/"
fig_dir = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/figures/Dis2Analysis/DESeq2/test/"
dir.create(fig_dir)
infopath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/pombeinfo/"
countpath = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/analysis/countData/"
DEseq_Table_path = "/Users/gregorybooth/Google\ Drive/SEAGATE_EXP/Fisher_collaboration/scripts/scriptsForPabitra/countData/DESeq_WGSnormed/"
filteredGL = read.table(paste(infopath, "SP_CompleteFilteredGenes.bed", sep = ""))

#Dis2 mutant DEseq data. 
FisherWT_dis2_11_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_FisherWT_dis2_11_gb_res_WGSnormed.txt", sep = ""), head = T)
FisherWT_dis2_11_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_FisherWT_dis2_11_pr_res_WGSnormed.txt", sep = ""), head = T)
FisherWT_dis2_11_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_FisherWT_dis2_11_term_res_WGSnormed.txt", sep = ""), head = T)
FisherWT_dis2Del_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_FisherWT_dis2Del_gb_res_WGSnormed.txt", sep = ""), head = T)
FisherWT_dis2Del_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_FisherWT_dis2Del_pr_res_WGSnormed.txt", sep = ""), head = T)
FisherWT_dis2Del_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_FisherWT_dis2Del_term_res_WGSnormed.txt", sep = ""), head = T)
FisherWT_dis2_T316A_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_dis2_T316A_gb_res_WGSnormed.txt", sep = ""), head = T)
FisherWT_dis2_T316A_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_dis2_T316A_pr_res_WGSnormed.txt", sep = ""), head = T)
FisherWT_dis2_T316A_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_dis2_T316A_term_res_WGSnormed.txt", sep = ""), head = T)
FisherWT_dis2_T316D_gb_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_dis2_T316D_gb_res_WGSnormed.txt", sep = ""), head = T)
FisherWT_dis2_T316D_pr_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_dis2_T316D_pr_res_WGSnormed.txt", sep = ""), head = T)
FisherWT_dis2_T316D_term_res = read.table(file = paste(DEseq_Table_path ,"DEseq2_AllGenes_dis2_T316D_term_res_WGSnormed.txt", sep = ""), head = T)
### 

plotFCranks = function(DESeqRes, fileName = "test.pdf"){
  sortDF = DESeqRes[order(DESeqRes[,2]),]
  sortDF$rank = 1:dim(sortDF)[1]
  sortDF$name = row.names(sortDF)
  p = ggplot(sortDF, aes(x=rank, y=log2FoldChange)) +
        geom_point() + 
        #geom_text(label=rownames(sortDF[1:10,])) # label bottom and top changed
        #geom_text(data=subset(sortDF, rank < 5 | rank > dim(sortDF)[1]-5),
        geom_text(data=subset(sortDF["SPBC776.02c",]),
            aes(label=name), size = 2, colour ="red")
   ggsave(p, filename = paste(fig_dir, fileName, sep = ""),device = pdf, width = 5, height = 5)
}
### make ranked plots of fold change
plotFCranks(FisherWT_dis2Del_gb_res, fileName = "dis2Del_rankedFC.pdf")
plotFCranks(FisherWT_dis2_11_gb_res, fileName = "dis2_11_rankedFC.pdf")
plotFCranks(FisherWT_dis2_T316A_gb_res, fileName = "dis2_T316A_rankedFC.pdf")
plotFCranks(FisherWT_dis2_T316D_gb_res, fileName = "dis2_T316D_rankedFC.pdf")
#################################################################################################
# function for plotting MA-plots by parsing out the sig. changed genes 
makeMA = function(DESeq_res, pval = 0.01, filename = 'Experiment.pdf', main = NULL, xlab = "Mean PRO-seq Counts (log10)", ylab = "fold change MB3-PP1/DMSO (log2)", geneList = NULL){
  if (is.null(geneList)){
    DESeq_res = DESeq_res
  }
  else{
    DESeq_res_filt = merge(DESeq_res, geneList, by.x = 0, by.y = 4)[,c(1:7)]
    rownames(DESeq_res_filt) = DESeq_res_filt[,1]
    DESeq_res = DESeq_res_filt[,c(2:7)]
  }
  unchangedGenes = DESeq_res[complete.cases(DESeq_res) & DESeq_res$padj>= pval, ]
  SigUpgenes = DESeq_res[complete.cases(DESeq_res) & DESeq_res$padj<= pval & DESeq_res$log2FoldChange > 0 , ]
  SigDownGenes = DESeq_res[complete.cases(DESeq_res) & DESeq_res$padj<= pval & DESeq_res$log2FoldChange < 0 , ]
  cat("Unchanged Genes = ", length(unchangedGenes[,1]), "\n", "Sig Up genes = ", length(SigUpgenes[,1]), "\n", "SigDownGenes = ", length(SigDownGenes[,1]), "\n" )
  pdf(file = paste(fig_dir, filename, sep = ""), width = 7, height = 7)
  plot(log10(unchangedGenes$baseMean), unchangedGenes$log2FoldChange, ylim = c(-3, 3), pch = 16, 
       main = main, xlab = xlab, ylab = ylab, panel.first = grid( nx = NULL, lty = 1, lwd = 2),
       cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, col = adjustcolor("black", alpha = 0.5))
  points(log10(SigUpgenes$baseMean), SigUpgenes$log2FoldChange, col = adjustcolor("steelblue", alpha = 0.5), pch = 16)
  points(log10(SigDownGenes$baseMean), SigDownGenes$log2FoldChange, col = adjustcolor("forestgreen",alpha = 0.5), pch = 16)
  legend("topleft", c(sprintf("Sig Up: N = %i", length(SigUpgenes[,1])), sprintf("Unchanged: N = %i", length(unchangedGenes[,1])), sprintf("Sig Down: N = %i", length(SigDownGenes[,1]))), text.col = c("steelblue", "black", "forestgreen"), cex = 1.5)
  abline(h=0, lty = 2, lwd = 2)
  dev.off()
}
## Apply the makeMA function to results dataframes generated above
# Comparisons between different Dis2 mutants:
# WT7 vs Dis2-11 
makeMA(FisherWT_dis2_11_gb_res, filename = "SP_FisherWT_dis2_11_gb_DESeqResults.pdf", main = "WT vs dis2-11 Change in Gene-body Counts", geneList = filteredGL)
makeMA(FisherWT_dis2_11_pr_res, filename = "SP_FisherWT_dis2_11_pr_DESeqResults.pdf", main = "WT vs dis2-11 Change in Promoter-prox Counts", geneList = filteredGL)
makeMA(FisherWT_dis2_11_term_res, filename = "SP_FisherWT_dis2_11_term_DESeqResults.pdf", main = "WT vs dis2-11 Change in Post-CPS Counts", geneList = filteredGL)
# WT7 vs Dis2Delete 
makeMA(FisherWT_dis2Del_gb_res, filename = "SP_FisherWT_dis2Del_gb_DESeqResults.pdf", main = "WT vs dis2Del Change in Gene-body Counts", geneList = filteredGL)
makeMA(FisherWT_dis2Del_pr_res, filename = "SP_FisherWT_dis2Del_pr_DESeqResults.pdf", main = "WT vs dis2Del Change in Promoter-prox Counts", geneList = filteredGL)
makeMA(FisherWT_dis2Del_term_res, filename = "SP_FisherWT_dis2Del_term_DESeqResults.pdf", main = "WT vs dis2Del Change in Post-CPS Counts", geneList = filteredGL)
# WT7 vs Dis2-T316A
makeMA(FisherWT_dis2_T316A_gb_res, filename = "SP_FisherWT_dis2_T316A_gb_DESeqResults.pdf", main = "WT vs dis2_T316A Change in Gene-body Counts", geneList = filteredGL)
makeMA(FisherWT_dis2_T316A_pr_res, filename = "SP_FisherWT_dis2_T316A_pr_DESeqResults.pdf", main = "WT vs dis2_T316A Change in Promoter-prox Counts", geneList = filteredGL)
makeMA(FisherWT_dis2_T316A_term_res, filename = "SP_FisherWT_dis2_T316A_term_DESeqResults.pdf", main = "WT vs dis2_T316A Change in Post-CPS Counts", geneList = filteredGL)
# WT7 vs Dis2-T316D
makeMA(FisherWT_dis2_T316D_gb_res, filename = "SP_FisherWT_dis2_T316D_gb_DESeqResults.pdf", main = "WT vs dis2_T316D Change in Gene-body Counts", geneList = filteredGL)
makeMA(FisherWT_dis2_T316D_pr_res, filename = "SP_FisherWT_dis2_T316D_pr_DESeqResults.pdf", main = "WT vs dis2_T316D Change in Promoter-prox Counts", geneList = filteredGL)
makeMA(FisherWT_dis2_T316D_term_res, filename = "SP_FisherWT_dis2_T316D_term_DESeqResults.pdf", main = "WT vs dis2_T316D Change in Post-CPS Counts", geneList = filteredGL)


