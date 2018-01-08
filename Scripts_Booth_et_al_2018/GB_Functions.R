library(bigWig)
library(lattice)
#############################################################################################################
## Functions for making scaled meta plots

# mround is used to round the length (in bases) to a number divisible by the 'base' argument
mround = function(x, base = 30){ 
  #base*ceiling(x/base)# base = 30 was chosen so that there will be a total of 60 windows (30 on each side of the midpoint of scaled window)
  base*floor(x/base)} # rather than rounding up, floor can be used to round down and prevent enchroaching on the unscaled region 

# scale_params assesses the length of your gene, returning a buffer length (1/2 the size of the region to be scaled), a step (window size for scaled region), and a scale facor (windowsize/10)
scale_params = function(start, end, step=10, buffer=300){
  length = end-start
  if (length <1000){
    buffer =  buffer
    step1  = step
    scaleFact = 1
  }
  else{
    length = ((end-start)-600)/2
    buffer = mround(length, base = 30)
    step1 = (buffer/30)
    scaleFact = step1/10 # scale factor is divided by your scaled window count to adjust as if it were a 10bp window (i.e. window size = 5bp will have a scale factor of 2)
  }
  return(c(buffer,step1,scaleFact))
}

collect.many.scaled = function (bed, bigWig.plus, bigWig.minus, flank=1000, FLbuffer = 300, step=10, 
                                do.sum = T) 
{
  #windowSize = ((2*flank) + (4*buffer))%/%step
  midPoint = (bed[, 2] + bed[, 3])/2
  TSS = bed[, 2]
  END = bed[,3]
  #neg.idxs = bed[, 6] == "-"
  #TSS[neg.idxs] = bed[neg.idxs, 3] - 1
  #END[neg.idxs] = bed[neg.idxs, 2] 
  
  prstart = (TSS - flank)
  prend = (TSS + FLbuffer)
  
  finstart = (END - FLbuffer)
  finend = (END + flank)
  #windowSize = (2*flank + 4*FLbuffer)/step
  windowSize = 320
  N = dim(bed)[1]
  result = matrix(nrow = N, ncol = windowSize)
  strands = NULL
  if (dim(bed)[2] >= 6) 
    strands = as.character(bed[, 6])
  else strands = rep("+", N)
  for (i in 1:N) {
    chrom = as.character(bed[i, 1])
    strand = strands[i]
    if (strand == "+") {
      bigWig = bigWig.plus
      midscaleParams = scale_params(bed[i,2], bed[i,3], step, buffer = 200)
      midbuf = midscaleParams[1]
      midstep = midscaleParams[2]
      midScaleFact = midscaleParams[3]
      midstart = (midPoint - midbuf)
      midend = (midPoint + midbuf)
      prrow = collect.counts(bigWig, chrom, prstart[i], prend[i], 
                             step, do.sum)
      midrow =  collect.counts(bigWig, chrom, midstart[i], midend[i], 
                               midstep, do.sum)/midScaleFact
      endrow = collect.counts(bigWig, chrom, finstart[i], finend[i], 
                              step, do.sum)
      row = c(prrow, midrow, endrow)
      if (length(row) == 320){
        result[i, ] = abs(row)}
    }
    else {
      bigWig = bigWig.minus
      if (is.null(bigWig.minus)) 
        bigWig = bigWig.plus
      midscaleParams = scale_params(bed[i,2], bed[i,3], step, buffer = 200)
      midbuf = midscaleParams[1]
      midstep = midscaleParams[2]
      midScaleFact = midscaleParams[3]
      midstart = (midPoint - midbuf)
      midend = (midPoint + midbuf)
      prrow = collect.counts(bigWig, chrom, prstart[i], prend[i], 
                             step, do.sum)
      midrow =  collect.counts(bigWig, chrom, midstart[i], midend[i], 
                               midstep, do.sum)/midScaleFact
      endrow = collect.counts(bigWig, chrom, finstart[i], finend[i], 
                              step, do.sum)
      row = c(prrow, midrow, endrow)
      if (length(row) == 320){
        result[i, ] = abs(rev(row))}
    }
  }
  rownames(result) = bed[,4]
  result[complete.cases(result), ]
}

meta.subsample.scaled = function (bed, bigWig.plus, bigWig.minus, flank=1000, buffer = 300, step=10, do.sum = T) 
{
  values = collect.many.scaled(bed, bigWig.plus, bigWig.minus, flank=flank, buffer, step, 
                               do.sum = do.sum)
  N = dim(values)[1]
  nPermut = 1000
  sampleFrac = 0.1
  #windowSize = ((prend-prstart)+(midend-midstart)+(finend-finstart))/step
  windowSize = 320
  result = matrix(nrow = nPermut, ncol = windowSize)
  M = as.integer(round(N * sampleFrac, 0))
  #values = collect.many.scaled(bed, bigWig.plus, bigWig.minus, flank=flank, buffer, step, 
  #do.sum = do.sum)
  for (i in 1:nPermut) {
    idx <- sample(N, size = M, replace = T)
    result[i, ] = colSums(values[idx, ])/M
  }
  ci9 = sapply(1:windowSize, function(idx) quantile(result[, 
                                                           idx], 0.875))
  ci1 = sapply(1:windowSize, function(idx) quantile(result[, 
                                                           idx], 0.125))
  ci5 = sapply(1:windowSize, function(idx) quantile(result[, 
                                                           idx], 0.5))
  return(list(result, ci9, ci1, ci5))
}

meta.plot.scaled.GROseq = function (result.plus, result.minus, step, xlab = "Scaled Gene Body", ylab = "Average Reads", col1 = 'red', col2 = 'blue', main = NULL, ylim = NULL, bothStrands = TRUE,  ...) {
  if(bothStrands){
    N = length(result.plus[[4]])
    x = ((1:N) - N/2) * step
    if (is.null(ylim)) 
      ylim = c(min(-result.minus[[4]])-.2, max(result.plus[[4]])+.2)
    plot(x, result.plus[[4]], type = "l", col = col1, lwd = 3, 
         xlab = xlab, ylab = ylab, ylim = ylim, main = main, xaxt = "n", cex.main=2, cex.axis = 1.5, cex.lab =1.5, ...)
    polygon(c(x, rev(x)), c(result.plus[[2]], rev(result.plus[[3]])), 
            col = colors()[339], border = NA)
    polygon(c(x, rev(x)), c(-result.minus[[3]], rev(-result.minus[[2]])), 
            col = colors()[407], border = NA)
    lines(x, result.plus[[4]], col = col1, lwd = 3)
    lines(x, -result.minus[[4]], col = col2, lwd = 3)
    lines(x, rep(0, N), lwd = 3, lty = 2)
    axis(side = 1, at = c(-1600,-600,-300,300,600,1600), labels = c(-1000,"TSS","+300","-300","CPA",1000), cex.axis = 1.5)
  }
  else{
    N = length(result.plus[[4]])
    x1 = 1:(N*step)
    x = ((1:N) - N/2) * step
    if (is.null(ylim)) 
      ylim = c(0, max(result.plus[[4]])+.4)
    plot(x, result.plus[[4]], type = "l", col = col1, lwd = 3, 
         xlab = xlab, ylab = ylab, ylim = ylim, main = main, xaxt = "n", cex.main=2, cex.axis = 1.5, cex.lab =1.5, ...)
    polygon(c(x, rev(x)), c(result.plus[[2]], rev(result.plus[[3]])), 
            col = colors()[339], border = NA)
    polygon(c(x, rev(x)), c(-result.minus[[3]], rev(-result.minus[[2]])), 
            col = colors()[407], border = NA)
    lines(x, result.plus[[4]], col = col1, lwd = 3)
    lines(x, -result.minus[[4]], col = col2, lwd = 3)
    lines(x, rep(0, N), lwd = 3, lty = 2)
    axis(side = 1, at = c(-1600,-600,-300,300,600,1600), labels = c(-1000,"TSS","+300","-300","CPA",1000), cex.axis = 1.5)
  }
  
}

meta.overlay = function (result,step, col1=rgb(0,0.5,0), col2= rgb(0,1,0,.2), strand = "")
{
  N = length(result[[4]])
  x = ((1:N) - N/2)* step
  
  if(strand == "+")
    #draw shade area
  {polygon(c(x, rev(x)), c(result[[2]], rev(result[[3]])), col=col2, border=NA)
   # redraw main Plot line on top 
   lines(x,result[[4]],col=col1, lwd=3)}
  else if(strand== "-")
    #draw shade area
  {polygon(c(x, rev(x)), c(-result[[2]], rev(-result[[3]])), col=col2, border=NA)
   # redraw main plot line on top
   lines(x,-result[[4]],col=col1, lwd=3)}
  else
  {print("ERROR: need strand")}
}

#############################################################################################################
## Functions for making Centered meta plots
total.reads = function(bigWig.plus, bigWig.minus){
  TotalReads <- bigWig.plus$mean * bigWig.plus$basesCovered + abs(bigWig.minus$mean * bigWig.minus$basesCovered)  
  return(TotalReads)
}

glist.2.bigWigFormat = function (geneList){
  newGL = geneList[,c(3,5,6,1)]
  newGL[,5]=0
  newGL[,6] =geneList[,4]
  return(newGL)
}

glist.for3pMeta = function (geneList){
  newGL = geneList[,c(3,6,5,1)]
  newGL[,5]=0
  newGL[,6] =geneList[,4]
  return(newGL)
}

meta.plot.GROseq = function (result.plus, result.minus, step, xlab = "Distance to center (bp)", ylab = "Average Reads", col1 = 'red', col2 = 'blue', main = NULL, ylim = NULL, bothStrands = TRUE, ...) {
  if(bothStrands){
    N = length(result.plus[[4]])
    x = ((1:N) - N/2) * step
    if (is.null(ylim)) 
      ylim = c(min(-result.minus[[4]]), max(result.plus[[4]]))
    plot(x, result.plus[[4]], type = "l", col = col1, lwd = 3, 
         xlab = xlab, ylab = ylab, ylim = ylim, main = main, ...)
    polygon(c(x, rev(x)), c(result.plus[[2]], rev(result.plus[[3]])), 
            col = colors()[339], border = NA)
    polygon(c(x, rev(x)), c(-result.minus[[3]], rev(-result.minus[[2]])), 
            col = colors()[407], border = NA)
    lines(x, result.plus[[4]], col = col1, lwd = 3)
    lines(x, -result.minus[[4]], col = col2, lwd = 3)
    lines(x, rep(0, N), lwd = 3, lty = 2)
  }
  else{
    N = length(result.plus[[4]])
    x = ((1:N) - N/2) * step
    if (is.null(ylim)) 
      ylim = c(min(result.plus[[4]]), max(result.plus[[4]]))
    plot(x, result.plus[[4]], type = "l", col = col1, lwd = 3, 
         xlab = xlab, ylab = ylab, ylim = ylim, main = main, ...)
    polygon(c(x, rev(x)), c(result.plus[[2]], rev(result.plus[[3]])), 
            col = colors()[339], border = NA)
    polygon(c(x, rev(x)), c(-result.minus[[3]], rev(-result.minus[[2]])), 
            col = colors()[407], border = NA)
    lines(x, result.plus[[4]], col = col1, lwd = 3)
    lines(x, -result.minus[[4]], col = col2, lwd = 3)
    lines(x, rep(0, N), lwd = 3, lty = 2)
  }
}

#############################################################################################################
## Functions for making Centered meta plots

rows.finite.df = function(df){
  for( i in c(1:ncol(df))){
    df = df[is.finite(df[,i]),]
  }
  return(df)
}

format.genelist = function(genes){
  GL = genes[,c(3,4,5,6,1)]
  return(GL) 
}

getCounts.sliding.pr <- function(wig.p, wig.m, map, genes, off1 = -50, off2 = 150, windowsize = 100) {
  N = dim(genes)[1]
  plusStrand = genes[,2] == '+'
  
  result = vector(mode="integer", length=N)
  mapresult = vector(mode="integer", length=N)
  scanlen = off2 - off1 - windowsize
  for (i in 1:N) {
    chrom = as.character(genes[i, 1])
    start = as.integer(genes[i, 3])
    end = as.integer(genes[i, 4])
    
    if (plusStrand[i]) {
      wig = wig.p
      
      scanResults = c()
      for (n in c(1:scanlen)){
        
        qStart = start + off1 + n
        qEnd = start + n + off1 + windowsize
        
        data = query.bigWig(wig, chrom, qStart, qEnd)
        scanResults = append(scanResults, abs(sum(data[,3])))
        
      }
      mapOffset = which.max(scanResults)  # grabs the index of the first instance of the maximum value ( = bp offset during scanning)
      result[i]  = max(scanResults)
      qMapStart = start + off1 + mapOffset + 36
      qMapEnd  = start + off1 + mapOffset + windowsize + 36
      mapdata = query.bigWig(map, chrom, qMapStart, qMapEnd)
      mapresult[i] = abs(sum(mapdata[,3]))
    } 
    else {
      wig = wig.m
      scanResults = c()
      if (is.null(wig.m))
        wig = wig.p
      
      for (n in c(1: scanlen)){
        qStart = end - off1 - n  - windowsize
        qEnd = end - off1 - n
        
        data = query.bigWig(wig, chrom, qStart, qEnd)
        scanResults = append(scanResults, abs(sum(data[,3])))
      }
      mapOffset = which.max(scanResults)  # grabs the index of the first instance of the maximum value ( = bp offset during scanning)
      result[i]  = max(scanResults)
      qMapStart = end - off1 - mapOffset - windowsize
      qMapEnd  = end - off1 - mapOffset 
      mapdata = query.bigWig(map, chrom, qMapStart, qMapEnd)
      mapresult[i] = abs(sum(mapdata[,3]))
    }
  }
  return(cbind(cbind(result), cbind(mapresult)))
}

getCounts.gb <- function(wig.p, wig.m, genes, off1) {
  N = dim(genes)[1]
  plusStrand = genes[,2] == '+'
  
  result = vector(mode="integer", length=N)
  
  for (i in 1:N) {
    chrom = as.character(genes[i, 1])
    start = as.integer(genes[i, 3])
    end = as.integer(genes[i, 4])
    
    if (plusStrand[i]) {
      wig = wig.p
      
      qStart = start + off1
      qEnd = end
      
      if (qStart > qEnd) {
        result[i] = NA
        next
      }
      
      data = query.bigWig(wig, chrom, qStart, qEnd)
      
      if (!is.null(data))
        result[i] = abs(sum(data[,3]))
    } else {
      wig = wig.m
      if (is.null(wig.m))
        wig = wig.p
      
      qStart = start
      qEnd = end - off1
      
      if (qStart > qEnd) {
        result[i] = NA
        next
      }
      
      data = query.bigWig(wig, chrom, qStart, qEnd)
      
      if (!is.null(data))
        result[i] = abs(sum(data[,3]))
    }
  }
  
  return(cbind(result))
}

getCounts.gb.mappability <- function(map, genes, off1) {
  N = dim(genes)[1]
  plusStrand = genes[,2] == '+'
  
  result = vector(mode="integer", length=N)
  
  for (i in 1:N) {
    chrom = as.character(genes[i, 1])
    start = as.integer(genes[i, 3])
    end = as.integer(genes[i, 4])
    
    if (plusStrand[i]) {
      
      qStart = start + off1 + 36
      qEnd = end + 36
      
      if (qStart > qEnd) {
        result[i] = NA
        next
      }
      
      data = query.bigWig(map, chrom, qStart, qEnd)
      
      if (!is.null(data))
        result[i] = abs(sum(data[,3]))
    } else {
      
      qStart = start
      qEnd = end - off1
      
      if (qStart > qEnd) {
        result[i] = NA
        next
      }
      
      data = query.bigWig(map, chrom, qStart, qEnd)
      
      if (!is.null(data))
        result[i] = abs(sum(data[,3]))
    }
  }
  
  return(cbind(result))
}

getCounts.sliding.term <- function(wig.p, wig.m, map, genes, off1 = 500, windowsize = 100) {
  N = dim(genes)[1]
  plusStrand = genes[,2] == '+'
  
  result = vector(mode="integer", length=N)
  mapresult = vector(mode="integer", length=N)
  scanlen = off1 - windowsize
  for (i in 1:N) {
    chrom = as.character(genes[i, 1])
    start = as.integer(genes[i, 3])
    end = as.integer(genes[i, 4])
    
    if (plusStrand[i]) {
      wig = wig.p
      
      scanResults = c()
      for (n in c(1:scanlen)){
        
        qStart = end + n
        qEnd = end + n + windowsize
        
        data = query.bigWig(wig, chrom, qStart, qEnd)
        scanResults = append(scanResults, abs(sum(data[,3])))
        
      }
      mapOffset = which.max(scanResults)  # grabs the index of the first instance of the maximum value ( = bp offset during scanning)
      result[i]  = max(scanResults)
      qMapStart = end + mapOffset + 36
      qMapEnd  = end + mapOffset + windowsize + 36
      mapdata = query.bigWig(map, chrom, qMapStart, qMapEnd)
      mapresult[i] = abs(sum(mapdata[,3]))
    } 
    else {
      wig = wig.m
      scanResults = c()
      if (is.null(wig.m))
        wig = wig.p
      
      for (n in c(1: scanlen)){
        qStart = start - n  - windowsize
        qEnd = start - n
        
        data = query.bigWig(wig, chrom, qStart, qEnd)
        scanResults = append(scanResults, abs(sum(data[,3])))
      }
      mapOffset = which.max(scanResults)  # grabs the index of the first instance of the maximum value ( = bp offset during scanning)
      result[i]  = max(scanResults)
      qMapStart = start - mapOffset - windowsize
      qMapEnd  = start - mapOffset 
      mapdata = query.bigWig(map, chrom, qMapStart, qMapEnd)
      mapresult[i] = abs(sum(mapdata[,3]))
    }
  }
  return(cbind(cbind(result), cbind(mapresult)))
}

CountTable <- function(wig.p, wig.m, map, genes) {
  MappedReads = total.reads(wig.p,wig.m)
  
  pr.data = getCounts.sliding.pr(wig.p, wig.m, map, genes)
  pr.counts = pr.data[,1]
  pr.map = pr.data[,2]
  pr.density = pr.counts / pr.map
  pr.density[pr.map == 0] = NA
  pr.densityRPKM = (pr.density*1000)/(MappedReads/1000000)
  
  gb.counts = getCounts.gb(wig.p, wig.m, genes, 200)
  gb.map = getCounts.gb.mappability(map, genes, 200)
  gb.density = gb.counts[,1] / gb.map[,1]
  gb.density[gb.map[,1] == 0] = NA
  gb.densityRPKM = (gb.density*1000)/(MappedReads/1000000)
  
  term.data = getCounts.sliding.term(wig.p, wig.m, map, genes)
  term.counts = term.data[,1]
  term.map = term.data[,2]
  term.density = term.counts / term.map
  term.density[term.map == 0] = NA
  term.densityRPKM = (term.density*1000)/(MappedReads/1000000)
  
  PI = pr.density/gb.density
  TI = term.density/gb.density
  
  geneNames = (genes[,5])
  chro = genes[,1]
  strand = genes[,2]
  start = as.numeric(genes[,3])
  end = as.numeric(genes[,4])
  result = cbind(data.frame(chro), geneNames, strand, start, end, pr.counts, pr.map, pr.density, pr.densityRPKM, gb.counts, gb.map, gb.density, gb.densityRPKM, term.counts, term.map, term.density, term.densityRPKM, PI, TI)
  colnames(result) <- c( "chr","genename",'strand','start','end', "pr_counts",  "pr_map", "pr_density", "pr_densityRPKM",  "gb_counts",  "gb_map", "gb_density", "gb_densityRPKM", "term_counts",  "term_map", "term_density", "term_densityRPKM", "PauseIndex", "TerminationIndex")
  #rownames(result)<-geneNames
  #result = as.data.frame(result)
  return(result)
}

fe.test= function(data, prReads, prmappL, gbReads, gbmappL, colName='')
{
  newtable = data.frame(data)
  rows = nrow(data)
  for(ii in 1:nrow(data))
  {
    fmat = matrix(c((data[ii, prReads]),(data[ii, prmappL]),(data[ii, gbReads]),(data[ii, gbmappL])), nrow=2) 
    
    f = fisher.test(fmat,alternative='greater')#it has to be a one-tailed fisher-exact test, or it will give you the de-enriched genes as well
    newtable[ii, colName] = p.adjust(f$p.value, method = 'bonferroni', n = rows)
  }
  return(newtable)
}

calcPI.sig= function(pauseIndexOutput){
  pauseIndexOutput = as.data.frame(pauseIndexOutput)
  
  pauseIndexOutput$EP=as.integer((pauseIndexOutput[,"pr_counts"]+pauseIndexOutput[,'pr_map'])*((pauseIndexOutput[,'gb_counts']/(pauseIndexOutput[,'gb_counts']+pauseIndexOutput[,'gb_map']))))
  
  pauseIndexOutput=pauseIndexOutput[!is.na(pauseIndexOutput$EP),] #remove NAs generated by previous calc (same used by fisher's exact test)
  
  pauseIndexOutput$EG=as.integer(pauseIndexOutput[,"pr_counts"]+pauseIndexOutput[,'pr_map']-pauseIndexOutput$EP)
  
  pauseIndexOutput=pauseIndexOutput[!is.na(pauseIndexOutput$EP),]
  
  new_data = fe.test(pauseIndexOutput, 'pr_counts', 'pr_map', 'gb_counts', 'gb_map', "Fexact_Pval")
  
  return(new_data)
}

#############################################################################################################
## Previous functions for calculating Pausing Index:  No sliding window.

getCounts.pr <- function(wig.p, wig.m, genes, off1, off2) {
  N = dim(genes)[1]
  plusStrand = genes[,2] == '+'
  
  result = vector(mode="integer", length=N)
  
  for (i in 1:N) {
    chrom = as.character(genes[i, 1])
    start = as.integer(genes[i, 3])
    end = as.integer(genes[i, 4])
    
    if (plusStrand[i]) {
      wig = wig.p
      
      qStart = start + off1
      qEnd = start + off2
      
      if (qStart > qEnd) {
        result[i] = NA
        next
      }
      
      data = query.bigWig(wig, chrom, qStart, qEnd)
      
      if (!is.null(data))
        result[i] = abs(sum(data[,3]))
    } else {
      wig = wig.m
      if (is.null(wig.m))
        wig = wig.p
      
      qStart = end - off2
      qEnd = end - off1
      
      if (qStart > qEnd) {
        result[i] = NA
        next
      }
      
      data = query.bigWig(wig, chrom, qStart, qEnd)
      
      if (!is.null(data))
        result[i] = abs(sum(data[,3]))
    }
  }
  
  return(cbind(result))
}

getCounts.pr.mappability <- function(map, genes, off1, off2) {
  N = dim(genes)[1]
  plusStrand = genes[,2] == '+'
  
  result = vector(mode="integer", length=N)
  
  for (i in 1:N) {
    chrom = as.character(genes[i, 1])
    start = as.integer(genes[i, 3])
    end = as.integer(genes[i, 4])
    
    if (plusStrand[i]) {
      
      qStart = start + off1 + 36
      qEnd = start + off2 + 36 
      
      if (qStart > qEnd) {
        result[i] = NA
        next
      }
      
      data = query.bigWig(map, chrom, qStart, qEnd)
      
      if (!is.null(data))
        result[i] = abs(sum(data[,3]))
    } else {
      
      qStart = end - off2
      qEnd = end - off1
      
      if (qStart > qEnd) {
        result[i] = NA
        next
      }
      
      data = query.bigWig(map, chrom, qStart, qEnd)
      
      if (!is.null(data))
        result[i] = abs(sum(data[,3]))
    }
  }
  
  return(cbind(result))
}

getCounts.term <- function(wig.p, wig.m, genes, off1, off2) {
  N = dim(genes)[1]
  plusStrand = genes[,2] == '+'
  
  result = vector(mode="integer", length=N)
  
  for (i in 1:N) {
    chrom = as.character(genes[i, 1])
    start = as.integer(genes[i, 3])
    end = as.integer(genes[i, 4])
    
    if (plusStrand[i]) {
      wig = wig.p
      
      qStart = end + off1
      qEnd = end + off2
      
      if (qStart > qEnd) {
        result[i] = NA
        next
      }
      
      data = query.bigWig(wig, chrom, qStart, qEnd)
      
      if (!is.null(data))
        result[i] = abs(sum(data[,3]))
    } else {
      wig = wig.m
      if (is.null(wig.m))
        wig = wig.p
      
      qStart = start - off2
      qEnd = start - off1
      
      if (qStart > qEnd) {
        result[i] = NA
        next
      }
      
      data = query.bigWig(wig, chrom, qStart, qEnd)
      
      if (!is.null(data))
        result[i] = abs(sum(data[,3]))
    }
  }
  
  return(cbind(result))
}

getCounts.term.mappability <- function(map, genes, off1, off2) {
  N = dim(genes)[1]
  plusStrand = genes[,2] == '+'
  
  result = vector(mode="integer", length=N)
  
  for (i in 1:N) {
    chrom = as.character(genes[i, 1])
    start = as.integer(genes[i, 3])
    end = as.integer(genes[i, 4])
    
    if (plusStrand[i]) {
      
      qStart = end + off1 + 36
      qEnd = end + off2 + 36 
      
      if (qStart > qEnd) {
        result[i] = NA
        next
      }
      
      data = query.bigWig(map, chrom, qStart, qEnd)
      
      if (!is.null(data))
        result[i] = abs(sum(data[,3]))
    } else {
      
      qStart = start - off2
      qEnd = start - off1
      
      if (qStart > qEnd) {
        result[i] = NA
        next
      }
      
      data = query.bigWig(map, chrom, qStart, qEnd)
      
      if (!is.null(data))
        result[i] = abs(sum(data[,3]))
    }
  }
  
  return(cbind(result))
}

approx.ratio.CI <- function(x1, x2, alpha=0.05) {
  t = qnorm(1 - alpha/2)
  n = x1 + x2
  zp = (t^2/2 + x1 + 1/2)^2 - ((n + t^2)/n) * (x1 + 1/2)^2
  zn = (t^2/2 + x1 - 1/2)^2 - ((n + t^2)/n) * (x1 - 1/2)^2
  
  a = (t^2/2 + x1 + 1/2 + sqrt(zp)) / (n + t^2/2 - x1 - 1/2 - sqrt(zp))
  b = (t^2/2 + x1 - 1/2 - sqrt(zn)) / (n + t^2/2 - x1 + 1/2 + sqrt(zn))
  
  return(c(b, a))
}

approx.ratios.CI <- function(num.counts, denom.counts, alpha=0.05) {
  stopifnot(length(num.counts) == length(denom.counts))
  N = length(num.counts)
  
  result = matrix(data=0, nrow=N, ncol=2)
  
  for (i in 1:N)
    result[i, ] = approx.ratio.CI(num.counts[i], denom.counts[i], alpha)
  
  return(result)
}

load.wigset <- function(row) {
  file = wig.table[row, 1]
  wig.p = NULL
  if (file != "")
    wig.p = load.bigWig(paste(basepath, file, sep=''))
  file = wig.table[row, 2]
  wig.m = NULL
  if (file != "")
    wig.m = load.bigWig(paste(basepath, file, sep=''))
  return(list(wig.p, wig.m, wig.table[row, 3]))
}

unload.wigset <- function(set) {
  if (!is.null(set[[1]]))
    unload.bigWig(set[[1]])
  if (!is.null(set[[2]]))
    unload.bigWig(set[[2]])
}

## function for calculating fisher's t-test pVal of pausing index 

DataFor.pause.index <- function(wig.p, wig.m, map, genes, name, alpha=0.05) {
  MappedReads = total.reads(wig.p,wig.m)
  pr.counts = getCounts.pr(wig.p, wig.m, genes, -50, 150)
  gb.counts = getCounts.gb(wig.p, wig.m, genes, 200)
  
  pr.map = getCounts.pr.mappability(map, genes, -50, 150)
  gb.map = getCounts.gb.mappability(map, genes, 200) 
  
  pr.density = pr.counts[,1] / pr.map[,1]
  pr.density[pr.map[,1] == 0] = NA
  gb.density = gb.counts[,1] / gb.map[,1]
  gb.density[gb.map[,1] == 0] = NA
  pr.densityRPKM = (pr.density*1000)/(MappedReads/1000000)
  gb.densityRPKM = (gb.density*1000)/(MappedReads/1000000)
  pause.index = pr.counts[,1] / gb.counts[,1]
  pause.index[gb.counts[,1] == 0] = NA
  
  cis = approx.ratios.CI(pr.counts[,1], gb.counts[,1], alpha)
  
  pause.index.scaled = pause.index * (gb.map[,1]/pr.map[,1])
  pause.index.scaled[gb.map[,1] == 0] = NA
  pause.index.scaled[pr.map[,1] == 0] = NA
  ci.low.scaled = cis[,1] * (gb.map[,1]/pr.map[,1])
  ci.low.scaled[gb.map[,1] == 0] = NA
  ci.low.scaled[pr.map[,1] == 0] = NA
  ci.high.scaled = cis[,2] * (gb.map[,1]/pr.map[,1])
  ci.high.scaled[gb.map[,1] == 0] = NA
  ci.high.scaled[pr.map[,1] == 0] = NA
  geneNames = as.character(genes[,5])
  result = cbind(as.numeric(pr.counts[,1]), gb.counts[,1], pr.map [,1], gb.map [,1], pr.density, pr.densityRPKM, gb.density, gb.densityRPKM, cis[,1], cis[,2], pause.index.scaled, ci.low.scaled, ci.high.scaled)
  colnames(result) <- c("pr_counts", "gb_counts", "pr_map", "gb_map", "pr_density","pr_densityRPKM", "gb_density", "gb_densityRPKM", "pause_index_CI_low", "pause_index_CI_high", "pause_index_scaled", "pause_index_CI_low_scaled", "pause_index_CI_high_scaled")
  rownames(result)<-geneNames
  return(result)
}

pause.index <- function(wig.p, wig.m, map, genes, name, alpha=0.05) {
  MappedReads = total.reads(wig.p,wig.m)
  pr.counts = getCounts.pr(wig.p, wig.m, genes, 0, 100) # because of remapped TSS, not taking any upstream regions (off1 = 0)
  gb.counts = getCounts.gb(wig.p, wig.m, genes, 200)
  term.counts = getCounts.term(wig.p, wig.m, genes, 0, 300)
  
  pr.map = getCounts.pr.mappability(map, genes, 0, 100)
  gb.map = getCounts.gb.mappability(map, genes, 200) 
  term.map = getCounts.term.mappability(map, genes, 0, 300)
  
  pr.density = pr.counts[,1] / pr.map[,1]
  pr.density[pr.map[,1] == 0] = NA
  gb.density = gb.counts[,1] / gb.map[,1]
  gb.density[gb.map[,1] == 0] = NA
  term.density = term.counts[,1] / term.map[,1]
  pr.densityRPKM = (pr.density*1000)/(MappedReads/1000000)
  gb.densityRPKM = (gb.density*1000)/(MappedReads/1000000)
  term.densityRPKM = (term.density*1000)/(MappedReads/1000000)
  pause.index = pr.counts[,1] / gb.counts[,1]
  pause.index[gb.counts[,1] == 0] = NA
  
  cis = approx.ratios.CI(pr.counts[,1], gb.counts[,1], alpha)
  
  pause.index.scaled = pause.index * (gb.map[,1]/pr.map[,1])
  pause.index.scaled[gb.map[,1] == 0] = NA
  pause.index.scaled[pr.map[,1] == 0] = NA
  ci.low.scaled = cis[,1] * (gb.map[,1]/pr.map[,1])
  ci.low.scaled[gb.map[,1] == 0] = NA
  ci.low.scaled[pr.map[,1] == 0] = NA
  ci.high.scaled = cis[,2] * (gb.map[,1]/pr.map[,1])
  ci.high.scaled[gb.map[,1] == 0] = NA
  ci.high.scaled[pr.map[,1] == 0] = NA
  geneNames = as.character(genes[,5])
  result = cbind(as.numeric(pr.counts[,1]), gb.counts[,1], term.counts[,1], pr.map [,1], gb.map [,1], term.map[,1], pr.density, pr.densityRPKM, gb.density, gb.densityRPKM, cis[,1], cis[,2], pause.index.scaled, ci.low.scaled, ci.high.scaled)
  colnames(result) <- c("pr_counts", "gb_counts", "term_counts", "pr_map", "gb_map", "term_map", "pr_density","pr_densityRPKM", "gb_density", "gb_densityRPKM", "pause_index_CI_low", "pause_index_CI_high", "pause_index_scaled", "pause_index_CI_low_scaled", "pause_index_CI_high_scaled")
  rownames(result)<-geneNames
  return(result)
}

#############################################################################################################
##  Heat map sorting 

downstream.Gene.sum = function(heatMatrix){
  firstbin = ncol(heatMatrix)/2 ## this grabs the middle bin (contains TSS)
  lastbin = ncol(heatMatrix)
  result = c()
  for( i in c(1: nrow(heatMatrix))){
    result[i]  = sum(heatMatrix[i,c(firstbin : lastbin)])
  }
  return(result)
}

SortExprQuantbyLength = function(genelist){
  genelist$genelength = genelist[,4]-genelist[,3]
  #genelist$score = 0
  sortedGL = genelist[order(genelist[,9], -genelist[,8]),] #sorting first by  expression quartile then increasing length
  return(sortedGL[,c(2,3,4,1,5,6,9,7,8)])
}

geneExpressionQuantiles = function(genelist, PI_outFinal, expThresh){
  expressedGenes = PI_outFinal[PI_outFinal$gb_densityRPKM > expThresh,]
  expressedGenes$ID = row.names(expressedGenes)
  expressedGL = genelist[genelist[,4]%in% row.names(expressedGenes),]
  gList_GBdens = merge(expressedGL, expressedGenes[,c(8,17)], by.x = 'V1', by.y = 'ID')
  #head(gList_GBdens)
  quants = summary(gList_GBdens[,7])
  quantVals = c()
  N = length(gList_GBdens[,7])
  for (i in 1:N){
    if (gList_GBdens[i,7]< quants[2]){
      quantVals = c(quantVals,1)
    }
    if (gList_GBdens[i,7]>= quants[2]){
      if (gList_GBdens[i,7]< quants[3]){
        quantVals = c(quantVals,2)
      }
    }
    if (gList_GBdens[i,7]>= quants[3]){
      if (gList_GBdens[i,7]< quants[5]){
        quantVals = c(quantVals,3)
      }
    }
    if (gList_GBdens[i,7]>= quants[5]){ 
      quantVals = c(quantVals,4)
    }
  }
  gList_GBdens$quantile = quantVals
  return(gList_GBdens)
}


#############################################################################################################
##  Count Table retrieving all reads within a set of coordinates

chrSwap = function(Glist){
  Glist[,1] = as.character(Glist[,1])
  Glist$V1[Glist$V1 == "chr01"] = "chrI"
  Glist$V1[Glist$V1 == "chr02"] = "chrII"
  Glist$V1[Glist$V1 == "chr03"] = "chrIII"
  return (Glist)  
}

getCounts <- function(wig.p, wig.m, genes) {
  N = dim(genes)[1]
  plusStrand = genes[,2] == '+'
  
  result = vector(mode="integer", length=N)
  
  for (i in 1:N) {
    chrom = as.character(genes[i, 1])
    start = as.integer(genes[i, 3])
    end = as.integer(genes[i, 4])
    
    if (plusStrand[i]) {
      wig = wig.p
      
      qStart = start
      qEnd = end
      
      if (qStart > qEnd) {
        result[i] = NA
        next
      }
      
      data = query.bigWig(wig, chrom, qStart, qEnd)
      
      if (!is.null(data))
        result[i] = abs(sum(data[,3]))
    } else {
      wig = wig.m
      if (is.null(wig.m))
        wig = wig.p
      
      qStart = start
      qEnd = end
      
      if (qStart > qEnd) {
        result[i] = NA
        next
      }
      
      data = query.bigWig(wig, chrom, qStart, qEnd)
      
      if (!is.null(data))
        result[i] = abs(sum(data[,3]))
    }
  }
  
  return(cbind(result))
}


getCounts.mappability <- function(map, genes) {
  N = dim(genes)[1]
  plusStrand = genes[,2] == '+'
  
  result = vector(mode="integer", length=N)
  
  for (i in 1:N) {
    chrom = as.character(genes[i, 1])
    start = as.integer(genes[i, 3])
    end = as.integer(genes[i, 4])
    
    if (plusStrand[i]) {
      
      qStart = start +36
      qEnd = end + 36
      
      if (qStart > qEnd) {
        result[i] = NA
        next
      }
      
      data = query.bigWig(map, chrom, qStart, qEnd)
      
      if (!is.null(data))
        result[i] = abs(sum(data[,3]))
    } else {
      
      qStart = start
      qEnd = end
      
      if (qStart > qEnd) {
        result[i] = NA
        next
      }
      
      data = query.bigWig(map, chrom, qStart, qEnd)
      
      if (!is.null(data))
        result[i] = abs(sum(data[,3]))
    }
  }
  
  return(cbind(result))
}


regionCountTable <- function(wig.p, wig.m, map, genes) {
  MappedReads = total.reads(wig.p,wig.m)
  region.counts = getCounts(wig.p, wig.m, genes)
  region.map = getCounts.mappability(map, genes)
  region.density = region.counts[,1] / region.map[,1]
  region.density[region.map[,1] == 0] = NA
  region.densityRPKM = (region.density*1000)/(MappedReads/1000000)
  geneNames = (genes[,5])
  chro = genes[,1]
  strand = genes[,2]
  start = as.numeric(genes[,3])
  end = as.numeric(genes[,4])
  result = cbind(data.frame(chro), geneNames, strand, start, end, region.counts[,1], region.map [,1], region.density, region.densityRPKM)
  colnames(result) <- c( "chr","genename",'strand','start','end', "region_counts",  "region_map", "region_density", "region_densityRPKM")
  #rownames(result)<-geneNames
  #result = as.data.frame(result)
  return(result)
}


############################################################################################################
# Functions for the "Area Under the Curve" analysis
# These functions scale all genes to 400 bins (~2bp bins; if unscaled), the first and last 200 bases of genes are unscaled, while everything inbetween is scaled to 200 bins
# are used to generate a cumulative density plot for each gene
# The area underneath the cumulative density plot describes the shape of the Pol II distribution on each gene. 
# you can then compare the shape (AUC) in WT and mutant strains. 

mround.AUC = function(x, base = 100){ 
  base*floor(x/base)}  # base = 30 was chosen so that there will be a total of 60 windows (30 on each side of the midpoint of scaled window)

# scale_params assesses the length of your gene, returning a buffer length (1/2 the size of the region to be scaled), a step (window size for scaled region), and a scale facor (windowsize/10)
scale_params.AUC = function(start, end, step=2, buffer=200){
  length = end-start
  if (length < 800){
    buffer =  buffer
    step  = step
    scaleFact = 1
  }
  else{
    length = ((end-start)-400)/2
    buffer = mround.AUC(length, base = 100)
    step = (buffer/100)
    scaleFact = step/2 # scale factor is divided by your scaled window count to adjust as if it were a 10bp window (i.e. window size = 5bp will have a scale factor of 2)
  }
  return(c(buffer,step,scaleFact))
}

collect.many.scaled.AUC = function (bed, bigWig.plus, bigWig.minus, FLbuffer = 200, step=2, 
                                do.sum = T) 
{
  #windowSize = ((2*flank) + (4*buffer))%/%step
  midPoint = (bed[, 2] + bed[, 3])/2
  TSS = bed[, 2]
  END = bed[,3]
  
  prstart = (TSS)
  prend = (TSS + FLbuffer)
  
  finstart = (END - FLbuffer)
  finend = (END)
  #windowSize = (2*flank + 4*FLbuffer)/step
  windowSize = 400
  N = dim(bed)[1]
  result = matrix(nrow = N, ncol = windowSize)
  strands = NULL
  if (dim(bed)[2] >= 6) 
    strands = as.character(bed[, 6])
  else strands = rep("+", N)
  for (i in 1:N) {
    chrom = as.character(bed[i, 1])
    strand = strands[i]
    if (strand == "+") {
      bigWig = bigWig.plus
      midscaleParams = scale_params.AUC(bed[i,2], bed[i,3], step, buffer = 200)
      midbuf = midscaleParams[1]
      midstep = midscaleParams[2]
      midScaleFact = midscaleParams[3]
      midstart = (midPoint - midbuf)
      midend = (midPoint + midbuf)
      prrow = collect.counts(bigWig, chrom, prstart[i], prend[i], 
                             step, do.sum)
      midrow =  collect.counts(bigWig, chrom, midstart[i], midend[i], 
                               midstep, do.sum)/midScaleFact
      endrow = collect.counts(bigWig, chrom, finstart[i], finend[i], 
                              step, do.sum)
      row = c(prrow, midrow, endrow)
      if (length(row) == 400){
        result[i, ] = abs(row)}
    }
    else {
      bigWig = bigWig.minus
      if (is.null(bigWig.minus)) 
        bigWig = bigWig.plus
      midscaleParams = scale_params.AUC(bed[i,2], bed[i,3], step, buffer = 200)
      midbuf = midscaleParams[1]
      midstep = midscaleParams[2]
      midScaleFact = midscaleParams[3]
      midstart = (midPoint - midbuf)
      midend = (midPoint + midbuf)
      prrow = collect.counts(bigWig, chrom, prstart[i], prend[i], 
                             step, do.sum)
      midrow =  collect.counts(bigWig, chrom, midstart[i], midend[i], 
                               midstep, do.sum)/midScaleFact
      endrow = collect.counts(bigWig, chrom, finstart[i], finend[i], 
                              step, do.sum)
      row = c(prrow, midrow, endrow)
      if (length(row) == 400){
        result[i, ] = abs(rev(row))}
    }
  }
  rownames(result) = bed[,4]
  result[complete.cases(result), ]
}


cum_sum_row = function(x){
  N = length(x)
  tot = sum(x)
  sumx = 0
  result = vector(length = N)
  for (i in 1:N){
    sumx = sumx + x[i]
    prop = sumx/tot
    result[i] = prop
  }
  return(result)
}

plot_shapes = function(df, startCol, endCol){
  for (i in startCol:endCol){ 
    plot(df[,i], type = "l", main = colnames(df[i]), 
         lty = 1, lwd = 3, col = 'darkgreen', 
         xlab = "Bin", ylab = "Cumulative PolII Density", 
         cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)}
}




############################################################################################################
# Functions for the "Area Under the Curve" analysis
# These functions scale all genes to 200 bins and are used to generate a cumulative density plot for each gene
# The area underneath the cumulative density plot describes the shape of the Pol II distribution on each gene. 
# you can then compare the shape (AUC) in WT and mutant strains. 
# The difference between these functions and the ones above isi that these scale the entire gene body.  


mround_full = function(x, base = 100){ 
  base*ceiling(x/base)}

scale_params_full = function(start, end, step=10){
  length = (end-start)/2
  buffer = mround_full(length, base = 100)
  step = (buffer/100) 
  return(c(buffer, step))
}

collect.many.scaled.full = function (bed, bigWig.plus, bigWig.minus, step=10, 
                                     do.sum = T) 
{
  midPoint = (bed[, 2] + bed[, 3])/2
  TSS = bed[, 2]
  END = bed[,3]
  #windowSize = (2*flank + 4*FLbuffer)/step
  windowSize = 200
  N = dim(bed)[1]
  result = matrix(nrow = N, ncol = windowSize)
  strands = NULL
  if (dim(bed)[2] >= 6) 
    strands = as.character(bed[, 6])
  else strands = rep("+", N)
  for (i in 1:N) {
    chrom = as.character(bed[i, 1])
    strand = strands[i]
    if (strand == "+") {
      bigWig = bigWig.plus
      midscaleParams = scale_params_full(bed[i,2], bed[i,3], step)
      midbuf = midscaleParams[1]  # buffer from midscaleParams
      midstep = midscaleParams[2] # step (note these will be of varying sizes depending on gene length, but for these purposes scaling is not necessary)
      #midScaleFact = midscaleParams[3]
      midstart = (midPoint - midbuf)
      midend = (midPoint + midbuf)
      midrow =  collect.counts(bigWig, chrom, midstart[i], midend[i], 
                               midstep, do.sum) #/midScaleFact
      row = midrow
      #if (length(row) == 200){
      result[i, ] = abs(row)}
    #}
    else {
      bigWig = bigWig.minus
      if (is.null(bigWig.minus)) 
        bigWig = bigWig.plus
      midscaleParams = scale_params_full(bed[i,2], bed[i,3], step)
      midbuf = midscaleParams[1]
      midstep = midscaleParams[2]
      #midScaleFact = midscaleParams[3]
      midstart = (midPoint - midbuf)
      midend = (midPoint + midbuf)
      midrow =  collect.counts(bigWig, chrom, midstart[i], midend[i], 
                               midstep, do.sum) #/midScaleFact
      row = midrow
      #if (length(row) == 200){
      result[i, ] = abs(rev(row))}
    #}
  }
  rownames(result) = bed[,4]
  result[complete.cases(result), ]
}


#####################################################################################
# Function designed to identify genes whose Pol II density profiles are influenced by upstream genes
# Purpose:  Calculate ratios of upstream region (-300 to TSS) and downstream region ((TSS+250) to (TSS+550))
# Function is to be applied to rows of a bed file containg gene information
# requires bigWig package


Upstream.Ratio = function(x, pbw, mbw, map){ # note X is the row (i.e. gene) to apply this function to
  chro = x[1]
  if (x[6] == "+"){
    bw = pbw
    upstart = as.numeric(x[2]) - 300
    upend = as.numeric(x[2])
    downstart = as.numeric(x[2]) + 250
    downend = as.numeric(x[2]) + 550
    upmap = sum(query.bigWig(map, chro, upstart, upend)[,3])
    upQuery = query.bigWig(bw, chro, upstart, upend)
    upDens = sum(upQuery[,3])/upmap
    downmap = sum(query.bigWig(map, chro, downstart, downend)[,3])
    downQuery = query.bigWig(bw, chro, downstart, downend)
    downDens = sum(downQuery[,3])/downmap
    updownRatio = (upDens/downDens)
  }
  else{
    bw = mbw
    upstart = as.numeric(x[3])
    upend = as.numeric(x[3]) + 300
    downstart = as.numeric(x[3]) - 550
    downend = as.numeric(x[3]) - 250
    upmap = sum(query.bigWig(map, chro, upstart, upend)[,3])
    upQuery = query.bigWig(bw, chro, upstart, upend)
    upDens = sum(upQuery[,3])/upmap
    downmap = sum(query.bigWig(map, chro, downstart, downend)[,3])
    downQuery = query.bigWig(bw, chro, downstart, downend)
    downDens = sum(downQuery[,3])/downmap
    updownRatio = (upDens/downDens)
  }
  return(updownRatio)
}


### example of use: 

#bed = read.table(file = "/Volumes/SEAGATE_EXP/Greg_YeastData/GeneLists_final/Pombe_08_23-15/observedTSS/pombe.ASM294v1.16.cleaned_sorted_PROcapObservedTSS_noOverlap.bed", sep = "\t")
#bed2 = read.table("/Volumes/SEAGATE_EXP/Greg_YeastData/GeneLists_final/Pombe_08_23-15/observedTSS/pombe.ASM294v1.16.cleaned_sorted.bed_PROcapObservedTSS_1000bpSep_coding.bed", stringsAsFactors = F)
#pbw = load.bigWig("/Volumes/SEAGATE_EXP/Greg_YeastData/PRO_Seq/Pombe/bigWig/3870_7157_12395_C53ARACXX_SP_WT_972h-COMBINED_REPS_GGCTAC_R1.fastq_sorted.bed_plus.bw")
#mbw = load.bigWig("/Volumes/SEAGATE_EXP/Greg_YeastData/PRO_Seq/Pombe/bigWig/3870_7157_12395_C53ARACXX_SP_WT_972h-COMBINED_REPS_GGCTAC_R1.fastq_sorted.bed_minus.bw")
#mappability = load.bigWig("/Volumes/SEAGATE_EXP/Greg_YeastData/PRO_Seq/Pombe/bedgraphs_Kbro/SPombeMappable_SacCerNonOverLap_expanded_ChrFix.bw")

#bed$testCalc = apply(bed, 1, Upstream.Ratio, pbw = pbw, mbw = mbw, map = mappability)

###########################################################################################################################
# use the average background region_density as lambda for a poisson distribution 
# Calculate P(at least GBdensity) for each gene based on the above poisson distribution
# Active genes have P < 0.01

# ActiveGeneProb should be applied to dataframe output from pause.index.  It will return the prob of observing a gene's GB density based on poisson distribution with a supplied lambda 
# note, lambda should be background rate/bp
ActiveGeneProb = function(row, lambda){
  GBdens = row["gb_density"]
  GBmap = row["gb_map"]
  P_active = ppois(GBdens, lambda * GBmap) ## calcs the prob of # of counts over a mappable region based on background rate 
  return(P_active)
}


############################################################################################################
# Functions for Shape analysis of promoters, gene bodies, and termination 
collect.many.pr = function (bed, bigWig.plus, bigWig.minus, halfWindow, step, 
                            do.sum = TRUE) 
{
  windowSize = (2 * halfWindow)%/%step
  midPoint = bed[, 2] + halfWindow
  neg.idxs = bed[, 6] == "-"
  midPoint[neg.idxs] = bed[neg.idxs, 3] - halfWindow - 1
  start = (midPoint - halfWindow)
  end = start + windowSize * step
  N = dim(bed)[1]
  result = matrix(nrow = N, ncol = windowSize)
  strands = NULL
  if (dim(bed)[2] >= 6) 
    strands = as.character(bed[, 6])
  else strands = rep("+", N)
  for (i in 1:N) {
    chrom = as.character(bed[i, 1])
    strand = strands[i]
    if (strand == "+") {
      bigWig = bigWig.plus
      row = collect.counts(bigWig, chrom, start[i], end[i], 
                           step, do.sum)
      result[i, ] = abs(row)
    }
    else {
      bigWig = bigWig.minus
      if (is.null(bigWig.minus)) 
        bigWig = bigWig.plus
      row = collect.counts(bigWig, chrom, start[i], end[i], 
                           step, do.sum)
      result[i, ] = abs(rev(row))
    }
  }
  rownames(result) = bed[,4]
  result[complete.cases(result), ]
}


collect.many.CPS = function (bed, bigWig.plus, bigWig.minus, halfWindow, step, 
                            do.sum = TRUE) 
{
  windowSize = (2 * halfWindow)%/%step
  midPoint = bed[, 3] + halfWindow
  neg.idxs = bed[, 6] == "-"
  midPoint[neg.idxs] = bed[neg.idxs, 2] - halfWindow - 1
  start = (midPoint - halfWindow)
  end = start + windowSize * step
  N = dim(bed)[1]
  result = matrix(nrow = N, ncol = windowSize)
  strands = NULL
  if (dim(bed)[2] >= 6) 
    strands = as.character(bed[, 6])
  else strands = rep("+", N)
  for (i in 1:N) {
    chrom = as.character(bed[i, 1])
    strand = strands[i]
    if (strand == "+") {
      bigWig = bigWig.plus
      row = collect.counts(bigWig, chrom, start[i], end[i], 
                           step, do.sum)
      result[i, ] = abs(row)
    }
    else {
      bigWig = bigWig.minus
      if (is.null(bigWig.minus)) 
        bigWig = bigWig.plus
      row = collect.counts(bigWig, chrom, start[i], end[i], 
                           step, do.sum)
      result[i, ] = abs(rev(row))
    }
  }
  rownames(result) = bed[,4]
  result[complete.cases(result), ]
}

##########################################################################################################
### ChIP-Chip functions
##########################################################################################################
## Functions for dealing with pos and neg values found in ChIP-chip bigWigs 
scale_params.ChIPchip = function(start, end, step=50, buffer=300){
  length = end-start
  if (length <1000){
    buffer =  buffer
    step  = step
    scaleFact = 1
  }
  else{
    length = ((end-start)-600)/2
    buffer = mround(length, base = step)
    step = (buffer/50)
    scaleFact = step/50 # scale factor is divided by your scaled window count to adjust as if it were a 10bp window (i.e. window size = 5bp will have a scale factor of 2)
  }
  return(c(buffer,step,scaleFact))
}

collect.many.ChIPchip <- function (bed, bigWig.plus, bigWig.minus, halfWindow, step, at.TSS = FALSE, 
                                   do.sum = FALSE) 
{
  windowSize = (2 * halfWindow)%/%step
  midPoint = (bed[, 2] + bed[, 3])/2
  if (at.TSS && dim(bed)[2] >= 6) {
    midPoint = bed[, 2]
    neg.idxs = bed[, 6] == "-"
    midPoint[neg.idxs] = bed[neg.idxs, 3] - 1
  }
  start = (midPoint - halfWindow)
  end = start + windowSize * step
  N = dim(bed)[1]
  result = matrix(nrow = N, ncol = windowSize)
  strands = NULL
  if (dim(bed)[2] >= 6) 
    strands = as.character(bed[, 6])
  else strands = rep("+", N)
  for (i in 1:N) {
    chrom = as.character(bed[i, 1])
    strand = strands[i]
    if (strand == "+") {
      bigWig = bigWig.plus
      row = collect.counts(bigWig, chrom, start[i], end[i], 
                           step, do.sum)
      result[i, ] = row
    }
    else {
      bigWig = bigWig.minus
      if (is.null(bigWig.minus)) 
        bigWig = bigWig.plus
      row = collect.counts(bigWig, chrom, start[i], end[i], 
                           step, do.sum)
      result[i, ] = rev(row)
    }
  }
  result
}

meta.subsample.ChIPchip = function (bed, bigWig.plus, bigWig.minus, halfWindow, step, at.TSS = FALSE, 
                                    do.sum = FALSE) 
{
  N = dim(bed)[1]
  nPermut = 1000
  sampleFrac = 0.1
  windowSize = (2 * halfWindow)%/%step
  result = matrix(nrow = nPermut, ncol = windowSize)
  M = as.integer(round(N * sampleFrac, 0))
  values = collect.many.ChIPchip(bed, bigWig.plus, bigWig.minus, halfWindow, 
                        step, at.TSS = at.TSS, do.sum = do.sum)
  for (i in 1:nPermut) {
    idx <- sample(N, size = M, replace = T)
    result[i, ] = colSums(values[idx, ])/M
  }
  ci9 = sapply(1:windowSize, function(idx) quantile(result[, 
                                                           idx], 0.875))
  ci1 = sapply(1:windowSize, function(idx) quantile(result[, 
                                                           idx], 0.125))
  ci5 = sapply(1:windowSize, function(idx) quantile(result[, 
                                                           idx], 0.5))
  return(list(result, ci9, ci1, ci5))
}



collect.many.scaled.ChIPchip = function (bed, bigWig.plus, bigWig.minus, flank=1000, FLbuffer = 300, step=50, 
                                         do.sum = T) 
{
  #windowSize = ((2*flank) + (4*buffer))%/%step
  midPoint = (bed[, 2] + bed[, 3])/2
  TSS = bed[, 2]
  END = bed[,3]
  #neg.idxs = bed[, 6] == "-"
  #TSS[neg.idxs] = bed[neg.idxs, 3] - 1
  #END[neg.idxs] = bed[neg.idxs, 2] 
  
  prstart = (TSS - flank)
  prend = (TSS + FLbuffer)
  
  finstart = (END - FLbuffer)
  finend = (END + flank)
  #windowSize = (2*flank + 4*FLbuffer)/step
  windowSize = 320
  N = dim(bed)[1]
  result = matrix(nrow = N, ncol = windowSize)
  strands = NULL
  if (dim(bed)[2] >= 6) 
    strands = as.character(bed[, 6])
  else strands = rep("+", N)
  for (i in 1:N) {
    chrom = as.character(bed[i, 1])
    strand = strands[i]
    if (strand == "+") {
      bigWig = bigWig.plus
      midscaleParams = scale_params(bed[i,2], bed[i,3], step, buffer = 200)
      midbuf = midscaleParams[1]
      midstep = midscaleParams[2]
      midScaleFact = midscaleParams[3]
      midstart = (midPoint - midbuf)
      midend = (midPoint + midbuf)
      prrow = collect.counts(bigWig, chrom, prstart[i], prend[i], 
                             step, do.sum)
      midrow =  collect.counts(bigWig, chrom, midstart[i], midend[i], 
                               midstep, do.sum)/midScaleFact
      endrow = collect.counts(bigWig, chrom, finstart[i], finend[i], 
                              step, do.sum)
      row = c(prrow, midrow, endrow)
      if (length(row) == 320){
        result[i, ] = (row)}
    }
    else {
      bigWig = bigWig.minus
      if (is.null(bigWig.minus)) 
        bigWig = bigWig.plus
      midscaleParams = scale_params(bed[i,2], bed[i,3], step, buffer = 200)
      midbuf = midscaleParams[1]
      midstep = midscaleParams[2]
      midScaleFact = midscaleParams[3]
      midstart = (midPoint - midbuf)
      midend = (midPoint + midbuf)
      prrow = collect.counts(bigWig, chrom, prstart[i], prend[i], 
                             step, do.sum)
      midrow =  collect.counts(bigWig, chrom, midstart[i], midend[i], 
                               midstep, do.sum)/midScaleFact
      endrow = collect.counts(bigWig, chrom, finstart[i], finend[i], 
                              step, do.sum)
      row = c(prrow, midrow, endrow)
      if (length(row) == 320){
        result[i, ] = (rev(row))}
    }
  }
  rownames(result) = bed[,4]
  result[complete.cases(result), ]
}

meta.subsample.scaled.ChIPchip = function (bed, bigWig.plus, bigWig.minus, flank=1000, buffer = 300, step=50, do.sum = T) 
{
  values = collect.many.scaled.ChIPchip(bed, bigWig.plus, bigWig.minus, flank=flank, buffer, step, 
                                        do.sum = do.sum)
  N = dim(values)[1]
  nPermut = 1000
  sampleFrac = 0.1
  #windowSize = ((prend-prstart)+(midend-midstart)+(finend-finstart))/step
  windowSize = 320
  result = matrix(nrow = nPermut, ncol = windowSize)
  M = as.integer(round(N * sampleFrac, 0))
  #values = collect.many.scaled(bed, bigWig.plus, bigWig.minus, flank=flank, buffer, step, 
  #do.sum = do.sum)
  for (i in 1:nPermut) {
    idx <- sample(N, size = M, replace = T)
    result[i, ] = colSums(values[idx, ])/M
  }
  ci9 = sapply(1:windowSize, function(idx) quantile(result[, 
                                                           idx], 0.875))
  ci1 = sapply(1:windowSize, function(idx) quantile(result[, 
                                                           idx], 0.125))
  ci5 = sapply(1:windowSize, function(idx) quantile(result[, 
                                                           idx], 0.5))
  return(list(result, ci9, ci1, ci5))
}

getCounts.ChIPchip <- function(wig.p, wig.m, genes) {
  N = dim(genes)[1]
  plusStrand = genes[,2] == '+'
  
  result = vector(mode="integer", length=N)
  
  for (i in 1:N) {
    chrom = as.character(genes[i, 1])
    start = as.integer(genes[i, 3])
    end = as.integer(genes[i, 4])
    
    if (plusStrand[i]) {
      wig = wig.p
      
      qStart = start
      qEnd = end
      
      if (qStart > qEnd) {
        result[i] = NA
        next
      }
      
      data = query.bigWig(wig, chrom, qStart, qEnd)
      
      if (!is.null(data))
        result[i] = sum(data[,3])
    } else {
      wig = wig.m
      if (is.null(wig.m))
        wig = wig.p
      
      qStart = start
      qEnd = end
      
      if (qStart > qEnd) {
        result[i] = NA
        next
      }
      
      data = query.bigWig(wig, chrom, qStart, qEnd)
      
      if (!is.null(data))
        result[i] = sum(data[,3])
    }
  }
  
  return(cbind(result))
}

#####################################################################################################################
# Heatmaps (for ploting two samples and their overlap)
#####################################################################################################################
## Want to be able to generate a figure with 3 heat maps
## the first and second will be absolute counts per bin for each sample. 
##The third will be log2 fold change per bin between the two samples. 
## heatmaps will be plotted from -250 to + maxLength.

compareHeats = function(genes, BWset1, BWset2, maxLength = 4000, step = 10, BWpath = bwpath){
  winNum = maxLength/step
  BW1pl = load.bigWig(paste(BWpath, BWset1[[1]], sep = ""))
  BW1mn = load.bigWig(paste(BWpath,BWset1[[2]], sep = ""))
  BW2pl = load.bigWig(paste(BWpath,BWset2[[1]], sep = ""))
  BW2mn = load.bigWig(paste(BWpath,BWset2[[2]], sep = ""))
  Untreated = collect.many(genes, BW1pl, BW1mn, halfWindow = maxLength, step = step, at.TSS = T, do.sum = T)[,c((winNum - (250/step)+1):(2*winNum))]
  Un_normed = Untreated*(1/as.numeric(BWset1[[4]]))
  Un_normed[Un_normed == 0] <- 0.1 # permits calculation of log2 FC for all bins 
  row.names(Un_normed) <- genes[,4]
  Treated = collect.many(genes, BW2pl, BW2mn, halfWindow = maxLength, step = step, at.TSS = T, do.sum = T)[,c((winNum - (250/step)+1):(2*winNum))]
  Tr_normed = Treated*(1/as.numeric(BWset2[[4]]))
  Tr_normed[Tr_normed == 0] <- 0.1
  row.names(Tr_normed) <- genes[,4]
  FC = log2(Tr_normed/Un_normed)
  row.names(FC) <- genes[,4]
  ## also prepare corresponding metaplot data for comparisons
  Un_meta = meta.subsample(genes, BW1pl, BW1mn, step = step, at.TSS = T, halfWindow = maxLength, do.sum = T)
  Un_metaNorm = meta.normalize(result = Un_meta, scaleFactor = 1/as.numeric(BWset1[[4]]))
  Un_metaDat = cbind(seq(-240, maxLength, step), Un_metaNorm[[4]][c((winNum - (250/step)+1):(2*winNum))], 
                     Un_metaNorm[[3]][c((winNum - (250/step)+1):(2*winNum))], Un_metaNorm[[2]][c((winNum - (250/step)+1):(2*winNum))], 1, 1)
  Tr_meta = meta.subsample(genes, BW2pl, BW2mn, step = step, at.TSS = T, halfWindow = maxLength, do.sum = T)
  Tr_metaNorm = meta.normalize(result = Tr_meta, scaleFactor = 1/as.numeric(BWset2[[4]]))
  Tr_metaDat = cbind(seq(-240, maxLength, step), Tr_metaNorm[[4]][c((winNum - (250/step)+1):(2*winNum))], 
                     Tr_metaNorm[[3]][c((winNum - (250/step)+1):(2*winNum))], Tr_metaNorm[[2]][c((winNum - (250/step)+1):(2*winNum))], 1, 2)
  metaData = rbind(Un_metaDat, Tr_metaDat)
  colnames(metaData) <- c("x", "mean", "lower", "upper", "sample", "treatment")
  #result = list(as.data.frame(Un_normed), as.data.frame(Tr_normed), as.data.frame(FC))
  result = list(Un_normed, Tr_normed, FC, metaData)
  return(result)
}



#Function below will resort the 3 heat maps generated by "compareHeats" based on the geneList provided
sortHeats = function(genes, heatList){
  sorted = list()
  for (i in 1:length(heatList)){
    heat = merge(genes, heatList[[i]], by.x = "V4", by.y = 0)
    heat = heat[,c(7:dim(heat)[2])]
    row.names(heat) = genes[,4]
    sorted[[i]] = heat
  }  
  return(sorted)
}

# Sort gene list by diff criteria: 
## sort filtered genes by length (short to long)
orderGenesbyLen = function(geneList, StoL = F){
  len = geneList[,3]-geneList[,2]
  geneList$length = len
  if (StoL){
    result = geneList[order(geneList$length),]
  }
  else{
    result = geneList[order(-geneList$length),]
  }
  return(result)
}

orderGenesbyPI = function(geneList, CountDat = UntreatedCountDat, HighToLow = F){
  GL = merge(geneList, CountDat, by.x = "V4", by.y = 0)[,c(2,3,4,1,5,6,18)] ## col 18 is "scaled" pause index
  if (HighToLow){
    result = GL[order(GL[,7]),]
  }
  else{
    result = GL[order(-GL[,7]),]
  }
  return(result)
}

orderGenesbyActivity = function(geneList, CountDat = UntreatedCountDat, HighToLow = T){
  GL = merge(geneList, CountDat, by.x = "V4", by.y = 0)[,c(2,3,4,1,5,6,15)] ## col 15 is gene body density (i.e. "activity" )
  if (HighToLow){
    result = GL[order(GL[,7]),]
  }
  else{
    result = GL[order(-GL[,7]),]
  }
  return(result)
}


# Plotting parameters #
library(RColorBrewer)
#### function for printing all three plots (untreated, treated, foldChange) on one fig
Plot3Heat = function(compareHeatRes, filename = "/test_LevPlot.pdf", allFC = F, distToTSS = 4000, step = 10){
  Col_breaksRaw = c(seq(0, 3, length.out = 100), 10000) ## this gives 101 total values for mapping colors to
  Col_breaksFC = c(seq(-10, 10, length.out = 100), 10000) 
  col.raw <- colorRampPalette(c("white", "#EECD86",	"#E18942",	"#B95835",	"#3D3242"))(100) #set 100 values to within this color range. 
  col.raw <- c(col.raw, 'black') # add black as the last value in the color mapping list. 
  col.fc <- colorRampPalette(c('blue', 'white', 'red'))(100)
  col.fc <- c(col.fc, 'black')
  ## use distToTSS and step to determine x axis ticks and labels
  if (distToTSS >= 10000){
    tcks = c(26)
    labs = c("TSS")
    for (i in 1: (distToTSS/5000)){
      tcks = c(tcks, i*(5000/step)+26)
      labs = c(labs, sprintf("%sKb", i*5))
    }
  }
  else {
    tcks = c(0, 26*(10/step), 126*(10/step), 226*(10/step), 326*(10/step))
    labs = c("-250", "TSS", "1Kb", "2Kb", "3Kb")
  }
  #generate each heat map
  if (distToTSS > 10000){
    if (allFC){
      Unplot = levelplot(t(compareHeatRes[[1]]), ylab = "genes log2(Fold Change)", xlab = "Dist. from TSS (bp)", main = "LacZ fold change (FP/DMSO)",
                         at = Col_breaksFC, col.regions = col.fc, colorkey = list(title = "LacZ (Fold Change)", at = seq(-10, 10, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks, labels = labs)), aspect = 0.33, useRaster = T)
      Trplot = levelplot(t(compareHeatRes[[2]]), ylab = "genes log2(Fold Change)", xlab = "Dist. from TSS (bp)", main = "KD fold change (FP/DMSO)",
                         at = Col_breaksFC, col.regions = col.fc, colorkey = list(title = "KD (Fold Change)", at = seq(-10, 10, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks, labels = labs)), aspect = 0.33, useRaster = T)
      FCplot = levelplot(t(compareHeatRes[[3]]), ylab = "genes (FC KD / FC LacZ)", xlab = "Dist. from TSS (bp)", main = "fold change in fold change",
                         at = Col_breaksFC, col.regions = col.fc, colorkey = list(title = "FC_KD / FC_LacZ", at = seq(-10, 10, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks, labels = labs)), aspect = 0.33, useRaster = T)
    }
    else{
      Unplot = levelplot(t(log10(compareHeatRes[[1]])), ylab = "genes log10(reads)", xlab = "Dist. from TSS (bp)", main = "PROseq DMSO",
                         at = Col_breaksRaw, col.regions = col.raw, colorkey = list(title = "log10(reads)", at = seq(0, 3, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks, labels = labs)), aspect = 0.33, useRaster = T)
      Trplot = levelplot(t(log10(compareHeatRes[[2]])), ylab = "genes log10(reads)", xlab = "Dist. from TSS (bp)", main = "PROseq FP (300nM)",
                         at = Col_breaksRaw, col.regions = col.raw, colorkey = list(title = "log10(reads)", at = seq(0, 3, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks, labels = labs)), aspect = 0.33, useRaster = T)
      FCplot = levelplot(t(compareHeatRes[[3]]), ylab = "genes log2(fold change)", xlab = "Dist. from TSS (bp)", main = "fold change (FP/DMSO)",
                         at = Col_breaksFC, col.regions = col.fc, colorkey = list(title = "log2(FC)", at = seq(-10, 10, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks , labels = labs)), aspect = 0.33, useRaster = T)
    }
    ##plotting to jpeg
    #jpeg(file = paste(fig_dir, filename, sep = ''), width = 7, height = 10, units = "in", res = 300)
    pdf(file = paste(fig_dir, filename, sep = ''), width = 7, height = 10)
    plot.new()
    print(Unplot,  split = c(1,1,1,3), more = T)
    print(Trplot,  split = c(1,2,1,3), more = T)
    print(FCplot,  split = c(1,3,1,3), more = F)
    dev.off()
  }
  else{
    if (allFC){
      Unplot = levelplot(t(compareHeatRes[[1]]), ylab = "genes log2(Fold Change)", xlab = "Dist. from TSS (bp)", main = "WT fold change (FP/DMSO)",
                         at = Col_breaksFC, col.regions = col.fc, colorkey = list(title = "LacZ (Fold Change)", at = seq(-10, 10, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks, labels = labs)), aspect = 3, useRaster = T)
      Trplot = levelplot(t(compareHeatRes[[2]]), ylab = "genes log2(Fold Change)", xlab = "Dist. from TSS (bp)", main = "Mut fold change (Treated/Untreated)",
                         at = Col_breaksFC, col.regions = col.fc, colorkey = list(title = "KD (Fold Change)", at = seq(-10, 10, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks, labels = labs)), aspect = 3, useRaster = T)
      FCplot = levelplot(t(compareHeatRes[[3]]), ylab = "genes (FC KD / FC LacZ)", xlab = "Dist. from TSS (bp)", main = "fold change in fold change",
                         at = Col_breaksFC, col.regions = col.fc, colorkey = list(title = "FC_KD / FC_LacZ", at = seq(-10, 10, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks, labels = labs)), aspect = 3, useRaster = T)
    }
    else{
      Unplot = levelplot(t(log10(compareHeatRes[[1]])), ylab = "genes log10(reads)", xlab = "Dist. from TSS (bp)", main = "PROseq Untreated",
                         at = Col_breaksRaw, col.regions = col.raw, colorkey = list(title = "log10(reads)", at = seq(0, 3, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks, labels = labs)), aspect = 3, useRaster = T)
      Trplot = levelplot(t(log10(compareHeatRes[[2]])), ylab = "genes log10(reads)", xlab = "Dist. from TSS (bp)", main = "PROseq Treated",
                         at = Col_breaksRaw, col.regions = col.raw, colorkey = list(title = "log10(reads)", at = seq(0, 3, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks, labels = labs)), aspect = 3, useRaster = T)
      FCplot = levelplot(t(compareHeatRes[[3]]), ylab = "genes log2(fold change)", xlab = "Dist. from TSS (bp)", main = "fold change (Treated/Untreated)",
                         at = Col_breaksFC, col.regions = col.fc, colorkey = list(title = "log2(FC)", at = seq(-10, 10, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks , labels = labs)), aspect = 3, useRaster = T)
    }
    ##plotting to jpeg
    #jpeg(file = paste(fig_dir, filename, sep = ''), width = 10, height = 7, units = "in", res = 300)
    pdf(file = paste(fig_dir, filename, sep = ''), width = 7, height = 10)
    plot.new()
    print(Unplot,  split = c(1,1,3,1), more = T)
    print(Trplot,  split = c(2,1,3,1), more = T)
    print(FCplot,  split = c(3,1,3,1), more = F)
    dev.off()
  }
}

#############################################################################################################
## functions below are modifications of the heatmap functions above.  
## These add trace lines for the gene boundaries. 
## NOTE: only use if you are anchoring on the TSS. 

compareHeats.traceGene = function(genes, BWset1, BWset2, maxLength = 4000, step = 10, BWpath = bwpath){
  winNum = maxLength/step
  WindowGeneEnd = 250/step + ceiling((genes[,3]- genes[,2])/step) # this defines the CPS of each gene
  BW1pl = load.bigWig(paste(BWpath, BWset1[[1]], sep = ""))
  BW1mn = load.bigWig(paste(BWpath,BWset1[[2]], sep = ""))
  BW2pl = load.bigWig(paste(BWpath,BWset2[[1]], sep = ""))
  BW2mn = load.bigWig(paste(BWpath,BWset2[[2]], sep = ""))
  Untreated = collect.many(genes, BW1pl, BW1mn, halfWindow = maxLength, step = step, at.TSS = T, do.sum = T)[,c((winNum - (250/step)+1):(2*winNum))]
  Un_normed = Untreated*(1/as.numeric(BWset1[[4]]))
  Un_normed[Un_normed == 0] <- 0.1 # permits calculation of log2 FC for all bins 
  row.names(Un_normed) <- genes[,4]
  Treated = collect.many(genes, BW2pl, BW2mn, halfWindow = maxLength, step = step, at.TSS = T, do.sum = T)[,c((winNum - (250/step)+1):(2*winNum))]
  Tr_normed = Treated*(1/as.numeric(BWset2[[4]]))
  Tr_normed[Tr_normed == 0] <- 0.1
  row.names(Tr_normed) <- genes[,4]
  FC = log2(Tr_normed/Un_normed)
  row.names(FC) <- genes[,4]
  # for each gene change the value of the cells for the start and end of each gene (change to 10000)
  for (i in 1:length(WindowGeneEnd)){
    if (WindowGeneEnd[i] < (winNum + (250/step))){
      Un_normed[i, WindowGeneEnd[i]] = 10000
      Tr_normed[i, WindowGeneEnd[i]] = 10000
      FC[i, WindowGeneEnd[i]] = 10000
    }
  }
  Un_normed[,(250/step)] = 10000
  Tr_normed[,(250/step)] = 10000
  FC[,(250/step)] = 10000
  ## also prepare corresponding metaplot data for comparisons
  Un_meta = meta.subsample(genes, BW1pl, BW1mn, step = step, at.TSS = T, halfWindow = maxLength, do.sum = T)
  Un_metaNorm = meta.normalize(result = Un_meta, scaleFactor = 1/as.numeric(BWset1[[4]]))
  Un_metaDat = cbind(seq(-240, maxLength, step), Un_metaNorm[[4]][c((winNum - (250/step)+1):(2*winNum))], 
                     Un_metaNorm[[3]][c((winNum - (250/step)+1):(2*winNum))], Un_metaNorm[[2]][c((winNum - (250/step)+1):(2*winNum))], 1, 1)
  Tr_meta = meta.subsample(genes, BW2pl, BW2mn, step = step, at.TSS = T, halfWindow = maxLength, do.sum = T)
  Tr_metaNorm = meta.normalize(result = Tr_meta, scaleFactor = 1/as.numeric(BWset2[[4]]))
  Tr_metaDat = cbind(seq(-240, maxLength, step), Tr_metaNorm[[4]][c((winNum - (250/step)+1):(2*winNum))], 
                     Tr_metaNorm[[3]][c((winNum - (250/step)+1):(2*winNum))], Tr_metaNorm[[2]][c((winNum - (250/step)+1):(2*winNum))], 1, 2)
  metaData = rbind(Un_metaDat, Tr_metaDat)
  colnames(metaData) <- c("x", "mean", "lower", "upper", "sample", "treatment")
  #result = list(as.data.frame(Un_normed), as.data.frame(Tr_normed), as.data.frame(FC))
  result = list(Un_normed, Tr_normed, FC, metaData)
  return(result)
}

Plot3Heat.traceGene = function(compareHeatRes, filename = "/test_LevPlot.jpeg", allFC = F, distToTSS = 4000, step = 10){
  Col_breaksRaw = c(seq(0, 3, length.out = 100), 10000) ## this gives 101 total values for mapping colors to
  Col_breaksFC = c(seq(-10, 10, length.out = 100), 10000) 
  col.raw <- colorRampPalette(c("white", "#EECD86",	"#E18942",	"#B95835",	"#3D3242"))(100) #set 100 values to within this color range. 
  col.raw <- c(col.raw, 'black') # add black as the last value in the color mapping list. 
  col.fc <- colorRampPalette(c('blue', 'white', 'red'))(100)
  col.fc <- c(col.fc, 'black')
  ## use distToTSS and step to determine x axis ticks and labels
  if (distToTSS >= 10000){
    tcks = c(26)
    labs = c("TSS")
    for (i in 1: (distToTSS/5000)){
      tcks = c(tcks, i*(5000/step)+26)
      labs = c(labs, sprintf("%sKb", i*5))
    }
  }
  else {
    tcks = c(0, 26*(10/step), 126*(10/step), 226*(10/step), 326*(10/step))
    labs = c("-250", "TSS", "1Kb", "2Kb", "3Kb")
  }
  #generate each heat map
  if (distToTSS > 10000){
    if (allFC){
      Unplot = levelplot(t(compareHeatRes[[1]]), ylab = "genes log2(Fold Change)", xlab = "Dist. from TSS (bp)", main = "LacZ fold change (FP/DMSO)",
                         at = Col_breaksFC, col.regions = col.fc, colorkey = list(title = "LacZ (Fold Change)", at = seq(-10, 10, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks, labels = labs)), aspect = 0.33, useRaster = T)
      Trplot = levelplot(t(compareHeatRes[[2]]), ylab = "genes log2(Fold Change)", xlab = "Dist. from TSS (bp)", main = "KD fold change (FP/DMSO)",
                         at = Col_breaksFC, col.regions = col.fc, colorkey = list(title = "KD (Fold Change)", at = seq(-10, 10, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks, labels = labs)), aspect = 0.33, useRaster = T)
      FCplot = levelplot(t(compareHeatRes[[3]]), ylab = "genes (FC KD / FC LacZ)", xlab = "Dist. from TSS (bp)", main = "fold change in fold change",
                         at = Col_breaksFC, col.regions = col.fc, colorkey = list(title = "FC_KD / FC_LacZ", at = seq(-10, 10, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks, labels = labs)), aspect = 0.33, useRaster = T)
    }
    else{
      Unplot = levelplot(t(log10(compareHeatRes[[1]])), ylab = "genes log10(reads)", xlab = "Dist. from TSS (bp)", main = "PROseq DMSO",
                         at = Col_breaksRaw, col.regions = col.raw, colorkey = list(title = "log10(reads)", at = seq(0, 3, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks, labels = labs)), aspect = 0.33, useRaster = T)
      Trplot = levelplot(t(log10(compareHeatRes[[2]])), ylab = "genes log10(reads)", xlab = "Dist. from TSS (bp)", main = "PROseq FP (300nM)",
                         at = Col_breaksRaw, col.regions = col.raw, colorkey = list(title = "log10(reads)", at = seq(0, 3, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks, labels = labs)), aspect = 0.33, useRaster = T)
      FCplot = levelplot(t(compareHeatRes[[3]]), ylab = "genes log2(fold change)", xlab = "Dist. from TSS (bp)", main = "fold change (FP/DMSO)",
                         at = Col_breaksFC, col.regions = col.fc, colorkey = list(title = "log2(FC)", at = seq(-10, 10, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks , labels = labs)), aspect = 0.33, useRaster = T)
    }
    ##plotting to jpeg
    #jpeg(file = paste(fig_dir, filename, sep = ''), width = 7, height = 10, units = "in", res = 300)
    pdf(file = paste(fig_dir, filename, sep = ''), width = 7, height = 10)
    plot.new()
    print(Unplot,  split = c(1,1,1,3), more = T)
    print(Trplot,  split = c(1,2,1,3), more = T)
    print(FCplot,  split = c(1,3,1,3), more = F)
    dev.off()
  }
  else{
    if (allFC){
      Unplot = levelplot(t(compareHeatRes[[1]]), ylab = "genes log2(Fold Change)", xlab = "Dist. from TSS (bp)", main = "WT fold change (FP/DMSO)",
                         at = Col_breaksFC, col.regions = col.fc, colorkey = list(title = "LacZ (Fold Change)", at = seq(-10, 10, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks, labels = labs)), aspect = 3, useRaster = T)
      Trplot = levelplot(t(compareHeatRes[[2]]), ylab = "genes log2(Fold Change)", xlab = "Dist. from TSS (bp)", main = "Mut fold change (Treated/Untreated)",
                         at = Col_breaksFC, col.regions = col.fc, colorkey = list(title = "KD (Fold Change)", at = seq(-10, 10, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks, labels = labs)), aspect = 3, useRaster = T)
      FCplot = levelplot(t(compareHeatRes[[3]]), ylab = "genes (FC KD / FC LacZ)", xlab = "Dist. from TSS (bp)", main = "fold change in fold change",
                         at = Col_breaksFC, col.regions = col.fc, colorkey = list(title = "FC_KD / FC_LacZ", at = seq(-10, 10, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks, labels = labs)), aspect = 3, useRaster = T)
    }
    else{
      Unplot = levelplot(t(log10(compareHeatRes[[1]])), ylab = "genes log10(reads)", xlab = "Dist. from TSS (bp)", main = "PROseq Untreated",
                         at = Col_breaksRaw, col.regions = col.raw, colorkey = list(title = "log10(reads)", at = seq(0, 3, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks, labels = labs)), aspect = 3, useRaster = T)
      Trplot = levelplot(t(log10(compareHeatRes[[2]])), ylab = "genes log10(reads)", xlab = "Dist. from TSS (bp)", main = "PROseq Treated",
                         at = Col_breaksRaw, col.regions = col.raw, colorkey = list(title = "log10(reads)", at = seq(0, 3, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks, labels = labs)), aspect = 3, useRaster = T)
      FCplot = levelplot(t(compareHeatRes[[3]]), ylab = "genes log2(fold change)", xlab = "Dist. from TSS (bp)", main = "fold change (Treated/Untreated)",
                         at = Col_breaksFC, col.regions = col.fc, colorkey = list(title = "log2(FC)", at = seq(-10, 10, length.out = 100)),
                         scales = list(y = list(draw = F), x = list(at = tcks , labels = labs)), aspect = 3, useRaster = T)
    }
    ##plotting to jpeg
    #jpeg(file = paste(fig_dir, filename, sep = ''), width = 10, height = 7, units = "in", res = 300)
    pdf(file = paste(fig_dir, filename, sep = ''), width = 7, height = 10)
    plot.new()
    print(Unplot,  split = c(1,1,3,1), more = T)
    print(Trplot,  split = c(2,1,3,1), more = T)
    print(FCplot,  split = c(3,1,3,1), more = F)
    dev.off()
  }
}
#############################################################################################################
