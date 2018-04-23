#! /usr/bin/python
# This script was written to take PRO-cap (+TAP) plus and minus bedgraphs, as well as a gene list, and find the base within 500bp around each annotated TSS that has the greatest read count.  This base is considered the "observed TSS". In cases where more than one base within the window has an equally high number of reads, the closest base to the annotation is selected.  For a base to be called the observed TSS it must have > 5 reads total.  This threshold should probably change depending on the read depth of an experiment, but most covered bases in my yeast data  have fewer than 5 reads.  This script outputs a .bed file in the same folder as the input genelist.  The bed file gives additional information about the observed TSS, including read counts, distance from the annotated TSS (absolute value (max 250)), etc.  

#earlier versions were only suitable for integer counts, however, because you will typically want to first normalize the plus and minus TAP libraries by their read depth (i.e. RPM), before subtracting, the script now allows for handling of floats.  

# becauase the absolute difference between two normalized samples varies from comparison to comparison, you must adjust the threshold of calling an observed TSS(so that after normalization and subtraction the difference is equivalent to 5 more in the Plus TAP sample). 

###

import time 
import sys 
import numpy as np

print 'usage:  python getObservedTSS.py plus.bedgraph minus.bedgraph geneList\n'

start_time  = time.time()

pbedgraph = sys.argv[1]
mbedgraph = sys.argv[2]
gList = sys.argv[3]

pbg = open(pbedgraph, "r") 
mbg = open(mbedgraph, "r") 
GL = open (gList, "r")

def MakeChrDict(bedgraph):
	chrDict = {}
	for line in bedgraph.readlines():
		wrds = line.strip().split()
		chro, start, count = wrds[0], wrds[1], wrds[3]
		if chro not in chrDict.keys():
			chrDict[chro] = {int(start):count}
		else:
			chrDict[chro][int(start)]=  float(count)
	return chrDict

print "making chromosome dictionary"
pchroDict = MakeChrDict(pbg)
mchroDict = MakeChrDict(mbg)

###  below is a function to find the best (highest read denstiy) 50mer window around the observed TSS
def best50bpWin(position, chro, bgDict): 
	winList = []
	sumList = []
	for n in range((position - 50), position): 
		winCounts = []
		start = n 
		end = n + 50
		for nj in range(start, end): 
			try: 
				winCounts += [abs(float(bgDict[chro][nj]))]
			except KeyError: 
				winCounts += [0]
		winList += [winCounts]
		sumList += [sum(winCounts)]
	bestIdx = sumList.index(max(sumList))
	return winList[bestIdx]

#### Below is a function to describe characteristics of the 50mer within which the TSS occurs.  defines the average base at which reads occur and then the StDev to get an idea of how dispersed reads are around the "observed TSS" 
def CountDispersion(countList): 
	#events = sum(countList)
	baseUse = []
	base = 0
	for count in countList: 
		base += 1
		baseUse += [base] * count #repeats base n = count times
	#mean = np.mean(baseUse)
	stDev = np.std(baseUse)
	return(stDev)


#chrSwap = {"chr01":"chrI", "chr02": "chrII", "chr03":"chrIII","chrMTR":"chrMTR"} 
def getObservedTSS(geneList, plusDict, minusDict): 
	outfile = open("%s_PROcapObservedTSS.bed" %(gList), 'w')
	outfile.write("chr\tObsTSS\tend\tgeneName\treadCount\tstrand\tDistFromAnnTSS\tTSS_50bp\n")
	for line in geneList.readlines(): 
		wrds = line.strip().split()
		name, chro, strand, start, stop = wrds[3], wrds[0], wrds[5], int(wrds[1]), int(wrds[2])
		if strand == "+": 
			CountList = []
			for ii in range(start-250, start+250): 
				try: 
					CountList += [abs(float(plusDict[chro][ii]))]
				except KeyError: 
					CountList +=[0] ## so all bases are considered for later indexing
			TSScount = max(abs(y) for y in CountList) #gets maximum read countaround TSS
			if TSScount > 0.39: ## sets a threshold for counts (i.e. 5) to be called "observed TSS"
				bestPositions = [i for i, x in enumerate(CountList) if x == TSScount] ## get all indices that have the maximum (i.e. ties for TSS) 
				distfromAnnTSS = [abs(250 - j) for j in bestPositions]
				TrueTSS = bestPositions[distfromAnnTSS.index(min(distfromAnnTSS))] ## gets the best position that is closest to annotated TSS if there are ties
				obsTSSpos = start - 250 + TrueTSS
			else: 
				obsTSSpos = start
				distfromAnnTSS  = [0]
				TSScount = CountList[250] ## count at annotated  TSS
			TSSregion = best50bpWin(obsTSSpos, chro, plusDict)
			#print TSSregion
			#dispersion = CountDispersion(TSSregion) 
			outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(chro, obsTSSpos, stop, name, TSScount, strand, min(distfromAnnTSS), TSSregion))
		else: 
			CountList = []
			for ii in range(stop-250, stop+250): 
				try: 
					CountList += [abs(float(minusDict[chro][ii]))]
				except KeyError: 
					CountList +=[0] ## so all bases are considered for later indexing     
			TSScount = max(abs(y) for y in CountList) #gets maximum read countaround TSS
			if TSScount > 0.39: 
				bestPositions = [i for i, x in enumerate(CountList) if x == TSScount] ## get all indices that have the maximum (i.e. ties for TSS) 
				distfromAnnTSS = [abs(250 - j) for j in bestPositions]                        
				TrueTSS = bestPositions[distfromAnnTSS.index(min(distfromAnnTSS))] ## gets the best position that is closest to annotated TSS if there are ties
				obsTSSpos = stop - 250 + TrueTSS +1 # must add the 1 back due to 0-based index coords (only used on minus strand)
			else: 
				obsTSSpos = stop
				distfromAnnTSS = [0]
				TSScount = CountList[250]
			TSSregion = best50bpWin(obsTSSpos, chro, minusDict)
			#print TSSregion
			#dispersion = CountDispersion(TSSregion)
			outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(chro, start, obsTSSpos, name, TSScount, strand, min(distfromAnnTSS), TSSregion))
	outfile.close()


getObservedTSS(GL, pchroDict, mchroDict) 

pbg.close()
mbg.close()
GL.close()

print "---------------- run-time: %s seconds  ----------------------" %(time.time() - start_time)


