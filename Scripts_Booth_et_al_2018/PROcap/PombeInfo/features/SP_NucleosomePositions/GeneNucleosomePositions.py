# -*- coding: utf-8 -*-
"""
Created on Wed Jan 28 13:09:01 2015

@author: Gregory Booth 
"""

### Script for taking gene Coordinates and Nucleosome coordinates and returning files listing coords for +1, +2, +3 nucleosomes 
## format of input files:  geneName genename2   chrom   strand  start   end

import sys 
import timeit 

startTime = timeit.default_timer()
#GL = sys.argv[1]
#NL = sys.argv[2]
GL = open("/Volumes/SEAGATE_EXP/Greg_YeastData/PRO_Seq/Pombe/CombinedWTdata/Kbro_info/pombeFull-lengthGenelistForGSP.txt_sorted.txt", 'r')
NL = open("/Volumes/SEAGATE_EXP/Greg_YeastData/sub_genelists/SP_NucleosomePositions/pombe_Nucl_centers_LFN_TF.bed.txt")
plus1Nucl = open("/Volumes/SEAGATE_EXP/Greg_YeastData/sub_genelists/SP_NucleosomePositions/SP_firstNuclcoords.txt", 'w')
plus2Nucl = open("/Volumes/SEAGATE_EXP/Greg_YeastData/sub_genelists/SP_NucleosomePositions/SP_secondNuclcoords.txt", 'w')
plus3Nucl = open("/Volumes/SEAGATE_EXP/Greg_YeastData/sub_genelists/SP_NucleosomePositions/SP_thirdNuclcoords.txt", 'w')
plus4Nucl = open("/Volumes/SEAGATE_EXP/Greg_YeastData/sub_genelists/SP_NucleosomePositions/SP_fourthNuclcoords.txt", 'w')

### split nucleosomes by chromosomal location.  

NuclDict = {}

for line in NL.readlines():
    wrds = line.strip().split()
    chro = wrds[1]
    center = wrds[2]
    if chro not in NuclDict.keys(): 
        NuclDict[chro] = [center]
    else: 
        NuclDict[chro]+= [center]
        
def orderGeneNucls(NuclDict, chro, strand, start, stop):
    NuclList = []    
    if chro in NuclDict.keys():
        if strand == '+':
            for nucl in NuclDict[chro]: 
                if int(nucl) >= int(start):
                    if int(nucl) <= int(stop):
                        NuclList += [nucl]
            NuclList = sorted(NuclList)
        else:
             for nucl in NuclDict[chro]: 
                if int(nucl) >= int(start):
                    if int(nucl) <= int(stop):
                        NuclList += [nucl]
             NuclList = sorted(NuclList, reverse = True)
    return NuclList



geneNuclDict = {}
for line in GL.readlines(): 
    wrds = line.strip().split() 
    Gene = wrds[1]
    chro = wrds[2]
    strand = wrds[3]
    start = wrds[4]
    stop = wrds[5]
    geneNuclDict[Gene] = orderGeneNucls(NuclDict, chro, strand, start, stop)

#print  geneNuclDict['YAL043C'], geneNuclDict['YAL041W'], geneNuclDict['YBL060W'],len(geneNuclDict['YAL043C']), len(geneNuclDict['YAL041W']), len(geneNuclDict['YBL060W']) 
GL.close()
GL1 = open("/Volumes/SEAGATE_EXP/Greg_YeastData/PRO_Seq/Pombe/CombinedWTdata/Kbro_info/pombeFull-lengthGenelistForGSP.txt_sorted.txt", 'r')


for line in GL1.readlines(): 
    	wrds = line.strip().split() 
    	gene = wrds[1]
    	chro = wrds[2]
    	strand = wrds[3]
    	start = wrds[4]
    	stop = wrds[5]
    	if strand =="+": #  
    		if len(geneNuclDict[gene])>0: #takes this first listed nucl center only if it's within 200bp from TSS
			if int(geneNuclDict[gene][0])< (200 + int(start)):
        			plus1Nucl.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(gene,"+1nucleosome", chro,strand, geneNuclDict[gene][0], int(geneNuclDict[gene][0])+1))
        	if len(geneNuclDict[gene])>1: 
            		plus2Nucl.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(gene,"+2nucleosome", chro,strand, geneNuclDict[gene][1], int(geneNuclDict[gene][1])+1))
            	if len(geneNuclDict[gene])>2: 
                	plus3Nucl.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(gene,"+3nucleosome", chro,strand, geneNuclDict[gene][2],int(geneNuclDict[gene][2])+1))
		if len(geneNuclDict[gene])>3:
			plus4Nucl.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(gene,"+4nucleosome", chro,strand, geneNuclDict[gene][3],int(geneNuclDict[gene][3])+1))
	else: 
		if len(geneNuclDict[gene])>0:
			if int(geneNuclDict[gene][0])> (int(stop)-200):
				plus1Nucl.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(gene,"+1nucleosome", chro,strand, geneNuclDict[gene][0], int(geneNuclDict[gene][0])+1))
		if len(geneNuclDict[gene])>1:
			plus2Nucl.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(gene,"+2nucleosome", chro,strand, geneNuclDict[gene][1], int(geneNuclDict[gene][1])+1))
		if len(geneNuclDict[gene])>2:
			plus3Nucl.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(gene,"+3nucleosome", chro,strand, geneNuclDict[gene][2],int(geneNuclDict[gene][2])+1))
		if len(geneNuclDict[gene])>3:
			plus4Nucl.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(gene,"+4nucleosome", chro,strand, geneNuclDict[gene][3],int(geneNuclDict[gene][3])+1))


plus1Nucl.close()
plus2Nucl.close()
plus3Nucl.close()
plus4Nucl.close()            
                        
stopTime = timeit.default_timer()
print stopTime-startTime, 'Seconds: runtime'
    
