#!/usr/bin/python

# 09-28-17 Edit
## Noticed that on the minus strand, the script was returning situations where the background sample had more signal.  This was because of a missing else statment (documented below), which automatically returned the negative value as if it were real (capped) signal.  This was a flaw that often lead to incorrect obsTSS calls for the minus strand. Do not use versions of this script older than this version. 

import os 
import sys 
import time

start_time = time.time()

pTAPbg = sys.argv[1]
mTAPbg = sys.argv[2]
strand = sys.argv[3]## plus or minus

pTAPdat = open(pTAPbg, 'r')
mTAPdat = open(mTAPbg, 'r') 

outbg = open("%s_BackgroundSubtracted.bedgraph"%(pTAPbg), 'w') 

def CreateReadDict(bg):
    nestedDict = {}	
    for line in bg.readlines():
        wrds = line.strip().split()
        chro, start, end, count = wrds[0],wrds[1],wrds[2], abs(float(wrds[3]))
        span = int(end) - int(start) 
        if chro not in nestedDict.keys():
            for i in range(span):
                nestedDict[chro] = {(int(start)+i): count}
            #nestedDict[chro] = {start: count}
        else: 
            for i in range(span):
                nestedDict[chro][int(start)+i]=count 
            #nestedDict[chro].update({start: count})
    return nestedDict


print "writing dictionary of covered positions in minus TAP sample" 
mTAPdict = CreateReadDict(mTAPdat)

print "subtracting background where data present (min = 0), else keeping +TAPcount ... be patient"
 
if strand == 'plus':
    for line in pTAPdat.readlines():
        wrds = line.strip().split()
        Pchro, Pstart, Pend, Pcount = wrds[0],wrds[1],wrds[2], abs(float(wrds[3]))
        #Pspan = int(Pend) - int(Pstart)
        #if Pchro != 'chrmt':
                #for i in range(Pspan): 
                        #if int(Pstart)+i in NETdict[Pchro].keys():
        try:
        	 signal = ((Pcount) - mTAPdict[Pchro][int(Pstart)])
		 if signal < 0: 
			signal = 0
			outbg.write("%s\t%s\t%s\t%s\n" %(Pchro, Pstart,Pend, signal)) 
        	 else:
			outbg.write("%s\t%s\t%s\t%s\n" %(Pchro, Pstart,Pend, signal))
       	except KeyError:
                 signal = Pcount
		 outbg.write("%s\t%s\t%s\t%s\n" %(Pchro, Pstart,Pend, signal))
                        
                                
elif strand =='minus':   
     for line in pTAPdat.readlines():
        wrds = line.strip().split()
        Pchro, Pstart, Pend, Pcount = wrds[0],wrds[1],wrds[2], abs(float(wrds[3]))
        #Pspan = int(Pend) - int(Pstart)
        #if Pchro != 'chrmt':
	try:
		signal = ((Pcount) - mTAPdict[Pchro][int(Pstart)])
		if signal < 0: 
			signal = 0 
			outbg.write("%s\t%s\t%s\t%s\n" %(Pchro, Pstart,Pend, signal))
		else: # this else statement was added 09-27-17 which corrected the issue 
			outbg.write("%s\t%s\t%s\t%s\n" %(Pchro, Pstart,Pend, -signal))
	except KeyError:
		signal = -Pcount
		outbg.write("%s\t%s\t%s\t%s\n" %(Pchro, Pstart,Pend, signal))

outbg.close()

#print "making bigwig file" 

#os.system("/home/macproadmin/users/GROseq_Utils/mkbigWig/bedGraphToBigWig ./PROseqNETseqRATIO_%s.begraph /media/macproadmin/SEAGATE_EXP/GROseq_Utils/annotations/sacCer3/SacCerChr.sizes_forPipeline.txt ./PROseqNETseqRATIO_%s.bw" %(strand, strand))


print("runtime: %s seconds" %(time.time() - start_time))

