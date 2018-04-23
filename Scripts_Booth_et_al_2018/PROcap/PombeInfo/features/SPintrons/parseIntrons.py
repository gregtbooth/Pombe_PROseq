#! /usr/bin/python 
## This is a script to split up introns based on whether they are a first, second, third, etc. intron.  each class will get return it's own bed file. 
import sys 
infile = sys.argv[1] 
intronList = open(infile, 'r')
out1 = open("1stIntrons.bed", "w") 
out2 = open("2ndIntrons.bed", "w") 
out3 = open("3rdIntrons.bed", "w") 
out4 = open("4thIntrons.bed", "w") 
out5 = open("5thIntrons.bed", "w") 
outLast =  open("LastIntrons.bed", "w") 

intronDict = {}

for line in intronList.readlines(): 
	wrds = line.strip().split()
	chro, start, stop, gene, strand = wrds[0], wrds[1], wrds[2], wrds[3], wrds[5]
	if gene not in intronDict.keys():
		intronDict[gene] = [[chro, start, stop, gene, strand]]
	else: 
		intronDict[gene] += [[chro, start, stop, gene, strand]]

for key in intronDict.keys():
	number = 0
	if intronDict[key][0][4] == "+": # if not plus strand, add introns in reverse order
		for intron in intronDict[key]:
			length = len(intronDict[key])
			number +=1
			if length > 1: 
				outLast.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(intronDict[key][length-1][0], intronDict[key][length-1][1], intronDict[key][length-1][2], intronDict[key][length-1][3], 0, intronDict[key][length-1][4]))

			if number == 1:
				out1.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(intronDict[key][0][0], intronDict[key][0][1], intronDict[key][0][2], intronDict[key][0][3], 0, intronDict[key][0][4]))
			elif number == 2:
				out2.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(intronDict[key][1][0], intronDict[key][1][1], intronDict[key][1][2], intronDict[key][1][3], 0, intronDict[key][1][4]))
			elif number == 3:
				out3.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(intronDict[key][2][0], intronDict[key][2][1], intronDict[key][2][2], intronDict[key][2][3], 0, intronDict[key][2][4]))
			elif number == 4:
				out4.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(intronDict[key][3][0], intronDict[key][3][1], intronDict[key][3][2], intronDict[key][3][3], 0, intronDict[key][3][4]))
			elif number == 5:
				out5.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(intronDict[key][4][0], intronDict[key][4][1], intronDict[key][4][2], intronDict[key][4][3], 0, intronDict[key][4][4]))
	else: 
		for intron in intronDict[key]:
			number +=1
			length = len(intronDict[key])
			if length > 1: 
				outLast.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(intronDict[key][-length][0], intronDict[key][-length][1], intronDict[key][-length][2], intronDict[key][-length][3], 0, intronDict[key][-length][4]))
			if number == 5:	
				out5.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(intronDict[key][-5][0], intronDict[key][-5][1], intronDict[key][-5][2], intronDict[key][-5][3], 0, intronDict[key][-5][4]))
			elif number == 4:
				out4.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(intronDict[key][-4][0], intronDict[key][-4][1], intronDict[key][-4][2], intronDict[key][-4][3], 0, intronDict[key][-4][4]))
			elif number == 3:
				out3.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(intronDict[key][-3][0], intronDict[key][-3][1], intronDict[key][-3][2], intronDict[key][-3][3], 0, intronDict[key][-3][4]))
			elif number == 2:
				out2.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(intronDict[key][-2][0], intronDict[key][-2][1], intronDict[key][-2][2], intronDict[key][-2][3], 0, intronDict[key][-2][4]))
			elif number == 1:
				out1.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(intronDict[key][-1][0], intronDict[key][-1][1], intronDict[key][-1][2], intronDict[key][-1][3], 0, intronDict[key][-1][4]))
	
intronList.close()
out1.close()
out2.close()
out3.close()
out4.close()
out5.close()
outLast.close()
