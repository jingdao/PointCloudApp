#!/usr/bin/python

import os
import glob
import sys

if len(sys.argv) >=2:
	inputDir = sys.argv[1]
else:
	inputDir = '.'

if len(sys.argv) >=5:
	start = int(sys.argv[2])
	inc = int(sys.argv[3])
	end = int(sys.argv[4])
else:
	start = 0
	inc = 1
	end = 0
		
if "-i" in sys.argv:
	ignoreZeros = True
else:
	ignoreZeros = False
if "-c" in sys.argv:
	saveCSV = True
else:
	saveCSV = False
truePositive = {}
falsePositive = {}
falseNegative = {}

numFiles = 0
numTrials = 0

def increment(tab,key):
	if key in tab:
		tab[key] += 1
	else:
		tab[key] = 1

categories = []
if os.path.isfile(inputDir+'/labelCategory.txt'):
	categoryFile = open(inputDir+'/labelCategory.txt')
	for l in categoryFile:
		categories.append(l.split()[1])
	categoryFile.close()	

indexFile = open(inputDir+'/svm_prediction.txt')
labelIndex = indexFile.readline().split()
indexFile.close()

statistics=[]

for i in range(start,end+1,inc):
	d = inputDir+'/clusters'+str(i).zfill(len(sys.argv[2]))
#	d = '.'
	if os.path.isdir(d) and os.path.isfile(d+'/labels.txt') and os.path.isfile(d+'/prediction.txt'):
		numFiles += 1
		if len(categories) == 0 and os.path.isfile(d+'/labelCategory.txt'):
			categoryFile = open(inputDir+'/labelCategory.txt')
			for l in categoryFile:
				categories.append(l.split()[1])
			categoryFile.close()	
		labels = open(d+'/labels.txt','r')
		prediction = open(d+'/prediction.txt','r')
		while True:
			try:
				l = int(labels.readline())
			except ValueError:
				break
			prob = prediction.readline().split()
			p = int(prob[0])
			if ignoreZeros and (l==0 or p==0):
				continue
			numTrials += 1
			if l == p:
				increment(truePositive,l)
			else:
				increment(falsePositive,p)
				increment(falseNegative,l)
			statLabel=[numTrials,categories[l],categories[p]]
			statProb = [0] * len(categories)
			for i in range(1,len(prob)):
				statProb[int(labelIndex[i])] = prob[i]
			statLabel.extend(statProb)
			statistics.append(statLabel)
		labels.close()
		prediction.close()

print "Processed " + str(numFiles) + " files (" + str(numTrials) + " trials)"
print "%2s %15s %4s %8s %10s %6s" % ("","Category","Num","Accuracy","Precision","Recall")
for i in range(len(categories)):
	TP = truePositive[i] if i in truePositive else 0
	FP = falsePositive[i] if i in falsePositive else 0
	FN = falseNegative[i] if i in falseNegative else 0
	num = TP + FN
	acc = (numTrials - FP - FN) * 100.0 / numTrials
	pre = TP * 100.0 / (TP + FP) if (TP + FP) != 0 else float('nan')
	rec = TP * 100.0 / (TP + FN) if (TP + FN) != 0 else float('nan')
	print "%2d %15s %4d %8.2f %10.2f %6.2f" % (i,categories[i],num,acc,pre,rec)

totalAcc = sum(truePositive.values()) * 100.0 / numTrials
totalPre = sum(truePositive.values()) * 100.0 / (sum(truePositive.values()) + sum(falsePositive.values()))
totalRec = sum(truePositive.values()) * 100.0 / (sum(truePositive.values()) + sum(falseNegative.values()))
print "%2s %15s %4d %8.2f %10.2f %6.2f" % ("","overall",numTrials,totalAcc,totalPre,totalRec)

if saveCSV:
	csvFile = open(inputDir+'/score.csv','w')
	csvFile.write("%6s,%15s,%15s," % ("Sample","Actual","Predicted"))
	for i in range(len(categories)):
		csvFile.write("%15s," % categories[i])
	csvFile.write('\n')
	for s in statistics:
		csvFile.write("%6d,%15s,%15s," % (s[0],s[1],s[2]))
		for i in range(len(categories)):
			csvFile.write("%15.4f," % float(s[3+i]))
		csvFile.write('\n')
	csvFile.close()
