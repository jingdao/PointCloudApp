#!/usr/bin/python

import os
import glob
import sys

if len(sys.argv) >=2:
	inputDir = sys.argv[1]
else:
	inputDir = '.'
	
if "-i" in sys.argv:
	ignoreZeros = True
else:
	ignoreZeros = False
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

for d in glob.glob(inputDir+'/*'):
	if os.path.isdir(d) and os.path.isfile(d+'/labels.txt'):
		numFiles += 1
		if len(categories) == 0 and os.path.isfile(d+'/labelCategory.txt'):
			categoryFile = open(inputDir+'/labelCategory.txt')
			for l in categoryFile:
				categories.append(l.split()[1])
			categoryFile.close()	
		labels = open(d+'/labels.txt','r')
		prediction = open(d+'/prediction.txt','r')
		if ignoreZeros:
			while True:
				try:
					l = int(labels.readline())
				except ValueError:
					break
				p = int(prediction.readline())
				if l==0 or p==0:
					continue
				numTrials += 1
				if l == p:
					increment(truePositive,l)
				else:
					increment(falsePositive,p)
					increment(falseNegative,l)
		else:
			while True:
				try:
					l = int(labels.readline())
				except ValueError:
					break
				p = int(prediction.readline())
				numTrials += 1
				if l == p:
					increment(truePositive,l)
				else:
					increment(falsePositive,p)
					increment(falseNegative,l)
		labels.close()
		prediction.close()

print "Processed " + str(numFiles) + " files (" + str(numTrials) + " trials)"
print "%2s %15s %4s %8s %10s %6s" % ("","Category","Num","Accuracy","Precision","Recall")
for i in truePositive:
	TP = truePositive[i] if i in truePositive else 0
	FP = falsePositive[i] if i in falsePositive else 0
	FN = falseNegative[i] if i in falseNegative else 0
	num = TP + FN
	acc = (numTrials - FP - FN) * 100.0 / numTrials
	pre = TP * 100.0 / (TP + FP)
	rec = TP * 100.0 / (TP + FN)
	print "%2d %15s %4d %8.2f %10.2f %6.2f" % (i,categories[i],num,acc,pre,rec)

totalAcc = sum(truePositive.values()) * 100.0 / numTrials
totalPre = sum(truePositive.values()) * 100.0 / (sum(truePositive.values()) + sum(falsePositive.values()))
totalRec = sum(truePositive.values()) * 100.0 / (sum(truePositive.values()) + sum(falseNegative.values()))
print "%2s %15s %4d %8.2f %10.2f %6.2f" % ("","overall",numTrials,totalAcc,totalPre,totalRec)
