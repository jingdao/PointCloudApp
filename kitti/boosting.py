#!/usr/bin/python

import numpy as np
from sklearn.svm import SVC, LinearSVC
import sys
import os
import random

if len(sys.argv) < 3:
	print './boosting.py trainDir/ testDir/'
	sys.exit(1)

trainDir = sys.argv[1]
testDir = sys.argv[2]
numParts = 5
random.seed(0)

class Node:
	def __init__(self,l,r,p,c,s,v):
		self.left=l
		self.right=r
		self.part=p
		self.label=c
		self.score=s
		self.svm=v

	def __str__(self):
		return "Node: part "+str(self.part)+" class "+str(self.label)+" score "+str(self.score)

def getLabels(labeldir):
	l=[]
	f=open(labeldir+'/labels.txt')
	for s in f:
		l.append(int(s))
	return l

def getFeatures(featuredir,sample_id,part_id):
	try:
		s=featuredir+'/'+str(sample_id)+'-cloud.pcd-'+str(part_id)+'-part.pcd-esf.pcd'
		f=open(s)
	except IOError:
		return None
	for i in range(11):
		f.readline()
	s=f.readline()
	l=[]
	for t in s.split():
		l.append(float(t))
	return l

def linearScaling(ft):
	min_acc=[]
	max_acc=[]
	for t in ft[0]:
		min_acc.append(t)
		max_acc.append(t)
	for i in range(1,len(ft)):
		for j in range(len(ft[i])):
			if ft[i][j] < min_acc[j]:
				min_acc[j] = ft[i][j]
			if ft[i][j] > max_acc[j]:
				max_acc[j] = ft[i][j]
	new_ft=[]
	for i in range(len(ft)):
		t=[]
		for j in range(len(ft[i])):
			if max_acc[j] == min_acc[j]:
				t.append(ft[i][j])
			else:
				t.append((ft[i][j] - min_acc[j]) / (max_acc[j] - min_acc[j]))
		new_ft.append(t)
	return new_ft

def gaussianScaling(ft):
	m = np.mean(ft,axis=0)
	s = np.std(ft,axis=0)
	new_ft=[]
	for i in range(len(ft)):
		t=[]
		for j in range(len(ft[i])):
			if s[j] > 0:
				t.append((ft[i][j]-m[j]) / s[j])
			else:
				t.append(ft[i][j])
		new_ft.append(t)
	return new_ft

categories = []
categoryFile = open(trainDir+'/../labelCategory.txt')
for l in categoryFile:
	categories.append(l.split()[1])
categoryFile.close()

trainLabels=getLabels(trainDir)
classes = set(trainLabels)
trainData={}
for p in range(numParts):
	trainData[p]=[]
	i=0
	while True:
		ft = getFeatures(trainDir,i,p)
		if not ft:
			break
		trainData[p].append(ft)
		i+=1
	trainData[p]=gaussianScaling(trainData[p])

testLabels=getLabels(testDir)
testData={}
for p in range(numParts):
	testData[p]=[]
	i=0
	while True:
		ft = getFeatures(testDir,i,p)
		if not ft:
			break
		testData[p].append(ft)
		i+=1
	testData[p]=gaussianScaling(testData[p])

def makeNode(trainData,trainLabels,visited):
	bestScore=0
	bestPart=None
	bestClass=None
	bestClassifier=None
	for p in range(numParts):
		for c in classes:
			if c in visited:
				continue
			subLabel = [l==c for l in trainLabels]
			while True:
				l1 = []
				f1 = []
				l2 = []
				f2 = []
				validation_index = random.sample(range(len(subLabel)),len(subLabel)/4)
				for i in range(len(subLabel)):
					if not i in validation_index:
						l1.append(subLabel[i])
						f1.append(trainData[p][i])
					else:
						l2.append(subLabel[i])
						f2.append(trainData[p][i])
				if any(l1) or not all(l1):
					break
			svc = LinearSVC(random_state=0)
			svc.fit(f1,l1)
			pd = svc.predict(f2)
			totalSamples = sum(l2)
			correctSamples = 0
			for i in range(len(pd)):
				if l2[i] and pd[i]:
					correctSamples += 1
#			score = 1.0 * correctSamples / totalSamples 
			score = correctSamples
			print p,c,score
			if score > bestScore:
				bestScore = score
				bestPart = p
				bestClass = c
				bestClassifier = svc
	return Node(None,None,bestPart,bestClass,bestScore,bestClassifier)

v=[]
nodes=[]
for i in range(len(classes)-1):
	currentData={}
	currentLabels=[]
	for p in range(numParts):
		currentData[p]=[]
	for j in range(len(trainLabels)):
		if not trainLabels[j] in v:
			for p in range(numParts):
				currentData[p].append(trainData[p][j])
			currentLabels.append(trainLabels[j])
	n = makeNode(currentData,currentLabels,v)
	nodes.append(n)
	v.append(n.label)


pd = nodes[0].svm.predict(testData[nodes[0].part])
totalRecall = 0
totalSamples = 0
for i in range(len(testLabels)):
	if testLabels[i] == n.label:
		totalSamples += 1
		if pd[i]:
			totalRecall += 1
print 1.0 * totalRecall / totalSamples
