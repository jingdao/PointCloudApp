#!/usr/bin/python

import sys
from sklearn.neighbors import NearestNeighbors
import numpy

if len(sys.argv) < 4:
	print "%s train/ test/ {fpfh,spin}\n" % sys.argv[0]
	sys.exit(0)

def parseLabelFile(filename):
	f = open(filename,'r')
	data = f.read()
	return numpy.array([int(l) for l in data.split()])

def parseDescFile(filename,n):
	f = open(filename,'r')
	numSamples = n
	while True:
		line = f.readline()
		if line.startswith('POINTS'):
			numPoints = int(line.split()[1])
			break
	f.readline()
	numSamples = min(numSamples,numPoints)
	samples = set(numpy.random.choice(range(numPoints),numSamples,False))
	features = []
	for i in range(numPoints):
		data = f.readline()
		if i in samples:
			data = [float(x) for x in data.split()]
			features.append(data)
	return features

trainLabels = parseLabelFile(sys.argv[1]+'/labels.txt')
testLabels = parseLabelFile(sys.argv[2]+'/labels.txt')
numSamples = 100
trainFeatures = []
for i in range(len(trainLabels)):
	f = parseDescFile('%s/%d-cloud.pcd-%s.pcd' % (sys.argv[1],i,sys.argv[3]), numSamples)
	trainFeatures.extend(f)
trainFeatures = numpy.array(trainFeatures)
print 'Parsed %d features from %d objects' % (len(trainFeatures),len(trainLabels))

categories = list(set(trainLabels))
TP = {c:0 for c in categories}
countMembers = {c:0 for c in categories}
nbrs = {}
for c in categories:
	mask = numpy.repeat(trainLabels==c,numSamples)
	X = trainFeatures[mask]
	nbrs[c] = NearestNeighbors(n_neighbors=1,algorithm='brute').fit(X)
	
for i in range(len(testLabels)):
	f = parseDescFile('%s/%d-cloud.pcd-%s.pcd' % (sys.argv[2],i,sys.argv[3]), numSamples)
	truth = testLabels[i]
	countMembers[truth] += 1
	prediction = None
	minDist = 0
	for c in categories:
		dist,index = nbrs[c].kneighbors(f)
		sumDist = sum(dist)
		if prediction is None or sumDist < minDist:
			prediction = c
			minDist = sumDist
	if prediction == truth:
		TP[truth] += 1

for i in categories:
	if not countMembers[i]==0:
		print "class %d (%2d samples): %2d (%.3f)" % (i,countMembers[i],TP[i],1.0 * TP[i] / countMembers[i])
totalMembers = sum(countMembers.values())
totalTP = sum(TP.values())
print "overall (%d samples): %2d (%.3f)" % (totalMembers,totalTP,1.0 * totalTP / totalMembers)
