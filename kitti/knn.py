#!/usr/bin/python
#perform KNN classification on data in libSVM format

from sklearn.neighbors import KNeighborsClassifier
import sys
import numpy as np
import os

inputFile = {'train':'./svm_train_scaled.txt', 'test':'./svm_test_scaled.txt'}
outputFile = './svm_prediction.txt'
numFeatures = 640
k = 1

if len(sys.argv) > 1:
	k = int(sys.argv[1])
if len(sys.argv) > 2:
	inputFile['train'] = sys.argv[2]
if len(sys.argv) > 3:
	inputFile['test'] = sys.argv[3]
if len(sys.argv) > 4:
	outputFile = sys.argv[4]
outputDir = os.path.dirname(outputFile)
	
#import data into matrix format
features = {'train':[], 'test':[]}
labels = {'train':[], 'test':[]}
for dataset in ['train', 'test']:
	file = open(inputFile[dataset])
	for line in file:
		tokens = line.split()
		labels[dataset].append(int(tokens[0]))
		f = [0] * numFeatures
		for i in range(1,len(tokens)):
			ind = int(tokens[i].split(':')[0])
			val = float(tokens[i].split(':')[1])
			f[ind - 1] = val
		features[dataset].append(f)

classes = set(labels['train'])
#print(features)
#print(labels)
	
#train classifier
neigh = KNeighborsClassifier(k,p=1,weights='distance')
neigh.fit(features['train'],labels['train'])

#test classifier
prediction = neigh.predict(features['test'])
proba = neigh.predict_proba(features['test'])
#for i in range(len(prediction)):
#	q = sorted(proba[i],reverse=True)
#	if q[0] - q[1] < 0.3:
#		prediction[i] = 0
print('Accuracy %.2f%%' % (neigh.score(features['test'],labels['test']) * 100.0))
#prediction=[]
#proba=[]
#dist, match = neigh.kneighbors(features['test'])
#numCategories = max(labels['train'])
#for i in range(len(labels['test'])):
#	score = {}
#	count = {}
#	p = []
#	for j in range(len(match[i])):
#		l = labels['train'][match[i][j]]
#		if l in score:
#			score[l] += 1/dist[i][j]
#			count[l] += 1
#		else:
#			score[l] = 1/dist[i][j]
#			count[l] = 1
#	for key in range(numCategories):
#		if key in score:
#			p.append(score[key] / count[key])
#		else:
#			p.append(0)
#	p = np.array(p)
#	proba.append(p / np.sum(p))
#	prediction.append(np.argmax(p)+1)

#find neighbors
neighborFile = open(outputDir+'/neighbors.txt','w')
for i in range(len(features['test'])):
	dist, match = neigh.kneighbors(features['test'][i])
	neighborFile.write(str(labels['test'][i])+' '+str(match[0][0])+' '+str(dist[0][0])+'\n')
neighborFile.close()

#output data
file = open(outputFile,'w')
file.write('labels ')
for c in classes:
	file.write(str(c)+' ')
file.write('\n')
for i in range(len(prediction)):
	l = prediction[i]
	file.write(str(l)+' ')
	for p in proba[i]:
		file.write(str(p)+' ')
	file.write('\n')
