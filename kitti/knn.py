#!/usr/bin/python
#perform KNN classification on data in libSVM format

from sklearn.neighbors import KNeighborsClassifier
import sys
import numpy as np

inputFile = {'train':'svm_train_scaled.txt', 'test':'svm_test_scaled.txt'}
outputFile = 'svm_prediction.txt'
numFeatures = 640
k = 3

if len(sys.argv) > 1:
	k = int(sys.argv[1])
if len(sys.argv) > 2:
	inputFile['train'] = sys.argv[2]
if len(sys.argv) > 3:
	inputFile['test'] = sys.argv[3]
if len(sys.argv) > 4:
	outputFile = sys.argv[4]
	
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

#print(features)
#print(labels)
	
#train classifier
neigh = KNeighborsClassifier(k,p=1)
neigh.fit(features['train'],labels['train'])

#test classifier
prediction = neigh.predict(features['test'])
print('Accuracy %.2f%%' % neigh.score(features['test'],labels['test']))

#output data
file = open(outputFile,'w')
file.write('labels\n')
for l in prediction:
	file.write(str(l)+'\n')
