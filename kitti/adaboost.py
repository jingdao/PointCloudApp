#!/usr/bin/python

import numpy as np
from sklearn.ensemble import AdaBoostClassifier
import sys
import os
import random

if len(sys.argv) > 1:
	dr = sys.argv[1]
else:
	dr = '.'
inputFile = {'train':dr+'/svm_train_scaled.txt', 'test':dr+'/svm_test_scaled.txt'}
outputFile = dr+'/svm_prediction.txt'

#import data into matrix format
features = {'train':[], 'test':[]}
sparse_features = {'train':[], 'test':[]}
labels = {'train':[], 'test':[]}
indices=[]
for dataset in ['train', 'test']:
	file = open(inputFile[dataset])
	for line in file:
		tokens = line.split()
		labels[dataset].append(int(tokens[0]))
		f = {}
		for i in range(1,len(tokens)):
			ind = int(tokens[i].split(':')[0])
			val = float(tokens[i].split(':')[1])
			f[ind] = val
		sparse_features[dataset].append(f)
		indices.extend(f.keys())

numFeatures = max(indices)
for dataset in ['train','test']:
	for sf in sparse_features[dataset]:
		f = [0] * numFeatures
		for ind in sf:
			f[ind-1] = sf[ind]
		features[dataset].append(f)

classes = set(labels['train'])

bdt = AdaBoostClassifier()
bdt.fit(features['train'],labels['train'])

#test classifier
prediction = bdt.predict(features['test'])
proba = bdt.predict_proba(features['test'])
print 'Accuracy %.2f%%' % (bdt.score(features['test'],labels['test'])*100)

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
