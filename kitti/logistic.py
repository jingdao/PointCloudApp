#!/usr/bin/python

import numpy as np
from sklearn import linear_model
import sys
import os

if len(sys.argv) > 1:
	dr = sys.argv[1]
else:
	dr = '.'
inputFile = {'train':dr+'/svm_train_scaled.txt', 'test':dr+'/svm_test_scaled.txt'}
outputFile = dr+'/svm_prediction.txt'
numFeatures = 640
useBinary = False

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

if useBinary:
	prediction = []
	proba=[]
	for c in classes:
		lr = linear_model.LogisticRegression(C=1,penalty='l1')
		lr.fit(features['train'],np.array(labels['train'])==c)
		p = np.array(lr.predict_proba(features['test']))
		proba.append(p[:,1])
	proba = np.transpose(np.array(proba))
	isValid = np.max(proba,axis=1) > 0.5
	prediction = (np.argmax(proba,axis=1)+1) * isValid.astype(np.int)
	for row in proba:
		row /= np.sum(row)
else:
	#train classifier
	lr = linear_model.LogisticRegression(C=1e6,penalty='l1')
	lr.fit(features['train'],labels['train'])
	#test classifier
	prediction = lr.predict(features['test'])
	proba = lr.predict_proba(features['test'])
	print('Accuracy %.2f%%' % lr.score(features['test'],labels['test']))

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
