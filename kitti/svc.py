#!/usr/bin/python

import numpy as np
from sklearn.svm import SVC, LinearSVC
import sys
import os
import random

if len(sys.argv) > 1:
	dr = sys.argv[1]
else:
	dr = '.'
inputFile = {'train':dr+'/svm_train_scaled.txt', 'test':dr+'/svm_test_scaled.txt'}
outputFile = dr+'/svm_prediction.txt'
numFeatures = 640
use_linear = True
plot_training_size = True

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

if plot_training_size:
	initial_size = 2
	step = 16
	db_f = features['train']
	db_l = labels['train']
	features['train'] = []
	labels['train'] = []
	for i in range(initial_size):
		j = random.randrange(len(db_f) / step)
		features['train'].extend(db_f[j*step:(j+1)*step])
		labels['train'].extend(db_l[j*step:(j+1)*step])
		del db_f[j*step:(j+1)*step]
		del db_l[j*step:(j+1)*step]
	while len(db_f) > 0:
		j = random.randrange(len(db_f) / step)
		features['train'].extend(db_f[j*step:(j+1)*step])
		labels['train'].extend(db_l[j*step:(j+1)*step])
		del db_f[j*step:(j+1)*step]
		del db_l[j*step:(j+1)*step]
		if use_linear:
			svc = LinearSVC(random_state=0)
		else:
			svc = SVC(C=2.0,gamma=0.0078,probability=True)
		svc.fit(features['train'],labels['train'])
		acc = svc.score(features['test'],labels['test'])
		print len(features['train']), '%.2f' % (acc*100)
else:
	#train classifier
	if use_linear:
		svc = LinearSVC(random_state=0)
	else:
		svc = SVC(C=2.0,gamma=0.0078,probability=True)
	svc.fit(features['train'],labels['train'])

#test classifier
prediction = svc.predict(features['test'])
if use_linear:
	proba = svc.decision_function(features['test'])
	for i in range(len(proba)):
		proba[i] = 1 / (1 + np.exp(-proba[i]))
		proba[i] = proba[i] / sum(proba[i])
#		if max(proba[i]) < 0.4:
#			prediction[i] = 0
else:
	proba = svc.predict_proba(features['test'])
print 'Accuracy %.2f%%' % (svc.score(features['test'],labels['test'])*100)

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
