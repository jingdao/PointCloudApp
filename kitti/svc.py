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
use_linear = True
plot_training_size = False
plot_weights = False

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

if plot_weights:
	colors = {0:'#ffffff',1:'#ff0000',2:'#00ff00',3:'#0000ff',4:'#ffff00',5:'#ff00ff',6:'#00ffff'}
	import matplotlib.pyplot as plt
	for i in range(len(svc.coef_)):
		w = svc.coef_[i]
		plt.plot(np.arange(len(w)),w,lw=2,color=colors[i])
	plt.legend(loc='upper left')
	plt.show()


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
