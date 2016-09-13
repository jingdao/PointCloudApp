#!/usr/bin/python
import sys
import os
import math

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

#train classifier
gm_mean = {}
gm_std = {}
gm_count = {}
for i in range(len(features['train'])):
	f = features['train'][i]
	l = labels['train'][i]
	if not l in gm_count:
		gm_count[l] = 0
		gm_mean[l] = [0] * numFeatures
		gm_std[l] = [0] * numFeatures
	gm_count[l] += 1
	for j in range(numFeatures):
		gm_mean[l][j] += f[j]
		gm_std[l][j] += f[j] * f[j]

for l in gm_count:
	for j in range(numFeatures):
		gm_mean[l][j] /= gm_count[l]
		gm_std[l][j] = math.sqrt(gm_std[l][j] / gm_count[l] - gm_mean[l][j] * gm_mean[l][j])

#test classifier
prediction = [0] * len(features['test'])
proba = [[1 for i in range(len(classes))] for j in range(len(features['test']))]
score = 0
for i in range(len(features['test'])):
	f = features['test'][i]
	l = labels['test'][i]
	max_proba = 0
	for k in gm_mean:
		for j in range(numFeatures):
			if not gm_std[k][j]==0:
				z = (f[j]- gm_mean[k][j]) / gm_std[k][j]
				if z > 0:
					p = 0.5 * math.erfc(z / math.sqrt(2))
				else:
					p = 0.5 + 0.5 * math.erf(z / math.sqrt(2))
				proba[i][k-1] *= p
			else:
				proba[i][k-1] *= 1 if f[j]==0 else 0
		if proba[i][k-1] > max_proba:
			max_proba = proba[i][k-1]
			prediction[i] = k
	if prediction[i] == l:
		score += 1
print 'Accuracy %.2f%%' % (1.0 * score / len(features['test'])*100)


#output data
file = open(outputFile,'w')
file.write('labels ')
for c in classes:
	file.write(str(c)+' ')
file.write('\n')
for i in range(len(prediction)):
	l = prediction[i]
	file.write(str(l))
	for p in proba[i]:
		file.write(' '+str(p))
	file.write('\n')
