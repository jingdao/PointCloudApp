#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

x={}
y={}
Y = []
mean = {}
stddev = {}
colors = {0:'white',1:'red',2:'green',3:'blue',4:'yellow',5:'purple',6:'orange'}
svmdata_file = None
args=[]
if len(sys.argv) > 1:
	if sys.argv[1].endswith('txt'):
		svmdata_file = sys.argv[1]
	else:
		for i in range(1,len(sys.argv)):
			args.append(int(sys.argv[i]))

categories = []
if os.path.isfile('labelCategory.txt'):
	categoryFile = open('labelCategory.txt')
else:
	categoryFile = open('../labelCategory.txt')
for l in categoryFile:
	categories.append(l)
categoryFile.close()
	
if svmdata_file:
	labels=[]
	file = open(svmdata_file)
	for line in file:
		tokens = line.split()
		l = int(tokens[0])
		labels.append(l)
		f = [0] * 640
		for i in range(1,len(tokens)):
			ind = int(tokens[i].split(':')[0])
			val = float(tokens[i].split(':')[1])
			f[ind - 1] = val
		if not l in x:
			x[l] = np.arange(640)
			y[l] = []
		y[l].append(f)
		Y.append(f)
else:
	labels = []
	labelFile = open('labels.txt')
	for l in labelFile:
		labels.append(int(l))
	labelFile.close()

	for i in range(len(labels)):
		f = open(str(i)+"-cloud.pcd-esf.pcd")
		for j in range(11):
			f.readline()
		line = f.readline().split()
		vals=[]
		for j in range(640):
			vals.append(float(line[j]))
		if not labels[i] in x:
			x[labels[i]] = np.arange(640)
			y[labels[i]] = []
		y[labels[i]].append(vals)
		Y.append(vals)

Y = np.array(Y)
mean = np.mean(Y,axis=0)
U,S,V = np.linalg.svd(Y)
v1 = V[0]
v2 = V[1]


if len(args) > 0:
	keys = args
else:
	keys = x.keys()
for key in keys:
	z1=[]
	z2=[]
	for i in range(len(y[key])):
		z1.append(np.dot(v1,np.array(y[key][i])-mean) / S[0])
		z2.append(np.dot(v2,np.array(y[key][i])-mean) / S[1])
	plt.plot(z1,z2,'o',label=categories[key],color=colors[key])

plt.legend(loc='upper left')
#plt.axis([0, 640, 0, 0.02])
plt.show()



