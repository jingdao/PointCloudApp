#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

x = {}
y = {}
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
		if not labels[i] in x:
			x[labels[i]] = np.arange(640)
			y[labels[i]] = []
		vals=[]
		for j in range(640):
			vals.append(float(line[j]))
		y[labels[i]].append(vals)

if len(args) > 0:
	keys = args
else:
	keys = x.keys()
for key in keys:
	mean[key] = np.mean(y[key],axis=0)
	stddev[key] = np.std(y[key],axis=0)
	plt.plot(x[key],mean[key],lw=2,label=categories[key],color=colors[key])
	plt.fill_between(x[key],mean[key]+stddev[key],mean[key]-stddev[key],facecolor=colors[key],alpha=0.5)

plt.legend(loc='upper left')
#plt.axis([0, 640, 0, 0.02])
plt.show()



