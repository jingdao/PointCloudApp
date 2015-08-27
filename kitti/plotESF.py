#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

args=[]
if len(sys.argv) > 1:
	for i in range(1,len(sys.argv)):
		args.append(int(sys.argv[i]))

x = {}
y = {}
mean = {}
stddev = {}
colors = {0:'white',1:'red',2:'green',3:'blue'}

labels = []
labelFile = open('labels.txt')
for l in labelFile:
	labels.append(int(l))
labelFile.close()

categories = []
if os.path.isfile('labelCategory.txt'):
	categoryFile = open('labelCategory.txt')
else:
	categoryFile = open('../labelCategory.txt')
for l in categoryFile:
	categories.append(l)
categoryFile.close()
	
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
plt.axis([0, 640, 0, 0.02])
plt.show()



