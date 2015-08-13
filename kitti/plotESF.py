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
		x[labels[i]] = []
		y[labels[i]] = []
	for j in range(640):
		x[labels[i]].append(j)
		y[labels[i]].append(float(line[j]))

if len(args) > 0:
	for key in args:
		plt.plot(x[key],y[key],'.',label=categories[key],markersize=5)
else:
	for key in x:
		plt.plot(x[key],y[key],'.',label=categories[key],markersize=5)

plt.legend(loc='upper left')
plt.axis([0, 640, 0, 0.02])
plt.show()



