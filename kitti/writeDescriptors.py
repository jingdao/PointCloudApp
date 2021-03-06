#!/usr/bin/python

import sys
import glob
import os
dir = sys.argv[1]

labels = None
if os.path.isfile(dir+'/labels.txt'):
	labels = []
	fd = open(dir+'/labels.txt')
	for line in fd:
		labels.append(int(line))
	print 'Read '+str(len(labels))+' labels'
	
outfile = open(dir+'/svmdata.txt','w')
if labels is None:
	maxN = 100
else:
	maxN = len(labels)
	
descriptorFile = open(dir+'/descriptor.pcd')
while not descriptorFile.readline().startswith("DATA ascii"):
    pass

for n in range(maxN):
	if labels is None:
		outfile.write('0 ')
	else:
		outfile.write(str(labels[n])+' ')
	features = descriptorFile.readline().split()
	n = 1
	for t in features:
		outfile.write(str(n)+':'+str(float(t))+' ')
		n += 1
	outfile.write('\n')
	
print 'Wrote labels to '+outfile.name
outfile.close()
	

