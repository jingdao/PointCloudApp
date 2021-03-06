#!/usr/bin/python

import sys
import glob
import os
dir = sys.argv[1]
desc = 'esf' if len(sys.argv) < 3 else sys.argv[2]

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
elif len(sys.argv) >= 4:
	maxN = int(sys.argv[3])
else:
	maxN = len(labels)
for n in range(maxN):
	f = dir + '/'+str(n)+'-cloud.pcd-'+desc+'.pcd'
	try:
		fd = open(f,'r')
	except IOError:
		break
	if labels is None:
		outfile.write('0 ')
	else:
		outfile.write(str(labels[n])+' ')
	for i in range(0,11):
		fd.readline()
	features = fd.readline().split()
	n = 1
	for t in features:
		outfile.write(str(n)+':'+str(float(t))+' ')
		n += 1
	outfile.write('\n')
	
print 'Wrote '+str(maxN)+' labels to '+outfile.name
outfile.close()
	

