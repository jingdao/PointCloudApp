#!/usr/bin/python

import sys

if len(sys.argv) < 3:
	print './concat_esf.py [input1.pcd input2.pcd ...] [output.pcd]'
	sys.exit(1)

infiles = []
for i in range(1,len(sys.argv)-1):
	infiles.append(open(sys.argv[i],'r'))
outfile = sys.argv[-1]
output=open(outfile,'w')

numPoints=[]
header=[]
features=""
readHeader=True
for f in infiles:
	while True:
		l = f.readline()
		if readHeader:
			header.append(l)
		if l.startswith("COUNT"):
			numPoints.append(int(l.split()[1]))
		elif l.startswith("DATA"):
			l = f.readline()
			features += l[:-1] + " "
			break
	readHeader=False

totalPoints=sum(numPoints)
for l in header:
	if l.startswith("COUNT"):
		output.write("COUNT "+str(totalPoints)+"\n")
	else:
		output.write(l)

output.write(features+"\n")
