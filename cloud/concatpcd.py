#!/usr/bin/python

import sys

if len(sys.argv) < 3:
	print './concatpcd.py [input1.pcd input2.pcd ...] [output.pcd]'
	sys.exit(1)

infiles = []
for i in range(1,len(sys.argv)-1):
	infiles.append(open(sys.argv[i],'r'))
outfile = sys.argv[-1]
output=open(outfile,'w')

numPoints=[]
header=[]
readHeader=True
for f in infiles:
	while True:
		l = f.readline()
		if readHeader:
			header.append(l)
		if l.startswith("POINTS"):
			numPoints.append(int(l.split()[1]))
		elif l.startswith("DATA"):
			break
	readHeader=False

totalPoints=sum(numPoints)
for l in header:
	if l.startswith("WIDTH"):
		output.write("WIDTH "+str(totalPoints)+"\n")
	elif l.startswith("POINTS"):
		output.write("POINTS "+str(totalPoints)+"\n")
	else:
		output.write(l)

for i in range(len(infiles)):
	for j in range(numPoints[i]):
		output.write(infiles[i].readline())
