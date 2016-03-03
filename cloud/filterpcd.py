#!/usr/bin/python

import sys

if len(sys.argv) < 9:
	print './filterpcd.py [input.pcd] [output.pcd] [minX] [minY] [minZ] [maxX] [maxY] [maxZ]'
	sys.exit(1)

infile = sys.argv[1]
outfile = sys.argv[2]
minX = float(sys.argv[3])
minY = float(sys.argv[4])
minZ = float(sys.argv[5])
maxX = float(sys.argv[6])
maxY = float(sys.argv[7])
maxZ = float(sys.argv[8])

input=open(infile,'r')
output=open(outfile,'w')

header=[]
while True:
	l = input.readline()
	header.append(l)
	if l.startswith('DATA'):
		break

data=[]
while True:
	l = input.readline().split()
	if len(l) < 3:
		break
	x = float(l[0])
	y = float(l[1])
	z = float(l[2])
	if x>=minX and x<=maxX and y>=minY and y<=maxY and z>=minZ and z<=maxZ:
		s=""
		vals = [x,y,z]
		vals.extend(l[3:])
		for v in vals:
			s+=str(v)+' '
		data.append(s)

for l in header:
	if l.startswith("WIDTH"):
		output.write("WIDTH "+str(len(data))+"\n")
	elif l.startswith("POINTS"):
		output.write("POINTS "+str(len(data))+"\n")
	else:
		output.write(l)

for l in data:
	output.write(l+'\n')
