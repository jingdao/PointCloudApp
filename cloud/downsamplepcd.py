#!/usr/bin/python

import sys
import math

if len(sys.argv) < 4:
	print './downsamplepcd.py [input.pcd] [output.pcd] [resolution]'
	sys.exit(1)

infile = sys.argv[1]
outfile = sys.argv[2]
resolution = float(sys.argv[3])

input=open(infile,'r')
output=open(outfile,'w')

header=[]
while True:
	l = input.readline()
	header.append(l)
	if l.startswith('DATA'):
		break

x = []
y = []
z = []
attr = []
while True:
	l = input.readline().split()
	if len(l) < 3:
		break
	x.append(float(l[0]))
	y.append(float(l[1]))
	z.append(float(l[2]))
	attr.append(l[3:])

xd = 0.5 * (min(x) + max(x))
yd = 0.5 * (min(y) + max(y))
zd = 0.5 * (min(z) + max(z))
pointset = set()
points = []

for i in range(len(x)):
	p = (int((x[i]-xd)/resolution),int((y[i]-yd)/resolution),int((z[i]-zd)/resolution))
	if not p in pointset:
		pointset.add(p)
		vals = [x[i],y[i],z[i]]
		vals.extend(attr[i])
		points.append(vals)

for l in header:
	if l.startswith("WIDTH"):
		output.write("WIDTH "+str(len(points))+"\n")
	elif l.startswith("POINTS"):
		output.write("POINTS "+str(len(points))+"\n")
	else:
		output.write(l)

for p in points:
	for v in p:
		output.write(str(v)+' ')
	output.write('\n')

