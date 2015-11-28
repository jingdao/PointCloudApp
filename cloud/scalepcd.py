#!/usr/bin/python

import sys

if len(sys.argv) < 4:
	print './scalepcd.py [input.pcd] [output.pcd] [scale_factor]'
	sys.exit(1)

infile = sys.argv[1]
outfile = sys.argv[2]
factor = float(sys.argv[3])

input=open(infile,'r')
output=open(outfile,'w')

while True:
	l = input.readline()
	output.write(l)
	if l.startswith('DATA'):
		break

while True:
	l = input.readline().split()
	if len(l) < 3:
		break
	x = float(l[0]) / factor
	y = float(l[1]) / factor
	z = float(l[2]) / factor
	vals = [x,y,z]
	vals.extend(l[3:])
	for v in vals:
		output.write(str(v)+' ')
	output.write('\n')
