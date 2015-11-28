#!/usr/bin/python

import sys

if len(sys.argv) < 4:
	print './translatepcd.py [input.pcd] [output.pcd] [x] [y] [z]'
	sys.exit(1)

infile = sys.argv[1]
outfile = sys.argv[2]
xt = float(sys.argv[3])
yt = float(sys.argv[4])
zt = float(sys.argv[5])

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
	x = float(l[0]) + xt
	y = float(l[1]) + yt
	z = float(l[2]) + zt
	vals = [x,y,z]
	vals.extend(l[3:])
	for v in vals:
		output.write(str(v)+' ')
	output.write('\n')
