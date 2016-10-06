#!/usr/bin/python

import sys
import random

if len(sys.argv) < 4:
	print './addnoise.py [input.pcd] [output.pcd] [std_dev]'
	sys.exit(1)

infile = sys.argv[1]
outfile = sys.argv[2]
std = float(sys.argv[3])

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
	x = float(l[0]) + random.gauss(0,std)
	y = float(l[1]) + random.gauss(0,std)
	z = float(l[2]) + random.gauss(0,std)
	vals = [x,y,z]
	vals.extend(l[3:])
	for v in vals:
		output.write(str(v)+' ')
	output.write('\n')
