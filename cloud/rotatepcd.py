#!/usr/bin/python

import sys
from math import cos,sin,pi

if len(sys.argv) < 4:
	print './rotatepcd.py [input.pcd] [output.pcd] [rx] [ry] [rz]'
	sys.exit(1)

infile = sys.argv[1]
outfile = sys.argv[2]
rx = float(sys.argv[3])
ry = float(sys.argv[4])
rz = float(sys.argv[5])
r1 = rx / 180 * pi
r2 = ry / 180 * pi
r3 = rz / 180 * pi
M = [0] * 9
M[0] = cos(r2)*cos(r3);
M[1] = - cos(r1) * sin(r3) + cos(r3) * sin(r1) * sin(r2);
M[2] = sin(r1) * sin(r3) + cos(r1) * cos(r3) * sin(r2);
M[3] = cos(r2) * sin(r3);
M[4] = cos(r1) * cos(r3) + sin(r1) * sin(r2) * sin(r3);
M[5] = - cos(r3) * sin(r1) + cos(r1) * sin(r2) * sin(r3);
M[6] = -sin(r2);
M[7] = cos(r2) * sin(r1);
M[8] = cos(r1) * cos(r2);

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
	xi = float(l[0])
	yi = float(l[1])
	zi = float(l[2])
	x = M[0]*xi + M[1]*yi + M[2]*zi
	y = M[3]*xi + M[4]*yi + M[5]*zi
	z = M[6]*xi + M[7]*yi + M[8]*zi
	vals = [x,y,z]
	vals.extend(l[3:])
	for v in vals:
		output.write(str(v)+' ')
	output.write('\n')
