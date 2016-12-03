#!/usr/bin/python

import sys
if len(sys.argv) < 3:
	print 'registerpcd.py matches.txt output.pcd'
	sys.exit(0)

import math
def transform(points,p):
	theta = p[0]
	tx = p[1]
	ty = p[2]
	tz = p[3]
	res = []
	for q in points:
		x = math.cos(theta) * q[0] -math.sin(theta) * q[1] + tx
		y = math.sin(theta) * q[0] +math.cos(theta) * q[1] + ty
		z = q[2] + tz
		res.append((x,y,z))
	return res

def getParam(q,p):
	t1 = math.atan2(p[1][0]-p[0][0],p[1][1]-p[0][1])
	t2 = math.atan2(q[1][0]-q[0][0],q[1][1]-q[0][1])
	theta = t2 - t1
	n = len(p)
	Tx = 0
	Ty = 0
	Tz = 0
	for i in range(n):
		Tx += q[i][0] - p[i][0]*math.cos(theta) + p[i][1]*math.sin(theta)
		Ty += q[i][1] - p[i][0]*math.sin(theta) - p[i][1]*math.cos(theta)
		Tz += q[i][2] - p[i][2]
	Tx /= n
	Ty /= n
	Tz /= n
	return (theta,Tx,Ty,Tz)

def savePCD(filename,points):
	f = open(filename,"w")
	l = len(points)
	header = """# .PCD v0.7 - Point Cloud Data file format
VERSION 0.7
FIELDS x y z
SIZE 4 4 4
TYPE F F F
COUNT 1 1 1
WIDTH %d
HEIGHT 1
VIEWPOINT 0 0 0 1 0 0 0
POINTS %d
DATA ascii
""" % (l,l)
	f.write(header)
	for p in points:
		f.write("%f %f %f\n"%(p[0],p[1],p[2]))
	f.close()
	print 'Saved %d points to %s' % (l,filename)

matches = open(sys.argv[1],'r')
source = []
param = []
previous = None
for l in matches:
	ll=l.split()
	source.append(ll[0])
	if previous is None:
		previous = []
		for p in ll[1:]:
			q = p.split(',')
			previous.append((float(q[0]),float(q[1]),float(q[2])))
	else:
		current = []
		for p in ll[1:]:
			if p=='null':
				current.append(None)
			else:
				q = p.split(',')
				current.append((float(q[0]),float(q[1]),float(q[2])))
		p = []
		q = []
		for i in range(len(previous)):
			if previous[i] is not None and current[i] is not None:
				p.append(previous[i])
				q.append(current[i])
		currentParam = getParam(p,q)
		param.append(currentParam)
		previous = current
matches.close()

print param

output = []
for i in range(len(source)):
	f=open(source[i],'r')
	for l in f:
		if l.startswith('DATA'):
			break
	points = []
	for l in f:
		q = l.split()
		points.append((float(q[0]),float(q[1]),float(q[2])))
	if not i==0:
		points = transform(points,param[i-1])
	output.extend(points)
	f.close()

savePCD(sys.argv[2],output)


