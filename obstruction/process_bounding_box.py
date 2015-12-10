#!/usr/bin/python
import sys
import numpy
#filter based on standard deviation
f = open('bounding_box.txt','r')
g = open('bounding_box_filtered.txt','w')

x=[]
y=[]
z=[]
lines=[]
for l in f:
	ll = l.split()
	x.append(float(ll[3]))
	y.append(float(ll[4]))
	z.append(float(ll[5]))
	lines.append(l)

sigma = 1.5
mx = numpy.mean(x)
my = numpy.mean(y)
mz = numpy.mean(z)
sx = numpy.std(x) * sigma
sy = numpy.std(y) * sigma
sz = numpy.std(z) * sigma

sys.stdout.write('filtered indices: ')
for i in range(len(x)):
	if abs(x[i]-mx)<sx and abs(y[i]-my)<sy and abs(z[i]-mz)<sz:
		g.write(lines[i])
	else:
		sys.stdout.write(str(i)+",")
sys.stdout.write('\n')
