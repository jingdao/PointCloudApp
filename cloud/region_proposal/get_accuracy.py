#!/usr/bin/python
import sys
if len(sys.argv) < 3:
	print './get_accuracy.py source/ target/'
	sys.exit(1)

def loadPCDCentroid(filename):
	try:
		pcd = open(filename,'r')
	except IOError:
		return None
	for l in pcd:
		if l.startswith('DATA'):
			break
	n = 0
	mx = 0
	my = 0
	mz = 0
	for l in pcd:
		ll = l.split()
		x = float(ll[0])
		y = float(ll[1])
		z = float(ll[2])
		mx += x
		my += y
		mz += z
		n += 1
	pcd.close()
	return (mx/n, my/n, mz/n)

import numpy

source = []
i=0
while True:
	m = loadPCDCentroid("%s/%d-cloud.pcd" % (sys.argv[1],i))
	if not m:
		break
	source.append(numpy.array(m))
	i += 1

target = []
i=0
while True:
	m = loadPCDCentroid("%s/%d-cloud.pcd" % (sys.argv[2],i))
	if not m:
		break
	target.append(numpy.array(m))
	i += 1

threshold = 2
correct = 0
RMSE = 0
for t in target:
	for s in source:
		d = numpy.linalg.norm(s-t)
		if d < threshold:
			correct += 1
			RMSE += d*d
			break
if correct > 0:
	RMSE = numpy.sqrt(RMSE / correct)
	print 'prec: %.2f recl: %.2f RMSE: %.2f' % (1.0*correct/len(source),1.0*correct/len(target),RMSE)
else:
	print 'prec: %.2f recl: %.2f RMSE: %.2f' % (0,0,0)
