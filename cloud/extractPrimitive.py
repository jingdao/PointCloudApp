#!/usr/bin/python
import sys
import numpy
import matplotlib.pyplot as plt

def findDeviations(A,limit):
	l = 10
	r = 1.5
	window = []
	for i in range(len(A)):
		window.append(A[max(0,i-l):min(len(A),i+l+1)])
	median = [sorted(w)[l] for w in window]
	peaks = []
	sequence = []
	for i in range(len(A)):
		if A[i] > median[i] * r:
			sequence.append(i)
		elif len(sequence) > 0:
			peaks.append(sequence)
			sequence = []
	if len(sequence) > 0:
		peaks.append(sequence)
	peaks = sorted(peaks,key = lambda x:max(A[x]),reverse=True)
	return peaks[:limit]

def cluster(points):
	S = set(points)
	centerX = []
	centerY = []
	for seed in points:
		if seed not in S:
			continue
		Q=[seed]
		C=[]
		while len(Q) > 0:
			q = Q.pop()
			if q not in S:
				continue
			C.append(q)
			S.remove(q)
			Q.append((q[0]-1,q[1]))
			Q.append((q[0],q[1]-1))
			Q.append((q[0]+1,q[1]))
			Q.append((q[0],q[1]+1))
		cx,cy = numpy.mean(C,0)
		centerX.append(cx)
		centerY.append(cy)
	return centerX,centerY

def PCD_getBounds(points,robust=False):
	if not robust:
		C = numpy.array(points)
		minX = numpy.min(C[:,0])
		maxX = numpy.max(C[:,0])
		minY = numpy.min(C[:,1])
		maxY = numpy.max(C[:,1])
		minZ = numpy.min(C[:,2])
		maxZ = numpy.max(C[:,2])
		return minX,minY,minZ,maxX,maxY,maxZ
	else:
		x = sorted([p[0] for p in points])
		y = sorted([p[1] for p in points])
		z = sorted([p[2] for p in points])
		l = len(points)
		low = l/4
		high = -l/4
		return x[low],y[low],min(z),x[high],y[high],max(z)

def PCD_getHash(points,resolution):
	x = [p[0] for p in points]
	y = [p[1] for p in points]
	z = [p[2] for p in points]
	xd = 0.5 * (min(x) + max(x))
	yd = 0.5 * (min(y) + max(y))
	zd = 0.5 * (min(z) + max(z))
	pointset = set()
	for i in range(len(x)):
		p = (int((x[i]-xd)/resolution),int((y[i]-yd)/resolution),int((z[i]-zd)/resolution))
		if not p in pointset:
			pointset.add(p)
	return pointset,(xd,yd,zd)

def PCD_xycluster(phash,box,points,resolution):
	clusters = []
	clusterID = {}
	S = set()
	for p in phash:
		S.add((p[0],p[1]))
	c = 0
	while len(S) > 0:
		seed = list(S)[0]
		Q = [seed]
		while len(Q) > 0:
			q = Q.pop()
			if q not in S:
				continue
			S.remove(q)
			clusterID[q] = c
			Q.append((q[0]-1,q[1]))
			Q.append((q[0],q[1]-1))
			Q.append((q[0]+1,q[1]))
			Q.append((q[0],q[1]+1))
		clusters.append([])
		c += 1
	for p in points:
		xi = int((p[0] - box[0]) / resolution)
		yi = int((p[1] - box[1]) / resolution)
		clusters[clusterID[(xi,yi)]].append(p)
	clusters = sorted(clusters,key=lambda x:len(x),reverse=True)
	return clusters

def getLines(points):
	minX,minY,minZ = numpy.min(points,0)
	maxX,maxY,maxZ = numpy.max(points,0)
	grid = 150
	threshold = 0.5
	dx = (maxX - minX) / grid
	dy = (maxY - minY) / grid
	resolution = min(dx,dy)
	gridX = (maxX - minX) / resolution
	gridY = (maxY - minY) / resolution
	acc = numpy.zeros((gridX+1,gridY+1))
	for p in points:
		idx = int((p[0]-minX) / resolution)
		idy = int((p[1]-minY) / resolution)
		acc[idx,idy] += 1
	mask = acc > 0
	sumX = numpy.sum(mask,0)
	sumY = numpy.sum(mask,1)
	dix = list(numpy.nonzero(sumX < max(sumX)*threshold)[0])
	diy = list(numpy.nonzero(sumY < max(sumY)*threshold)[0])
	for i in diy:
		for j in dix:
			acc[i,j] = 0
	dix = list(numpy.nonzero(sumX >= max(sumX)*threshold)[0])
	diy = list(numpy.nonzero(sumY >= max(sumY)*threshold)[0])
	plt.imshow(mask.transpose(),cmap='Greys')
	plt.figure()
	plt.imshow(acc.transpose() > 0,cmap='Greys')
	plt.show()
	nh = 1
	nv = 1
	mapx = {diy[0]:0}
	mapy = {dix[0]:0}
	for i in range(1,len(diy)):
		if diy[i] > diy[i-1] + 1:
			nv += 1
		mapx[diy[i]] = nv-1
	for i in range(1,len(dix)):
		if dix[i] > dix[i-1] + 1:
			nh += 1
		mapy[dix[i]] = nh-1
	horizontal = [[] for i in range(nh)]
	vertical = [[] for i in range(nv)]
	for p in points:
		idx = int((p[0]-minX) / resolution)
		idy = int((p[1]-minY) / resolution)
		if idx in mapx:
			vertical[mapx[idx]].append(p)
		if idy in mapy:
			horizontal[mapy[idy]].append(p)
	return horizontal,vertical

def loadPCD(filename):
	pcd = open(filename,'r')
	for l in pcd:
		if l.startswith('DATA'):
			break
	points = []
	for l in pcd:
		ll = l.split()
		x = float(ll[0])
		y = float(ll[1])
		z = float(ll[2])
		points.append([x,y,z])
	pcd.close()
	points = numpy.array(points)
	return points

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

def saveOBJ(filename,dimensions):
	f = open(filename,"w")
	f.write("v %f %f %f\n" % (dimensions[0],dimensions[1],dimensions[2]))
	f.write("v %f %f %f\n" % (dimensions[3],dimensions[1],dimensions[2]))
	f.write("v %f %f %f\n" % (dimensions[0],dimensions[4],dimensions[2]))
	f.write("v %f %f %f\n" % (dimensions[3],dimensions[4],dimensions[2]))
	f.write("v %f %f %f\n" % (dimensions[0],dimensions[1],dimensions[5]))
	f.write("v %f %f %f\n" % (dimensions[3],dimensions[1],dimensions[5]))
	f.write("v %f %f %f\n" % (dimensions[0],dimensions[4],dimensions[5]))
	f.write("v %f %f %f\n" % (dimensions[3],dimensions[4],dimensions[5]))
	f.write("f 1 2 4 3\n")
	f.write("f 1 2 6 5\n")
	f.write("f 1 3 7 5\n")
	f.write("f 2 4 8 6\n")
	f.write("f 3 4 8 7\n")
	f.write("f 5 6 8 7\n")
	f.close()

if len(sys.argv) < 3:
	print 'Usage: extractPrimitive.py input.pcd outDir/'
	sys.exit(1)

if 'level' in sys.argv[1]:
	target = loadPCD(sys.argv[1])
	h,v = getLines(target)
	b = 0
	for i in range(len(h)):
		bounds = PCD_getBounds(h[i])
		savePCD("%s/%d-cloud.pcd" % (sys.argv[2],b),h[i])
		saveOBJ("%s/beam%d.obj" % (sys.argv[2],b), bounds)
		b += 1
	for i in range(len(v)):
		bounds = PCD_getBounds(v[i])
		savePCD("%s/%d-cloud.pcd" % (sys.argv[2],b),v[i])
		saveOBJ("%s/beam%d.obj" % (sys.argv[2],b), bounds)
		b += 1
	print 'Found %d horizontal %d vertical' % (len(h),len(v))
elif 'storey' in sys.argv[1]:
	target = loadPCD(sys.argv[1])
	r = 0.2
	phash,box = PCD_getHash(target,r)
	clusters = PCD_xycluster(phash,box,target,r)
	b = 0
	targetBounds = PCD_getBounds(target)
	for i in range(min(len(clusters),20)):
		savePCD("%s/%d-cloud.pcd" % (sys.argv[2],i), clusters[i])
		bounds = PCD_getBounds(clusters[i],True)
		if bounds[3]-bounds[0] < 3 and bounds[4]-bounds[1] < 3 and bounds[5]-bounds[2] > targetBounds[5]-targetBounds[2]-0.5:
			saveOBJ("%s/column%d.obj" % (sys.argv[2],b), bounds)
			b += 1
	print 'Found %d columns' % b
else:
	points = loadPCD(sys.argv[1])
	minX,minY,minZ = numpy.min(points,0)
	maxX,maxY,maxZ = numpy.max(points,0)
	print 'Loaded %s (%f,%f,%f,%f,%f,%f)' % (sys.argv[1],minX,minY,minZ,maxX,maxY,maxZ)
	horiz_grid = 100
	vert_grid = 100
	dx = (maxX - minX) / horiz_grid
	dy = (maxY - minY) / horiz_grid
	dz = (maxZ - minZ) / vert_grid
	acc = numpy.zeros((horiz_grid+1,horiz_grid+1,vert_grid+1))
	for p in points:
		idx = int((p[0]-minX) / dx)
		idy = int((p[1]-minY) / dy)
		idz = int((p[2]-minZ) / dz)
		acc[idx,idy,idz] = 1

	levels = numpy.sum(acc,(0,1))
	peaks = findDeviations(levels,6)
	print 'Found levels',peaks
#	peakVal = [i for j in peaks for i in j]
#	plt.plot(levels,'.',mew=1,markersize=10)
#	plt.plot(peakVal,levels[peakVal],'x',mew=2,markersize=10)
#	plt.xlabel('Z-coordinate bins')
#	plt.ylabel('Point density')
#	plt.show()

	peaks = sorted(peaks,key=lambda x:x[0])
	levels = [[] for l in range(len(peaks))]
	storeys = [[] for l in range(len(peaks)-1)]
	for p in points:
		idz = int((p[2]-minZ) / dz)
		for l in range(len(peaks)):
			if idz >= peaks[l][0] and idz <= peaks[l][-1]:
				levels[l].append((p[0],p[1],p[2]))
				break
			elif l < len(peaks)-1 and idz < peaks[l+1][0] and idz > peaks[l][-1]:
				storeys[l].append((p[0],p[1],p[2]))
				break
	for l in range(len(levels)):
		savePCD("%s/level%d.pcd"%(sys.argv[2],l),levels[l])
	for l in range(len(storeys)):
		savePCD("%s/storey%d.pcd"%(sys.argv[2],l),storeys[l])

