#!/usr/bin/python
import sys
import numpy
import matplotlib.pyplot as plt

def diff(a,b):
	return 1.0*(a-b)/b

def findpeaks(A,limit):
	peaks=[]
	if A[0] > A[1]:
		peaks.append((0,diff(A[0],A[1])))
	for i in range(1,len(A)-1):
		if A[i] > A[i-1] and A[i] > A[i+1]:
			peaks.append((i,max(diff(A[i],A[i-1]),diff(A[i],A[i+1]))))
	if A[-1] > A[-2]:
		peaks.append((len(A)-1,diff(A[-1],A[-2])))
	median = sorted(A)[(len(A)-1)/2]
	peaks = filter(lambda x:A[x[0]]>median,peaks)
	peaks = sorted(peaks,key=lambda x:x[1],reverse=True)
	return [x[0] for x in peaks[:limit]]

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

if len(sys.argv) < 3:
	print 'Usage: extractPrimitive.py input.pcd outDir/'
	sys.exit(1)

pcd = open(sys.argv[1],'r')
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
minX,minY,minZ = numpy.min(points,0)
maxX,maxY,maxZ = numpy.max(points,0)

print 'Loaded %s (%f,%f,%f,%f,%f,%f)' % (sys.argv[1],minX,minY,minZ,maxX,maxY,maxZ)

horiz_grid = 100
vert_grid = 30
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
peaks = findpeaks(levels,5)
print 'Found levels',peaks
plt.plot(levels,'.')
plt.show()

peaks = sorted(peaks)
for l in range(len(peaks)-1):
	im = numpy.sum(acc[:,:,peaks[l]+1:peaks[l+1]],2)
	low = numpy.min(im)
	high = numpy.max(im)
	mask = im > high - 5
	x,y = numpy.nonzero(mask)
	x,y = cluster(zip(x,y))
	plt.imshow(im)
	plt.scatter(y,x,marker='o',edgecolors='w',facecolors='none',s=50)
	plt.title('Storey %d (%d columns)' % (l,len(x)))
	plt.show()
	plt.imshow(acc[:,:,peaks[l+1]])
	plt.title('Level %d' % l)
	plt.show()

#levels = [[] for l in range(len(peaks)-1)]
#for p in points:
#	idz = int((p[2]-minZ) / dz)
#	for l in range(len(peaks)-1):
#		if idz == peaks[l+1]:
#			levels[l].append((p[0],p[1],p[2]))
#			break
#for l in range(len(levels)):
#	savePCD("%s/level%d.pcd"%(sys.argv[2],l),levels[l])

