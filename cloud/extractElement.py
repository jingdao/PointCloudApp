#!/usr/bin/python
import sys
import numpy
import matplotlib.pyplot as plt

if len(sys.argv) < 4:
	print 'Usage: extractElement bim.obj input.pcd outDir/'
	sys.exit(1)

class Cuboid:
	def __init__(self,corner):
		self.corner = corner
		self.length = [0]*3
		self.axes = [numpy.zeros(3) for i in range(3)]

def convexHull(P):
	n = len(P)
	k = 0
	H = [None] * (2*n)
	if n==0:
		return []
	P = sorted(P,key=lambda x:x[1])
	P = sorted(P,key=lambda x:x[0])

	for i in range(n):
		while k >= 2 and numpy.cross(H[k-1]-H[k-2], P[i]-H[k-2]) <= 0:
			k -= 1
		H[k] = P[i]
		k += 1

	t = k+1
	for i in range(n-2,-1,-1):
		while k >= t and numpy.cross(H[k-1]-H[k-2], P[i]-H[k-2]) <= 0:
			k -= 1
		H[k] = P[i]
		k += 1

	return H[:k-1]

def getCuboidFromCorner(v):
	c = Cuboid(v[0])
	norm = [0]*8
	for i in range(1,8):
		v[i] -= c.corner
		norm[i]  = numpy.sqrt(v[i].dot(v[i]))
	norm_sorted = sorted([(i,norm[i]) for i in range(1,8)], key = lambda x:x[1])
	threshold = numpy.sqrt(norm_sorted[0][1]*norm_sorted[0][1]+norm_sorted[1][1]*norm_sorted[1][1])
	if abs(norm_sorted[2][1] - threshold) > abs(norm_sorted[3][1] - threshold):
		selector = [0,1,2]
	else:
		selector = [0,1,3]
	for j in range(3):
		i = norm_sorted[selector[j]][0]
		c.length[j] = norm[i]
		c.axes[j] = v[i] / norm[i]
	return c

def getCuboidFromPoints(v):
	z = v[:,2]
	zmin = min(z)
	zmax = max(z)
	v = v[(z==zmin) | (z==zmax)]
	hull = numpy.array(convexHull(v[:,:2]))
	edgeAngles = []
	for i in range(len(hull) - 1):
		theta = numpy.arctan2(hull[i+1][1] - hull[i][1] , hull[i+1][0] - hull[i][0])
		edgeAngles.append(theta)
	minArea = float('inf')
	for theta in edgeAngles:
		R = numpy.array([[numpy.cos(theta),numpy.sin(theta)],[-numpy.sin(theta),numpy.cos(theta)]])
		vr = R.dot(v[:,:2].transpose())
		xmin = min(vr[0])
		xmax = max(vr[0])
		ymin = min(vr[1])
		ymax = max(vr[1])
		A = (xmax-xmin) * (ymax-ymin)
		if A < minArea:
			minArea = A
			c = Cuboid(numpy.hstack((R.transpose().dot([xmin,ymin]),zmin)))
			c.length[0] = xmax-xmin
			c.length[1] = ymax-ymin
			c.length[2] = zmax-zmin
			c.axes[0] = numpy.hstack((R[:,0],0))
			c.axes[1] = numpy.hstack((R[:,1],0))
			c.axes[2] = numpy.array([0,0,1])
	return c

v=set()
elements=[]
obj = open(sys.argv[1],'r')
for l in obj:
	if l[0:2]=='v ':
		ll = l.split()
		x = float(ll[1])
		y = float(ll[2])
		z = float(ll[3])
		v.add((x,y,z))
	elif l.startswith('usemtl '):
#		if l.startswith('usemtl IfcWallStandardCase'):
		v = numpy.array([numpy.array(i) for i in v])
		if len(v) == 8:
			c = getCuboidFromCorner(v)
		else:
			c = getCuboidFromPoints(v)
#		print c.length,c.axes
		elements.append(c)
		v=set()
obj.close()
print 'Found %d elements' % len(elements)
if len(elements)==0:
	sys.exit(0)

pcd = open(sys.argv[2],'r')
for l in pcd:
	if l.startswith('DATA'):
		break
points = []
for l in pcd:
	ll = l.split()
	x = float(ll[0])
	y = float(ll[1])
	z = float(ll[2])
	points.append(numpy.array([x,y,z]))
pcd.close()

j=0
margin = 0.1
for i in range(len(elements)):
	c = elements[i]
	subset = []
#	for p in points:
#		q = p - c.corner
#		q0 = q.dot(c.axes[0])
#		q1 = q.dot(c.axes[1])
#		q2 = q.dot(c.axes[2])
#		if q0 >= -margin and q0 <= c.length[0]+margin and \
#			q1 >= -margin and q1 <= c.length[1]+margin and \
#			q2 >= -margin and q2 <= c.length[2]+margin:
#			subset.append(p)
	p = numpy.array(points)
	q = p - c.corner
	q0 = q.dot(c.axes[0])
	q1 = q.dot(c.axes[1])
	q2 = q.dot(c.axes[2])
	valid = (q0 >= -margin) & (q0 <= c.length[0]+margin) & \
			(q1 >= -margin) & (q1 <= c.length[1]+margin) & \
			(q2 >= -margin) & (q2 <= c.length[2]+margin)
	subset = p[valid]
	if len(subset) < 1:
		print c.axes
		continue
	filename = sys.argv[3]+'/'+str(j)+'-cloud.pcd'
	j += 1
	f=open(filename,'w')
	f.write("# .PCD v0.7 - Point Cloud Data file format\n")
	f.write("VERSION 0.7\n")
	f.write("FIELDS x y z\n")
	f.write("SIZE 4 4 4\n")
	f.write("TYPE F F F\n")
	f.write("COUNT 1 1 1\n")
	f.write("WIDTH "+str(len(subset))+"\n")
	f.write("HEIGHT 1\n")
	f.write("VIEWPOINT 0 0 0 1 0 0 0\n")
	f.write("POINTS "+str(len(subset))+"\n")
	f.write("DATA ascii\n")
	for p in subset:
		f.write("%f %f %f\n" % (p[0],p[1],p[2]))
	f.close()
	print 'Wrote '+str(len(subset))+' points to '+filename
