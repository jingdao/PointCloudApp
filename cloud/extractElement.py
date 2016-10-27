#!/usr/bin/python
import sys
import numpy

if len(sys.argv) < 4:
	print 'Usage: extractElement bim.obj input.pcd outDir/'
	sys.exit(1)

class Cuboid:
	def __init__(self,corner):
		self.corner = corner
		self.length = [0]*3
		self.axes = [numpy.zeros(3) for i in range(3)]

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
	elif len(v) > 0:
		v = [numpy.array(i) for i in v]
		if len(v)==8:
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
#			print c.corner,c.length,c.axes
			elements.append(c)
		v=set()
obj.close()

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

for i in range(len(elements)):
	filename = sys.argv[3]+'/'+str(i)+'-cloud.pcd'
	f=open(filename,'w')
	c = elements[i]
	subset = []
	for p in points:
		q = p - c.corner
		q0 = q.dot(c.axes[0])
		q1 = q.dot(c.axes[1])
		q2 = q.dot(c.axes[2])
		if q0 >= 0 and q0 <= c.length[0] and \
			q1 >= 0 and q1 <= c.length[1] and \
			q2 >= 0 and q2 <= c.length[2]:
			subset.append(p)
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
