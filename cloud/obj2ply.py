#!/usr/bin/python
import sys
import re

if len(sys.argv) < 3:
	print 'Usage: obj2ply.py bim.obj outDir/'
	sys.exit(1)

vertices=[]
faces=[]
name=""
material=""
elements=[]
offset=0

index=0
labels=open(sys.argv[2]+'/labels.txt','a')
while True:
	try:
		f=open('%s/%d.ply'%(sys.argv[2],index),'r')
	except IOError:
		break
	index+=1
	f.close()

def writePLY():
	global index
	filename = '%s/%d.ply' % (sys.argv[2],index)
	index += 1
	f=open(filename,'w')
	header="""ply
format ascii 1.0
element vertex %d
property float x
property float y
property float z
element face %d
property list uchar int vertex_indices
end_header
""" % (len(vertices),len(faces))
	f.write(header)
	for v in vertices:
		f.write(v)
	for v in faces:
		f.write(v)
	f.close()
	labels.write("%s %s\n"%(name,material))
	print 'Wrote %d points to %s'%(len(vertices),filename)

obj = open(sys.argv[1],'r')
for l in obj:
	if l[0:2]=='g ':
		if len(vertices)>0:
			writePLY()
			offset += len(vertices)
			vertices=[]
			faces=[]
		name = l[2:].strip()
	elif l[0:2]=='v ':
		vertices.append(l[2:])
	elif l.startswith('usemtl '):
		material=l[7:].strip()
	elif l[0:2]=='f ':
		ll=re.split(' |//',l[2:])
		vid = [str(int(i)-1-offset) for i in ll[::2]]
		faces.append(str(len(vid))+' '+' '.join(vid)+'\n')
writePLY()
labels.close()

