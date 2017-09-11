#!/usr/bin/python
import sys
import re
import os

if len(sys.argv) < 3:
	print 'Usage: obj2ply.py bim.obj [outDir/ out.ply]'
	sys.exit(1)

vertices=[]
faces=[]
name=""
material=""
elements=[]
offset=0
labels=None

if os.path.isdir(sys.argv[2]):
	index=0
	labels=open(sys.argv[2]+'/labels.txt','a')
	while True:
		try:
			f=open('%s/%d.ply'%(sys.argv[2],index),'r')
		except IOError:
			break
		index+=1
		f.close()
	filename = '%s/%d.ply' % (sys.argv[2],index)
else:
	filename = sys.argv[2]

def writePLY(filename):
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
	if labels is not None:
		labels.write("%s %s\n"%(name,material))
	print 'Wrote %d vertices %d faces to %s'%(len(vertices),len(faces),filename)

obj = open(sys.argv[1],'r')
for l in obj:
	if l[0:2]=='g ' and labels is not None:
		if len(vertices)>0:
			writePLY(filename)
			offset += len(vertices)
			vertices=[]
			faces=[]
		name = l[2:].strip()
	elif l[0:2]=='v ':
		vertices.append(l[2:])
	elif l.startswith('usemtl '):
		material=l[7:].strip()
	elif l[0:2]=='f ':
		if '/' in l:
			ll=re.split(' |//',l[2:])
			vid = [str(int(i)-1-offset) for i in ll[::2]]
		else:
			ll=l[2:].split()
			vid = [str(int(i)-1-offset) for i in ll]
		faces.append(str(len(vid))+' '+' '.join(vid)+'\n')
writePLY(filename)
if labels is not None:
	labels.close()

