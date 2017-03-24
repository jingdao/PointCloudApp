#!/usr/bin/python
import subprocess
import re
d={}
f = open('labels.txt','r')
i = 0
for l in f:
	ll = re.split(':| ',l)
	cls = ll[0]
	if not cls[0].isdigit():
		cls = cls.replace('(','_').replace(')','_')
		if cls not in d:
			d[cls] = 1
		else:
			d[cls] += 1
	i += 1
f.close()

dd = sorted(d.items(),key=lambda x:x[1],reverse=True)
for s in dd:
	print s

f = open('labels.txt','r')
i = 0
d2 = {}
for l in f:
	ll = re.split(':| ',l)
	cls = ll[0]
	cls = cls.replace('(','_').replace(')','_')
	if not cls[0].isdigit() and d[cls] > 2 and d[cls] < 40:
		if cls not in d2:
			print cls
			subprocess.call(['mkdir',cls])
			d2[cls] = 1
		subprocess.call(['cp','%d.pcd'%i,cls])
	i += 1
f.close()
