#!/usr/bin/python

import numpy
import matplotlib.pyplot as plt

before=[]
after=[]

f = open('meanshift.txt','r')
for l in f:
	s = l.split()
	if s[0]=='before':
		before.append((float(s[1]),float(s[2])))
	elif s[0]=='after':
		after.append((float(s[1]),float(s[2])))
f.close()

before = numpy.array(before)
after = numpy.array(after)
plt.plot(before[:,0],before[:,1],'o',label='original points',ms=5)
plt.plot(after[:,0],after[:,1],'x',label='cluster centers',ms=20,mew=10)
plt.legend()
plt.xlabel('x')
plt.ylabel('y')
plt.show()
