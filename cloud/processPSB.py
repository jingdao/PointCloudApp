#!/usr/bin/python

import os
import time

psb_dir='/home/jd/Documents/benchmark'
script_dir='/home/jd/Documents/PointCloudApp/cloud'

categoryFile = open('labelCategory.txt','r')
i=0
categories={}
for l in categoryFile:
	categories[l.split()[1]] = i
	i += 1
categoryFile.close()

models={}
classFile = open(psb_dir+'/classification/v1/coarse1/coarse1Train.cla')
while True:
	l = classFile.readline()
	if len(l) <= 0:
		break
	l = l.split()
	if not len(l) == 3:
		continue
	categoryName = l[0]
	if categoryName in categories:
		categorySize = int(l[2])
		models[categoryName] = []
		for i in range(categorySize):
			models[categoryName].append(int(classFile.readline()))
classFile.close()

labelFile = open('labels.txt','w')
index=0
for modelName in models:
	for m in models[modelName]:
		offFile = psb_dir+'/db/'+str(m/100)+'/m'+str(m)+'/m'+str(m)+'.off'
		plyFile = str(index)+'.ply'
		pcdFile = str(index)+'-cloud.pcd'
		os.system(script_dir+'/off2ply '+offFile+' '+plyFile)
		os.system('pcl_mesh2pcd -level 1 '+plyFile+' '+pcdFile+' > /dev/null &')
		time.sleep(5)
		os.system('xdotool key --delay 1000 q q')
		labelFile.write(str(categories[modelName])+'\n')
		index += 1
	print 'Processed category '+modelName+': '+str(len(models[modelName]))+' files'
labelFile.close()

		





