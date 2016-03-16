#!/usr/bin/python

import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def findSingleBlock(arr):
	countBlocks=0
	i=0
	while i<len(arr):
		if countBlocks==0 and arr[i]>0:
			countBlocks += 1
		elif countBlocks > 0:
			if arr[i] > 0:
				countBlocks += 1
			else:
				break
		i += 1
	while i<len(arr):
		if arr[i] > 0:
			return 0
		i += 1
	return countBlocks

def findDoubleBlock(arr):
	countBlock1=0
	countBlock2=0
	countGap=0
	i=0
	while i<len(arr):
		if countBlock1==0 and arr[i]>0:
			countBlock1 += 1
		elif countBlock1 > 0:
			if arr[i] > 0:
				countBlock1 += 1
			else:
				break
		i += 1
	while i<len(arr):
		if countBlock2==0:
			if arr[i] > 0:
				countBlock2 += 1
			else:
				countGap += 1
		elif countBlock2 > 0:
			if arr[i] > 0:
				countBlock2 += 1
			else:
				break
		i += 1
	while i<len(arr):
		if arr[i] > 0:
			return 0
		i += 1
	if countBlock1 > 0 and countBlock2 > 0:
		return countGap
	else:
		return 0

def find2DBlock(arr):
	i=0
	blocks=[]
	while i<len(arr):
		j=findSingleBlock(arr[i])
		if len(blocks)==0 and j>0:
			blocks.append(j)
		elif len(blocks) > 0:
			if j > 0:
				blocks.append(j)
			else:
				break
		i += 1
	while i<len(arr):
		j=findSingleBlock(arr[i])
		if j > 0:
			return 0
		i += 1
	if len(blocks)==0:
		return 0
	else:
		return max(blocks)

def detectBoom(grid,length,width,height):
	numBoomElement=0
	maxBoomElement=0
	x=0
	boom_threshold = min(width,height) / 2
	while x<length:
		j = find2DBlock(grid[x])
		if numBoomElement==0 and j>0 and j<boom_threshold:
			numBoomElement += 1
		elif numBoomElement > 0:
			if j>0 and j<boom_threshold:
				numBoomElement += 1
				maxBoomElement = max(maxBoomElement,numBoomElement)
			else:
				numBoomElement = 0
		x += 1
	return maxBoomElement

def sideToMid(grid,length,width,height):
	leftCount=0
	midCount=0
	rightCount=0
	leftCount2=0
	midCount2=0
	rightCount2=0
	for x in range(length/2):
		for z in range(height):
			for y in range(width):
				if grid[x][z][y] > 0:
					if y < width/3:
						leftCount += 1
					elif y < 2*width/3:
						midCount += 1
					else:
						rightCount += 1
	for x in range(length/2,length):
		for z in range(height):
			for y in range(width):
				if grid[x][z][y] > 0:
					if y < width/3:
						leftCount2 += 1
					elif y < 2*width/3:
						midCount2 += 1
					else:
						rightCount2 += 1
	m1 = 1.0 * (leftCount + rightCount) / midCount
	m2 = 1.0 * (leftCount2 + rightCount2) / midCount2
	return min(m1,m2) / max(m1,m2)

def frontToBack(grid,length,width,height):
	frontCount=0
	backCount=0
	for x in range(length/2):
		for z in range(height):
			for y in range(width):
				if grid[x][z][y] > 0:
					backCount += 1
	for x in range(length/2,length):
		for z in range(height):
			for y in range(width):
				if grid[x][z][y] > 0:
					frontCount += 1
	return 1.0 * min(frontCount,backCount) /max(frontCount,backCount) 

def detectBucket(grid,length,width,height):
	numBucketElement=0
	for x in range(length):
		for z in range(height-1):
			if findSingleBlock(grid[x][z+1])>3 and findDoubleBlock(grid[x][z])>3 and findDoubleBlock(grid[x][z])<5:
				numBucketElement += 1
				break
	return numBucketElement > 0

def topToBottom(grid,length,width,height):
	topCount1=0
	bottomCount1=0
	topCount2=0
	bottomCount2=0
	topCount3=0
	bottomCount3=0
	for x in range(length):
		for z in range(height):
			for y in range(width):
				if grid[x][z][y] > 0:
					if z < height/2:
						if x < length/3:
							bottomCount1+=1
						elif x < 2*length/3:
							bottomCount2+=1
						else:
							bottomCount3+=1
					else:
						if x < length/3:
							topCount1+=1
						elif x < 2*length/3:
							topCount2+=1
						else:
							topCount3+=1
	t1 = topCount1 / (bottomCount1+1) * (bottomCount2+bottomCount3) / (topCount2 + topCount3)
	t2 = topCount3 / (bottomCount3+1) * (bottomCount1+bottomCount2) / (topCount1 + topCount2)
#	return 1.0 * max(t1,t2) / (min(t1,t2) + 1)
	return min(t1,t2) / (max(t1,t2) + 1)

def midToFrontBack(grid,length,width,height):
	frontCount=0
	midCount=0
	backCount=0
	for x in range(length/3):
		for z in range(height):
			for y in range(width):
				if grid[x][z][y] > 0:
					backCount += 1
	for x in range(length/3,length*2/3):
		for z in range(height):
			for y in range(width):
				if grid[x][z][y] > 0:
					midCount += 1
	for x in range(length*2/3,length):
		for z in range(height):
			for y in range(width):
				if grid[x][z][y] > 0:
					frontCount += 1
	return 1.0 * midCount / (frontCount + backCount)
	
def gradToPeak(grid,length,width,height):
	maxheight=[0]*length
	for x in range(length):
		for z in range(height):
			for y in range(width):
				if grid[x][z][y] > 0:
					maxheight[x] = max(maxheight[x],z)
	peak_id = np.argmax(maxheight)
	g=[0,0]
	if peak_id > 0:
		g[0]=np.polyfit(range(peak_id+1),maxheight[:peak_id+1],1)[0]
	if peak_id < len(maxheight) - 1:
		g[1]=np.polyfit(range(len(maxheight),peak_id,-1),maxheight[peak_id:],1)[0]
	return min(g[0],g[1])/max(g[0],g[1])

def distToPeak(grid,length,width,height):
	maxheight=[0]*length
	for x in range(length):
		for z in range(height):
			for y in range(width):
				if grid[x][z][y] > 0:
					maxheight[x] = max(maxheight[x],z)
	peak_id = np.argmax(maxheight)
	p1 = 1.0 * peak_id / length
	p2 = 1.0 * (length-peak_id) / length
	return min(p1,p2)/max(p1,p2)

def detectBackhoe(grid,length,width,height):
	countTip=0
	for x in range(length):
		for y in range(width):
			if grid[x][height-1][y] > 0:
				countTip += 1
	return countTip

def detectBlade(grid,length,width,height):
	numBladeElement=0
	for x in range(length):
		countVertical=0
		for z in range(height):
			if findSingleBlock(grid[x][z])>2 and findSingleBlock(grid[x][z])<5:
				countVertical += 1
		if countVertical>2 and countVertical<5:
			numBladeElement += 1
	return numBladeElement >= 2

def getFeatures(grid,length,width,height):
	criteria0 = 1.0 * length / min(width,height)
	criteria1 = 1.0 * height / width
	criteria2 = frontToBack(grid,length,width,height)
	criteria3 = midToFrontBack(grid,length,width,height)
	criteria4 = sideToMid(grid,length,width,height)
	criteria5 = topToBottom(grid,length,width,height)
	criteria6 = gradToPeak(grid,length,width,height)
	criteria7 = distToPeak(grid,length,width,height)

#	criteria8 = detectBoom(grid,length,width,height)
	return [
		criteria0,
		criteria1,
		criteria2,
		criteria3,
		criteria4,
		criteria5,
		criteria6,
		criteria7,
		criteria8,
#		criteria9,
	]

def plotColoredRegion(x,y,colors,categories):
	stepSize = 0.05
	for c in x.keys():
		for n in range(len(x[c])):
			x[c][n] += np.random.randn() * stepSize 
			y[c][n] += np.random.randn() * stepSize
	xv = [i for c in x for i in x[c]]
	yv = [i for c in y for i in y[c]]
	xmin=np.min(xv) - 1
	xmax=np.max(xv) + 1
	ymin=np.min(yv) - 1
	ymax=np.max(yv) + 1
	fillRegion=False
	if fillRegion:
		colormap={}
		for c in colors.keys():
			hexvalue = int(colors[c][1:],16)
			red = (hexvalue>>16)&255
			green = (hexvalue>>8)&255
			blue = hexvalue&255
			scale=0.5 / 255
			colormap[c] = [scale*red,scale*green,scale*blue]
		xdim = int((xmax-xmin) / stepSize)
		ydim = int((ymax-ymin) / stepSize)
		palette = np.zeros((ydim,xdim,3))
		fixed = np.zeros((ydim,xdim),dtype=bool)
		for c in x.keys():
			for n in range(len(x[c])):
				xi = int((x[c][n]-xmin)/stepSize)
				yi = int((y[c][n]-ymin)/stepSize)
				palette[yi][xi] = colormap[c]
				fixed[yi][xi] = True
		updated=True
		while updated:
			updated=False
			visited = np.zeros((ydim,xdim),dtype=bool)
			for i in range(1,ydim-1):
				for j in range(1,xdim-1):
					if fixed[i][j]:
						if not fixed[i-1][j]:
							palette[i-1][j] = palette[i][j]
							visited[i-1][j] = True
							updated=True
						if not fixed[i+1][j]:
							palette[i+1][j] = palette[i][j]
							visited[i+1][j] = True
							updated=True
						if not fixed[i][j-1]:
							palette[i][j-1] = palette[i][j]
							visited[i][j-1] = True
							updated=True
						if not fixed[i][j+1]:
							palette[i][j+1] = palette[i][j]
							visited[i][j+1] = True
							updated=True
			fixed += visited
		palette[0][0] = palette[0][1]
		palette[0][xdim-1] = palette[1][xdim-1]
		palette[ydim-1][0] = palette[ydim-1][1]
		palette[ydim-1][xdim-1] = palette[ydim-1][xdim-2]
		plt.figure()
		plt.imshow(palette,interpolation='nearest')
		for key in x:
			plt.plot((x[key]-xmin)/stepSize,(y[key]-ymin)/stepSize,'o',markersize=10,label=categories[key],color=colors[key])
	else:
		plt.figure()
		for key in x:
			plt.plot(x[key],y[key],'.',markersize=10,label=categories[key],color=colors[key])
		plt.axis([xmin,xmax,ymin,ymax])
	plt.legend(loc='upper left')
	matplotlib.rcParams.update({'font.size': 7})
	plt.show()

if __name__ == "__main__":
	if len(sys.argv) < 2:
		print './part_feature.py -v [labels.txt ../labelCategory.txt *-cloud.pcd.og]'
		print './part_feature.py 0.og prediction.txt'
		print './part_feature.py 0.og 0-esf.pcd'
		sys.exit(1)

	debug=sys.argv[1]=='-v'
	if debug:
		labels = open('labels.txt','r')
		categories=[]
		categoryFile = open('../labelCategory.txt')
		for l in categoryFile:
			categories.append(l)
		categoryFile.close()
		c1={}
		c2={}
		i=0
		while True:
			try:
				f=open(str(i)+"-cloud.pcd.og",'r')
			except IOError:
				break
			l=f.readline().split()
			length = int(l[0])
			width = int(l[1])
			height = int(l[2])
			l=f.readline().split()
			grid=[[[float(l[x+y*length+z*length*width]) for y in range(width)] for z in range(height)] for x in range(length)]
			f.close()
			lb = int(labels.readline())
#			sys.stdout.write(str(lb)+' ')
			criteria = getFeatures(grid,length,width,height)
			if not lb in c1.keys():
				c1[lb]=[]
				c2[lb]=[]
			c1[lb].append(criteria[5])
			c2[lb].append(criteria[6])
#			if not lb==4 and not lb==3 and not lb==2:
#			print i,lb,criteria
			i += 1
		colors = {0:'#ffffff',1:'#ff0000',2:'#00ff00',3:'#0000ff',4:'#ffff00',5:'#ff00ff',6:'#00ffff'}
		plotColoredRegion(c1,c2,colors,categories)
	else:
		f=open(sys.argv[1],'r')
		l=f.readline().split()
		length = int(l[0])
		width = int(l[1])
		height = int(l[2])

		l=f.readline().split()
		grid=[[[float(l[x+y*length+z*length*width]) for y in range(width)] for z in range(height)] for x in range(length)]
		f.close()
		criteria = getFeatures(grid,length,width,height)
		if sys.argv[2].endswith('-esf.pcd'):
			esf = open(sys.argv[2],'w')
			esf.write("# .PCD v0.7 - Point Cloud Data file format\n")
			esf.write("VERSION 0.7\n")
			esf.write("FIELDS esf\n")
			esf.write("SIZE 4\n")
			esf.write("TYPE F\n")
			esf.write("COUNT "+str(len(criteria))+"\n")
			esf.write("WIDTH 1\n")
			esf.write("HEIGHT 1\n")
			esf.write("VIEWPOINT 0 0 0 1 0 0 0\n")
			esf.write("POINTS 1\n")
			esf.write("DATA ascii\n");
			for c in criteria:
				esf.write(str(c)+' ')
			esf.write("\n")
			esf.close()
		else:
			prediction = open(sys.argv[2],'a')
			if criteria[2] + criteria[3] + criteria[4] >= 2 and criteria[0] > 2.2:
				prediction.write('4\n')
			elif criteria[0] > 3.5:
				prediction.write('3\n')
			elif criteria[0] < 2.2:
				prediction.write('2\n')
			elif criteria[2] and criteria[6] < 5:
				prediction.write('1\n')
			else:
				prediction.write('5\n')
			prediction.close()
	
