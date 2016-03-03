#!/usr/bin/python

import sys

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
	while x<length:
		j = find2DBlock(grid[x])
		if numBoomElement==0 and j>0 and j<4:
			numBoomElement += 1
		elif numBoomElement > 0:
			if j>0 and j<4:
				numBoomElement += 1
				maxBoomElement = max(maxBoomElement,numBoomElement)
			else:
				numBoomElement = 0
		x += 1
	return maxBoomElement >= 3

def detectBoom2(grid,length,width,height):
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
	return midCount*(leftCount2+rightCount2) < 0.4 * midCount2*(leftCount+rightCount)


def detectBoom3(grid,length,width,height):
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
	return backCount > 1.2 * frontCount

def detectBucket(grid,length,width,height):
	numBucketElement=0
	for x in range(length):
		for z in range(height-1):
			if findSingleBlock(grid[x][z+1])>3 and findDoubleBlock(grid[x][z])>3 and findDoubleBlock(grid[x][z])<5:
				numBucketElement += 1
				break
	return numBucketElement > 0

def detectBucket2(grid,length,width,height):
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
	return topCount1*(bottomCount2+bottomCount3) < 0.5 * bottomCount1*(topCount2+topCount3) or \
			topCount3*(bottomCount1+bottomCount2) < 0.5 * bottomCount3*(topCount1+topCount2)

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

if __name__ == "__main__":
	if len(sys.argv) < 3:
		print './part_feature.py 0.og prediction.txt'
		sys.exit(1)

	debug=True
	if debug:
		labels = open('labels.txt','r')
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
			criteria1 = detectBoom(grid,length,width,height)
			criteria2 = detectBoom2(grid,length,width,height)
			criteria3 = detectBoom3(grid,length,width,height)
			criteria4 = 1.0 * length / min(width,height)
			criteria5 = detectBucket2(grid,length,width,height)
			criteria6 = detectBackhoe(grid,length,width,height)
#			if not lb==4 and not lb==3 and not lb==2:
			print i,lb,criteria1,criteria2,criteria3,criteria4,criteria5,criteria6
			i += 1
	else:
		f=open(sys.argv[1],'r')
		l=f.readline().split()
		length = int(l[0])
		width = int(l[1])
		height = int(l[2])

		l=f.readline().split()
		grid=[[[float(l[x+y*length+z*length*width]) for y in range(width)] for z in range(height)] for x in range(length)]
		f.close()
		prediction = open(sys.argv[2],'a')
		criteria1 = detectBoom(grid,length,width,height)
		criteria2 = detectBoom2(grid,length,width,height)
		criteria3 = detectBoom3(grid,length,width,height)
		criteria4 = 1.0 * length / min(width,height)
		criteria5 = detectBucket2(grid,length,width,height)
		criteria6 = detectBackhoe(grid,length,width,height)
		if criteria1 + criteria2 + criteria3 >= 2 and criteria4 > 2.2:
			prediction.write('4\n')
		elif criteria4 > 3.5:
			prediction.write('3\n')
		elif criteria4 < 2.2:
			prediction.write('2\n')
		elif criteria1 and criteria6 < 5:
			prediction.write('1\n')
		else:
			prediction.write('5\n')
		prediction.close()
	
