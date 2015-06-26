#include "kdtree.h"

KdTree::KdTree(PCD* pointcloud) {
	if (pointcloud->numPoints == 0)
		return;
	kdtreeLeaves = (KdTreeNode*) calloc(pointcloud->numPoints, sizeof(KdTreeNode));
	int* pointList = (int*) malloc(pointcloud->numPoints * sizeof(int));
	for (int i=0;i<pointcloud->numPoints;i++) {
		pointList[i] = i;
		kdtreeLeaves[i].index = i;
		kdtreeLeaves[i].isLeaf = true;
	}
	this->kdtreeDepth = 0;
	this->float_data = pointcloud->float_data;
	this->cloud = pointcloud;
	kdtreeRoot = buildKdTreeNode(pointList,pointcloud->numPoints,0);
	free(pointList);
}

KdTree::~KdTree() {
	if (kdtreeLeaves) free(kdtreeLeaves);
	for (size_t i=0;i<kdtreeBranches.size();i++)
		free(kdtreeBranches[i]);
}

KdTree::KdTreeNode* KdTree::buildKdTreeNode(int* pointList, int n, int depth) {
	if (depth > kdtreeDepth)
		kdtreeDepth = depth;
	if (n == 1) //only one point
		return kdtreeLeaves + *pointList; 
	int dimension = depth % 3;

	//find median element
	int i1 = rand() % n, i2 = rand() % n, i3 = rand() % n;
	int medianIndex;
	float f1 = float_data[pointList[i1] * 4 + dimension];
	float f2 = float_data[pointList[i2] * 4 + dimension];
	float f3 = float_data[pointList[i3] * 4 + dimension];
	if (f1 >= f2 && f1 <= f3) medianIndex = i1;
	else if (f2 >= f1 && f2 <= f3) medianIndex = i2;
	else medianIndex = i3;
//	int medianIndex = rand() % n;
	int median = pointList[medianIndex];
	float medianValue = float_data[median * 4 + dimension];
	KdTreeNode* currentNode = (KdTreeNode*) calloc(1,sizeof(KdTreeNode));
	currentNode->value = medianValue;
	kdtreeBranches.push_back(currentNode);

	//partition around median
	int leftSize, *leftPoints = pointList;
	int rightSize, *rightPoints;
	int tmp;
	pointList[medianIndex] = pointList[n-1];
	int i=-1,j=-1;
	do {
		j++;
		if (float_data[pointList[j] * 4 + dimension] <= medianValue) {
			i++;
			tmp = pointList[i];
			pointList[i] = pointList[j];
			pointList[j] = tmp;
		}
	} while (j < n - 2);
	pointList[n-1] = pointList[i+1];
	pointList[i+1] = median;
	if (i + 1 == n - 1) //median is last element
		leftSize = i + 1;
	else 
		leftSize = i + 2;
	rightSize = n - leftSize;
	rightPoints = pointList + leftSize;
//	printf("n %d i,j %d %d left %d right %d median %d\n",n,i,j,leftSize,rightSize,medianIndex);

	//get left and right children
	if (leftSize > 0)
		currentNode->left = buildKdTreeNode(leftPoints,leftSize,depth+1);
	if (rightSize > 0)
		currentNode->right = buildKdTreeNode(rightPoints,rightSize,depth+1);

//	printf("depth %d left %d right %d median %d %f\n",depth,leftSize,rightSize,medianIndex,medianValue);
	return currentNode;
}

void KdTree::search(std::vector<int> *indices,float x,float y,float z,float distance) {
	Cube range = {x-distance,y-distance,z-distance,x+distance,y+distance,z+distance};
	Cube boundingBox = getBoundingBox();
	searchNode(kdtreeRoot,indices,&range,boundingBox,0);
}

void KdTree::searchNode(KdTreeNode* node,std::vector<int> *indices,Cube* range,Cube cell,int depth){
	if (!node) return; //nothing here
	if (node->isLeaf) {
		float x = float_data[node->index * 4];
		float y = float_data[node->index * 4 + 1];
		float z = float_data[node->index * 4 + 2];
		if (x >= range->x1 && x <= range->x2 && 
			y >= range->y1 && y <= range->y2 && 
			z >= range->z1 && z <= range->z2)
			indices->push_back(node->index);
	} else {
		if (cell.x1 >= range->x1 && cell.x2 <= range->x2 &&
			cell.y1 >= range->y1 && cell.y2 <= range->y2 &&
			cell.z1 >= range->z1 && cell.z2 <= range->z2) {
			//region fully contained in range
			getIndexFromSubTree(node->left,indices);
			getIndexFromSubTree(node->right,indices);
		} else if ((cell.x2 >= range->x1 && cell.x1 <= range->x2) &&
			(cell.y2 >= range->y1 && cell.y1 <= range->y2) &&
			(cell.z2 >= range->z1 && cell.z1 <= range->z2)) {
			//region intersects range
			int dimension = depth % 3;
			Cube leftCell = cell;
			Cube rightCell = cell;
			//trim cells according to splitting plane
			switch (dimension) {
				case 0:
					leftCell.x2 = node->value;
					rightCell.x1 = node->value;
					break;
				case 1:
					leftCell.y2 = node->value;
					rightCell.y1 = node->value;
					break;
				case 2:
					leftCell.z2 = node->value;
					rightCell.z1 = node->value;
					break;
			}
			searchNode(node->left,indices,range,leftCell,depth+1);
			searchNode(node->right,indices,range,rightCell,depth+1);
		}
//		else {
//			printf("rejected\n");
//			printf("cell: %f %f %f %f %f %f\n",cell.x1,cell.y1,cell.z1,cell.x2,cell.y2,cell.z2);
//			printf("range: %f %f %f %f %f %f\n",range->x1,range->y1,range->z1,range->x2,range->y2,range->z2);
//		}
	}
}

void KdTree::getIndexFromSubTree(KdTreeNode* node, std::vector<int> *indices) {
	if (!node) return;
	if (node->isLeaf)
		indices->push_back(node->index);
	else {
		getIndexFromSubTree(node->left,indices);
		getIndexFromSubTree(node->right,indices);
	}
}

KdTree::Cube KdTree::getBoundingBox() {
	float minX=float_data[0],maxX=float_data[0];
	float minY=float_data[1],maxY=float_data[1];
	float minZ=float_data[2],maxZ=float_data[2];
	for (int i=1;i<cloud->numPoints;i++) {
		if (float_data[i * 4] < minX) minX = float_data[i * 4];
		else if (float_data[i * 4] > maxX) maxX = float_data[i * 4];
		if (float_data[i * 4 + 1] < minY) minY = float_data[i * 4 + 1];
		else if (float_data[i * 4 + 1] > maxY) maxY = float_data[i * 4 + 1];
		if (float_data[i * 4 + 2] < minZ) minZ = float_data[i * 4 + 2];
		else if (float_data[i * 4 + 2] > maxZ) maxZ = float_data[i * 4 + 2];
	}
	Cube res = {minX, minY, minZ, maxX, maxY, maxZ};
	return res;
}
