#include "kdtree.h"

KdTree::KdTree(PCD* pointcloud) {
	if (!pointcloud || pointcloud->numPoints == 0)
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
	findMedian(pointList,n,n/2,dimension); //partitions around median
	int medianIndex = n / 2;
	int median = pointList[medianIndex];
	float medianValue = float_data[median * 4 + dimension];
	KdTreeNode* currentNode = (KdTreeNode*) calloc(1,sizeof(KdTreeNode));
	currentNode->value = medianValue;
	kdtreeBranches.push_back(currentNode);

	//partition around median
	int leftSize = medianIndex, *leftPoints = pointList;
	int rightSize = n - leftSize, *rightPoints = pointList + leftSize;
//	printf("n %d i,j %d %d left %d right %d median %d\n",n,i,j,leftSize,rightSize,medianIndex);

	//get left and right children
	if (leftSize > 0)
		currentNode->left = buildKdTreeNode(leftPoints,leftSize,depth+1);
	if (rightSize > 0)
		currentNode->right = buildKdTreeNode(rightPoints,rightSize,depth+1);

//	printf("depth %d left %d right %d median %d %f\n",depth,leftSize,rightSize,medianIndex,medianValue);
	return currentNode;
}

void KdTree::findMedian(int* pointList, int n, int target, int dimension) {
	if (n < 2) return;
	int pivot_index = rand() % n;
	int pivot = pointList[pivot_index];
	float pivotValue = float_data[pivot * 4 + dimension];
	pointList[pivot_index] = pointList[n-1];
	int i=-1,j=-1,tmp;
	do {
		j++;
		if (float_data[pointList[j] * 4 + dimension] < pivotValue || 
			(float_data[pointList[j] * 4 + dimension] == pivotValue && rand()%2 == 0)) {
			i++;
			tmp = pointList[i];
			pointList[i] = pointList[j];
			pointList[j] = tmp;
		}
	} while (j < n - 2);
	i++;
	pointList[n-1] = pointList[i];
	pointList[i] = pivot;
	if (i == target) return;
	else if (i > target) findMedian(pointList, i , target, dimension);
	else findMedian(pointList + i + 1, n - i - 1, target - i - 1, dimension);
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
