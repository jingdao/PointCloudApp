#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include "pcd.h"

class PCD;

class KdTree {
	public:	
		struct KdTreeNode {
			bool isLeaf;
			union {
				int index;
				float value;
			};
			KdTreeNode* left;
			KdTreeNode* right;
		};
		struct Cube {float x1,y1,z1,x2,y2,z2;};

		KdTreeNode* kdtreeLeaves;
		std::vector<KdTreeNode*> kdtreeBranches;
		KdTreeNode* kdtreeRoot;
		int kdtreeDepth;
		float* float_data;
		PCD* cloud;

		KdTree(PCD* pointcloud);
		~KdTree();
		KdTreeNode* buildKdTreeNode(int* pointList, int n, int depth);
		void search(std::vector<int> *indices,float x,float y,float z,float distance);
		void searchNode(KdTreeNode* node,std::vector<int> *indices,Cube* range,Cube cell,int depth);
		void getIndexFromSubTree(KdTreeNode* node, std::vector<int> *indices);
		Cube getBoundingBox();
};
