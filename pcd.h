#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dirent.h>
#include <vector>
#define DEFAULT_FOCAL_LENGTH 1.0f

class PCD {
	public:
		enum PCD_field_type{
			X,
			Y,
			Z,
			RGB,
			NORMAL,
			MOMENT_INVARIANT
		};

		enum PCD_data_storage {
			ASCII,
			BINARY,
			NONE
		};

		enum PCD_file_format{
			VERSION,
			FIELDS,
			SIZE,
			TYPE,
			COUNT,
			WIDTH,
			HEIGHT,
			VIEWPOINT,
			POINTS,
			DATA,
			POINT_CLOUD
		};

		struct PCD_Field{
			PCD_field_type type;
			int size;
			int count;
			char data_type;
		};

		struct Quaternion{
			float r;
			float i;
			float j;
			float k;
		} ;

		struct Camera{
			Quaternion position;
			Quaternion rotation;
			float focal_length;
		};

		struct Plane {float a,b,c,d;};
		struct Cube {float x1,y1,z1,x2,y2,z2;};

		struct KdTreeNode {
			bool isLeaf;
			union {
				int index;
				float value;
			};
			KdTreeNode* left;
			KdTreeNode* right;
		};

		float viewpoint[7];
		Camera camera;
		int numFields;
		PCD_Field* fields;
		char* data;
		float* float_data;
		int numPoints;
		int capacity;
		PCD_data_storage data_storage;
		Cube boundingBox;
		KdTreeNode* kdtreeLeaves;
		std::vector<KdTreeNode*> kdtreeBranches;
		KdTreeNode* kdtreeRoot;
		int kdtreeDepth;

		PCD(const char* fileName);
		PCD(int size);
		~PCD();
		bool expand(int size);
		void drawLine(float x1,float y1,float z1,float x2,float y2,float z2,float rgb,float interval);
		void drawAxis(float x,float y,float z,float* rotation,float length,float interval);
		void loadDescriptor(const char* filename);
		void writeToPCD(const char* filename);
		void writeToPLY(const char* filename);
		void writeClustersToPCD(std::vector<std::vector<int>> *indices,const char* dir);

		Plane segmentPlane(int iter,float threshold,float inlierRatio);
		void filterPlane(std::vector<int> *ind,Plane p,float threshold,bool positiveOrNegative);
		PCD* extractIndices(std::vector<int> *ind);

		Cube getBoundingBox();
		KdTreeNode* buildKdTreeNode(int* pointList, int n, int depth);
		void buildKdTree();
		void searchKdTree(std::vector<int> *indices,float x,float y,float z,float distance);
		void searchKdTreeNode(KdTreeNode* node,std::vector<int> *indices,Cube* range,Cube cell,int depth);
		void getIndexFromSubTree(KdTreeNode* node, std::vector<int> *indices);
		void euclideanCluster(std::vector<std::vector<int>> *indices,float distance,size_t minSize,size_t maxSize,size_t maxClusters);

		static PCD* LoadFromMatrix(const char* fileName);
		static PCD* LoadFromBundler(const char* fileName);
		static PCD* LoadFromKITTI(const char* fileName);
		static PCD* LoadFromCameraParameter(const char* dir);
		static PCD* LoadFromPTS(const char* fileName);
		static float colormap(float f);
};
