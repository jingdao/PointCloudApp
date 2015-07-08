#pragma once
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dirent.h>
#include <vector>
#include "hashtable.h"
#include "kdtree.h"
#include "descriptor.h"
#define DEFAULT_FOCAL_LENGTH 1.0f

class KdTree;

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

		float viewpoint[7];
		Camera camera;
		int numFields;
		PCD_Field* fields;
		char* data;
		float* float_data;
		int numPoints;
		int capacity;
		PCD_data_storage data_storage;
		KdTree* kdtree;

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
		void writeToOFF(const char* filename);

		static PCD* LoadFromPLY(const char* fileName);
		static PCD* LoadFromMatrix(const char* fileName);
		static PCD* LoadFromBundler(const char* fileName);
		static PCD* LoadFromKITTI(const char* fileName, const char* oxts);
		static PCD* LoadFromCameraParameter(const char* dir);
		static PCD* LoadFromPTS(const char* fileName);
		static PCD* LoadFromCluster(const char* dir);
		static float colormap(float f);

		Quaternion quaternionFromAngle(float rx,float ry, float rz);
		Quaternion quaternionMult(Quaternion qa, Quaternion qb);
		Quaternion quaternionInv(Quaternion q);
		void rotate(Quaternion q);
		void translate(float x, float y, float z);
		Quaternion getCentroid(std::vector<int> *indices);
		int get3DProjection(unsigned int* pixels,int screenWidth, int screenHeight, int pointSize, bool zoomToFit);

		Plane segmentPlane(int iter,float threshold,float inlierRatio);
		void filterPlane(std::vector<int> *ind,Plane p,float threshold,bool positiveOrNegative);
		PCD* extractIndices(std::vector<int> *ind);
		void euclideanClustering(std::vector<std::vector<int>> *indices,float distance,size_t minSize,size_t maxSize,size_t maxClusters);
		void euclideanClusteringKDTree(std::vector<std::vector<int>> *indices,float distance,size_t minSize,size_t maxSize,size_t maxClusters);

};
