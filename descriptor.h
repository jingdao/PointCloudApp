#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <math.h>
#include "pcd.h"
#include "normal.h"

class PCD;

class Descriptor {
	public:
		int numPoints;
		int dimension;
		float* data;
		PCD* pointcloud;

		double principalLengths[3];
		double principalAxes[3][3];
		double bbCenter[3];

		Descriptor(const char* filename);
		Descriptor(PCD* p);
		~Descriptor();

		void getPCA(double* lambda_real,double* v);
		void getPCA_XY(double* lambda_real,double* v);
		void setPCA(double* lambda_real,double* v);
		void kMeansClustering(std::vector<std::vector<int>> *indices);
};
