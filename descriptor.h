#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <math.h>

class Descriptor {
	public:
		int numPoints;
		int dimension;
		float* data;

		Descriptor(const char* filename);
		~Descriptor();

		void kMeansClustering(std::vector<std::vector<int>> *indices);
};
