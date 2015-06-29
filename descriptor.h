#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

class Descriptor {
	public:
		int numPoints;
		int dimension;
		float* data;

		Descriptor(const char* filename);
		~Descriptor();
};
