#include "descriptor.h"

Descriptor::Descriptor(const char* filename) {
	FILE* f = fopen(filename,"r");
	if (!f) {
		printf("File not found: %s\n", filename);
		return;
	}
	char buf[256];
	while (fgets(buf, 256, f)) {
		if (sscanf(buf, "POINTS %d", &numPoints) == 1) {
		} else if (sscanf(buf,"COUNT %d",&dimension)==1) {
		} else if (strncmp(buf,"DATA ascii",10)==0) 
			break;
	}
	float val;
	char* c;
	data = new float[numPoints * dimension];
	for (int i=0;i<numPoints;i++) {
		fgets(buf,256,f);
		c = buf;
		for (int j=0;j<dimension;j++) {
			val = strtod(c,&c);
			data[i * dimension + j] = val;
		}
	}
	fclose(f);
}

Descriptor::~Descriptor() {
	if (data) delete[] data;
}

void Descriptor::kMeansClustering(std::vector<std::vector<int>> *indices) {
	int k = 10;
//	int k = (int) sqrt(numPoints / 2);
	//initialize
	for (int i=0;i<k;i++) {
		std::vector<int> v;
		indices->push_back(v);
	}
	int* associations = new int[numPoints]();
	float* means = new float[k * dimension];
	int iterations = 0;

	//choose initial means
	for (int i=0;i<k;i++) {
		int index = rand() % numPoints;
		memcpy(means + i * dimension, data + index * dimension, dimension * sizeof(float));
	}

	bool update = true;
	while (update) {
		update = false;
		for (int i=0;i<k;i++)
			(*indices)[i].clear();

		//assignment step
		for (int i=0;i<numPoints;i++) {
			//find nearest mean
			int minIndex;
			float minDistance;
			for (int j=0;j<k;j++) {
				float distance = 0;
				for (int l=0;l<dimension;l++)
					distance += (data[i*dimension+l] - means[j*dimension+l]) * (data[i*dimension+l] - means[j*dimension+l]);
				if (j==0 || distance < minDistance) {
					minIndex = j;
					minDistance = distance;
				}
			}
			if (minIndex != associations[i])
				update = true;
			associations[i] = minIndex;
			(*indices)[minIndex].push_back(i);
		}

		//update step
		for (int i=0;i<k;i++) {
			for (int j=0;j<dimension;j++) {
				means[i * dimension + j] = 0;
				size_t clusterSize = (*indices)[i].size();
				for (size_t l=0;l<clusterSize;l++) {
					means[i * dimension + j] += data[(*indices)[i][l] * dimension + j];
				}
				means[i * dimension + j] /= clusterSize;
			}
		}
		iterations++;
	}

//	printf("K-Means: Found %lu clusters (%d iterations)\n",indices->size(),iterations);

	delete[] associations;
	delete[] means;
}
