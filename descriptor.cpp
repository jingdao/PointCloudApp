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
}

Descriptor::~Descriptor() {
	if (data) delete[] data;
}
