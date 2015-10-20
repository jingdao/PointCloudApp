#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <math.h>
#define NUM_BINS 100

struct PCD {
	int numPoints;
	float* float_data;
};
enum PCD_data_storage {
	ASCII,
	BINARY,
	NONE
};
struct Quaternion{
	float r;
	float i;
	float j;
	float k;
};

PCD* NewPCD(const char* fileName) {
	PCD* pcd = new PCD();
	PCD_data_storage data_storage = NONE;
	FILE* f = fopen(fileName, "r");
	if (!f) {
		printf("File not found: %s\n", fileName);
		return NULL;
	}
	char buf[256];
	int pointsParsed = 0;
	while (fgets(buf, 256, f)) {
		if (sscanf(buf, "POINTS %d", &pcd->numPoints) == 1) {
			pcd->float_data = (float*)malloc(4 * pcd->numPoints * sizeof(float));
		} else if (strncmp(buf,"DATA ascii",10)==0) {
			data_storage = ASCII;
		} else if (strncmp(buf,"DATA binary_compressed",23)==0) {
			data_storage = BINARY;
			fread(pcd->float_data,sizeof(float),pcd->numPoints*4,f);
			break;
		}
		else if (data_storage == ASCII) {
			if (sscanf(buf, "%f %f %f %f", pcd->float_data+pointsParsed * 4, pcd->float_data+pointsParsed * 4 + 1,
				pcd->float_data+pointsParsed * 4 + 2, pcd->float_data+pointsParsed * 4 + 3) >= 3) {
				pointsParsed++;
			}
		}
	}
	fclose(f);
	return pcd;
}

Quaternion getCentroid(PCD* pcd) {
	float sumX=0, sumY=0, sumZ=0;
	for (int i=0;i<pcd->numPoints;i++) {
		sumX += pcd->float_data[i * 4];
		sumY += pcd->float_data[i * 4 + 1];
		sumZ += pcd->float_data[i * 4 + 2];
	}
	Quaternion q = {0, sumX / pcd->numPoints, sumY / pcd->numPoints, sumZ / pcd->numPoints};
	return q;
}

int main(int argc,char* argv[]) {
	if (argc < 3) {
		printf("Usage: ./rd_desc input.pcd output-cloud.pcd-esf.pcd");
		return 1;
	}

	PCD* cloud = NewPCD(argv[1]);
	if (!cloud)
		return 1;

	Quaternion q = getCentroid(cloud);
	float minX=cloud->float_data[0],maxX=cloud->float_data[0];
	float minY=cloud->float_data[1],maxY=cloud->float_data[1];
	float minZ=cloud->float_data[2],maxZ=cloud->float_data[2];
	for (int i=1;i<cloud->numPoints;i++) {
		if (cloud->float_data[i * 4] < minX) minX = cloud->float_data[i * 4];
		else if (cloud->float_data[i * 4] > maxX) maxX = cloud->float_data[i * 4];
		if (cloud->float_data[i * 4 + 1] < minY) minY = cloud->float_data[i * 4 + 1];
		else if (cloud->float_data[i * 4 + 1] > maxY) maxY = cloud->float_data[i * 4 + 1];
		if (cloud->float_data[i * 4 + 2] < minZ) minZ = cloud->float_data[i * 4 + 2];
		else if (cloud->float_data[i * 4 + 2] > maxZ) maxZ = cloud->float_data[i * 4 + 2];
	}
//	printf("Centroid: %f %f %f\n",q.i,q.j,q.k);
//	printf("Bounding box: x:(%.2f %.2f) y:(%.2f %.2f) z:(%.2f %.2f)\n",minX,maxX,minY,maxY,minZ,maxZ);

	float *r = new float[cloud->numPoints];
	float minR,maxR;
	for (int i=0;i<cloud->numPoints;i++) {
		r[i] = 0;
		r[i] += (cloud->float_data[i*4] - q.i) * (cloud->float_data[i*4] - q.i);
		r[i] += (cloud->float_data[i*4+1] - q.j) * (cloud->float_data[i*4+1] - q.j);
		r[i] += (cloud->float_data[i*4+2] - q.k) * (cloud->float_data[i*4+2] - q.k);
		r[i] = sqrt(r[i]);
		if (i==0 || r[i] < minR)
			minR = r[i];
		if (i==0 || r[i] > maxR)
			maxR = r[i];
	}

	int *count = new int[NUM_BINS]();
	for (int i=0;i<cloud->numPoints;i++) {
		float ratio = (r[i] - minR) / (maxR - minR);
		int idx = (int) (ratio * NUM_BINS);
		count[idx>=NUM_BINS ? NUM_BINS-1 : idx]++;
	}

	FILE* f = fopen(argv[2],"w");
	fprintf(f,"# .PCD v0.7 - Point Cloud Data file format\n"
	"VERSION 0.7\n"
	"FIELDS esf\n"
	"SIZE 4\n"
	"TYPE F\n"
	"COUNT %d\n"
	"WIDTH 1\n"
	"HEIGHT 1\n"
	"VIEWPOINT 0 0 0 1 0 0 0\n"
	"POINTS 1\n"
	"DATA ascii\n",NUM_BINS);
	float scale = 1.0 / cloud->numPoints;
	for (int i=0;i<NUM_BINS;i++) {
		fprintf(f,"%f ",scale * count[i]);
	}
	fclose(f);
//	printf("Saved to %s\n",argv[2]);
	delete[] count;
	delete[] r;
	return 0;
}
