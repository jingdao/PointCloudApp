#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <math.h>
#define BASE_HASH_CONSTANT 0.618033988
#define STEP_HASH_CONSTANT 0.707106781
#define STRING_HASH_CONSTANT 5381
#define CANDIDATE_PLANES 8

struct HPoint {
	int x,y,z;
	unsigned char r,g,b;
	int label;
	int curvature;
};
struct Point {
	float x,y,z;
};
struct HPCD {
	int numGrid;
	int numPoints;
	int maxSize;
	float leafSize;
	float minX,minY,minZ,maxX,maxY,maxZ;
	HPoint** data;
	bool hasColor;
};
struct Plane {
	double a,b,c,d;
};
enum PCD_data_storage {
	ASCII,
	BINARY,
	NONE
};

inline int baseHash(int size, int hashKey) {
	return (int)(size*((BASE_HASH_CONSTANT*hashKey) - (int)(BASE_HASH_CONSTANT*hashKey)));
}

inline int stepHash(int size, int hashKey) {
	int res = (int)(size*((STEP_HASH_CONSTANT*hashKey) - (int)(STEP_HASH_CONSTANT*hashKey)));
	//make step size odd since table size is power of 2
	return res % 2 ? res : res + 1; 
}

inline int getIntKey(int x,int y,int z) {
	int h = STRING_HASH_CONSTANT;
	h = (h << 5) + h + x;
	h = (h << 5) + h + y;
	h = (h << 5) + h + z;
	if (h < 0)
		return -h;
	else return h;
}

HPCD* HPCD_Init(const char* inFile, float resolution) {
	HPCD* res = new HPCD;
	FILE* f = fopen(inFile, "rb");
	if (!f) {
		printf("%s not found\n", inFile);
		return NULL;
	}
	char buf[256];
	std::vector<float> float_data;
	std::vector<unsigned char> color_data;
	std::vector<int> label_data;
	PCD_data_storage data_storage = NONE;
	int totalPoints = 0;
	while (fgets(buf, 256, f)) {
		if (sscanf(buf, "POINTS %d", &totalPoints) == 1) {
		}
		else if (strncmp(buf, "FIELDS x y z rgb", 16) == 0) {
			res->hasColor = true;
		}
		else if (strncmp(buf, "FIELDS x y z", 12) == 0) {
			res->hasColor = false;
		}
		else if (strncmp(buf, "DATA ascii", 10) == 0) {
			data_storage = ASCII;
			break;
		}
	}
	if (data_storage != ASCII)
		return NULL;

	long int pos = ftell(f);
	fseek(f, 0, SEEK_END);
	long int end = ftell(f);
	char* c = new char[end - pos];
	fseek(f, pos, SEEK_SET);
	fread(c, 1, end - pos, f);

	float x, y, z;
	unsigned char r, g, b;
	x = strtod(c, &c);
	y = strtod(c, &c);
	z = strtod(c, &c);
	if (res->hasColor) {
		strtol(c, &c, 10);
	}
	res->minX = res->maxX = x;
	res->minY = res->maxY = y;
	res->minZ = res->maxZ = z;
	
	//printf("start parsing %s\n", inFile);
	if (res->hasColor) {
		for (int i = 1; i<totalPoints; i++) {
			int rgb = 0;
			x = strtod(c, &c);
			y = strtod(c, &c);
			z = strtod(c, &c);
			rgb = strtol(c, &c, 10);
			r = (rgb >> 16) & 0xFF;
			g = (rgb >> 8) & 0xFF;
			b = rgb & 0xFF;
			if (x < res->minX) res->minX = x;
			else if (x > res->maxX) res->maxX = x;
			if (y < res->minY) res->minY = y;
			else if (y > res->maxY) res->maxY = y;
			if (z < res->minZ) res->minZ = z;
			else if (z > res->maxZ) res->maxZ = z;
			float_data.push_back(x);
			float_data.push_back(y);
			float_data.push_back(z);
			color_data.push_back(r);
			color_data.push_back(g);
			color_data.push_back(b);
		}
	}
	else {
		for (int i = 1; i<totalPoints; i++) {
			x = strtod(c, &c);
			y = strtod(c, &c);
			z = strtod(c, &c);
			if (x < res->minX) res->minX = x;
			else if (x > res->maxX) res->maxX = x;
			if (y < res->minY) res->minY = y;
			else if (y > res->maxY) res->maxY = y;
			if (z < res->minZ) res->minZ = z;
			else if (z > res->maxZ) res->maxZ = z;
			float_data.push_back(x);
			float_data.push_back(y);
			float_data.push_back(z);
		}
	}
	fclose(f);
	//printf("end parsing %s\n", inFile);
	totalPoints = float_data.size() / 3;
	float minDist = res->maxX - res->minX;
	if (res->maxY - res->minY < minDist)
		minDist = res->maxY - res->minY;
	if (res->maxZ - res->minZ < minDist)
		minDist = res->maxZ - res->minZ;
	res->leafSize = resolution;
	res->numGrid = minDist / res->leafSize;
	//res->leafSize = minDist / res->numGrid;
	res->maxSize = 8;
	while (res->maxSize < 4 * totalPoints)
		res->maxSize *= 2;
	res->data = new HPoint*[res->maxSize]();
	int i, j = 0, k, l = 0;
	res->numPoints = 0;
	if (res->hasColor){
		for (i = 0; i<totalPoints; i++) {
			float x = float_data[j++];
			float y = float_data[j++];
			float z = float_data[j++];
			unsigned char r = color_data[l++];
			unsigned char g = color_data[l++];
			unsigned char b = color_data[l++];
			int xi = (int)((x - res->minX) / res->leafSize);
			int yi = (int)((y - res->minY) / res->leafSize);
			int zi = (int)((z - res->minZ) / res->leafSize);
			int ikey = getIntKey(xi, yi, zi);
			int key = baseHash(res->maxSize, ikey);
			int step = stepHash(res->maxSize, ikey);
			for (k = 0; k < res->maxSize; k++) {
				HPoint* h = res->data[key];
				if (!h) {
					HPoint* p = new HPoint;
					p->x = xi;
					p->y = yi;
					p->z = zi;
					p->r = r;
					p->g = g;
					p->b = b;
					res->data[key] = p;
					res->numPoints++;
					break;
				}
				else if (h->x == xi && h->y == yi && h->z == zi){
					break;
				}
				else {
					key += step;
					key %= res->maxSize;
				}
			}
		}
	}
	else {
		for (i = 0; i<totalPoints; i++) {
			float x = float_data[j++];
			float y = float_data[j++];
			float z = float_data[j++];
			int xi = (int)((x - res->minX) / res->leafSize);
			int yi = (int)((y - res->minY) / res->leafSize);
			int zi = (int)((z - res->minZ) / res->leafSize);
			int ikey = getIntKey(xi, yi, zi);
			int key = baseHash(res->maxSize, ikey);
			int step = stepHash(res->maxSize, ikey);
			for (k = 0; k < res->maxSize; k++) {
				HPoint* h = res->data[key];
				if (!h) {
					HPoint* p = new HPoint;
					p->x = xi;
					p->y = yi;
					p->z = zi;
					p->r = 255;
					p->g = 255;
					p->b = 255;
					res->data[key] = p;
					res->numPoints++;
					break;
				}
				else if (h->x == xi && h->y == yi && h->z == zi){
					break;
				}
				else {
					key += step;
					key %= res->maxSize;
				}
			}
		}
	}
	printf("Processed point cloud (numPoints:%d maxSize:%d leafSize:%f)\n", res->numPoints, res->maxSize, res->leafSize);
	printf("Bounding box: x:(%.2f %.2f) y:(%.2f %.2f) z:(%.2f %.2f)\n", res->minX, res->maxX, res->minY, res->maxY, res->minZ, res->maxZ);
	return res;
}

void HPCD_resize(HPCD* res) { //to save memory
	int newSize = 8; 
	while (newSize < 4 * res->numPoints)
		newSize *= 2;
	HPoint** newdata = new HPoint*[newSize]();
	int i,k;
	for (i=0;i<res->maxSize;i++) {
		HPoint* h = res->data[i];
		if (!h)
			continue;
		int ikey = getIntKey(h->x,h->y,h->z);
		int key = baseHash(newSize,ikey);
		int step = stepHash(newSize,ikey);
		for (k = 0; k < newSize; k++) {
			if (!newdata[key]) {
				newdata[key] = h;
				break;
			} else {
				key += step;
				key %= newSize;
			}
		}
	}
	delete[] res->data;
	res->data = newdata;
	res->maxSize = newSize;
	printf("Processed point cloud (numPoints:%d maxSize:%d leafSize:%f)\n",res->numPoints,res->maxSize,res->leafSize);
}

void HPCD_write(const char* filename, HPCD* pointcloud) {
	FILE* f = fopen(filename, "w");
	if (!f) {
		printf("Cannot write to file: %s\n", filename);
		return;
	}
	else if (pointcloud->hasColor) {
		fprintf(f, "# .PCD v0.7 - Point Cloud Data file format\n"
			"VERSION 0.7\n"
			"FIELDS x y z rgb\n"
			"SIZE 4 4 4 4\n"
			"TYPE F F F I\n"
			"COUNT 1 1 1 1\n"
			"WIDTH %d\n"
			"HEIGHT 1\n"
			"VIEWPOINT 0 0 0 1 0 0 0\n"
			"POINTS %d\n"
			"DATA ascii\n", pointcloud->numPoints, pointcloud->numPoints);
		for (int i = 0; i<pointcloud->maxSize; i++) {
			HPoint *h = pointcloud->data[i];
			if (h) {
				int rgb = (h->r << 16) | (h->g << 8) | h->b;
				fprintf(f, "%f %f %f %d\n",
					pointcloud->minX + h->x * pointcloud->leafSize,
					pointcloud->minY + h->y * pointcloud->leafSize,
					pointcloud->minZ + h->z * pointcloud->leafSize,
					rgb);
			}
		}
	}
	else {
		fprintf(f, "# .PCD v0.7 - Point Cloud Data file format\n"
			"VERSION 0.7\n"
			"FIELDS x y z\n"
			"SIZE 4 4 4\n"
			"TYPE F F F\n"
			"COUNT 1 1 1\n"
			"WIDTH %d\n"
			"HEIGHT 1\n"
			"VIEWPOINT 0 0 0 1 0 0 0\n"
			"POINTS %d\n"
			"DATA ascii\n", pointcloud->numPoints, pointcloud->numPoints);
		for (int i = 0; i<pointcloud->maxSize; i++) {
			HPoint *h = pointcloud->data[i];
			if (h) {
				fprintf(f, "%f %f %f\n",
					pointcloud->minX + h->x * pointcloud->leafSize,
					pointcloud->minY + h->y * pointcloud->leafSize,
					pointcloud->minZ + h->z * pointcloud->leafSize);
			}
		}
	}
	fclose(f);
	printf("Wrote %d points to %s\n", pointcloud->numPoints, filename);
}

void HPCD_writePoints(char* filename,std::vector<Point> *points) {
	FILE* f = fopen(filename, "w");
	if (!f) {
		printf("Cannot write to file: %s\n", filename);
		return;
	}
	fprintf(f,"# .PCD v0.7 - Point Cloud Data file format\n"
	"VERSION 0.7\n"
	"FIELDS x y z\n"
	"SIZE 4 4 4\n"
	"TYPE F F F\n"
	"COUNT 1 1 1\n"
	"WIDTH %lu\n"
	"HEIGHT 1\n"
	"VIEWPOINT 0 0 0 1 0 0 0\n"
	"POINTS %lu\n"
	"DATA ascii\n",points->size(),points->size());
	for (size_t i=0;i<points->size();i++) {
		fprintf(f,"%f %f %f\n",points->at(i).x,points->at(i).y,points->at(i).z);
	}
	fclose(f);
	printf("Wrote %lu points to %s\n",points->size(),filename);
}

Plane segmentPlane(HPCD* cloud,std::vector<Point> *filtered, int iter, float segThreshold, float inlierRatio) {
	int maxInliers = 0;
	int optimumInliers = (int) (inlierRatio * cloud->numPoints);
	int distanceThreshold;
	int previousNumPoints = cloud->numPoints;
	Plane bestPlane = { 0, 0, 1, 0 };
	int i, k = 0;
	int* pointdata = new int[cloud->numPoints * 3];
	for (int j = 0; j<cloud->maxSize; j++) {
		if (cloud->data[j]) {
			pointdata[k++] = cloud->data[j]->x;
			pointdata[k++] = cloud->data[j]->y;
			pointdata[k++] = cloud->data[j]->z;
		}
	}
	for (i = 0; i<iter; i++) {
		//Pick 3 points
		int *p0 = pointdata + (rand() % cloud->numPoints * 3);
		int *p1 = pointdata + (rand() % cloud->numPoints * 3);
		int *p2 = pointdata + (rand() % cloud->numPoints * 3);
		int p0x = p0[0], p0y = p0[1], p0z = p0[2];
		int p1x = p1[0], p1y = p1[1], p1z = p1[2];
		int p2x = p2[0], p2y = p2[1], p2z = p2[2];
		Plane currentPlane = {
			(p1y - p0y)*(p2z - p0z) - (p2y - p0y)*(p1z - p0z),
			(p1z - p0z)*(p2x - p0x) - (p2z - p0z)*(p1x - p0x),
			(p1x - p0x)*(p2y - p0y) - (p2x - p0x)*(p1y - p0y),
			0
		};
		currentPlane.d = -(currentPlane.a * p0x + currentPlane.b * p0y + currentPlane.c * p0z);
		if (currentPlane.a == 0 && currentPlane.b == 0 && currentPlane.c == 0)
			continue; //picked collinear points
		distanceThreshold = sqrt(
			currentPlane.a * currentPlane.a +
			currentPlane.b * currentPlane.b +
			currentPlane.c * currentPlane.c) * segThreshold / cloud->leafSize;
		int numInliers = 0;
		for (int j = 0; j<cloud->numPoints; j++) {
			if (fabs(currentPlane.a * pointdata[j * 3] +
				currentPlane.b * pointdata[j * 3 + 1] +
				currentPlane.c * pointdata[j * 3 + 2] +
				currentPlane.d)
				< distanceThreshold)
				numInliers++;
		}
		if (numInliers > maxInliers) {
			maxInliers = numInliers;
			bestPlane = currentPlane;
			if (maxInliers > optimumInliers)
				break;
		}
	}
	double norm = sqrt(
		bestPlane.a * bestPlane.a +
		bestPlane.b * bestPlane.b +
		bestPlane.c * bestPlane.c);
	distanceThreshold = norm * segThreshold / cloud->leafSize;
	double RMSE = 0;
	for (int j = 0; j<cloud->maxSize; j++) {
		HPoint* h = cloud->data[j];
		if (h && fabs(bestPlane.a * h->x +
			bestPlane.b * h->y +
			bestPlane.c * h->z +
			bestPlane.d)
			< distanceThreshold) {
			double d = 1.0 * (bestPlane.a * h->x + bestPlane.b * h->y + bestPlane.c * h->z + bestPlane.d) / distanceThreshold;
			RMSE += d * d;
			//delete cloud->data[j];
			Point p = {
				cloud->minX + h->x * cloud->leafSize,
				cloud->minY + h->y * cloud->leafSize,
				cloud->minZ + h->z * cloud->leafSize
			};
			filtered->push_back(p);
			cloud->data[j] = NULL;
			cloud->numPoints--;
		}
	}
	bestPlane.a /= norm;
	bestPlane.b /= norm;
	bestPlane.c /= norm;
	bestPlane.d /= norm;
	printf("RANSAC: %d points %d iters (%f,%f,%f,%f) (RMSE %f)\n", cloud->numPoints, i, bestPlane.a, bestPlane.b, bestPlane.c, bestPlane.d, sqrt(RMSE / maxInliers));
	delete[] pointdata;
	HPCD_resize(cloud);
	return bestPlane;
}

int main(int argc, char* argv[]) {
	if (argc < 4) {
		printf("Usage: %s source.pcd target.pcd output.pcd\n",argv[0]);
		exit(0);
	}

	HPCD* source = HPCD_Init(argv[1],0.01);
	HPCD* target = HPCD_Init(argv[2],0.01);
	Plane sourcePlane[CANDIDATE_PLANES];
	Plane targetPlane[CANDIDATE_PLANES];
	std::vector<Point> sourceFiltered;
	std::vector<Point> targetFiltered;
	for (int i=0;i<CANDIDATE_PLANES;i++) {
		sourcePlane[i] = segmentPlane(source,&sourceFiltered,5000,0.1,0.2);
		targetPlane[i] = segmentPlane(target,&targetFiltered,5000,0.1,0.2);
		char buffer[128];
		sprintf(buffer,"test/%d-cloud.pcd",2*i);
		HPCD_writePoints(buffer,&sourceFiltered);
		sprintf(buffer,"test/%d-cloud.pcd",2*i+1);
		HPCD_writePoints(buffer,&targetFiltered);
		sourceFiltered.clear();
		targetFiltered.clear();
	}

	HPCD_write("source.pcd",source);
	HPCD_write("target.pcd",target);

}
