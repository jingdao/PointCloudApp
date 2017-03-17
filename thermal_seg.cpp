#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <math.h>
#include <unistd.h>
#include <SDL/SDL.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <ft2build.h>
#include FT_FREETYPE_H
#define BASE_HASH_CONSTANT 0.618033988
#define STEP_HASH_CONSTANT 0.707106781
#define STRING_HASH_CONSTANT 5381
#define USE_THERMAL 1
#define DISPLAY_LABEL 1

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
struct Cluster {
	int index, size;
};
struct Feature {
	float length,width,height,elevation;
};

FT_Library ft;
FT_Face face;
FT_GlyphSlot glyph;
GLuint textures;
const int labelWidth=80;
const int labelHeight=30,fontpixels=20;
unsigned char raster[labelWidth * labelHeight * 3];
double cameraX=20,cameraY=20,cameraZ=10;
double centerX=0,centerY=0,centerZ=0;
double upX=0,upY=0,upZ=1;
int mouseIndex = 0;
int previousX,previousY;
double scrollSpeed = 1.01;
HPCD* source;
std::vector< std::vector<float> > boxes;
std::vector<Feature> features;
std::vector<char*> classes;
char* className[4] = {NULL,"person","display","light"};
int classCount[4];

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

HPCD* HPCD_copy(HPCD* cloud, bool deep) {
	HPCD* h = new HPCD;
	memcpy(h, cloud, sizeof(HPCD));
	h->data = new HPoint*[h->maxSize];
	if (deep) {
		for (int i = 0; i < cloud->maxSize; i++) {
			HPoint* p = cloud->data[i];
			if (p) {
				h->data[i] = new HPoint;
				*(h->data[i]) = *p;
			}
			else {
				h->data[i] = NULL;
			}
		}
	}
	else {
		memcpy(h->data, cloud->data, h->maxSize * sizeof(HPoint*));
	}
	return h;
}

int HPCD_find(HPCD* cloud, int x, int y, int z) {
	int ikey = getIntKey(x, y, z);
	int j = baseHash(cloud->maxSize, ikey);
	int step = stepHash(cloud->maxSize, ikey);
	for (int k = 0; k<cloud->maxSize; k++) {
		HPoint* h = cloud->data[j];
		if (!h) {
			return -1;
		}
		else if (h->x == x && h->y == y && h->z == z){
			return j;
		}
		else {
			j += step;
			j %= cloud->maxSize;
		}
	}
	return -1;
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
	filtered->clear();
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

int compareCluster(const void* v1, const void* v2) {
	Cluster* c1 = (Cluster*)v1;
	Cluster* c2 = (Cluster*)v2;
	return c2->size - c1->size;
}

int euclideanClustering(HPCD* cloud, int clusterThreshold, int totalClusters, int minsize) {
	bool* visited = new bool[cloud->maxSize]();
	std::vector<Cluster> clusters;
	int numClusters = 0;
	int numCombinations = 1;
	for (int i = 0; i < clusterThreshold; i++)
		numCombinations *= 6;
	for (int i = 0; i<cloud->maxSize; i++) {
		HPoint* seed = cloud->data[i];
		if (!seed || visited[i] || seed->label < 0) continue;
		Cluster c = { numClusters, 0 };
		std::vector<int> Q;
		Q.push_back(i);
		visited[i] = true;
		while (Q.size() > 0) {
			int p = Q[Q.size() - 1];
			Q.pop_back();
			HPoint* h = cloud->data[p];
			h->label = numClusters;
			c.size++;
			
			for (int k = 0; k < numCombinations; k++) {
				int l = k;
				int xi = h->x;
				int yi = h->y;
				int zi = h->z;
				for (int m = 0; m < clusterThreshold;m++) {
					switch (l % 6) {
					case 0: if (xi >= h->x) xi += 1; break;
					case 1: if (xi <= h->x) xi -= 1; break;
					case 2: if (yi >= h->y) yi += 1; break;
					case 3: if (yi <= h->y) yi -= 1; break;
					case 4: if (zi >= h->z) zi += 1; break;
					case 5: if (zi <= h->z) zi -= 1; break;
					}
					int j = HPCD_find(cloud, xi, yi, zi);
					HPoint* h2 = cloud->data[j];
					if (j >= 0 && !visited[j] && h2->label >= 0) {
						Q.push_back(j);
						visited[j] = true;
					}
					l /= 6;
				}
			}	
		}
		//		printf("Cluster %d: %d\n",numClusters,c.size);
		clusters.push_back(c);
		numClusters++;
	}
	Cluster *sortedClusters = new Cluster[numClusters];
	memcpy(sortedClusters, clusters.data(), numClusters*sizeof(Cluster));
	qsort(sortedClusters, numClusters, sizeof(Cluster), compareCluster);
	int *idx = new int[numClusters];
	for (int i = 0; i<numClusters; i++)
		idx[i] = -1;
	int k;
	for (k = 0; k<totalClusters && k<numClusters && sortedClusters[k].size > minsize; k++)
		idx[sortedClusters[k].index] = k;
	for (int i = 0; i<cloud->maxSize; i++) {
		HPoint* h = cloud->data[i];
		if (h && h->label>=0)
			h->label = idx[h->label];
	}
	delete[] idx;
	delete[] visited;
	delete[] sortedClusters;
	return k;
}

void getPCA(std::vector<Point> *cloud, std::vector<float> *box) {
	double cov[4] = {}; //column major
	Point center = {};
	for (int i = 0; i < cloud->size(); i++) {
		center.x += cloud->at(i).x;
		center.y += cloud->at(i).y;
		center.z += cloud->at(i).z;
	}
	center.x /= cloud->size(); 
	center.y /= cloud->size();
	center.z /= cloud->size();
	for (int j = 0; j<cloud->size(); j++) {
		float deltaP[2] = {
			cloud->at(j).x- center.x,
			cloud->at(j).y- center.y,
		};
		cov[0] += deltaP[0] * deltaP[0];
		cov[1] += deltaP[1] * deltaP[0];
		cov[2] += deltaP[0] * deltaP[1];
		cov[3] += deltaP[1] * deltaP[1];
	}
	cov[0] /= cloud->size() * cloud->size();
	cov[1] /= cloud->size() * cloud->size();
	cov[2] /= cloud->size() * cloud->size();
	cov[3] /= cloud->size() * cloud->size();
	float trace = cov[0] + cov[3];
	float det = cov[0] * cov[3] - cov[1] * cov[2];
	float L1 = trace / 2 + sqrt(trace*trace / 4 - det);
	float L2 = trace / 2 - sqrt(trace*trace / 4 - det);
	float minScale[3], maxScale[3];
	float v[9] = {
		0,0,0,
		0,0,0,
		0,0,1
	};
	if (cov[2] != 0) {
		v[0] = L1 - cov[3];
		v[1] = L2 - cov[3];
		v[3] = v[4] = cov[2];
	}
	else if (cov[1] != 0) {
		v[0] = v[1] = cov[1];
		v[3] = L1 - cov[0];
		v[4] = L2 - cov[0];
	}
	else {
		v[0] = v[4] = 1;
	}
	float m1 = sqrt(v[0] * v[0] + v[3] * v[3]);
	float m2 = sqrt(v[1] * v[1] + v[4] * v[4]);
	v[0] /= m1;
	v[3] /= m1;
	v[1] /= m2;
	v[4] /= m2;
	for (int j = 0; j<cloud->size(); j++) {
		for (int i = 0; i<3; i++) {
			float dotProduct =
				cloud->at(j).x * v[i * 3] +
				cloud->at(j).y * v[i * 3 + 1] +
				cloud->at(j).z * v[i * 3 + 2];
			if (j == 0 || dotProduct < minScale[i])
				minScale[i] = dotProduct;
			if (j == 0 || dotProduct > maxScale[i])
				maxScale[i] = dotProduct;
		}
	}
	float bbCenter[3] = {0,0,0};
	for (int i = 0; i<3; i++) {
		bbCenter[0] += (minScale[i] + maxScale[i]) / 2 * v[i * 3];
		bbCenter[1] += (minScale[i] + maxScale[i]) / 2 * v[i * 3 + 1];
		bbCenter[2] += (minScale[i] + maxScale[i]) / 2 * v[i * 3 + 2];
	}
	for (int i = 0; i<8; i++) {
		float coords[3];
		for (int j = 0; j<3; j++) {
			coords[j] = bbCenter[j];
			for (int axis = 0; axis<3; axis++) {
				float sign = (i & 1 << axis) ? 1 : -1;
				coords[j] += sign * (maxScale[axis]-minScale[axis]) / 2 * v[axis * 3 + j];
			}
			box->push_back(coords[j]);
		}
	}
	Feature f = {
		maxScale[0] - minScale[0],
		maxScale[1] - minScale[1],
		maxScale[2] - minScale[2],
		bbCenter[2]
	};
//	printf("%f %f %f %f\n",f.length,f.width,f.height,f.elevation);
	int classMember=0;
	if (f.length > 0.21 && f.length < 2.36 && f.width > 0.14 && f.width < 1.34 && f.height > 0.73 && f.height < 3.03 && f.elevation < 0.5)
		classMember=1;
//	if (f.length > 0.3 && f.length < 2 && f.width > 0.3 && f.height > 0.5 && f.elevation < 0.5)
//		classMember=1;
	else if (f.length > 1 && f.length < 5 && f.height > 0.5 && f.elevation < 1.5)
		classMember=2;
	else if (f.length < 0.4 && f.width < 0.3 && f.elevation > 0.7)
		classMember=3;
	classes.push_back(className[classMember]);
	classCount[classMember]++;
	features.push_back(f);
}

void getBoundingBox(HPCD* h) {
	boxes.clear();
	int xrange = (h->maxX - h->minX) / h->leafSize + 1;
	int yrange = (h->maxY - h->minY) / h->leafSize + 1;
	std::vector< std::vector<Point> > points;
	int numClusters = 0;
	for (int i = 0; i < h->maxSize; i++) {
		HPoint* p = h->data[i];
		if (p && p->label >= numClusters)
			numClusters = p->label + 1;
	}
	boxes.resize(numClusters);
	points.resize(numClusters);
	for (int i = 0; i < h->maxSize; i++) {
		HPoint* p = h->data[i];
		if (p && p->label >= 0) {
//			Point q = { (float)(p->x - xrange / 2), (float)(p->y - yrange / 2), (float)p->z};
			Point q = {
				h->minX + p->x * h->leafSize,
				h->minY + p->y * h->leafSize,
				h->minZ + p->z * h->leafSize,
			};
			points[p->label].push_back(q);
		}
	}
	for (int i = 0; i < numClusters; i++) {
		getPCA(&points[i], &boxes[i]);
	}
}

void drawLine(int i, int j, int k) {
	glBegin(GL_LINES);
	glVertex3d(boxes[i][j * 3], boxes[i][j * 3 + 1], boxes[i][j * 3 + 2]);
	glVertex3d(boxes[i][k * 3], boxes[i][k * 3 + 1], boxes[i][k * 3 + 2]);
	glEnd();
}

void render_text(const char *text, unsigned char *data) {
	memset(data,0,labelWidth*labelHeight*3);
	const char *p;
	unsigned char color[] = {255,255,255};
	int x = 0;
	for(p = text; *p; p++) {
		if(FT_Load_Char(face, *p, FT_LOAD_RENDER))
			continue;
		unsigned char *src = glyph->bitmap.buffer;
		int k=0;
		for (int i=0;i<glyph->bitmap.rows;i++) {
			unsigned char *dest = data + ((glyph->bitmap_top - i + fontpixels/2) * labelWidth + x + glyph->bitmap_left)* 3;
			for (int j=0;j<glyph->bitmap.width;j++) {
				memset(dest,*src,3); // draw in grayscale
				//*dest = *src; //draw in red
				src++;
				dest+=3;
			}
		}
		x += glyph->bitmap_left + glyph->bitmap.width;
//		x += glyph->bitmap.width;
		if (x >= labelWidth)
			break;
	}
}

void drawText(int i) {
	glRasterPos3f((boxes[i][12] + boxes[i][15] + boxes[i][18] + boxes[i][21]) / 4, (boxes[i][13] + boxes[i][16] + boxes[i][19] + boxes[i][22]) / 4, boxes[i][14]);
	render_text(classes[i],raster);
	glDrawPixels(labelWidth,labelHeight,GL_RGB,GL_UNSIGNED_BYTE,raster);
}

void draw() {
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glPushMatrix();
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(cameraX,cameraY,cameraZ,centerX,centerY,centerZ,upX,upY,upZ);

	glLineWidth(1.0);
	glColor3ub(255, 255, 255);
	for (int i = 0; i<boxes.size(); i++) {
#if DISPLAY_LABEL
		if (classes[i]) {
#endif
			drawLine(i, 0, 1);
			drawLine(i, 0, 2);
			drawLine(i, 1, 3);
			drawLine(i, 2, 3);
			drawLine(i, 0, 4);
			drawLine(i, 1, 5);
			drawLine(i, 2, 6);
			drawLine(i, 3, 7);
			drawLine(i, 4, 5);
			drawLine(i, 4, 6);
			drawLine(i, 5, 7);
			drawLine(i, 6, 7);
#if DISPLAY_LABEL
			drawText(i);
		}
#endif
	}

	glPointSize(1.0);
	glBegin(GL_POINTS);
	for (int i=0;i<source->maxSize;i++) {
		HPoint* h = source->data[i];
		if (h) {
			glColor3ub(h->r,h->g,h->b);
			glVertex3d(
				source->minX+h->x*source->leafSize,
				source->minY+h->y*source->leafSize,
				source->minZ+h->z*source->leafSize);
		}
	}
	glEnd();

	glFlush();
	SDL_GL_SwapBuffers();
	glPopMatrix();
	glPopAttrib();
}

int main(int argc, char* argv[]) {
	if (argc < 2) {
		printf("Usage: %s source.pcd\n",argv[0]);
		exit(0);
	}

	srand(0);
	FT_Init_FreeType(&ft);
	FT_New_Face(ft,"/usr/share/fonts/truetype/freefont/FreeSans.ttf",0,&face);
	FT_Set_Pixel_Sizes(face,fontpixels,fontpixels);
	glyph = face->glyph;

	source = HPCD_Init(argv[1],0.1);
	HPCD* model = HPCD_Init(argv[1],0.1);
#if USE_THERMAL
	int minIntensity=-1,maxIntensity=-1;
	//absolute threshold
	for (int i=0;i<model->maxSize;i++) {
		HPoint* h = model->data[i];
		if (h) {
			int intensity = (h->r << 16) | (h->g << 8) | h->b;
			if (minIntensity < 0 || intensity < minIntensity)
				minIntensity = intensity;
			if (maxIntensity < 0 || intensity > maxIntensity)
				maxIntensity = intensity;
			if (intensity > 60000) {
				h->r = 255; h->g = 0; h->b = 0;
				h->label = 0;
			} else {
				h->r = 100; h->g = 100; h->b = 100;
				h->label = -1;
			}
		}
	}
#else
	std::vector<Point> filtered;
	for (int i=0;i<10;i++)
		segmentPlane(model,&filtered,1000,0.1,0.5);
	for (int i=0;i<model->maxSize;i++) {
		HPoint* h = model->data[i];
		if (h)
			h->label = 0;
	}
#endif
	int n = euclideanClustering(model,3,100,1);
	printf("%d clusters found\n",n);
	memset(classCount,0,4*sizeof(int));
	getBoundingBox(model);
	for (int i=0;i<4;i++)
		printf("%s: %d\n",i==0 ? "None" : className[i], classCount[i]);
#if !USE_THERMAL
	for (int i=0;i<n;i++) {
		int r = rand()%255;
		int g = rand()%255;
		int b = rand()%255;
		for (int j=0;j<source->maxSize;j++) {
			HPoint* h = source->data[j];
			if (h) {
				int id = HPCD_find(model,
					(int)(h->x*source->leafSize/model->leafSize),
					(int)(h->y*source->leafSize/model->leafSize),
					(int)(h->z*source->leafSize/model->leafSize));
				if (id >= 0) {
					if (model->data[id]->label == i) {
#if DISPLAY_LABEL
						if (classes[i]) {
#endif
						h->r = r;
						h->g = g;
						h->b = b;
#if DISPLAY_LABEL
						} else {
							h->r = 100;
							h->g = 100;
							h->b = 100;
						}
#endif
					}
				} else {
					h->r = 100;
					h->g = 100;
					h->b = 100;
				}
			}
		}
	}
#endif

	SDL_Init(SDL_INIT_VIDEO);
	SDL_WM_SetCaption("Point Cloud", NULL);
	SDL_SetVideoMode(1800,1000, 32, SDL_OPENGL);
    glEnable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	glDisable(GL_CULL_FACE);
	glPolygonMode(GL_FRONT, GL_FILL);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(70,(double)640/480,1,1000);

	int interval = 10000;
	SDL_Event event;
	while (SDL_PollEvent(&event)); //clear event buffer
	draw();
	while (true) {
		while (SDL_PollEvent(&event)) {
			switch(event.type){
				case SDL_KEYDOWN:
					switch( event.key.keysym.sym ){
						case SDLK_LEFT:
						cameraX -= 1;
						break;
						case SDLK_RIGHT:
						cameraX += 1;
						break;
						case SDLK_UP:
						cameraZ += 1;
						break;
						case SDLK_DOWN:
						cameraZ -= 1;
						break;
						default:
						break;
					}
					draw();
					break;
				case SDL_MOUSEBUTTONDOWN:
					if (event.button.button == SDL_BUTTON_WHEELUP) {
						cameraX /= scrollSpeed;
						cameraY /= scrollSpeed;
						cameraZ /= scrollSpeed;
						draw();
					} else if (event.button.button == SDL_BUTTON_WHEELDOWN) {
						cameraX *= scrollSpeed;
						cameraY *= scrollSpeed;
						cameraZ *= scrollSpeed;
						draw();
					} else {
						mouseIndex = event.button.button == SDL_BUTTON_LEFT ? 1 : 2;
						previousX = event.button.x;
						previousY = event.button.y;
					}
					break;
				case SDL_MOUSEBUTTONUP:
					mouseIndex = 0;
					break;
				case SDL_MOUSEMOTION:
					if (mouseIndex == 1) {
						double rho = sqrt(cameraX*cameraX+cameraY*cameraY);
						double xstep = cameraY / rho;
						double ystep = -cameraX / rho;
						cameraX += 0.05 * (event.motion.x-previousX) * xstep;
						cameraY += 0.05 * (event.motion.x-previousX) * ystep;
						cameraZ += 0.05 * (event.motion.y-previousY);
						previousX = event.motion.x;
						previousY = event.motion.y;
						draw();
					}
					break;
				case SDL_QUIT:
					exit(0);
					break;
			}
		}
		usleep(interval);
	}
}
