#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include <math.h>
#include <time.h>
#include "lapack.h"
#define BASE_HASH_CONSTANT 0.618033988
#define STEP_HASH_CONSTANT 0.707106781
#define STRING_HASH_CONSTANT 5381
#define PROFILE 0
#define SAVE_ASSEMBLY 1
#define VERBOSE 0
#define NUM_PARTS 1
#define REDUCED_FEATURE 1
#define USE_DENSITY 1
#define GLOBAL_NORM 0
#define DENSITY_GRID 9

enum PCD_data_storage {
	ASCII,
	BINARY,
	NONE
};

struct HPoint {
	int x,y,z;
	unsigned char r,g,b;
	int label;
};

struct Box {
	float minX,minY,minZ,maxX,maxY,maxZ;
};

struct HPCD {
	int numGrid;
	int numPoints;
	int maxSize;
	float leafSize;
	float minX,minY,minZ,maxX,maxY,maxZ;
	HPoint** data;
	bool hasColor;
	Box oriented_box;
};

struct Dimensions {
	int length,width,height;
};

inline float distance(HPoint a,HPoint b) {
	return sqrt((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) + (a.z-b.z)*(a.z-b.z));
}

inline float angle(float a,float b,float c) {
	return fabs(acos((a*a + b*b - c*c) / (2 * a * b)));
}

inline int baseHash(int size, int hashKey) {
	return (int)(size*((BASE_HASH_CONSTANT*hashKey) - (int)(BASE_HASH_CONSTANT*hashKey)));
}

inline int stepHash(int size, int hashKey) {
	int res = (int)(size*((STEP_HASH_CONSTANT*hashKey) - (int)(STEP_HASH_CONSTANT*hashKey)));
	//make step size odd since table size is power of 2
	return res % 2 ? res : res + 1; 
}

void getCentroid(HPCD* pcd,double* center) {
	double sumX=0, sumY=0, sumZ=0;
	for (int i=0;i<pcd->maxSize;i++) {
		HPoint* p = pcd->data[i];
		if (p) {
			sumX += p->x;
			sumY += p->y;
			sumZ += p->z;
		}
	}
	center[0] = sumX / pcd->numPoints;
	center[1] = sumY / pcd->numPoints;
	center[2] = sumZ / pcd->numPoints;
}

void eigenvalue(int N,double* A,double* lambda_real,double* lambda_imag,double* v) {

	int info,ldvl=1,ldvr=N,lwork=15*N;	
	double *work = new double[lwork]();
	char jobvl = 'N', jobvr = 'V';
	dgeev_(&jobvl,&jobvr, &N, A, &N, lambda_real, lambda_imag,
	    NULL,&ldvl, v, &ldvr ,work, &lwork, &info);
//	printf("info: %d\n",info);
//	printf("optimal: %f\n",work[0]);
	if (info!=0) {
		printf("Error in subroutine dgeev_ (info=%d)\n",info);
	}
	delete[] work;
}

Box getBoundingBox(HPCD* pointcloud) {
	Box b;
	bool init = true;
	for (int i=0;i<pointcloud->maxSize;i++) {
		HPoint* p = pointcloud->data[i];
		if (p) {
			if (init) {
				b.minX = p->x;
				b.maxX = p->x;
				b.minY = p->y;
				b.maxY = p->y;
				b.minZ = p->z;
				b.maxZ = p->z;
				init = false;
			} else {
				if (p->x < b.minX) b.minX = p->x;
				if (p->x > b.maxX) b.maxX = p->x;
				if (p->y < b.minY) b.minY = p->y;
				if (p->y > b.maxY) b.maxY = p->y;
				if (p->z < b.minZ) b.minZ = p->z;
				if (p->z > b.maxZ) b.maxZ = p->z;
			}
		}
	}
	return b;
}

void getPCA_XY(HPCD* pointcloud) {
	double cov[4] = {}; //column major
	double center[3];
	getCentroid(pointcloud,center);
	for (int j=0;j<pointcloud->maxSize;j++) {
		HPoint *p = pointcloud->data[j];
		if (p) {
			double deltaP[2] = {
				p->x - center[0],
				p->y - center[1],
			};
			cov[0] += deltaP[0] * deltaP[0];
			cov[1] += deltaP[1] * deltaP[0];
			cov[2] += deltaP[0] * deltaP[1];
			cov[3] += deltaP[1] * deltaP[1];
		}
	}
	//compute PCA
	double lambda_imag[2];
	double lambda_real[2];
	double v[4];
	eigenvalue(2,cov,lambda_real,lambda_imag,v);
	//sort by eigenvalue
	if (lambda_real[0] < lambda_real[1]) {
		double tmpL = lambda_real[1];
		lambda_real[1] = lambda_real[0];
		lambda_real[0] = tmpL;
		double tmpV[2];
		memcpy(tmpV,v+2,2*sizeof(double));
		memcpy(v+2,v,2*sizeof(double));
		memcpy(v,tmpV,2*sizeof(double));
	}
	//project to principal axes
	for (int j=0;j<pointcloud->maxSize;j++) {
		HPoint *p = pointcloud->data[j];
		if (p) {
			float dotProductX = (p->x - center[0]) * v[0] + (p->y -center[1]) * v[1];
			float dotProductY = (p->x - center[0]) * v[2] + (p->y - center[1]) * v[3];
			p->x = (int) (dotProductX);
			p->y = (int) (dotProductY);
			p->z -= center[2];
		}
	}
	//fix direction of major axis
	int frontMinX=pointcloud->numGrid;
	int frontMaxX=0;
	int frontMinY=pointcloud->numGrid;
	int frontMaxY=-pointcloud->numGrid;
	int frontMinZ=pointcloud->numGrid;
	int frontMaxZ=-pointcloud->numGrid;
	int backMinX=pointcloud->numGrid;
	int backMaxX=0;
	int backMinY=pointcloud->numGrid;
	int backMaxY=-pointcloud->numGrid;
	int backMinZ=pointcloud->numGrid;
	int backMaxZ=-pointcloud->numGrid;
	int numFront = 0;
	int numBack = 0;
	for (int j=0;j<pointcloud->maxSize;j++) {
		HPoint *p = pointcloud->data[j];
		if (p) {
			if (p->x > 0) {
				numFront++;
				if (abs(p->x) > frontMaxX) frontMaxX = abs(p->x);
				if (abs(p->x) < frontMinX) frontMinX = abs(p->x);
				if (p->y > frontMaxY) frontMaxY = p->y;
				if (p->y < frontMinY) frontMinY = p->y;
				if (p->z > frontMaxZ) frontMaxZ = p->z;
				if (p->z < frontMinZ) frontMinZ = p->z;
			} else {
				numBack++;
				if (abs(p->x) > backMaxX) backMaxX = abs(p->x);
				if (abs(p->x) < backMinX) backMinX = abs(p->x);
				if (p->y > backMaxY) backMaxY = p->y;
				if (p->y < backMinY) backMinY = p->y;
				if (p->z > backMaxZ) backMaxZ = p->z;
				if (p->z < backMinZ) backMinZ = p->z;
			}
		}
	}
	int frontVolume = (frontMaxX-frontMinX) * (frontMaxY-frontMinY) * (frontMaxZ-frontMinZ);
	int backVolume = (backMaxX-backMinX) * (backMaxY-backMinY) * (backMaxZ-backMinZ);
	if (frontVolume > backVolume /*&& numFront > numBack*/) {
		for (int j=0;j<pointcloud->maxSize;j++) {
			HPoint *p = pointcloud->data[j];
			if (p) {
				p->x = -p->x;
			}
		}
	}
	pointcloud->oriented_box = getBoundingBox(pointcloud);
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

inline int HPCD_find(HPCD* cloud,int x,int y,int z) {
	int ikey = getIntKey(x,y,z);
	int j = baseHash(cloud->maxSize,ikey);
	int step = stepHash(cloud->maxSize,ikey);
	for (int k=0;k<cloud->maxSize;k++) {
		HPoint* h = cloud->data[j];
		if (!h) {
			return -1;
		} else if (h->x == x && h->y == y && h->z == z){
			return j;
		} else {
			j += step;
			j %= cloud->maxSize;
		}
	}
	return -1;
}

HPCD* HPCD_Init(char* inFile,int numGrid,Box* box,std::vector<float> *float_data) {
	HPCD* res = new HPCD();
	FILE* f = fopen(inFile,"r");
	if (!f) {
		printf("%s not found\n",inFile);
		return NULL;
	}
	char buf[256];
	std::vector<unsigned char> color_data;
	PCD_data_storage data_storage = NONE;
	int totalPoints = 0;
	while (fgets(buf, 256, f)) {
		if (sscanf(buf, "POINTS %d", &totalPoints) == 1) {
		} else if (strncmp(buf,"FIELDS x y z rgb",16)==0) {
			res->hasColor = true;
		} else if (strncmp(buf,"FIELDS x y z",12)==0) {
			res->hasColor = false;
		} else if (strncmp(buf,"DATA ascii",10)==0) {
			data_storage = ASCII;
			break;
		}
	}
	if (data_storage != ASCII)
		return NULL;
	if (box) {
		res->minX = box->minX;
		res->maxX = box->maxX;
		res->minY = box->minY;
		res->maxY = box->maxY;
		res->minZ = box->minZ;
		res->maxZ = box->maxZ;
		float x,y,z;
		unsigned char r,g,b;
		int rgb;
		if (res->hasColor) {
			for (int i=0;i<totalPoints;i++) {
				if (!fgets(buf,256,f))
					break;
				if (!sscanf(buf, "%f %f %f %d",&x,&y,&z,&rgb) == 6)
					break;
				if (x >= res->minX && x <= res->maxX &&
					y >= res->minY && y <= res->maxY &&
					z >= res->minZ && z <= res->maxZ) {
					float_data->push_back(x);
					float_data->push_back(y);
					float_data->push_back(z);
					r = (rgb >> 16) & 0xFF;
					g = (rgb >> 8) & 0xFF;
					b = rgb & 0xFF;
					color_data.push_back(r);
					color_data.push_back(g);
					color_data.push_back(b);
				}
			}
		} else {
			for (int i=0;i<totalPoints;i++) {
				if (!fgets(buf,256,f))
					break;
				if (!sscanf(buf, "%f %f %f",&x,&y,&z) == 3)
					break;
				if (x >= res->minX && x <= res->maxX &&
					y >= res->minY && y <= res->maxY &&
					z >= res->minZ && z <= res->maxZ) {
					float_data->push_back(x);
					float_data->push_back(y);
					float_data->push_back(z);
				}
			}
		}
	} else {
		float x,y,z;
		unsigned char r,g,b;
		int rgb;
		if (!fgets(buf,256,f))
			return NULL;
		sscanf(buf, "%f %f %f",&x,&y,&z);
		res->minX = res->maxX = x;
		res->minY = res->maxY = y;
		res->minZ = res->maxZ = z;
		if (res->hasColor) {
			for (int i=1;i<totalPoints;i++) {
				if (!fgets(buf,256,f))
					break;
				if (!sscanf(buf, "%f %f %f %d",&x,&y,&z,&rgb) == 6)
					break;
				if (x < res->minX) res->minX = x;
				else if (x > res->maxX) res->maxX = x;
				if (y < res->minY) res->minY = y;
				else if (y > res->maxY) res->maxY = y;
				if (z < res->minZ) res->minZ = z;
				else if (z > res->maxZ) res->maxZ = z;
				float_data->push_back(x);
				float_data->push_back(y);
				float_data->push_back(z);
				r = (rgb >> 16) & 0xFF;
				g = (rgb >> 8) & 0xFF;
				b = rgb & 0xFF;
				color_data.push_back(r);
				color_data.push_back(g);
				color_data.push_back(b);
			}
		} else {
			for (int i=1;i<totalPoints;i++) {
				if (!fgets(buf,256,f))
					break;
				if (!sscanf(buf, "%f %f %f",&x,&y,&z) == 3)
					break;
				if (x < res->minX) res->minX = x;
				else if (x > res->maxX) res->maxX = x;
				if (y < res->minY) res->minY = y;
				else if (y > res->maxY) res->maxY = y;
				if (z < res->minZ) res->minZ = z;
				else if (z > res->maxZ) res->maxZ = z;
				float_data->push_back(x);
				float_data->push_back(y);
				float_data->push_back(z);
			}
		}
	}
	fclose(f);
#if PROFILE
	clock_gettime(CLOCK_MONOTONIC,&toc);
	printf("Profile (Initialization): %f\n",toc.tv_sec - tic.tv_sec + 0.000000001 * toc.tv_nsec - 0.000000001 * tic.tv_nsec);
	tic = toc;
#endif
	totalPoints = float_data->size() / 3;
	float minDist = res->maxX - res->minX;
	if (res->maxY - res->minY < minDist)
		minDist = res->maxY - res->minY;
	if (res->maxZ - res->minZ < minDist)
		minDist = res->maxZ - res->minZ;
	res->numGrid = numGrid;
	res->leafSize = minDist / res->numGrid;
	res->maxSize = 8; 
	while (res->maxSize < 4 * totalPoints)
		res->maxSize *= 2;
	res->data = new HPoint*[res->maxSize]();
	int i,j=0,k,l=0;
	res->numPoints = 0;
	if (res->hasColor){
		for (i=0;i<totalPoints;i++) {
			float x = (*float_data)[j++];
			float y = (*float_data)[j++];
			float z = (*float_data)[j++];
			unsigned char r = color_data[l++];
			unsigned char g = color_data[l++];
			unsigned char b = color_data[l++];
			int xi = (int) ((x-res->minX)/res->leafSize);
			int yi = (int) ((y-res->minY)/res->leafSize);
			int zi = (int) ((z-res->minZ)/res->leafSize);
			int ikey = getIntKey(xi,yi,zi);
			int key = baseHash(res->maxSize,ikey);
			int step = stepHash(res->maxSize,ikey);
			for (k = 0; k < res->maxSize; k++) {
				HPoint* h = res->data[key];
				if (!h) {
					HPoint* p = new HPoint();
					p->x = xi;
					p->y = yi;
					p->z = zi;
					p->r = r;
					p->g = g;
					p->b = b;
					res->data[key] = p;
					res->numPoints++;
					break;
				} else if (h->x == xi && h->y == yi && h->z == zi){
					break;
				} else {
					key += step;
					key %= res->maxSize;
				}
			}
		}
	} else {
		for (i=0;i<totalPoints;i++) {
			float x = (*float_data)[j++];
			float y = (*float_data)[j++];
			float z = (*float_data)[j++];
			int xi = (int) ((x-res->minX)/res->leafSize);
			int yi = (int) ((y-res->minY)/res->leafSize);
			int zi = (int) ((z-res->minZ)/res->leafSize);
			int ikey = getIntKey(xi,yi,zi);
			int key = baseHash(res->maxSize,ikey);
			int step = stepHash(res->maxSize,ikey);
			for (k = 0; k < res->maxSize; k++) {
				HPoint* h = res->data[key];
				if (!h) {
					HPoint* p = new HPoint();
					p->x = xi;
					p->y = yi;
					p->z = zi;
					p->r = p->g = p->b = 0;
					res->data[key] = p;
					res->numPoints++;
					break;
				} else if (h->x == xi && h->y == yi && h->z == zi){
					break;
				} else {
					key += step;
					key %= res->maxSize;
				}
			}
		}
	}
//	printf("Processed point cloud (numPoints:%d maxSize:%d leafSize:%f)\n",res->numPoints,res->maxSize,res->leafSize);
//	printf("Bounding box: x:(%.2f %.2f) y:(%.2f %.2f) z:(%.2f %.2f)\n",res->minX,res->maxX,res->minY,res->maxY,res->minZ,res->maxZ);
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

void HPCD_write(char* filename,HPCD* pointcloud) {
	FILE* f = fopen(filename, "w");
	if (!f) {
		printf("Cannot write to file: %s\n", filename);
		return;
	}
	if (pointcloud->hasColor) {
		fprintf(f,"# .PCD v0.7 - Point Cloud Data file format\n"
		"VERSION 0.7\n"
		"FIELDS x y z rgb\n"
		"SIZE 4 4 4 4\n"
		"TYPE F F F I\n"
		"COUNT 1 1 1 1\n"
		"WIDTH %d\n"
		"HEIGHT 1\n"
		"VIEWPOINT 0 0 0 1 0 0 0\n"
		"POINTS %d\n"
		"DATA ascii\n",pointcloud->numPoints,pointcloud->numPoints);
		for (int i=0;i<pointcloud->maxSize;i++) {
			HPoint *h = pointcloud->data[i];
			if (h) {
				int rgb = (h->r << 16) | (h->g << 8) | h->b;
				fprintf(f,"%f %f %f %d\n",
					pointcloud->minX + h->x * pointcloud->leafSize,
					pointcloud->minY + h->y * pointcloud->leafSize,
					pointcloud->minZ + h->z * pointcloud->leafSize,
					rgb);
			}
		}
	} else {
		fprintf(f,"# .PCD v0.7 - Point Cloud Data file format\n"
		"VERSION 0.7\n"
		"FIELDS x y z\n"
		"SIZE 4 4 4\n"
		"TYPE F F F\n"
		"COUNT 1 1 1\n"
		"WIDTH %d\n"
		"HEIGHT 1\n"
		"VIEWPOINT 0 0 0 1 0 0 0\n"
		"POINTS %d\n"
		"DATA ascii\n",pointcloud->numPoints,pointcloud->numPoints);
		for (int i=0;i<pointcloud->maxSize;i++) {
			HPoint *h = pointcloud->data[i];
			if (h) {
				fprintf(f,"%f %f %f\n",
					pointcloud->minX + h->x * pointcloud->leafSize,
					pointcloud->minY + h->y * pointcloud->leafSize,
					pointcloud->minZ + h->z * pointcloud->leafSize);
			}
		}
	}
	fclose(f);
#if VERBOSE
	printf("Wrote %d points to %s\n",pointcloud->numPoints,filename);
#endif
}

void HPCD_delete(HPCD* cloud) {
	for (int i=0;i<cloud->maxSize;i++)
		if (cloud->data[i])
			delete cloud->data[i];
	delete[] cloud->data;
	delete cloud;
}

void equalize(std::vector<float> *data, int numBins, float* hist, float mx) {
	if (mx < 0) {
		mx = data->at(0);
		for (unsigned int i=1;i<data->size();i++) {
			if (data->at(i) > mx) mx = data->at(i);
		}
	}
	float scale = 1.0 * numBins / mx;
	for (unsigned int i=0;i<data->size();i++) {
		hist[(int)(data->at(i) * scale)] ++;
	}
	for (int i=0;i<numBins;i++) {
		hist[i] /= data->size();
	}
}

std::vector<float> getHistogram(HPCD* cloud) {
	std::vector<float> hist;
	int numBins = 64;
	int numSamples = 20000;
	int* pointdata = new int[cloud->numPoints*3];
	int k=0;
	for (int j=0;j<cloud->maxSize;j++) {
		if (cloud->data[j]) {
			pointdata[k++] = cloud->data[j]->x;
			pointdata[k++] = cloud->data[j]->y;
			pointdata[k++] = cloud->data[j]->z;
		}
	}
	std::vector<float> distance_data;
	std::vector<float> area_data;
	std::vector<float> angle_data;
	for (int n=0;n<numSamples;n++) {
		int id0,id1,id2;
		do {
			id0 = (rand() % cloud->numPoints) * 3;
			id1 = (rand() % cloud->numPoints) * 3;
			id2 = (rand() % cloud->numPoints) * 3;
		} while (id0==id1 || id0==id2 || id1==id2);
		HPoint p0 = {pointdata[id0],pointdata[id0+1],pointdata[id0+2],0,0,0,0};
		HPoint p1 = {pointdata[id1],pointdata[id1+1],pointdata[id1+2],0,0,0,0};
		HPoint p2 = {pointdata[id2],pointdata[id2+1],pointdata[id2+2],0,0,0,0};
		float a = distance(p0,p1);
		float b = distance(p0,p2);
		float c = distance(p1,p2);
		float s = (a + b + c) / 2;
		float area = sqrt(s*(s-a)*(s-b)*(s-c));
		if (isnan(area) || area < 1) { //collinear
			n--; continue;
		}
		distance_data.push_back(a);
		distance_data.push_back(b);
		distance_data.push_back(c);
		angle_data.push_back(angle(a,b,c));
		angle_data.push_back(angle(b,c,a));
		angle_data.push_back(angle(c,a,b));
		area_data.push_back(area);
	}
	float dx = cloud->oriented_box.maxX - cloud->oriented_box.minX;
	float dy = cloud->oriented_box.maxY - cloud->oriented_box.minY;
	float dz = cloud->oriented_box.maxZ - cloud->oriented_box.minZ;
	std::vector<float> dv;
	dv.push_back(dx);
	dv.push_back(dy);
	dv.push_back(dz);
	std::sort(dv.begin(),dv.end());
#if REDUCED_FEATURE
	hist.resize(numBins,0);
#if GLOBAL_NORM
	float global_distance = sqrt(dx*dx + dy*dy + dz*dz);
	equalize(&distance_data,numBins,&hist[0],global_distance);
#else
	equalize(&distance_data,numBins,&hist[0],-1);
#endif
#else
	hist.resize(numBins*3,0);
#if GLOBAL_NORM
	float global_distance = sqrt(dx*dx + dy*dy + dz*dz);
	float global_area = 0.5 * dv[2] * sqrt(dv[0]*dv[0] + dv[1]*dv[1]); 
	equalize(&angle_data,numBins,&hist[0],M_PI);
	equalize(&area_data,numBins,&hist[numBins],global_area);
	equalize(&distance_data,numBins,&hist[2*numBins],global_distance);
#else
	equalize(&angle_data,numBins,&hist[0],-1);
	equalize(&area_data,numBins,&hist[numBins],-1);
	equalize(&distance_data,numBins,&hist[2*numBins],-1);
#endif
#endif
	delete[] pointdata;
	return hist;
}

std::vector<float> getDensity(HPCD* cloud,Dimensions *dim) {
	std::vector<int> density;
	std::vector<float> hist;
	Box b = cloud->oriented_box;
	int gridscale[3];
	float gridLength[3] = {
		b.maxX - b.minX,
		b.maxY - b.minY,
		b.maxZ - b.minZ
	};
	int minIndex = 0;
	float minDimension = (b.maxX - b.minX);
	for (int i=1;i<3;i++) {
		if (gridLength[i] < minDimension) {
			minDimension = gridLength[i];
			minIndex = i;
		}
	}
	for (int i=0;i<3;i++) {
		if (i==minIndex) {
			gridscale[i] = DENSITY_GRID;
		} else {
			gridscale[i] = DENSITY_GRID * (gridLength[i] / minDimension);
		}
	}
	dim->length = gridscale[0];
	dim->width = gridscale[1];
	dim->height = gridscale[2];
	density.resize(dim->length*dim->width*dim->height);
	float scaleX = 1.0 * (gridscale[0] - 1) / gridLength[0];
	float scaleY = 1.0 * (gridscale[1] - 1) / gridLength[1];
	float scaleZ = 1.0 * (gridscale[2] - 1) / gridLength[2];
	for (int i=0;i<cloud->maxSize;i++) {
		HPoint* p = cloud->data[i];
		if (p) {
			int xi = (int) round((p->x - b.minX) * scaleX);
			int yi = (int) round((p->y - b.minY) * scaleY);
			int zi = (int) round((p->z - b.minZ) * scaleZ);
			int pi = xi + yi*dim->length + zi*dim->width*dim->length;
			density[pi]++;
		}
	}
	for (size_t i=0;i<density.size();i++) {
		hist.push_back(1.0 * density[i] / cloud->numPoints);
	}
	return hist;
}

std::vector<float> getAssemblyHistogram(std::vector<HPCD> assembly) {
	std::vector<float> histN;
	for (size_t i=0;i<assembly.size();i++) {
		std::vector<float> hist = getHistogram(&assembly[i]);
		histN.insert(histN.end(),hist.begin(),hist.end());
	}
	return histN;
}

void writeHistogram(char* fileName,std::vector<float> hist) {
	FILE* f = fopen(fileName,"w");
	if (!f)
		return;
	fprintf(f,"# .PCD v0.7 - Point Cloud Data file format\n"
		"VERSION 0.7\n"
		"FIELDS esf\n"
		"SIZE 4\n"
		"TYPE F\n"
		"COUNT %lu\n"
		"WIDTH 1\n"
		"HEIGHT 1\n"
		"VIEWPOINT 0 0 0 1 0 0 0\n"
		"POINTS 1\n"
		"DATA ascii\n",hist.size());
	for (unsigned int i=0;i<hist.size();i++) {
		fprintf(f,"%f ",hist[i]);
	}
	fprintf(f,"\n");
#if VERBOSE
	printf("Wrote %lu-bin histogram to %s\n",hist.size(),fileName);
#endif
	fclose(f);
}

void writeOccupancyGrid(char* fileName,std::vector<float> hist,Dimensions dim) {
	FILE* f = fopen(fileName,"w");
	if (!f)
		return;
	fprintf(f,"%d %d %d\n",dim.length,dim.width,dim.height);
	for (size_t i=0;i<hist.size();i++) {
		fprintf(f,"%f ",hist[i]);
	}
	fprintf(f,"\n");
	fclose(f);
}

std::vector<HPCD> make_part(HPCD* cloud) {
	std::vector<HPCD> assembly;
	std::vector< std::vector<HPoint*> > points;
	points.resize(NUM_PARTS);
	for (int i=0;i<cloud->maxSize;i++) {
		HPoint *h = cloud->data[i];
		if (h) {
			bool fb = h->x >= 0;
			bool lr = h->y >= 0;
			bool ud = h->z >= 0;
			if (points.size() >= 8)
				points[(ud<<2)|(lr<<1)|fb].push_back(h);
			else if (points.size() >= 4)
				points[(ud<<1)|fb].push_back(h);
			else if (points.size() >= 2)
				points[fb].push_back(h);
			if (points.size() % 2 == 1)
				points[points.size() - 1].push_back(h);
		}
	}
	for (unsigned int i=0;i<points.size();i++) {
		HPCD part = *cloud;
		part.maxSize = points[i].size();
		part.numPoints = points[i].size();
		part.data = new HPoint*[part.numPoints];
		memcpy(part.data,points[i].data(),points[i].size()*sizeof(HPoint*));
		assembly.push_back(part);
	}
	return assembly;
}

int main(int argc, char* argv[]) {
	if (argc<2) {
		printf("./esfc 0-cloud.pcd\n");
		return 1;
	}

//	srand(time(NULL));
	srand(0);
	int numGrid = 200;
	std::vector<float> float_data;
	HPCD* cloud = HPCD_Init(argv[1],numGrid,NULL,&float_data);
	if (!cloud) {
		return 1;
	}
//	HPCD_resize(cloud);
	getPCA_XY(cloud);
	char buffer[128];
	std::vector<HPCD> assembly = make_part(cloud);
#if SAVE_ASSEMBLY
	for (size_t i=0;i<assembly.size();i++) {
		sprintf(buffer,"%s-%lu-part.pcd",argv[1],i);
		HPCD_write(buffer,&assembly[i]);
	}
#endif

#if USE_DENSITY
	Dimensions dim;
	std::vector<float> density = getDensity(cloud,&dim);
	sprintf(buffer,"%s.og",argv[1]);
	writeOccupancyGrid(buffer,density,dim);
#else
	sprintf(buffer,"%s-esf.pcd",argv[1]);
	writeHistogram(buffer,getAssemblyHistogram(assembly));
#endif
	HPCD_delete(cloud);
	for (size_t i=0;i<assembly.size();i++) {
		delete[] assembly[i].data;
	}
}
