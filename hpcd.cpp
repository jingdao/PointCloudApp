#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <math.h>
#include <time.h>
#define BASE_HASH_CONSTANT 0.618033988
#define STEP_HASH_CONSTANT 0.707106781
#define STRING_HASH_CONSTANT 5381
#define PROFILE 0
#define USE_LEAF 0

#if PROFILE
	struct timespec start,tic,toc;
#endif

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

struct Box {
	float minX,minY,minZ,maxX,maxY,maxZ;
};

struct Plane {
	int a,b,c,d;
};

struct Cluster {
	int index,size;
};

struct Coordinate {
	int x,y;
};

enum PCD_data_storage {
	ASCII,
	BINARY,
	NONE
};

void colormap (float f,unsigned char *r,unsigned char *g,unsigned char *b) {
	*r=0;
	*g=0;
	*b=0;
	if (f<=0) {
		*b = 128;
	} else if (f <= 0.25) {
		*g = (unsigned char) f / 0.25 * 255;
		*b = (unsigned char) 128 * (1 - f / 0.25);
	} else if (f <= 0.5) {
		*g = 255;
		*r = (unsigned char) (f - 0.25) / 0.25 * 255;
	} else if (f <= 0.75) {
		*r = 255;
		*g = (unsigned char) 255 + (0.5 - f) / 0.25 * 127;
	} else if (f <= 1) {
		*r = 255;
		*g = (unsigned char) 128 * (1 - f) / 0.25;
	} else {
		*r = 255;
	}
}

int compareCluster(const void* v1,const void* v2) {
	Cluster* c1 = (Cluster*) v1;
	Cluster* c2 = (Cluster*) v2;
	return c2->size - c1->size;
}

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

inline HPoint* HPCD_add(HPCD* cloud,int x,int y,int z) {
	if (x<0 || y<0 || z<0)
		return NULL;
	int xrange = (cloud->maxX - cloud->minX) / cloud->leafSize + 1;
	int yrange = (cloud->maxY - cloud->minY) / cloud->leafSize + 1;
	int zrange = (cloud->maxZ - cloud->minZ) / cloud->leafSize + 1;
	if (x>=xrange || y>= yrange || z >=zrange)
		return NULL;
	int ikey = getIntKey(x,y,z);
	int j = baseHash(cloud->maxSize,ikey);
	int step = stepHash(cloud->maxSize,ikey);
	HPoint* added = NULL;
	for (int k=0;k<cloud->maxSize;k++) {
		HPoint* h = cloud->data[j];
		if (!h) {
			HPoint* p = new HPoint();
			p->x = x;
			p->y = y;
			p->z = z;
			cloud->numPoints++;
			cloud->data[j] = p;
			added = p;
			break;
		} else if (h->x == x && h->y == y && h->z == z){
			break;
		} else {
			j += step;
			j %= cloud->maxSize;
		}
	}
	if (added && cloud->maxSize < cloud->numPoints * 4)
		HPCD_resize(cloud);
	return added;
}

HPCD* HPCD_Init(char* inFile,int numGrid,Box* box,std::vector<float> *float_data) {
	HPCD* res = new HPCD;
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
		if (res->hasColor) {
			for (int i=0;i<totalPoints;i++) {
				if (!fgets(buf,256,f))
					break;
//				if (!sscanf(buf, "%f %f %f %hhu %hhu %hhu",&x,&y,&z,&r,&g,&b) == 6)
//					break;
				int rgb=0;
				if (!sscanf(buf,"%f %f %f %d",&x,&y,&z,&rgb) == 6)
					break;
				r = (rgb >> 16) & 0xFF;
				g = (rgb >> 8) & 0xFF;
				b = rgb & 0xFF;
				if (x >= res->minX && x <= res->maxX &&
					y >= res->minY && y <= res->maxY &&
					z >= res->minZ && z <= res->maxZ) {
					float_data->push_back(x);
					float_data->push_back(y);
					float_data->push_back(z);
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
//				if (!sscanf(buf, "%f %f %f %hhu %hhu %hhu",&x,&y,&z,&r,&g,&b) == 6)
//					break;
				int rgb=0;
				if (!sscanf(buf,"%f %f %f %d",&x,&y,&z,&rgb) == 6)
					break;
				r = (rgb >> 16) & 0xFF;
				g = (rgb >> 8) & 0xFF;
				b = rgb & 0xFF;
				if (x < res->minX) res->minX = x;
				else if (x > res->maxX) res->maxX = x;
				if (y < res->minY) res->minY = y;
				else if (y > res->maxY) res->maxY = y;
				if (z < res->minZ) res->minZ = z;
				else if (z > res->maxZ) res->maxZ = z;
				float_data->push_back(x);
				float_data->push_back(y);
				float_data->push_back(z);
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
					HPoint* p = new HPoint;
					p->x = xi;
					p->y = yi;
					p->z = zi;
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
	printf("Processed point cloud (numPoints:%d maxSize:%d leafSize:%f)\n",res->numPoints,res->maxSize,res->leafSize);
	printf("Bounding box: x:(%.2f %.2f) y:(%.2f %.2f) z:(%.2f %.2f)\n",res->minX,res->maxX,res->minY,res->maxY,res->minZ,res->maxZ);
	return res;
}

HPCD* HPCD_Init_KITTI(char* inFile,int numGrid,Box* box,std::vector<float> *float_data) {
	HPCD* res = new HPCD;
	res->hasColor = false;
	FILE* f = fopen(inFile,"r");
	if (!f) {
		printf("%s not found\n",inFile);
		return NULL;
	}
	fseek (f, 0, SEEK_END);   // non-portable
	int size_bytes=ftell(f);
	int totalPoints = size_bytes / 4 / sizeof(float);
	fseek(f,0,SEEK_SET);
	if (box) {
		res->minX = box->minX;
		res->maxX = box->maxX;
		res->minY = box->minY;
		res->maxY = box->maxY;
		res->minZ = box->minZ;
		res->maxZ = box->maxZ;
		float x,y,z,r;
		for (int i=0;i<totalPoints;i++) {
			if (fread(&x,sizeof(float),1,f) != 1)
				break;
			if (fread(&y,sizeof(float),1,f) != 1)
				break;
			if (fread(&z,sizeof(float),1,f) != 1)
				break;
			if (fread(&r,sizeof(float),1,f) != 1)
				break;
			if (x >= res->minX && x <= res->maxX &&
			    y >= res->minY && y <= res->maxY &&
			    z >= res->minZ && z <= res->maxZ) {
				float_data->push_back(x);
				float_data->push_back(y);
				float_data->push_back(z);
			}
		}
	} else {
		float x,y,z,r;
		if (fread(&x,sizeof(float),1,f) != 1)
			return NULL;
		if (fread(&y,sizeof(float),1,f) != 1)
			return NULL;
		if (fread(&z,sizeof(float),1,f) != 1)
			return NULL;
		if (fread(&r,sizeof(float),1,f) != 1)
			return NULL;
		res->minX = res->maxX = x;
		res->minY = res->maxY = y;
		res->minZ = res->maxZ = z;
		for (int i=1;i<totalPoints;i++) {
			if (fread(&x,sizeof(float),1,f) != 1)
				break;
			if (fread(&y,sizeof(float),1,f) != 1)
				break;
			if (fread(&z,sizeof(float),1,f) != 1)
				break;
			if (fread(&r,sizeof(float),1,f) != 1)
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
	int i,j=0,k;
	res->numPoints = 0;
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
				HPoint* p = new HPoint;
				p->x = xi;
				p->y = yi;
				p->z = zi;
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
	printf("Processed point cloud (numPoints:%d maxSize:%d leafSize:%f)\n",res->numPoints,res->maxSize,res->leafSize);
	printf("Bounding box: x:(%.2f %.2f) y:(%.2f %.2f) z:(%.2f %.2f)\n",res->minX,res->maxX,res->minY,res->maxY,res->minZ,res->maxZ);
	return res;
}

HPCD* HPCD_Init_PTS(char* inFile,int numGrid,Box* box,std::vector<float> *float_data) {
	std::vector<unsigned char> color_data;
	HPCD* res = new HPCD;
	res->hasColor = true;
	FILE* f = fopen(inFile,"r");
	if (!f) {
		printf("%s not found\n",inFile);
		return NULL;
	}
	char buf[256];
	int totalPoints = 0;
	if (fgets(buf,256,f))
		sscanf(buf,"%d",&totalPoints);
	else
		return NULL;
	if (box) {
		res->minX = box->minX;
		res->maxX = box->maxX;
		res->minY = box->minY;
		res->maxY = box->maxY;
		res->minZ = box->minZ;
		res->maxZ = box->maxZ;
		float x,y,z,d;
		unsigned char r,g,b;
		for (int i=0;i<totalPoints;i++) {
			if (!fgets(buf,256,f))
				break;
			if (!sscanf(buf, "%f %f %f %f %hhu %hhu %hhu",&x,&y,&z,&d,&r,&g,&b) == 7)
				break;
			if (x >= res->minX && x <= res->maxX &&
			    y >= res->minY && y <= res->maxY &&
			    z >= res->minZ && z <= res->maxZ) {
				float_data->push_back(x);
				float_data->push_back(y);
				float_data->push_back(z);
				color_data.push_back(r);
				color_data.push_back(g);
				color_data.push_back(b);
			}
		}
	} else {
		float x,y,z,d;
		unsigned char r,g,b;
		if (!fgets(buf,256,f))
			return NULL;
		sscanf(buf, "%f %f %f %f %hhu %hhu %hhu",&x,&y,&z,&d,&r,&g,&b);
		res->minX = res->maxX = x;
		res->minY = res->maxY = y;
		res->minZ = res->maxZ = z;
		for (int i=1;i<totalPoints;i++) {
			if (!fgets(buf,256,f))
				break;
			if (!sscanf(buf, "%f %f %f %f %hhu %hhu %hhu",&x,&y,&z,&d,&r,&g,&b) == 7)
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
			color_data.push_back(r);
			color_data.push_back(g);
			color_data.push_back(b);
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
			} else if (h->x == xi && h->y == yi && h->z == zi){
				break;
			} else {
				key += step;
				key %= res->maxSize;
			}
		}
	}
	printf("Processed point cloud (numPoints:%d maxSize:%d leafSize:%f)\n",res->numPoints,res->maxSize,res->leafSize);
	printf("Bounding box: x:(%.2f %.2f) y:(%.2f %.2f) z:(%.2f %.2f)\n",res->minX,res->maxX,res->minY,res->maxY,res->minZ,res->maxZ);
	return res;
}

void HPCD_delete(HPCD* cloud) {
	for (int i=0;i<cloud->maxSize;i++)
		if (cloud->data[i])
			delete cloud->data[i];
	delete[] cloud->data;
	delete cloud;
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
	printf("Wrote %d points to %s\n",pointcloud->numPoints,filename);
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

void HPCD_writeLeafClusters(char* outDir,HPCD* cloud, std::vector<std::vector<int> > *indices) {
	char buffer[128];
	for (size_t i=0;i<indices->size();i++) {
		std::vector<int> idx = (*indices)[i];
		snprintf(buffer,128,"%s/%lu-cloud.pcd",outDir,i);
		FILE* f = fopen(buffer, "w");
		if (!f) {
			continue;
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
		"DATA ascii\n",idx.size(),idx.size());
		for (size_t i=0;i<idx.size();i++) {
			HPoint *h = cloud->data[idx[i]];
			fprintf(f,"%f %f %f\n",
				cloud->minX + h->x * cloud->leafSize,
				cloud->minY + h->y * cloud->leafSize,
				cloud->minZ + h->z * cloud->leafSize);
		}
		fclose(f);
	}
	printf("Saved %lu point clouds to %s\n",indices->size(),outDir);
}

void HPCD_writeClusters(char* outDir,HPCD* cloud, std::vector<float> *float_data,int numClusters) {
	char buffer[128];
	size_t numPoints = float_data->size()/3;
	std::vector<float> *clusters = new std::vector<float>[numClusters];
	std::vector<float> outliers;
	int j=0;
	for (size_t i=0;i<numPoints;i++) {
		float x = (*float_data)[j++]; 
		float y = (*float_data)[j++]; 
		float z = (*float_data)[j++];
		int xi = (int) ((x-cloud->minX)/cloud->leafSize);
		int yi = (int) ((y-cloud->minY)/cloud->leafSize);
		int zi = (int) ((z-cloud->minZ)/cloud->leafSize);
		int ikey = getIntKey(xi,yi,zi);
		int key = baseHash(cloud->maxSize,ikey);
		int step = stepHash(cloud->maxSize,ikey);
		int label = -1;
		for (int k=0;k<cloud->maxSize;k++) {
			HPoint* h = cloud->data[key];
			if (!h)
				break;
			else if (h->x == xi && h->y == yi && h->z == zi) {
				label = h->label;
				break;
			} else {
				key += step;
				key %= cloud->maxSize;
			}
		}
		if (label >= 0) {
			clusters[label].push_back(x);
			clusters[label].push_back(y);
			clusters[label].push_back(z);
		} else {
			outliers.push_back(x);
			outliers.push_back(y);
			outliers.push_back(z);
		}
	}
	for (int i=0;i<numClusters;i++) {
		size_t clusterSize = clusters[i].size() / 3;
		if (clusterSize == 0)
			break;
		sprintf(buffer,"%s/%d-cloud.pcd",outDir,i);
		FILE* f = fopen(buffer, "w");
		if (!f) {
			printf("Cannot write to file: %s\n", buffer);
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
		"DATA ascii\n",clusterSize,clusterSize);
		for (size_t l=0;l<clusterSize;l++)
			fprintf(f,"%f %f %f\n",clusters[i][l*3],clusters[i][l*3+1],clusters[i][l*3+2]);
		fclose(f);
	}
	size_t clusterSize = outliers.size() / 3;
	sprintf(buffer,"%s/outlier.pcd",outDir);
	FILE* f = fopen(buffer, "w");
	if (!f) {
		printf("Cannot write to file: %s\n", buffer);
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
	"DATA ascii\n",clusterSize,clusterSize);
	for (size_t i=0;i<clusterSize;i++)
		fprintf(f,"%f %f %f\n",outliers[i*3],outliers[i*3+1],outliers[i*3+2]);
	fclose(f);
	delete[] clusters;
	printf("Saved %d point clouds to %s\n",numClusters,outDir);
}

HPCD* HPCD_combine(std::vector<HPCD*> *region) {
	HPCD* res = new HPCD;
	memcpy(res,region->at(0),sizeof(HPCD));
	res->data = new HPoint*[res->maxSize]();
	res->numPoints = 0;
	for (size_t i=0;i<region->size();i++) {
		HPCD* h = region->at(i);
		for (int j=0;j<h->maxSize;j++) {
			if (h->data[j]) {
				res->data[j] = h->data[j];
				res->numPoints++;
			}
		}
	}
	return res;
}

void HPCD_divide(HPCD* cloud, std::vector<HPCD*> *region, int numRegions) {
	for (int i=0;i<numRegions;i++) {
		for (int j=0;j<numRegions;j++) {
			HPCD* h = new HPCD;
			memcpy(h,cloud,sizeof(HPCD));
			h->numPoints = 0;
			h->data = new HPoint*[h->maxSize]();
			region->push_back(h);
		}
	}
	int xrange = (cloud->maxX - cloud->minX) / cloud->leafSize + 1;
	int yrange = (cloud->maxY - cloud->minY) / cloud->leafSize + 1;
	for (int k=0;k<cloud->maxSize;k++) {
		HPoint* p = cloud->data[k];
		if (p) {
			int i = p->x * numRegions / xrange;
			int j = p->y * numRegions / yrange;
			HPCD* h = region->at(i * numRegions + j);
			h->data[k] = p;
			h->numPoints++;
		}
	}
}

void HPCD_dilate(HPCD* cloud) {
	for (int i=0;i<cloud->maxSize;i++) {
		HPoint* h = cloud->data[i];
		if (h) {
			HPoint* p = HPCD_add(cloud,h->x+1,h->y,h->z);
			if (p) {p->r=h->r;p->g=h->g;p->b=h->b;}
			p = HPCD_add(cloud,h->x-1,h->y,h->z);
			if (p) {p->r=h->r;p->g=h->g;p->b=h->b;}
			p = HPCD_add(cloud,h->x,h->y+1,h->z);
			if (p) {p->r=h->r;p->g=h->g;p->b=h->b;}
			p = HPCD_add(cloud,h->x,h->y-1,h->z);
			if (p) {p->r=h->r;p->g=h->g;p->b=h->b;}
		}
	}
}

void segmentPlane(HPCD* cloud, int iter,float inlierRatio) {
	int maxInliers = 0;
	int optimumInliers = inlierRatio * cloud->numPoints;
	int distanceThreshold;
	Plane bestPlane = {0,0,1,0};
	int i,k=0;
	int* pointdata = new int[cloud->numPoints*3];
	for (int j=0;j<cloud->maxSize;j++) {
		HPoint* h = cloud->data[j];
		if (h) {
			pointdata[k++] = h->x;
			pointdata[k++] = h->y;
			pointdata[k++] = h->z;
		}
	}
	for (i=0;i<iter;i++) {
		//Pick 3 points
		int *p0 = pointdata + (rand() % cloud->numPoints * 3);
		int *p1 = pointdata + (rand() % cloud->numPoints * 3);
		int *p2 = pointdata + (rand() % cloud->numPoints * 3);
		int p0x=p0[0],p0y=p0[1],p0z=p0[2];
		int p1x=p1[0],p1y=p1[1],p1z=p1[2];
		int p2x=p2[0],p2y=p2[1],p2z=p2[2];
		Plane currentPlane = {
			(p1y-p0y)*(p2z-p0z) - (p2y-p0y)*(p1z-p0z),
			(p1z-p0z)*(p2x-p0x) - (p2z-p0z)*(p1x-p0x),
			(p1x-p0x)*(p2y-p0y) - (p2x-p0x)*(p1y-p0y),
			0
		};
		currentPlane.d = -(currentPlane.a * p0x + currentPlane.b * p0y + currentPlane.c * p0z);
		if (currentPlane.a == 0 && currentPlane.b == 0 && currentPlane.c ==0 )
			continue; //picked collinear points
		distanceThreshold = sqrt(
			currentPlane.a * currentPlane.a + 
			currentPlane.b * currentPlane.b + 
			currentPlane.c * currentPlane.c);
		int numInliers = 0;
		for (int j=0;j<cloud->numPoints;j++) {
			if ( abs( currentPlane.a * pointdata[j*3] +
				 currentPlane.b * pointdata[j*3+1] +
				 currentPlane.c * pointdata[j*3+2] +
				 currentPlane.d )
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
	distanceThreshold = sqrt(
		bestPlane.a * bestPlane.a + 
		bestPlane.b * bestPlane.b + 
		bestPlane.c * bestPlane.c);
	double RMSE = 0;
	for (int j=0;j<cloud->maxSize;j++) {
		HPoint* h = cloud->data[j];
		if ( h && abs( bestPlane.a * h->x +
			 bestPlane.b * h->y +
			 bestPlane.c * h->z +
			 bestPlane.d )
			 < distanceThreshold) {
			double d = 1.0 * (bestPlane.a * h->x + bestPlane.b * h->y + bestPlane.c * h->z + bestPlane.d) / distanceThreshold;
			RMSE += d * d;
			delete cloud->data[j];
			cloud->data[j] = NULL;
			cloud->numPoints--;
		}
	}
	printf("RANSAC: %d points %d iters (%d,%d,%d,%d) (RMSE %f)\n",cloud->numPoints,i,bestPlane.a,bestPlane.b,bestPlane.c,bestPlane.d,sqrt(RMSE / maxInliers));
	delete[] pointdata;
}

void extractPlane(HPCD* cloud, float segThreshold, int iter,float inlierRatio, std::vector<Point> *points) {
	int maxInliers = 0;
	int optimumInliers = inlierRatio * cloud->numPoints;
	int distanceThreshold;
	Plane bestPlane = {0,0,1,0};
	int i,k=0;
	int* pointdata = new int[cloud->numPoints*3];
	for (int j=0;j<cloud->maxSize;j++) {
		HPoint* h = cloud->data[j];
		if (h) {
			pointdata[k++] = h->x;
			pointdata[k++] = h->y;
			pointdata[k++] = h->z;
		}
	}
	for (i=0;i<iter;i++) {
		//Pick 3 points
		int *p0 = pointdata + (rand() % cloud->numPoints * 3);
		int *p1 = pointdata + (rand() % cloud->numPoints * 3);
		int *p2 = pointdata + (rand() % cloud->numPoints * 3);
		int p0x=p0[0],p0y=p0[1],p0z=p0[2];
		int p1x=p1[0],p1y=p1[1],p1z=p1[2];
		int p2x=p2[0],p2y=p2[1],p2z=p2[2];
		Plane currentPlane = {
			(p1y-p0y)*(p2z-p0z) - (p2y-p0y)*(p1z-p0z),
			(p1z-p0z)*(p2x-p0x) - (p2z-p0z)*(p1x-p0x),
			(p1x-p0x)*(p2y-p0y) - (p2x-p0x)*(p1y-p0y),
			0
		};
		currentPlane.d = -(currentPlane.a * p0x + currentPlane.b * p0y + currentPlane.c * p0z);
		if (currentPlane.a == 0 && currentPlane.b == 0 && currentPlane.c ==0 )
			continue; //picked collinear points
		distanceThreshold = sqrt(
			currentPlane.a * currentPlane.a + 
			currentPlane.b * currentPlane.b + 
			currentPlane.c * currentPlane.c) * segThreshold / cloud->leafSize;
		int numInliers = 0;
		for (int j=0;j<cloud->numPoints;j++) {
			if ( abs( currentPlane.a * pointdata[j*3] +
				 currentPlane.b * pointdata[j*3+1] +
				 currentPlane.c * pointdata[j*3+2] +
				 currentPlane.d )
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
	distanceThreshold = sqrt(
		bestPlane.a * bestPlane.a + 
		bestPlane.b * bestPlane.b + 
		bestPlane.c * bestPlane.c) * segThreshold / cloud->leafSize;
	double RMSE = 0;
	for (int j=0;j<cloud->maxSize;j++) {
		HPoint* h = cloud->data[j];
		if ( h && abs( bestPlane.a * h->x +
			 bestPlane.b * h->y +
			 bestPlane.c * h->z +
			 bestPlane.d )
			 < distanceThreshold) {
			double d = 1.0 * (bestPlane.a * h->x + bestPlane.b * h->y + bestPlane.c * h->z + bestPlane.d) / distanceThreshold;
			RMSE += d * d;
			Point p = {
				cloud->minX + h->x * cloud->leafSize,
				cloud->minY + h->y * cloud->leafSize,
				cloud->minZ + h->z * cloud->leafSize
			};
			points->push_back(p);
			delete cloud->data[j];
			cloud->data[j] = NULL;
			cloud->numPoints--;
		}
	}
	printf("RANSAC: %d points %d iters (%d,%d,%d,%d) (RMSE %f)\n",cloud->numPoints,i,bestPlane.a,bestPlane.b,bestPlane.c,bestPlane.d,sqrt(RMSE / maxInliers));
	delete[] pointdata;
}

void segmentLowest(HPCD* cloud) {
	int xrange = (cloud->maxX - cloud->minX) / cloud->leafSize + 1;
	int yrange = (cloud->maxY - cloud->minY) / cloud->leafSize + 1;
	int* lowest = new int[xrange * yrange];
	for (int i=0;i<xrange*yrange;i++)
		lowest[i] = -1;
	for (int i=0;i<cloud->maxSize;i++) {
		HPoint* h = cloud->data[i];
		if (h) {
			h->label = 0;
			int id = h->x * yrange + h->y;
			if (lowest[id]==-1 || cloud->data[lowest[id]]->z > h->z)
				lowest[id] = i;
		}
	}
	for (int i=0;i<xrange*yrange;i++) {
		if (lowest[i]!=-1) {
			cloud->data[lowest[i]]->label = 1;
		}
	}
	for (int i=0;i<cloud->maxSize;i++) {
		HPoint* h = cloud->data[i];
		if (h && h->label) {
			cloud->data[i] = NULL;
			cloud->numPoints--;
		}
	}
	delete[] lowest;
}

void calculateCurvature(HPCD* cloud,int colorscheme) { //0:none 1:grayscale 2:color
	int min_curvature=-1,max_curvature=-1;
	int r=3;
	for (int i=0;i<cloud->maxSize;i++) {
		HPoint* h = cloud->data[i];
		if (h) {
			int sdx=0,sdy=0,sdz=0;
			for (int dx=-r;dx<=r;dx++) {
				for (int dy=-r;dy<=r;dy++) {
					for (int dz=-r;dz<=r;dz++) {
						if (HPCD_find(cloud,h->x+dx,h->y+dy,h->z+dz) >= 0) {
							sdx += dx;
							sdy += dy;
							sdz += dz;
						}
					}
				}
			}
			h->curvature = abs(sdx) + abs(sdy) + abs(sdz);
			if (min_curvature < 0 || h->curvature < min_curvature)
				min_curvature = h->curvature;
			if (max_curvature < 0 || h->curvature > max_curvature)
				max_curvature = h->curvature;
		}
	}
	if (colorscheme) {
		float factor = 1.0 / (max_curvature - min_curvature);
		for (int i=0;i<cloud->maxSize;i++) {
			HPoint* h = cloud->data[i];
			if (h) {
				if (colorscheme==1) {
					h->r = 255 * factor * (h->curvature-min_curvature);
					h->g = 255 * factor * (h->curvature-min_curvature);
					h->b = 255 * factor * (h->curvature-min_curvature);
				} else
					colormap(factor*(h->curvature-min_curvature) , &h->r, &h->g, &h->b);
			}
		}
	}
}

void euclideanLeafClustering(HPCD* cloud, std::vector<std::vector<int> > *indices,size_t minSize,size_t maxSize,size_t maxClusters) {
	bool* visited = new bool[cloud->maxSize]();
	for (int i=0;i<cloud->maxSize;i++) {
		if (!cloud->data[i] || visited[i]) continue;
		std::vector<int> P,Q;
		Q.push_back(i);
		visited[i] = true;
		while (Q.size() > 0) {
			int p = Q[Q.size()-1];
			Q.pop_back();
			P.push_back(p);
			HPoint* h = cloud->data[p];
//			for (int dx=-distance;dx<=distance;dx++) {
//				for (int dy=-distance;dy<=distance;dy++) {
//					for (int dz=-distance;dz<=distance;dz++) {
//						int j = HPCD_find(cloud,h->x+dx,h->y+dy,h->z+dz);
//						if (j>=0 && ! visited[j]) {
//							Q.push_back(j);
//							visited[j] = true;
//						}
//					} 
//				}
//			}
			int j;
			j = HPCD_find(cloud,h->x-1,h->y,h->z);
			if (j>=0 && ! visited[j]) {
				Q.push_back(j);
				visited[j] = true;
			}
			j = HPCD_find(cloud,h->x+1,h->y,h->z);
			if (j>=0 && ! visited[j]) {
				Q.push_back(j);
				visited[j] = true;
			}
			j = HPCD_find(cloud,h->x,h->y-1,h->z);
			if (j>=0 && ! visited[j]) {
				Q.push_back(j);
				visited[j] = true;
			}
			j = HPCD_find(cloud,h->x,h->y+1,h->z);
			if (j>=0 && ! visited[j]) {
				Q.push_back(j);
				visited[j] = true;
			}
			j = HPCD_find(cloud,h->x,h->y,h->z-1);
			if (j>=0 && ! visited[j]) {
				Q.push_back(j);
				visited[j] = true;
			}
			j = HPCD_find(cloud,h->x,h->y,h->z+1);
			if (j>=0 && ! visited[j]) {
				Q.push_back(j);
				visited[j] = true;
			}
		}
		if (P.size() >= minSize && P.size() <= maxSize)
			indices->push_back(P);
		if (indices->size() >= maxClusters)
			break;
	}
	delete[] visited;
}

void euclideanClustering(HPCD* cloud,int totalClusters) {
	bool* visited = new bool[cloud->maxSize]();
	std::vector<Cluster> clusters;
	int numClusters = 0;
	for (int i=0;i<cloud->maxSize;i++) {
		HPoint* seed = cloud->data[i];
		if (!seed || visited[i]) continue;
		Cluster c = {numClusters,0};
		std::vector<int> Q;
		Q.push_back(i);
		visited[i] = true;
		while (Q.size() > 0) {
			int p = Q[Q.size()-1];
			Q.pop_back();
			HPoint* h = cloud->data[p];
			h->label = numClusters;
			c.size++;
			int j;
			j = HPCD_find(cloud,h->x-1,h->y,h->z);
			if (j>=0 && !visited[j]) {
				Q.push_back(j);
				visited[j] = true;
			}
			j = HPCD_find(cloud,h->x+1,h->y,h->z);
			if (j>=0 && !visited[j]) {
				Q.push_back(j);
				visited[j] = true;
			}
			j = HPCD_find(cloud,h->x,h->y-1,h->z);
			if (j>=0 && !visited[j]) {
				Q.push_back(j);
				visited[j] = true;
			}
			j = HPCD_find(cloud,h->x,h->y+1,h->z);
			if (j>=0 && !visited[j]) {
				Q.push_back(j);
				visited[j] = true;
			}
			j = HPCD_find(cloud,h->x,h->y,h->z-1);
			if (j>=0 && !visited[j]) {
				Q.push_back(j);
				visited[j] = true;
			}
			j = HPCD_find(cloud,h->x,h->y,h->z+1);
			if (j>=0 && !visited[j]) {
				Q.push_back(j);
				visited[j] = true;
			}
		}
//		printf("Cluster %d: %d\n",numClusters,c.size);
		clusters.push_back(c);
		numClusters++;
	}
	Cluster *sortedClusters = new Cluster[numClusters];
	memcpy(sortedClusters,clusters.data(),numClusters*sizeof(Cluster));
	qsort(sortedClusters,numClusters,sizeof(Cluster),compareCluster);
	int *idx = new int[numClusters];
	for (int i=0;i<numClusters;i++)
		idx[i] = -1;
	for (int i=0;i<totalClusters && i<numClusters;i++)
		idx[sortedClusters[i].index] = i;
	for (int i=0;i<cloud->maxSize;i++) {
		HPoint* h = cloud->data[i];
		if (h)
			h->label = idx[h->label];
	}
	delete[] idx;
	delete[] visited;
	delete[] sortedClusters;
}

int euclideanClusteringXY(HPCD* cloud, int totalClusters) {
	int clusterThreshold = 1;
	int minsize = 100;
	std::vector<Cluster> clusters;
	int numClusters = 0;
	int xrange = (cloud->maxX - cloud->minX) / cloud->leafSize + 1;
	int yrange = (cloud->maxY - cloud->minY) / cloud->leafSize + 1;
	bool* valid = new bool[xrange * yrange]();
	int* grid = new int[xrange * yrange];
	for (int i = 0; i<xrange*yrange; i++)
		grid[i] = -1;
	for (int i = 0; i < cloud->maxSize; i++) {
		HPoint* h = cloud->data[i];
		if (h)
			valid[h->x*yrange + h->y] = true;
	}
	Coordinate p;
	for (p.x = 0; p.x < xrange; p.x++) {
		for (p.y = 0; p.y < yrange; p.y++) {
			int id = p.x * yrange + p.y;
			if (!valid[id] || grid[id]>=0)
				continue;
			std::vector<Coordinate> Q;
			Q.push_back(p);
			Cluster c = { numClusters, 0 };
			while (Q.size() > 0) {
				Coordinate q = Q[Q.size() - 1];
				Q.pop_back();
				int id = q.x * yrange + q.y;
				if (!valid[id] || grid[id]>=0)
					continue;
				c.size++;
				grid[id] = numClusters;
				for (int xi = q.x - clusterThreshold; xi <= q.x + clusterThreshold; xi++) {
					for (int yi = q.y - clusterThreshold; yi <= q.y + clusterThreshold; yi++) {
						if (xi >= 0 && yi >= 0 && xi < xrange && yi < yrange) {
							Coordinate r = { xi, yi };
							Q.push_back(r);
						}
					}
				}
			}
			clusters.push_back(c);
			numClusters++;
		}
	}
	for (int i = 0; i<cloud->maxSize; i++) {
		HPoint* h = cloud->data[i];
		if (h)
			h->label = grid[h->x * yrange + h->y];
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
		if (h)
			h->label = idx[h->label];
	}
	delete[] valid;
	delete[] idx;
	delete[] sortedClusters;
	return k;
}

void rgb2lab(unsigned char R,unsigned char G,unsigned char B,float *_L, float *_a, float *_b) {
	float r = R / 255.0;
	float g = G / 255.0;
	float b = B / 255.0;
	r = (r > 0.04045 ? pow((r + 0.055) / 1.055, 2.4) : r / 12.92) * 100.0;
	g = (g > 0.04045 ? pow((g + 0.055) / 1.055, 2.4) : g / 12.92) * 100.0;
	b = (b > 0.04045 ? pow((b + 0.055) / 1.055, 2.4) : b / 12.92) * 100.0;
	float X = r * 0.4124 + g * 0.3576 + b * 0.1805;
	float Y = r * 0.2126 + g * 0.7152 + b * 0.0722;
	float Z = r * 0.0193 + g * 0.1192 + b * 0.9505;

	float var_X = X / 95.047;
	float var_Y = Y / 100;
	float var_Z = Z / 108.883;

	if ( var_X > 0.008856 ) var_X = pow(var_X,( 1.0/3 ));
	else var_X = (903.3*var_X + 16) / 116;
	if ( var_Y > 0.008856 ) var_Y = pow(var_Y,( 1.0/3 ));
	else var_Y = (903.3*var_Y + 16) / 116;
	if ( var_Z > 0.008856 ) var_Z = pow(var_Z,( 1.0/3 ));
	else var_Z = (903.3*var_Z + 16) / 116;

	*_L = ( 116 * var_Y ) - 16;
	if (*_L < 0)
		*_L = 0;
	*_a = 500 * ( var_X - var_Y );
	*_b = 200 * ( var_Y - var_Z );
}

bool validNeighbor(HPoint* p1, HPoint* p2) {
	float colorThreshold = 35000;
	int curvatureThreshold = 200;
	float colorDiff = 0;
	colorDiff += (p1->r - p2->r) * (p1->r - p2->r);
	colorDiff += (p1->g - p2->g) * (p1->g - p2->g);
	colorDiff += (p1->b - p2->b) * (p1->b - p2->b);
//	float L1,a1,b1,L2,a2,b2;
//	rgb2lab(p1->r,p1->g,p1->b,&L1,&a1,&b1);
//	rgb2lab(p2->r,p2->g,p2->b,&L2,&a2,&b2);
//	colorDiff += (L1 - L2) * (L1 - L2);
//	colorDiff += (a1 - a2) * (a1 - a2);
//	colorDiff += (b1 - b2) * (b1 - b2);
	if (colorDiff > colorThreshold)
		return false;
//	if (p2->curvature > curvatureThreshold)
//		return false;
	if (abs(p1->curvature - p2->curvature) > curvatureThreshold)
		return false;
	return true;
}

void HPCD_convertColor(HPCD* cloud) {
	std::vector<float> listL;
	std::vector<float> listA;
	std::vector<float> listB;
	float minL=9999,minA=9999,minB=9999;
	float maxL=-1,maxA=-1,maxB=-1;
	for (int i=0;i<cloud->maxSize;i++) {
		HPoint* p = cloud->data[i];
		if (p) {
			float l,a,b;
			rgb2lab(p->r,p->g,p->b,&l,&a,&b);
			listL.push_back(l);
			listA.push_back(a);
			listB.push_back(b);
			if (l > maxL) maxL = l;
			if (l < minL) minL = l;
			if (a > maxA) maxA = a;
			if (a < minA) minA = a;
			if (b > maxB) maxB = b;
			if (b < minB) minB = b;
		}
	}
	int j=0;
	for (int i=0;i<cloud->maxSize;i++) {
		HPoint* p = cloud->data[i];
		if (p) {
			p->r = (unsigned char) ((listL[j] - minL) / (maxL - minL) * 255);
			p->g = (unsigned char) ((listA[j] - minA) / (maxA - minA) * 255);
			p->b = (unsigned char) ((listB[j] - minB) / (maxB - minB) * 255);
			j++;
		}
	}
}

void colorClustering(HPCD* cloud,int totalClusters) {
	bool* visited = new bool[cloud->maxSize]();
	std::vector<Cluster> clusters;
	int numClusters = 0;
	for (int i=0;i<cloud->maxSize;i++) {
		HPoint* seed = cloud->data[i];
		if (!seed || visited[i]) continue;
		Cluster c = {numClusters,0};
		std::vector<int> Q;
		Q.push_back(i);
		visited[i] = true;
		while (Q.size() > 0) {
			int p = Q[Q.size()-1];
			Q.pop_back();
			HPoint* h = cloud->data[p];
			h->label = numClusters;
			c.size++;
			if (h->curvature > 200)
				continue;
			int j;
			j = HPCD_find(cloud,h->x-1,h->y,h->z);
			if (j>=0 && !visited[j] && validNeighbor(h,cloud->data[j])) {
				Q.push_back(j);
				visited[j] = true;
			}
			j = HPCD_find(cloud,h->x+1,h->y,h->z);
			if (j>=0 && !visited[j] && validNeighbor(h,cloud->data[j])) {
				Q.push_back(j);
				visited[j] = true;
			}
			j = HPCD_find(cloud,h->x,h->y-1,h->z);
			if (j>=0 && !visited[j] && validNeighbor(h,cloud->data[j])) {
				Q.push_back(j);
				visited[j] = true;
			}
			j = HPCD_find(cloud,h->x,h->y+1,h->z);
			if (j>=0 && !visited[j] && validNeighbor(h,cloud->data[j])) {
				Q.push_back(j);
				visited[j] = true;
			}
			j = HPCD_find(cloud,h->x,h->y,h->z-1);
			if (j>=0 && !visited[j] && validNeighbor(h,cloud->data[j])) {
				Q.push_back(j);
				visited[j] = true;
			}
			j = HPCD_find(cloud,h->x,h->y,h->z+1);
			if (j>=0 && !visited[j] && validNeighbor(h,cloud->data[j])) {
				Q.push_back(j);
				visited[j] = true;
			}
		}
//		printf("Cluster %d: %d\n",numClusters,c.size);
		clusters.push_back(c);
		numClusters++;
	}
	Cluster *sortedClusters = new Cluster[numClusters];
	memcpy(sortedClusters,clusters.data(),numClusters*sizeof(Cluster));
	qsort(sortedClusters,numClusters,sizeof(Cluster),compareCluster);
	int *idx = new int[numClusters];
	for (int i=0;i<numClusters;i++)
		idx[i] = -1;
	for (int i=0;i<totalClusters && i<numClusters;i++)
		idx[sortedClusters[i].index] = i;
	for (int i=0;i<cloud->maxSize;i++) {
		HPoint* h = cloud->data[i];
		if (h)
			h->label = idx[h->label];
	}
	delete[] idx;
	delete[] visited;
	delete[] sortedClusters;
}

void writeParam(char* filename,HPCD* cloud,float segment,float cluster,int mincluster,int maxcluster) {
	FILE* f = fopen(filename,"w");
	if (!f)
		return;
	fprintf(f,"leaf: %f %f %f\n",cloud->leafSize,cloud->leafSize,cloud->leafSize);
	fprintf(f,"segmentation: %f\n",segment);
	fprintf(f,"clustering: %d-%d %f\n",mincluster,maxcluster,cluster);
	fprintf(f,"x range: %f %f\n",cloud->minX,cloud->maxX);
	fprintf(f,"y range: %f %f\n",cloud->minY,cloud->maxY);
	fprintf(f,"z range: %f %f\n",cloud->minZ,cloud->maxZ);
	fclose(f);
}

int main(int argc,char* argv[]) {
	if (argc < 3) {
		printf("Usage: ./hpcd in.pcd out/ [-n numGrid] [-b x1 y1 z1 x2 y2 z2]\n");
		return 1;
	}
	
	char* inFile = argv[1];
	char* outDir = argv[2];
//	Box box = {140,-224,164,240,-137,178};
//	Box box = {-513,-212,167,-452,-147,174};
	Box *box = NULL;
	int numGrid = 100;
	std::vector<float> float_data;
//	srand(time(NULL));
	srand(0);
	for (int i=3;i<argc;i++) {
		if (strcmp(argv[i],"-n")==0 && i+1<argc)
			numGrid = atoi(argv[++i]);
		else if (strcmp(argv[i],"-b")==0 && i+6<argc) {
			box = new Box;
			box->minX = strtod(argv[++i],NULL);
			box->minY = strtod(argv[++i],NULL);
			box->minZ = strtod(argv[++i],NULL);
			box->maxX = strtod(argv[++i],NULL);
			box->maxY = strtod(argv[++i],NULL);
			box->maxZ = strtod(argv[++i],NULL);
		}
	}

#if PROFILE
	clock_gettime(CLOCK_MONOTONIC,&start);
	tic = start;
#endif
	HPCD* cloud = NULL;
	if (strcmp(inFile+strlen(inFile)-4,".pcd")==0)
		cloud = HPCD_Init(inFile,numGrid,box,&float_data);
	else if (strcmp(inFile+strlen(inFile)-4,".bin")==0)
		cloud = HPCD_Init_KITTI(inFile,numGrid,box,&float_data);
	else if (strcmp(inFile+strlen(inFile)-4,".pts")==0)
		cloud = HPCD_Init_PTS(inFile,numGrid,box,&float_data);
	if (!cloud || cloud->numPoints <= 0)
		return 1;

	char buf[128];
	sprintf(buf,"%s/filtered.pcd",outDir);
	HPCD_resize(cloud);
#if PROFILE
	clock_gettime(CLOCK_MONOTONIC,&toc);
	printf("Profile (Downsampling): %f\n",toc.tv_sec - tic.tv_sec + 0.000000001 * toc.tv_nsec - 0.000000001 * tic.tv_nsec);
	tic = toc;
#else
	HPCD_write(buf,cloud);
	sprintf(buf,"%s/param.txt",outDir);
	writeParam(buf,cloud,cloud->leafSize,cloud->leafSize,cloud->numPoints/1000,cloud->numPoints/2);
#endif
	segmentLowest(cloud);
	segmentLowest(cloud);
	HPCD_resize(cloud);
#if PROFILE
	clock_gettime(CLOCK_MONOTONIC,&toc);
	printf("Profile (Ground segmentation): %f\n",toc.tv_sec - tic.tv_sec + 0.000000001 * toc.tv_nsec - 0.000000001 * tic.tv_nsec);
	tic = toc;
#else
	sprintf(buf,"%s/original.pcd",outDir);
	HPCD_write(buf,cloud);
#endif
#if USE_LEAF
	std::vector<std::vector<int> > indices;
	euclideanLeafClustering(cloud,&indices,cloud->numPoints/1000,cloud->numPoints/2,200);
#else
	int numClusters = 200;
//	euclideanClustering(cloud,numClusters);
	euclideanClusteringXY(cloud,numClusters);
//	colorClustering(cloud,numClusters);
#endif
#if PROFILE
	clock_gettime(CLOCK_MONOTONIC,&toc);
	printf("Profile (Clustering): %f\n",toc.tv_sec - tic.tv_sec + 0.000000001 * toc.tv_nsec - 0.000000001 * tic.tv_nsec);
	printf("Profile (Total): %f\n",toc.tv_sec - start.tv_sec + 0.000000001 * toc.tv_nsec - 0.000000001 * start.tv_nsec);
#endif
#if USE_LEAF
	HPCD_writeLeafClusters(outDir,cloud,&indices);
#else
	HPCD_writeClusters(outDir,cloud,&float_data,numClusters);
#endif
	HPCD_delete(cloud);
	if (box)
		delete box;

	return 0;
}
