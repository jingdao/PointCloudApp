#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <math.h>
#include <time.h>
#define BASE_HASH_CONSTANT 0.618033988
#define STEP_HASH_CONSTANT 0.707106781
#define STRING_HASH_CONSTANT 5381
#define PROFILE 1
#define USE_LEAF 0

#if PROFILE
	struct timespec start,tic,toc;
#endif

struct HPoint {
	int x,y,z;
	int label;
};

struct HPCD {
	int numGrid;
	int numPoints;
	int maxSize;
	float leafSize;
	float minX,minY,minZ,maxX,maxY,maxZ;
	HPoint** data;
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

enum PCD_data_storage {
	ASCII,
	BINARY,
	NONE
};

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

HPCD* HPCD_Init(char* inFile,int numGrid,Box* box,std::vector<float> *float_data) {
	HPCD* res = new HPCD;
	FILE* f = fopen(inFile,"r");
	if (!f) {
		printf("%s not found\n",inFile);
		return NULL;
	}
	char buf[256];
	PCD_data_storage data_storage = NONE;
	int totalPoints = 0;
	while (fgets(buf, 256, f)) {
		if (sscanf(buf, "POINTS %d", &totalPoints) == 1) {
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
	} else {
		float x,y,z;
		if (!fgets(buf,256,f))
			return NULL;
		sscanf(buf, "%f %f %f",&x,&y,&z);
		res->minX = res->maxX = x;
		res->minY = res->maxY = y;
		res->minZ = res->maxZ = z;
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

HPCD* HPCD_Init_KITTI(char* inFile,int numGrid,Box* box,std::vector<float> *float_data) {
	HPCD* res = new HPCD;
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
	fclose(f);
	printf("Wrote %d points to %s\n",pointcloud->numPoints,filename);
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

void segmentPlane(HPCD* cloud, int iter,float inlierRatio) {
	int maxInliers = 0;
	int optimumInliers = inlierRatio * cloud->numPoints;
	int distanceThreshold;
	Plane bestPlane = {0,0,1,0};
	int i,k=0;
	int* pointdata = new int[cloud->numPoints*3];
	for (int j=0;j<cloud->maxSize;j++) {
		if (cloud->data[j]) {
			pointdata[k++] = cloud->data[j]->x;
			pointdata[k++] = cloud->data[j]->y;
			pointdata[k++] = cloud->data[j]->z;
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
	for (int j=0;j<cloud->maxSize;j++) {
		HPoint* h = cloud->data[j];
		if ( h && abs( bestPlane.a * h->x +
			 bestPlane.b * h->y +
			 bestPlane.c * h->z +
			 bestPlane.d )
			 < distanceThreshold) {
			delete cloud->data[j];
			cloud->data[j] = NULL;
			cloud->numPoints--;
		}
	}
	printf("RANSAC: %d points %d iters (%d,%d,%d,%d)\n",cloud->numPoints,i,bestPlane.a,bestPlane.b,bestPlane.c,bestPlane.d);
	delete[] pointdata;
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
	segmentPlane(cloud,10000,0.5);
	HPCD_resize(cloud);
	segmentPlane(cloud,10000,0.5);
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
	int numClusters = 50;
	euclideanClustering(cloud,numClusters);
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
