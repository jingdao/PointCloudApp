#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <math.h>
#include <time.h>
#include "hpcd.h"
#include "lapack.h"
#define CLOUD_LEAF_SIZE 0.1
#define CENTROID_LEAF_SIZE 0.1
#define VOTE_THRESHOLD 10
#define RADIUS 1.0

struct PointNormal {
	float x,y,z,nx,ny,nz;
};

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

void getNormalVector(std::vector<Point> *neighbors,PointNormal *pn) {
	Point center = {};
	for (size_t i=0;i<neighbors->size();i++) {
		center.x += neighbors->at(i).x;
		center.y += neighbors->at(i).y;
		center.z += neighbors->at(i).z;
	}
	center.x /= neighbors->size();
	center.y /= neighbors->size();
	center.z /= neighbors->size();
	double cov[9] = {};
	for (size_t j=0;j<neighbors->size();j++) {
		float deltaP[3] = {
			neighbors->at(j).x - center.x,
			neighbors->at(j).y - center.y,
			neighbors->at(j).z - center.z
		};
		cov[0] += deltaP[0] * deltaP[0];
		cov[1] += deltaP[1] * deltaP[0];
		cov[2] += deltaP[2] * deltaP[0];
		cov[3] += deltaP[0] * deltaP[1];
		cov[4] += deltaP[1] * deltaP[1];
		cov[5] += deltaP[2] * deltaP[1];
		cov[6] += deltaP[0] * deltaP[2];
		cov[7] += deltaP[1] * deltaP[2];
		cov[8] += deltaP[2] * deltaP[2];
	}
	//compute PCA
	double lambda_real[3], lambda_imag[3], v[9];
	eigenvalue(3,cov,lambda_real,lambda_imag,v);
	//get normal and curvature
	double minLambda;
	int minIndex;
	for (int j=0;j<3;j++) {
		if (j==0 || lambda_real[j] < minLambda) {
			minIndex = j;
			minLambda = lambda_real[j];
		}
	}
	pn->nx = v[minIndex * 3];
	pn->ny = v[minIndex * 3 + 1];
	pn->nz = v[minIndex * 3 + 2];
	//normal points to origin?
//	if (pn->x * pn->nx + pn->y * pn->ny + pn->z * pn->nz > 0) {
//		pn->nx = -pn->nx;
//		pn->ny = -pn->ny;
//		pn->nz = -pn->nz;
//	}
}

void getPointNormals(HPCD* cloud, std::vector<PointNormal> *norm, int clusterThreshold) {
	int numCombinations = 1;
	for (int i = 0; i < clusterThreshold; i++)
		numCombinations *= 6;
	int minNeighbors=cloud->numPoints;
	int maxNeighbors=0;
	float meanNeighbors=0;
	for (int i=0;i<cloud->maxSize;i++) {
		HPoint* h = cloud->data[i];
		if (!h)
			continue;
		PointNormal pn = {
			cloud->minX + h->x * cloud->leafSize,
			cloud->minY + h->y * cloud->leafSize,
			cloud->minZ + h->z * cloud->leafSize,
			0,0,0};
		std::vector<Point> neighbors;
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
				if (j >= 0) {
					Point p = {
						cloud->minX + cloud->data[j]->x * cloud->leafSize,
						cloud->minY + cloud->data[j]->y * cloud->leafSize,
						cloud->minZ + cloud->data[j]->z * cloud->leafSize,
						0,0,0
					};
					neighbors.push_back(p);
				}
				l /= 6;
			}
		}	
		if (neighbors.size() < minNeighbors) minNeighbors = neighbors.size();
		if (neighbors.size() > maxNeighbors) maxNeighbors = neighbors.size();
		meanNeighbors += neighbors.size();
		getNormalVector(&neighbors,&pn);
		norm->push_back(pn);
	}
	printf("normal neighborhood: min %d max %d mean %.2f\n",minNeighbors,maxNeighbors,meanNeighbors/cloud->numPoints);
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

int main(int argc,char* argv[]) {
	if (argc < 3) {
		printf("Usage: %s in.pcd out.pcd\n",argv[0]);
		return 1;
	}
	
	char* inFile = argv[1];
	char* outFile = argv[2];

	HPCD* cloud = HPCD_Init(inFile,CLOUD_LEAF_SIZE);
	std::vector<PointNormal> norm;
	getPointNormals(cloud,&norm,2);

	std::vector<Point> candidates;
	std::vector<Point> spheres;

	HPCD* centroids = new HPCD;
	centroids->numPoints = 0;
	centroids->maxSize = 8;
	centroids->leafSize = CENTROID_LEAF_SIZE;
	centroids->minX = cloud->minX;
	centroids->minY = cloud->minY;
	centroids->minZ = cloud->minZ;
	centroids->data = new HPoint*[8]();
	float r = RADIUS;
	for (size_t i=0;i<norm.size();i++) {
		float x = norm[i].x + r * norm[i].nx;
		float y = norm[i].y + r * norm[i].ny;
		float z = norm[i].z + r * norm[i].nz;
		int xi = (x - centroids->minX) / centroids->leafSize;
		int yi = (y - centroids->minY) / centroids->leafSize;
		int zi = (z - centroids->minZ) / centroids->leafSize;
		HPoint* h = HPCD_add(centroids,xi,yi,zi);
		if (h) {
			h->label++;
			if (h->label == VOTE_THRESHOLD) {
				Point p = {
					xi * centroids->leafSize + centroids->minX,
					yi * centroids->leafSize + centroids->minY,
					zi * centroids->leafSize + centroids->minZ,
					0,0,0
				};
				candidates.push_back(p);
			}
		}
		x = norm[i].x - r * norm[i].nx;
		y = norm[i].y - r * norm[i].ny;
		z = norm[i].z - r * norm[i].nz;
		xi = (x - centroids->minX) / centroids->leafSize;
		yi = (y - centroids->minY) / centroids->leafSize;
		zi = (z - centroids->minZ) / centroids->leafSize;
		h = HPCD_add(centroids,xi,yi,zi);
		if (h) {
			h->label++;
			if (h->label == VOTE_THRESHOLD) {
				Point p = {
					xi * centroids->leafSize + centroids->minX,
					yi * centroids->leafSize + centroids->minY,
					zi * centroids->leafSize + centroids->minZ,
					0,0,0
				};
				candidates.push_back(p);
			}
		}
	}
	printf("%lu sphere candidates identified\n",candidates.size());

	for (int i=0;i<cloud->maxSize;i++) {
		HPoint* h = cloud->data[i];
		if (h) {
			Point q = {
				h->x * cloud->leafSize + cloud->minX,
				h->y * cloud->leafSize + cloud->minY,
				h->z * cloud->leafSize + cloud->minZ,
				0,0,0
			};
			for (size_t j=0;j<candidates.size();j++) {
				float dx = q.x - candidates[j].x;
				float dy = q.y - candidates[j].y;
				float dz = q.z - candidates[j].z;
				if (fabs(dx*dx + dy*dy + dz*dz - r*r) < CLOUD_LEAF_SIZE) {
					spheres.push_back(q);
					break;
				}
			}
		}
	}	

	HPCD_writePoints(outFile,&spheres);
	HPCD_delete(cloud);
	HPCD_delete(centroids);

}
