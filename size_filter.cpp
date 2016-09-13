#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <math.h>
#include <unistd.h>
#include <dirent.h>
#include "lapack.h"
#define MARGIN 0.3

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

void getPCA_XY(PCD* pointcloud,double* lambda_real,double* v) {
	double cov[4] = {}; //column major
	Quaternion center = getCentroid(pointcloud);
	for (int j=0;j<pointcloud->numPoints;j++) {
		float deltaP[2] = {
			pointcloud->float_data[j * 4] - center.i,
			pointcloud->float_data[j * 4 + 1] - center.j,
		};
		cov[0] += deltaP[0] * deltaP[0];
		cov[1] += deltaP[1] * deltaP[0];
		cov[2] += deltaP[0] * deltaP[1];
		cov[3] += deltaP[1] * deltaP[1];
	}
	//compute PCA
	double lambda_imag[2];
	eigenvalue(2,cov,lambda_real,lambda_imag,v);
	//complete 3D structure
	lambda_real[2] = -1;
	v[4] = v[3];
	v[3] = v[2];
	v[2] = v[5] = v[6] = v[7] = 0;
	v[8] = 1;
}

void setPCA(PCD* pointcloud,double* lambda_real, double* v,double* principalLengths, double* principalAxes, double* bbCenter) {
	//sort by eigenvalue
	double tmpL,tmpV[3];
	if (lambda_real[0] < lambda_real[1]) {
		tmpL = lambda_real[1];
		lambda_real[1] = lambda_real[0];
		lambda_real[0] = tmpL;
		memcpy(tmpV, v+3, 3 * sizeof(double));
		memcpy(v+3, v, 3 * sizeof(double));
		memcpy(v, tmpV, 3 * sizeof(double));
	}
	if (lambda_real[1] < lambda_real[2]) {
		tmpL = lambda_real[2];
		lambda_real[2] = lambda_real[1];
		lambda_real[1] = tmpL;
		memcpy(tmpV, v+6, 3 * sizeof(double));
		memcpy(v+6, v+3, 3 * sizeof(double));
		memcpy(v+3, tmpV, 3 * sizeof(double));
	}
	if (lambda_real[0] < lambda_real[2]) {
		tmpL = lambda_real[2];
		lambda_real[2] = lambda_real[0];
		lambda_real[0] = tmpL;
		memcpy(tmpV, v+6, 3 * sizeof(double));
		memcpy(v+6, v, 3 * sizeof(double));
		memcpy(v, tmpV, 3 * sizeof(double));
	}
//	printf("lambda: %f %f %f\n",lambda_real[0],lambda_real[1],lambda_real[2]);
//	printf("v1: %f %f %f\n",v[0],v[1],v[2]);
//	printf("v2: %f %f %f\n",v[3],v[4],v[5]);
//	printf("v3: %f %f %f\n",v[6],v[7],v[8]);
	//compute projections
	float minScale[3], maxScale[3];
	for (int j=0;j<pointcloud->numPoints;j++) {
		for (int i=0;i<3;i++) {
			float dotProduct =
				pointcloud->float_data[j*4] * v[i*3] +
				pointcloud->float_data[j*4+1] * v[i*3+1] +
				pointcloud->float_data[j*4+2] * v[i*3+2];
			if (j==0 || dotProduct < minScale[i])
				minScale[i] = dotProduct;
			if (j==0 || dotProduct > maxScale[i])
				maxScale[i] = dotProduct;				
		}
	}
	float cx=0,cy=0,cz=0;
	for (int i=0;i<3;i++) {
		cx+=(minScale[i]+maxScale[i])/2 * v[i*3];
		cy+=(minScale[i]+maxScale[i])/2 * v[i*3+1];
		cz+=(minScale[i]+maxScale[i])/2 * v[i*3+2];
	}
	//set member variables
	bbCenter[0] = cx;
	bbCenter[1] = cy;
	bbCenter[2] = cz;
	for (int i=0;i<3;i++) {
		principalLengths[i] = maxScale[i] - minScale[i];
		for (int j=0;j<3;j++)
			principalAxes[i*3 + j] = v[i * 3 + j];
	}
}

void addPointsToVector(PCD* pointcloud,std::vector<float> *v) {
	for (int i=0;i<pointcloud->numPoints;i++) {
		v->push_back(pointcloud->float_data[i*4]);
		v->push_back(pointcloud->float_data[i*4+1]);
		v->push_back(pointcloud->float_data[i*4+2]);
	}
}

float getVolume(PCD* pointcloud) {
	float minX,maxX,minY,maxY,minZ,maxZ;
	for (int i=0;i<pointcloud->numPoints;i++) {
		float X = pointcloud->float_data[i*4];
		float Y = pointcloud->float_data[i*4+1];
		float Z = pointcloud->float_data[i*4+2];
		if (i==0 || X < minX)
			minX = X;
		if (i==0 || X > maxX)
			maxX = X;
		if (i==0 || Y < minY)
			minY = Y;
		if (i==0 || Y > maxY)
			maxY = Y;
		if (i==0 || Z < minZ)
			minZ = Z;
		if (i==0 || Z > maxZ)
			maxZ = Z;
	}
	float leafsize = (maxZ - minZ) / 20;
	int xrange = (maxX-minX) / leafsize + 1;
	int yrange = (maxY-minY) / leafsize + 1;
	int zrange = (maxZ-minZ) / leafsize + 1;
	float volume = 0;
	bool* occupied = new bool[xrange*yrange*zrange]();
	for (int i=0;i<pointcloud->numPoints;i++) {
		float X = pointcloud->float_data[i*4];
		float Y = pointcloud->float_data[i*4+1];
		float Z = pointcloud->float_data[i*4+2];
		int xi = (X - minX) / leafsize;
		int yi = (Y - minY) / leafsize;
		int zi = (Z - minZ) / leafsize;
		int index = (zi*yrange + yi) * xrange + xi;
		if (!occupied[index]) {
			occupied[index] = true;
			volume += 1;
		}
	}
//	for (int j=0;j<yrange;j++) {
//		for (int k=0;k<zrange;k++) {
//			int i=0;
//			int index = (k*yrange+j) * xrange + i;
//			volume += xrange;
//			while (i<xrange && !occupied[index]) {
//				volume--;
//				i++;
//				index++;
//			}
//			i=xrange-1;
//			index = (k*yrange+j) * xrange + i;
//			while (i>=0 && !occupied[index]) {
//				volume--;
//				i--;
//				index--;
//			}
//		}
//	}
	delete[] occupied;
	return volume / (xrange * yrange * zrange);
}

float getCenterOfMassHeight(PCD* pointcloud) {
	float h = 0;
	for (int i=0;i<pointcloud->numPoints;i++) {
		h += pointcloud->float_data[i*4+2];
	}
	return h / pointcloud->numPoints;
}

int main(int argc,char* argv[]) {
	if (argc < 4) {
		printf("./size_filter trainFolder testFolder outputFolder\n");
		return 1;
	}

	DIR *dir = opendir(argv[1]);
	if (!dir) {
		printf("Cannot open %s\n",argv[1]);
		return 1;
	}
	struct dirent *dp;
	char * file_name;
	char buffer[128];
	char command[512];
	float minX=-1,maxX=-1;
	float minY=-1,maxY=-1;
	float minZ=-1,maxZ=-1;
	float minH=-1,maxH=-1;
	int num_train = 0, num_test=0, inliers=0, outliers=0;
	while ((dp=readdir(dir)) != NULL) {
		char* file_name = dp->d_name;            
		int l = strlen(file_name);
		if (strncmp(file_name+l-10,"-cloud.pcd",10)==0) {
			sprintf(buffer,"%s/%s",argv[1],file_name);
			PCD* p = NewPCD(buffer);
			if (p) {
				double lambda[3], v[9];
				double principalLengths[3], principalAxes[9], bbCenter[3];
				getPCA_XY(p,lambda,v);
				setPCA(p,lambda,v,principalLengths,principalAxes,bbCenter);
				float comHeight = getCenterOfMassHeight(p);
//				printf("Loaded %13s %6.2f %6.2f %6.2f\n",
//					file_name,principalLengths[0],principalLengths[1],principalLengths[2]);
				if (principalLengths[2] > 50) {
					principalLengths[0] *= 0.0254;
					principalLengths[1] *= 0.0254;
					principalLengths[2] *= 0.0254;
					comHeight *= 0.0254;
				}
				if (minX < 0 || principalLengths[0] < minX)
					minX = principalLengths[0];
				if (maxX < 0 || principalLengths[0] > maxX)
					maxX = principalLengths[0];
				if (minY < 0 || principalLengths[1] < minY)
					minY = principalLengths[1];
				if (maxY < 0 || principalLengths[1] > maxY)
					maxY = principalLengths[1];
				if (minZ < 0 || principalLengths[2] < minZ)
					minZ = principalLengths[2];
				if (maxZ < 0 || principalLengths[2] > maxZ)
					maxZ = principalLengths[2];
				if (minH < 0 || comHeight < minH)
					minH = comHeight;
				if (maxH < 0 || comHeight > maxH)
					maxH = comHeight;
				num_train++;
				delete[] p->float_data;
				delete p;
			}
		}
	}
	minX *= (1 - MARGIN);
	minY *= (1 - MARGIN);
	minZ *= (1 - MARGIN);
	minH *= (1 - MARGIN);
	maxX *= (1 + MARGIN);
	maxY *= (1 + MARGIN);
	maxZ *= (1 + MARGIN);
	maxH *= (1 + MARGIN);
	printf("Trained %d samples (%5.2f-%5.2f) (%5.2f-%5.2f) (%5.2f-%5.2f) (%5.2f,%5.2f)\n",num_train,minX,maxX,minY,maxY,minZ,maxZ,minH,maxH);

	std::vector<float> outlier_points;
	sprintf(buffer,"%s/outlier.pcd",argv[2]);
	PCD* background = NewPCD(buffer);
	if (background)
		addPointsToVector(background,&outlier_points);
	dir = opendir(argv[2]);
	if (!dir) {
		printf("Cannot open %s\n",argv[2]);
		return 1;
	}
	float baseHeight=-1;
	while ((dp=readdir(dir)) != NULL) {
		char* file_name = dp->d_name;            
		int l = strlen(file_name);
		if (strncmp(file_name+l-10,"-cloud.pcd",10)==0) {
			sprintf(buffer,"%s/%s",argv[2],file_name);
			PCD* p = NewPCD(buffer);
			for (int i=0;i<p->numPoints;i++) {
				float z = p->float_data[i*4+2];
				if (baseHeight < 0 || z < baseHeight)
					baseHeight = z;
			}
			delete[] p->float_data;
			delete p;
		}
	}
	dir = opendir(argv[2]);
	while ((dp=readdir(dir)) != NULL) {
		char* file_name = dp->d_name;            
		int l = strlen(file_name);
		if (strncmp(file_name+l-10,"-cloud.pcd",10)==0) {
			sprintf(buffer,"%s/%s",argv[2],file_name);
			PCD* p = NewPCD(buffer);
			if (p) {
				double lambda[3], v[9];
				double principalLengths[3], principalAxes[9], bbCenter[3];
				getPCA_XY(p,lambda,v);
				setPCA(p,lambda,v,principalLengths,principalAxes,bbCenter);
				float comHeight = getCenterOfMassHeight(p) - baseHeight;
//				printf("Loaded %13s %6.2f %6.2f %6.2f\n",file_name,principalLengths[0],principalLengths[1],principalLengths[2]);
				if (principalLengths[2] > 50) {
					principalLengths[0] *= 0.0254;
					principalLengths[1] *= 0.0254;
					principalLengths[2] *= 0.0254;
					comHeight *= 0.0254;
				}
				num_test++;
				if (principalLengths[0]>minX && principalLengths[0]<maxX &&
					principalLengths[1]>minY && principalLengths[1]<maxY &&
					principalLengths[2]>minZ && principalLengths[2]<maxZ
					/*&& comHeight>minH && comHeight<maxH*/) {
					sprintf(command,"cp %s %s/%d-cloud.pcd",buffer,argv[3],inliers);
					printf("%13s %5.2f %5.2f %5.2f %5.2f \033[32m/\033[0m\n",file_name,
						principalLengths[0],principalLengths[1],principalLengths[2], comHeight);
					system(command);
					inliers++;
				} else {
					printf("%13s %5.2f %5.2f %5.2f %5.2f \033[31mx\033[0m\n",file_name,
						principalLengths[0],principalLengths[1],principalLengths[2], comHeight);
					outliers++;
					if (background)
						addPointsToVector(p,&outlier_points);
				}
				delete[] p->float_data;
				delete p;
			}
		}
	}
	printf("Tested %d samples (%d inliers) (%d outliers)\n",num_test,inliers,outliers);
	sprintf(buffer,"%s/prediction.txt",argv[3]);
	FILE* labels = fopen(buffer,"w");
	for (int i=0;i<inliers;i++)
		fwrite("0\n",1,2,labels);
	fclose(labels);
	if (background) {
		int outlier_numPoints = outlier_points.size()/3;
		sprintf(buffer,"%s/outlier.pcd",argv[3]);
		FILE* f = fopen(buffer,"w");
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
		"DATA ascii\n",outlier_numPoints,outlier_numPoints);
		for (int i=0;i<outlier_numPoints;i++)
			fprintf(f,"%f %f %f\n",outlier_points[i*3],outlier_points[i*3+1],outlier_points[i*3+2]);
		fclose(f);
	}
}

