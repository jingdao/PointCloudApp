#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <math.h>
#include <time.h>
#include "hpcd.h"
#define CLOUD_LEAF_SIZE 0.1
#define FEATURE_SIZE 64

float maxScore = 0;
float corr_threshold = 0.8;
float merge_threshold = 0.5;

void normalize(float* f, int k) {
	float s = 0;
	for (int i = 0; i < k; i++)
		s += f[i] * f[i];
	s = sqrt(s);
	for (int i = 0; i < k; i++)
		f[i] /= s;
}

void lbpHist(HPCD* cloud, float* feature) {
	for (int i=0;i<cloud->maxSize;i++) {
		HPoint* h = cloud->data[i];
		if (h) {
			int id=0;
			if (HPCD_find(cloud,h->x-1,h->y,h->z) >= 0) id+=1;
			if (HPCD_find(cloud,h->x,h->y-1,h->z) >= 0) id+=2;
			if (HPCD_find(cloud,h->x-1,h->y,h->z-1) >= 0) id+=4;
			if (HPCD_find(cloud,h->x+1,h->y,h->z) >= 0) id+=8;
			if (HPCD_find(cloud,h->x,h->y+1,h->z) >= 0) id+=16;
			if (HPCD_find(cloud,h->x,h->y,h->z+1) >= 0) id+=32;
			feature[id] += 1;
		}
	}
	normalize(feature,FEATURE_SIZE);
}

void readPoints(char* filename,std::vector<Point> *points) {
	bool header = true;
	FILE* f = fopen(filename,"r");
	if (!f) {
		printf("Cannot read from file: %s\n", filename);
		return;
	}
	char buf[256];
	while (fgets(buf, 256, f)) {
		if (header) {
			Point p;
			if (sscanf(buf,"%f %f %f",&p.x,&p.y,&p.z)==3)
				points->push_back(p);
		} else if (strncmp(buf,"DATA ascii",10)==0) {
			header = false;
		} 
	}
	fclose(f);
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

void saveLBP(char* filename, HPCD* cloud, float* feature) {
	FILE* f = fopen(filename,"w");
	for (int i=0;i<FEATURE_SIZE;i++)
		fprintf(f,"%f ",feature[i]);
	fprintf(f,"\n%f %f %f\n",cloud->maxX-cloud->minX,cloud->maxY-cloud->minY,cloud->maxZ-cloud->minZ);
	fclose(f);
}

void saveDetections(char* filename,HPCD* scene,HPCD* model,std::vector<Point> *modelPoints,std::vector<Point> *detection) {
	std::vector<Point> output;
	for (size_t i=0;i<modelPoints->size();i++) {
		Point h = modelPoints->at(i);
		for (size_t j=0;j<detection->size();j++) {
			Point p = {
				scene->minX + detection->at(j).x * scene->leafSize + h.x - model->minX,
				scene->minY + detection->at(j).y * scene->leafSize + h.y - model->minY,
				scene->minZ + detection->at(j).z * scene->leafSize + h.z - model->minZ,
				0,0,0};
			output.push_back(p);
		}
	}
	HPCD_writePoints(filename,&output);
}

void saveDetectionsToFolder(char* foldername,HPCD* scene,HPCD* model,std::vector<Point> *modelPoints,std::vector<Point> *detection) {
	char filename[128];
	for (size_t j=0;j<detection->size();j++) {
		std::vector<Point> output;
		sprintf(filename,"%s/%lu-cloud.pcd",foldername,j);
		for (size_t i=0;i<modelPoints->size();i++) {
			Point h = modelPoints->at(i);
			Point p = {
				scene->minX + detection->at(j).x * scene->leafSize + h.x - model->minX,
				scene->minY + detection->at(j).y * scene->leafSize + h.y - model->minY,
				scene->minZ + detection->at(j).z * scene->leafSize + h.z - model->minZ,
				0,0,0};
			output.push_back(p);
		}
		HPCD_writePoints(filename,&output);
	}
}

void meanShift(std::vector<Point> *points,std::vector<Point> *mu,float R) {
	bool* visited = new bool[points->size()]();

	while (true) {
		size_t n=0;
		while (n < points->size() && visited[n])
			n++;
		if (n >= points->size())
			break;

		Point seed = points->at(n);
		std::vector<Point> Q;
		std::vector<Point> S;
		Q.push_back(seed);
		visited[n] = true;

		while (Q.size() > 0) {
			Point p = Q[Q.size() - 1];
			Q.pop_back();
			S.push_back(p);
			for (size_t i=0;i<points->size();i++) {
				if (!visited[i]) {
					float dist = 0;
					dist += (points->at(i).x - p.x) * (points->at(i).x - p.x);
					dist += (points->at(i).y - p.y) * (points->at(i).y - p.y);
					dist += (points->at(i).z - p.z) * (points->at(i).z - p.z);
					if (dist < R * R) {
						Q.push_back(points->at(i));
						visited[i] = true;
					}
				}
			}
		}

		if (S.size() < 3)
			continue;
		Point c = {0,0,0,0,0,0};
		for (size_t i=0;i<S.size();i++) {
			c.x += S[i].x;
			c.y += S[i].y;
			c.z += S[i].z;
		}
		c.x /= S.size();
		c.y /= S.size();
		c.z /= S.size();
		mu->push_back(c);
	}
	printf("mean shift: before %lu after %lu\n",points->size(),mu->size());
//	for (size_t i=0;i<points->size();i++) {
//		printf("before %f %f %f\n",points->at(i).x,points->at(i).y,points->at(i).z);
//	}
//	for (size_t i=0;i<mu->size();i++) {
//		printf("after %f %f %f\n",mu->at(i).x,mu->at(i).y,mu->at(i).z);
//	}
	delete[] visited;
}

void anms(HPCD* scene,std::vector<int> *candidates,std::vector<Point> *mu, int clusterThreshold) {
	int numCombinations = 1;
	for (int i = 0; i < clusterThreshold; i++)
		numCombinations *= 6;
	for (size_t i=0;i<candidates->size();i++) {
		HPoint* hp = scene->data[candidates->at(i)];
		if (hp->label==0)
			continue;
		std::vector<HPoint*> Q;
		Point centroid = {};
		int count=0;
		Q.push_back(hp);
		while (Q.size() > 0) {
			HPoint* h = Q[Q.size()-1];
			Q.pop_back();
			if (h->label==0)
				continue;
			h->label=0;
			count++;
			centroid.x += h->x;
			centroid.y += h->y;
			centroid.z += h->z;
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
					int j = HPCD_find(scene, xi, yi, zi);
					if (j >= 0) {
						Q.push_back(scene->data[j]);
					}
					l /= 6;
				}
			}	
		}
		centroid.x /= count;
		centroid.y /= count;
		centroid.z /= count;
		mu->push_back(centroid);
	}
	printf("anms: before %lu after %lu (max score %f)\n",candidates->size(),mu->size(),maxScore);
}

void windowCorrelation(HPCD* scene,HPCD* model, std::vector<Point> *detection) {
	std::vector<int> offsetX;
	std::vector<int> offsetY;
	std::vector<int> offsetZ;
	std::vector<Point> candidate_points;
	std::vector<int> candidates;
	for (int i=0;i<model->maxSize;i++) {
		HPoint* h = model->data[i];
		if (h) {
			offsetX.push_back(h->x);
			offsetY.push_back(h->y);
			offsetZ.push_back(h->z);
		}
	}
	for (int i=0;i<scene->maxSize;i++) {
		HPoint* h = scene->data[i];
		if (h) {
			h->label=0;
			int count = 0;
			for (size_t k=0;k<offsetX.size();k++) {
				if (HPCD_find(scene,h->x+offsetX[k],h->y+offsetY[k],h->z+offsetZ[k]) >= 0)
					count++;
			}
			float score = 1.0 * count / model->numPoints;
			if (score > maxScore)
				maxScore = score;
			if (score > corr_threshold) {
				Point p = {h->x,h->y,h->z,0,0,0};
				candidate_points.push_back(p);
				candidates.push_back(i);
				h->label=1;
//				printf("%f %d %d %d\n",score,h->x,h->y,h->z);
			}
		}
	}
	//non-maximal suppression
	meanShift(&candidate_points,detection,merge_threshold / CLOUD_LEAF_SIZE);
//	anms(scene,&candidates,detection,3);
}

int main(int argc,char* argv[]) {
	if (argc < 4) {
//		printf("Usage: %s in.pcd out.lbp3d\n",argv[0]);
//		printf("Usage: %s scene.pcd model.pcd detection.pcd\n",argv[0]);
		printf("Usage: %s scene.pcd model.pcd detection [-c corr_threshold -m merge_threshold/\n",argv[0]);
		return 1;
	}
	
	char* sceneFile = argv[1];
	char* modelFile = argv[2];
	for (int i=4;i<argc;i++) {
		if (strcmp(argv[i],"-c")==0 && i+1<argc)
			corr_threshold = atof(argv[++i]);
		else if (strcmp(argv[i],"-m")==0 && i+1<argc)
			merge_threshold = atof(argv[++i]);
	}

//	HPCD* cloud = HPCD_Init(inFile,CLOUD_LEAF_SIZE);
//	float* feature = new float[FEATURE_SIZE]();
//	lbpHist(cloud,feature);
//	saveLBP(outFile,cloud,feature);
//	delete[] feature;
//	HPCD_delete(cloud);

	std::vector<Point> modelPoints;
	readPoints(modelFile,&modelPoints);
	HPCD* scene = HPCD_Init(sceneFile,CLOUD_LEAF_SIZE);
	HPCD* model = HPCD_Init(modelFile,CLOUD_LEAF_SIZE);
//	printf("%d %d\n",scene->numPoints,model->numPoints);
	std::vector<Point> detection;
	windowCorrelation(scene,model,&detection);
//	saveDetections(argv[3],scene,model,&modelPoints,&detection);
	saveDetectionsToFolder(argv[3],scene,model,&modelPoints,&detection);
	HPCD_delete(scene);
	HPCD_delete(model);
}
