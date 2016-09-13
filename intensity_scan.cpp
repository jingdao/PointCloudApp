#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>

typedef struct {
	float x,y,z,intensity;
} Point;

std::vector<Point> LoadFromPTS(char* fileName) {
	std::vector<Point> pcd;
	FILE* f = fopen(fileName, "r");
	if (!f) {
		printf("File not found: %s\n", fileName);
		return pcd;
	}
	char buf[256];
	int pointsParsed = 0;
	float x,y,z,its;
	int r,g,b;
	float minIts = 99999;
	float maxIts = -99999;
	while (fgets(buf, 256, f)) {
		int n = atoi(buf);
		for (int i=0;i<n;i++) {
			fgets(buf,256,f);
			if (sscanf(buf, "%f %f %f %f %d %d %d",&x,&y,&z,&its,&r,&g,&b) == 7) {
				Point p = {x,y,z,its};
				pcd.push_back(p);
				if (its < minIts) minIts = its;
				if (its > maxIts) maxIts = its;
				pointsParsed++;
			} else {
				printf("Error parsing %s\n",fileName);
				printf("Line %d: %s\n",pointsParsed,buf);
				break;
			}
		}
	}
	fclose(f);
	printf("Parsed %s (%d points) min %f max %f\n",fileName,pointsParsed,minIts,maxIts);
	for (int i=0;i<pointsParsed;i++)
		pcd[i].intensity = (pcd[i].intensity - minIts) / (maxIts - minIts);
	return pcd;
}

//void writeToPCD(std::vector<Point> *p, char* fileName) {
//	FILE* f = fopen(fileName, "w");
//	fprintf(f,"# .PCD v0.7 - Point Cloud Data file format\n"
//	"VERSION 0.7\n"
//	"FIELDS x y z intensity\n"
//	"SIZE 4 4 4 4\n"
//	"TYPE F F F F\n"
//	"COUNT 1 1 1 1\n"
//	"WIDTH %lu\n"
//	"HEIGHT 1\n"
//	"VIEWPOINT 0 0 0 1 0 0 0\n"
//	"POINTS %lu\n"
//	"DATA ascii\n",p->size(),p->size());
//	for (size_t i=0;i<p->size();i++) {
//		Point q = p->at(i);
//		fprintf(f,"%f %f %f %f\n",q.x,q.y,q.z,q.intensity);
//	}
//	fclose(f);
//}

void writeToPCD(std::vector<Point> *p, char* fileName) {
	FILE* f = fopen(fileName, "w");
	fprintf(f,"# .PCD v0.7 - Point Cloud Data file format\n"
	"VERSION 0.7\n"
	"FIELDS x y z rgb\n"
	"SIZE 4 4 4 4\n"
	"TYPE F F F I\n"
	"COUNT 1 1 1 1\n"
	"WIDTH %lu\n"
	"HEIGHT 1\n"
	"VIEWPOINT 0 0 0 1 0 0 0\n"
	"POINTS %lu\n"
	"DATA ascii\n",p->size(),p->size());
	for (size_t i=0;i<p->size();i++) {
		Point q = p->at(i);
		int r = q.intensity * 255;
		int rgb = (r<<16) | (r<<8) | r;
		fprintf(f,"%f %f %f %d\n",q.x,q.y,q.z,rgb);
	}
	fclose(f);
}

int main(int argc, char* argv[]) {

	if (argc < 2) {
		printf("%s cloud.pts/cloud.pcd [out.pcd]\n",argv[0]);
		return 1;
	}

	std::vector<Point> p = LoadFromPTS(argv[1]);
	writeToPCD(&p,argv[2]);

}
