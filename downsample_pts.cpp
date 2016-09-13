#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#define BASE_HASH_CONSTANT 0.618033988
#define STEP_HASH_CONSTANT 0.707106781
#define STRING_HASH_CONSTANT 5381

typedef struct {
	int x,y,z;
} HPoint;

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

int main(int argc,char* argv[]) {
	if (argc < 4) {
		printf("Usage: %s in.pts out.pcd resolution\n",argv[0]);
		return 1;
	}

	FILE* f = fopen(argv[1], "r");
	if (!f) {
		printf("File not found: %s\n", argv[1]);
		return 1;
	}

	char buf[256];
	float resolution = atof(argv[3]);
	float x,y,z,its;
	int r,g,b,rgb;
	float minIts,maxIts;
	float minX,maxX,minY,maxY,minZ,maxZ; 
	bool use_intensity = false;
	int pointsParsed = 0;
	while (fgets(buf, 256, f)) {
		int n = atoi(buf);
		fgets(buf,256,f);
		sscanf(buf, "%f %f %f %f %d %d %d",&x,&y,&z,&its,&r,&g,&b);
		pointsParsed++;
		minIts = maxIts = its;
		minX = maxX = x;
		minY = maxY = y;
		minZ = maxZ = z;
		for (int i=1;i<n;i++) {
			fgets(buf,256,f);
			int scanned = sscanf(buf, "%f %f %f %f %d %d %d",&x,&y,&z,&its,&r,&g,&b);
			pointsParsed++;
			if (scanned == 7 || scanned == 4) {
				if (scanned == 4)
					use_intensity = true;
				if (its < minIts) minIts = its;
				if (its > maxIts) maxIts = its;
				if (x < minX) minX = x;
				else if (x > maxX) maxX = x;
				if (y < minY) minY = y;
				else if (y > maxY) maxY = y;
				if (z < minZ) minZ = z;
				else if (z > maxZ) maxZ = z;
			} else {
				printf("Error parsing %s\n",argv[1]);
				printf("Line %d: %s\n",pointsParsed,buf);
				return 1;
			}
		}
	}
	fclose(f);

	FILE* out = fopen(argv[2],"w");
	int checkpoint1 = fprintf(out,"# .PCD v0.7 - Point Cloud Data file format\n"
	"VERSION 0.7\n"
	"FIELDS x y z rgb\n"
	"SIZE 4 4 4 4\n"
	"TYPE F F F I\n"
	"COUNT 1 1 1 1\n"
	"WIDTH ");
	int checkpoint2 = fprintf(out,
	"%10d\n"
	"HEIGHT 1\n"
	"VIEWPOINT 0 0 0 1 0 0 0\n"
	"POINTS ",0);
	fprintf(out,"%10d\n"
	"DATA ascii\n",0);
	f = fopen(argv[1],"r");
	int numPoints = 0;
	HPoint** table = new HPoint*[pointsParsed];
	while (fgets(buf,256,f)) {
		int n = atoi(buf);
		for (int i=0;i<n;i++) {
			fgets(buf,256,f);
			sscanf(buf, "%f %f %f %f %d %d %d",&x,&y,&z,&its,&r,&g,&b);
			int xi = (int) ((x-minX)/resolution);
			int yi = (int) ((y-minY)/resolution);
			int zi = (int) ((z-minZ)/resolution);
			int ikey = getIntKey(xi,yi,zi);
			int key = baseHash(pointsParsed,ikey);
			int step = stepHash(pointsParsed,ikey);
			for (int k = 0; k < pointsParsed; k++) {
				HPoint* h = table[key];
				if (!h) {
					HPoint* p = new HPoint;
					p->x = xi;
					p->y = yi;
					p->z = zi;
					table[key] = p;
					numPoints++;
					if (use_intensity) {
						r = (its-minIts) / (maxIts-minIts) * 255;
						rgb = (r<<16) | (r<<8) | r;
					} else
						rgb = (r<<16) | (g<<8) | b;
					fprintf(out,"%f %f %f %d\n",x,y,z,rgb);
					break;
				} else if (h->x == xi && h->y == yi && h->z == zi){
					break;
				} else {
					key += step;
					key %= pointsParsed;
				}
			}
		}
	}
	fclose(f);
	fseek(out,checkpoint1,0);
	fprintf(out,"%10d",numPoints);
	fseek(out,checkpoint1+checkpoint2,0);
	fprintf(out,"%10d",numPoints);
	fclose(out);
}
	
