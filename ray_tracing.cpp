#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#define USE_Y_VERTICAL 0
#define RADIAL_SAMPLING 1
#define ZERO_THRESHOLD 0.01

//http://www.scratchapixel.com/old/lessons/3d-basic-lessons/lesson-10-polygonal-objects/

struct Point {
	double x,y,z;
};

struct Vector {
	double x,y,z;
};

struct Plane {
	double a,b,c,d;
};

struct Triangle {
	size_t id1,id2,id3;
};

double magnitude(Vector v) {
	return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

Vector normalize(Vector v) {
	double m = magnitude(v);
	Vector w = {v.x/m, v.y/m, v.z/m};
	return w;
}

void readPLY(char* filename, std::vector<Point> *vertices, std::vector<Triangle> *faces) {
	FILE* f = fopen(filename, "r");
	if (!f) {
		printf("File not found: %s\n", filename);
		return;
	}
	char buf[256];
	int numVertex,numFace;
	while (fgets(buf, 256, f)) {
		if (sscanf(buf,"element vertex %d",&numVertex)==1) {
		} else if (sscanf(buf,"element face %d",&numFace)==1) {
		} else if (strncmp(buf,"end_header",10)==0) {
			for (int i=0;i<numVertex;i++) {
				fgets(buf,256,f);
				Point p;
				if (sscanf(buf, "%lf %lf %lf",&(p.x),&(p.y),&(p.z)) == 3) {
					vertices->push_back(p);
				} else {
					printf("Error parsing %s\n",filename);
					printf("Line %d: %s\n",i,buf);
					break;
				}
			}
			for (int i=0;i<numFace;i++) {
				fgets(buf,256,f);
				Triangle t;
				if (sscanf(buf, "3 %lu %lu %lu",&(t.id1),&(t.id2),&(t.id3)) == 3) {
					faces->push_back(t);
				} else {
					printf("Error parsing %s\n",filename);
					printf("Line %d: %s\n",i,buf);
					break;
				}
			}
			break;
		}
	}
	fclose(f);
}

void writeToPCD(char* filename,std::vector<Point> *pointcloud) {
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
	"DATA ascii\n",pointcloud->size(),pointcloud->size());
	for (size_t i=0;i<pointcloud->size();i++) {
		fprintf(f,"%f %f %f\n",(*pointcloud)[i].x,(*pointcloud)[i].y,(*pointcloud)[i].z);
	}
	fclose(f);
	printf("Wrote %lu points to %s\n",pointcloud->size(),filename);
}

Plane getPlane(std::vector<Point> *vertices, Triangle tri) {
	Point p1 = (*vertices)[tri.id1];
	Point p2 = (*vertices)[tri.id2];
	Point p3 = (*vertices)[tri.id3];
	Vector v1 = {p1.x - p2.x, p1.y - p2.y, p1.z - p2.z};
	Vector v2 = {p1.x - p3.x, p1.y - p3.y, p1.z - p3.z};
	Vector crossProduct = {
		v1.y * v2.z - v1.z * v2.y,
		v1.z * v2.x - v1.x * v2.z,
		v1.x * v2.y - v1.y * v2.x
	};
	crossProduct = normalize(crossProduct);
	Plane plane = {
		crossProduct.x,
		crossProduct.y,
		crossProduct.z,
		p1.x*crossProduct.x + p1.y*crossProduct.y + p1.z*crossProduct.z
	};
	return plane;
}

bool intersects(Point rayOrigin,Vector rayDirection,Plane plane,double *distance) {
	//if ray perpendicular to triangle
	if (fabs(rayDirection.x*plane.a + rayDirection.y*plane.b + rayDirection.z*plane.c) < ZERO_THRESHOLD)
		return false;
	//let intersection P = rayOrigin + distance * rayDirection
	//					P lies in plane described by ax+by+cz=d
	//calculate distance
	*distance = -(plane.a*rayOrigin.x + plane.b*rayOrigin.y + plane.c*rayOrigin.z - plane.d) / 
				(plane.a*rayDirection.x + plane.b*rayDirection.y + plane.c*rayDirection.z);
	if (*distance <= 0)
		return false;
	return true;
}

bool triangleContains(std::vector<Point> *vertices,Triangle tri,Plane plane,Point target) {
	Point p1 = (*vertices)[tri.id1];
	Point p2 = (*vertices)[tri.id2];
	Point p3 = (*vertices)[tri.id3];
	Vector v1 = {p2.x - p1.x, p2.y - p1.y, p2.z - p1.z};
	Vector v2 = {p3.x - p2.x, p3.y - p2.y, p3.z - p2.z};
	Vector v3 = {p1.x - p3.x, p1.y - p3.y, p1.z - p3.z};
	Vector c1 = {target.x - p1.x, target.y - p1.y, target.z - p1.z};
	Vector c2 = {target.x - p2.x, target.y - p2.y, target.z - p2.z};
	Vector c3 = {target.x - p3.x, target.y - p3.y, target.z - p3.z};
	double tripleproduct1 =
		plane.a * (v1.y * c1.z - v1.z * c1.y) + 
		plane.b * (v1.z * c1.x - v1.x * c1.z) + 
		plane.c * (v1.x * c1.y - v1.y * c1.x);
	double tripleproduct2 =
		plane.a * (v2.y * c2.z - v2.z * c2.y) + 
		plane.b * (v2.z * c2.x - v2.x * c2.z) + 
		plane.c * (v2.x * c2.y - v2.y * c2.x);
	double tripleproduct3 =
		plane.a * (v3.y * c3.z - v3.z * c3.y) + 
		plane.b * (v3.z * c3.x - v3.x * c3.z) + 
		plane.c * (v3.x * c3.y - v3.y * c3.x);

	return tripleproduct1 > 0 && tripleproduct2 > 0 && tripleproduct3 > 0;
}

int main(int argc, char* argv[]) {
	if (argc < 3) {
		printf("Usage: ./ray_tracing input.ply outputFolder\n");
		return 1;
	}

	std::vector<Point> vertices;
	std::vector<Triangle> faces;
	std::vector<Plane> planes;
	std::vector<Point> pointcloud;

	readPLY(argv[1],&vertices,&faces);

	//get bounding box
	double minX=vertices[0].x,maxX=vertices[0].x;
	double minY=vertices[0].y,maxY=vertices[0].y;
	double minZ=vertices[0].z,maxZ=vertices[0].z;
	for (size_t i=1;i<vertices.size();i++) {
		if (vertices[i].x < minX) minX = vertices[i].x;
		else if (vertices[i].x > maxX) maxX = vertices[i].x;
		if (vertices[i].y < minY) minY = vertices[i].y;
		else if (vertices[i].y > maxY) maxY = vertices[i].y;
		if (vertices[i].z < minZ) minZ = vertices[i].z;
		else if (vertices[i].z > maxZ) maxZ = vertices[i].z;
	}
	printf("Bounding box: x:(%.2f %.2f) y:(%.2f %.2f) z:(%.2f %.2f)\n",minX,maxX,minY,maxY,minZ,maxZ);
	Point centroid = {
		(minX + maxX) / 2,
		(minY + maxY) / 2,
		(minZ + maxZ) / 2
	};

	//get normals
	for (size_t i=0;i<faces.size();i++) {
		Plane v = getPlane(&vertices,faces[i]);
		planes.push_back(v);
	}

	double resolution = 0.001; //radians
	int numCameras=16;
	char buffer[128];
#if USE_Y_VERTICAL
	double radius = (maxX-minX) > (maxZ-minZ) ? (maxX-minX)*2 : (maxZ-minZ)*2; 
#else
	double radius = (maxX-minX) > (maxY-minY) ? (maxX-minX)*2 : (maxY-minY)*2; 
#endif
	double noise_sigma = 0 * radius;
	double alpha=M_PI/2/numCameras;
	for (int k=0;k<numCameras;k++) {
		pointcloud.clear();
#if USE_Y_VERTICAL
		Point cameraOrigin = {
			centroid.x + radius * sin(alpha) + noise_sigma * rand() / RAND_MAX,
			centroid.y + noise_sigma * rand() / RAND_MAX,
			centroid.z + radius * cos(alpha) + noise_sigma * rand() / RAND_MAX
		};
#else
		Point cameraOrigin = {
			centroid.x + radius * cos(alpha) + noise_sigma * rand() / RAND_MAX,
			centroid.y + radius * sin(alpha) + noise_sigma * rand() / RAND_MAX,
			centroid.z + noise_sigma * rand() / RAND_MAX
		};
#endif
		Vector principalDirection = {
			centroid.x - cameraOrigin.x,
			centroid.y - cameraOrigin.y,
			centroid.z - cameraOrigin.z
		};
		principalDirection = normalize(principalDirection);
#if RADIAL_SAMPLING
		for (double theta=-M_PI/4;theta<M_PI/4;theta+=resolution) {
			 for (double phi=-M_PI/4;phi<M_PI/4;phi+=resolution) {
#if USE_Y_VERTICAL
				Vector rayDirection = {
					principalDirection.z * sin(theta) * cos(phi) + principalDirection.x * cos(theta) * cos(phi),
					sin(phi),
					principalDirection.z * cos(theta) * cos(phi) - principalDirection.x * sin(theta) * cos(phi)
				};
#else
				Vector rayDirection = {
					principalDirection.x * cos(theta) * cos(phi) - principalDirection.y * sin(theta) * cos(phi),
					principalDirection.x * sin(theta) * cos(phi) + principalDirection.y * cos(theta) * cos(phi),
					sin(phi)
				};
#endif
				Point rayOrigin = cameraOrigin;
#else
		for (double theta=minZ;theta<maxZ;theta+=(maxZ-minZ)*resolution) {
			for (double phi=-radius/4;phi<radius/4;phi+=(radius/2)*resolution) {
				Vector rayDirection = principalDirection;
#if USE_Y_VERTICAL
				Point rayOrigin = {
					cameraOrigin.x + phi * principalDirection.z,
					theta,
					cameraOrigin.z - phi * principalDirection.x
				};
#else
				Point rayOrigin = {
					cameraOrigin.x + phi * principalDirection.y,
					cameraOrigin.y - phi * principalDirection.x,
					theta
				};
#endif
#endif
				bool isValid = false;
				Point closestPoint;
				double minDistance = DBL_MAX;
				//find intersection for each triangle
				for (size_t i=0;i<faces.size();i++) {
					double distance;
					if (intersects(rayOrigin,rayDirection,planes[i],&distance)) {
						if (distance < minDistance) {
							Point intersection = {
								rayOrigin.x + rayDirection.x * distance,
								rayOrigin.y + rayDirection.y * distance,
								rayOrigin.z + rayDirection.z * distance
							};
							if (triangleContains(&vertices,faces[i],planes[i],intersection)) {
								isValid = true;
								closestPoint = intersection;
								minDistance = distance;
							}
						}
					}
				}
				if (isValid) {
					pointcloud.push_back(closestPoint);
				}
			 }
		}

		alpha += 2*M_PI/numCameras;
		if (pointcloud.size() > 0) {
			sprintf(buffer,"%s/%d-cloud.pcd",argv[2],k);
			writeToPCD(buffer,&pointcloud);
		}
	}
	return 0;
}
