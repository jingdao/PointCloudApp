#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

struct Point {
	float x,y,z;
};

struct Vector {
	float x,y,z;
};

struct Plane {
	float a,b,c,d;
};

struct Triangle {
	size_t id1,id2,id3;
};

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
				if (sscanf(buf, "%f %f %f",&(p.x),&(p.y),&(p.z)) == 3) {
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
	Plane plane = {
		crossProduct.x,
		crossProduct.y,
		crossProduct.z,
		p1.x*crossProduct.x + p1.y*crossProduct.y + p1.z*crossProduct.z
	};
	return plane;
}

bool intersects(Point rayOrigin,Vector rayDirection,Plane plane,std::vector<Point> *vertices,Triangle tri,float *distance) {
	//if ray perpendicular to triangle
	if (fabs(rayDirection.x*plane.a + rayDirection.y*plane.b + rayDirection.z*plane.c) < 0.1)
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


int main(int argc, char* argv[]) {
	if (argc < 3) {
		printf("./ray_tracing input.ply output.pcd\n");
		return 1;
	}

	std::vector<Point> vertices;
	std::vector<Triangle> faces;
	std::vector<Plane> planes;

	readPLY(argv[1],&vertices,&faces);

	//get bounding box
	float minX=vertices[0].x,maxX=vertices[0].x;
	float minY=vertices[0].y,maxY=vertices[0].y;
	float minZ=vertices[0].z,maxZ=vertices[0].z;
	for (size_t i=1;i<vertices.size();i++) {
		if (vertices[i].x < minX) minX = vertices[i].x;
		else if (vertices[i].x > maxX) maxX = vertices[i].x;
		if (vertices[i].y < minY) minY = vertices[i].y;
		else if (vertices[i].y > maxY) maxY = vertices[i].y;
		if (vertices[i].z < minZ) minZ = vertices[i].z;
		else if (vertices[i].z > maxZ) maxZ = vertices[i].z;
	}
	printf("Bounding box: x:(%.1f %.1f) y:(%.1f %.1f) z:(%.1f %.1f)\n",minX,maxX,minY,maxY,minZ,maxZ);

	//get normals
	for (size_t i=0;i<faces.size();i++) {
		Plane v = getPlane(&vertices,faces[i]);
		planes.push_back(v);
	}

	Point rayOrigin = {0,0,-1};
	Vector rayDirection = {0,0,1};
	std::vector<Point> v;
	Point p1 = {0,0,0};
	Point p2 = {0,1,0};
	Point p3 = {1,0,0};
	Triangle t = {0,1,2};
	v.push_back(p1);
	v.push_back(p2);
	v.push_back(p3);
	Plane pl = getPlane(&v,t);
	float distance;
	intersects(rayOrigin,rayDirection,pl,&v,t,&distance);


	return 0;
}
