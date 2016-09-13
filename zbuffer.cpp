#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <dirent.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/osmesa.h>
#define NUM_CAMERAS 4
#define RESOLUTION 600

//http://www.scratchapixel.com/old/lessons/3d-basic-lessons/lesson-10-polygonal-objects/

struct Point {
	double x,y,z;
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
				Point p;
				if (fgets(buf,256,f) && sscanf(buf, "%lf %lf %lf",&(p.x),&(p.y),&(p.z)) == 3) {
					vertices->push_back(p);
				} else {
					printf("Error parsing %s\n",filename);
					printf("Line %d: %s\n",i,buf);
					break;
				}
			}
			for (int i=0;i<numFace;i++) {
				Triangle t;
				if (fgets(buf,256,f) && sscanf(buf, "3 %lu %lu %lu",&(t.id1),&(t.id2),&(t.id3)) == 3) {
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

int main(int argc, char* argv[]) {
	if (argc < 3) {
		printf("Usage: ./zbuffer input.ply [output.pcd,outputFolder]\n");
		return 1;
	}

	std::vector<Point> vertices;
	std::vector<Triangle> faces;
	std::vector<Point> pointcloud;
	char buffer[128];
	bool merge = opendir(argv[2]) == NULL;

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
	for (size_t i = 0; i < vertices.size(); i++) {
		vertices[i].x -= centroid.x;
		vertices[i].y -= centroid.y;
		vertices[i].z -= centroid.z;
	}

	int width = RESOLUTION;
	int height = RESOLUTION;
	OSMesaContext ctx;
	ctx = OSMesaCreateContextExt(OSMESA_RGB, 32, 0, 0, NULL );
	unsigned char * pbuffer = new unsigned char [3 * width * height];
	// Bind the buffer to the context and make it current
	if (!OSMesaMakeCurrent(ctx, (void*)pbuffer, GL_UNSIGNED_BYTE, width, height))
		printf("fail to create MESA context\n");
	OSMesaPixelStore(OSMESA_Y_UP, 0);

    glEnable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	glDisable(GL_CULL_FACE);
	glPolygonMode(GL_FRONT, GL_FILL);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	float fov = 70;
	float fov_scale = 2 * tan(fov / 2 / 180 * M_PI);
	float zfar = 100000;
	gluPerspective(fov,1,1,zfar);
	glViewport(0, 0, width, height);
	float cameraX = maxX - minX;
	float cameraY = maxY - minY;
	float cameraZ = maxZ - minZ;
	float cx = 0.5 * (width + 1);
	float cy = 0.5 * (height + 1);
	unsigned int* depth = new unsigned int[width * height];
	float rho = sqrt(cameraX*cameraX + cameraY*cameraY);
	float theta = atan2(cameraY, cameraX);
	int depthBits=0;
	glGetIntegerv(GL_DEPTH_BITS, &depthBits);
	printf("depth buffer bits %d\n",depthBits);

	for (int k = 0; k < NUM_CAMERAS; k++) {
//		if (!merge)
//			pointcloud.clear();
		float rx = rho * cos(theta);
		float ry = rho * sin(theta);
		theta += 2 * 3.14159265 / NUM_CAMERAS;
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		gluLookAt(rx,ry, cameraZ, 0,0,0, 0,0,1);
		GLfloat R[16] =  {};
		glGetFloatv(GL_MODELVIEW_MATRIX, R);
//		printf("%f %f %f %f\n%f %f %f %f\n%f %f %f %f\n%f %f %f %f\n",R[0],R[1],R[2],R[3],R[4],R[5],R[6],R[7],R[8],R[9],R[10],R[11],R[12],R[13],R[14],R[15]);
//		printf("camera: %f %f %f\n",rx,ry,cameraZ);
		glBegin(GL_TRIANGLES);
		glColor3ub(150, 150, 150);
		for (size_t i = 0; i < faces.size(); i++) {
			Point p1 = vertices[faces[i].id1];
			Point p2 = vertices[faces[i].id2];
			Point p3 = vertices[faces[i].id3];
			glVertex3f(p1.x, p1.y, p1.z);
			glVertex3f(p2.x, p2.y, p2.z);
			glVertex3f(p3.x, p3.y, p3.z);
		}
		glEnd();
		glFinish(); // done rendering
		GLint outWidth, outHeight, bitPerDepth;
		GLboolean ret = OSMesaGetDepthBuffer(ctx, &outWidth, &outHeight, &bitPerDepth, (void**)&depth);
		for (int i = 0; i < width; i++) {
			for (int j = 0; j < height; j++) {
				float buffer_z = (float) depth[j * width + i] / 0xFFFFFFFF;
				if (buffer_z > 0 && buffer_z < 1) {
//					float z = -1 / (1 - buffer_z);
					float z = 1 / (buffer_z - 1 - buffer_z / zfar);
					float x = (i - cx) * -z / width * fov_scale;
					float y = (j - cy) * -z / height * fov_scale;
					x -= R[12];
					y -= R[13];
					z -= R[14];
					Point p = {
						R[0] * x + R[1] * y + R[2] * z,
						R[4] * x + R[5] * y + R[6] * z,
						R[8] * x + R[9] * y + R[10] * z
					};
					pointcloud.push_back(p);
				}
			}
		}
//		printf("pointcloud: %lu\n",pointcloud.size());
		if (!merge && pointcloud.size() > 0) {
			int n=0;
			while (true) {
				sprintf(buffer,"%s/%d-cloud.pcd",argv[2],n);
				FILE* f = fopen(buffer,"r");
				if (!f) {
					writeToPCD(buffer,&pointcloud);
					break;
				}
				fclose(f);
				n++;
			}
		}
	}
	if (merge && pointcloud.size() > 0) {
		sprintf(buffer,"%s",argv[2]);
		writeToPCD(buffer,&pointcloud);
	}

//	sprintf(buffer,"%s/vertex.pcd",argv[2]);
//	writeToPCD(buffer,&vertices);

	OSMesaDestroyContext(ctx);
	return 0;
}
