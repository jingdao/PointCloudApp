#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <math.h>
#include <unistd.h>
#include <SDL/SDL.h>
#include <GL/gl.h>
#include <GL/glu.h>

struct PCD {
	int numPoints;
	float* float_data;
};
enum PCD_data_storage {
	ASCII,
	BINARY,
	NONE
};
struct Point {
	float x,y,z;
};
struct Triangle {
	int id1, id2, id3;
};
double cameraX=20,cameraY=20,cameraZ=10;
double centerX=0,centerY=0,centerZ=0;
double upX=0,upY=0,upZ=1;
int mouseIndex = 0;
int previousX,previousY;
double scrollSpeed = 1.01;

std::vector<PCD*> cloud;
PCD* building = NULL;
std::vector< std::vector<Point> > boxes;
std::vector<Point> vertices;
std::vector<Triangle> faces;

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

void getPCA(PCD *cloud, std::vector<Point> *box) {
	double cov[4] = {}; //column major
	Point center = {};
	for (int i = 0; i < cloud->numPoints; i++) {
		center.x += cloud->float_data[i*4];
		center.y += cloud->float_data[i*4+1];
		center.z += cloud->float_data[i*4+2];
	}
	center.x /= cloud->numPoints; 
	center.y /= cloud->numPoints;
	center.z /= cloud->numPoints;
	for (int j = 0; j<cloud->numPoints; j++) {
		float deltaP[2] = {
			cloud->float_data[j*4]- center.x,
			cloud->float_data[j*4+1]- center.y,
		};
		cov[0] += deltaP[0] * deltaP[0];
		cov[1] += deltaP[1] * deltaP[0];
		cov[2] += deltaP[0] * deltaP[1];
		cov[3] += deltaP[1] * deltaP[1];
	}
	cov[0] /= cloud->numPoints * cloud->numPoints;
	cov[1] /= cloud->numPoints * cloud->numPoints;
	cov[2] /= cloud->numPoints * cloud->numPoints;
	cov[3] /= cloud->numPoints * cloud->numPoints;
	float trace = cov[0] + cov[3];
	float det = cov[0] * cov[3] - cov[1] * cov[2];
	float L1 = trace / 2 + sqrt(trace*trace / 4 - det);
	float L2 = trace / 2 - sqrt(trace*trace / 4 - det);
	float minScale[3], maxScale[3];
	float v[9] = {
		0,0,0,
		0,0,0,
		0,0,1
	};
	if (cov[2] != 0) {
		v[0] = L1 - cov[3];
		v[1] = L2 - cov[3];
		v[3] = v[4] = cov[2];
	}
	else if (cov[1] != 0) {
		v[0] = v[1] = cov[1];
		v[3] = L1 - cov[0];
		v[4] = L2 - cov[0];
	}
	else {
		v[0] = v[4] = 1;
	}
	float m1 = sqrt(v[0] * v[0] + v[3] * v[3]);
	float m2 = sqrt(v[1] * v[1] + v[4] * v[4]);
	v[0] /= m1;
	v[3] /= m1;
	v[1] /= m2;
	v[4] /= m2;
	for (int j = 0; j<cloud->numPoints; j++) {
		for (int i = 0; i<3; i++) {
			float dotProduct =
				cloud->float_data[j*4] * v[i * 3] +
				cloud->float_data[j*4+1] * v[i * 3 + 1] +
				cloud->float_data[j*4+2] * v[i * 3 + 2];
			if (j == 0 || dotProduct < minScale[i])
				minScale[i] = dotProduct;
			if (j == 0 || dotProduct > maxScale[i])
				maxScale[i] = dotProduct;
		}
	}
	float bbCenter[3] = {0,0,0};
	for (int i = 0; i<3; i++) {
		bbCenter[0] += (minScale[i] + maxScale[i]) / 2 * v[i * 3];
		bbCenter[1] += (minScale[i] + maxScale[i]) / 2 * v[i * 3 + 1];
		bbCenter[2] += (minScale[i] + maxScale[i]) / 2 * v[i * 3 + 2];
	}
	for (int i = 0; i<8; i++) {
		float coords[3];
		for (int j = 0; j<3; j++) {
			coords[j] = bbCenter[j];
			for (int axis = 0; axis<3; axis++) {
				float sign = (i & 1 << axis) ? 1 : -1;
				coords[j] += sign * (maxScale[axis]-minScale[axis]) / 2 * v[axis * 3 + j];
			}
		}
		Point p = {coords[0],coords[1],coords[2]};
		box->push_back(p);
	}
}

void loadOBJ(char* fileName,std::vector<Point> *vt, std::vector<Triangle> *fc) {
	FILE* f = fopen(fileName,"r");
	char buf[256];
	while (fgets(buf,256,f)) {
		if (buf[0]=='v' && buf[1]==' ') {
			char *c = buf + 2;
			Point p;
			p.x = strtod(c,&c);
			p.y = strtod(c,&c);
			p.z = strtod(c,&c);
			vt->push_back(p);
		} else if (buf[0]=='f' && buf[1]==' ') {
			char* c = buf + 2;
			Triangle t;
			t.id1 = strtol(c,&c,10) - 1;
			c = strchr(c,' ');
			t.id2 = strtol(c+1,&c,10) - 1;
			c = strchr(c,' ');
			t.id3 = strtol(c+1,&c,10) - 1;
			fc->push_back(t);
		}
	}
	fclose(f);
}

void centerMap() {
	int p=building->numPoints;
	float cx=0,cy=0,cz=0;
	for (int n=0;n<building->numPoints;n++) {
		cx += building->float_data[n*4];
		cy += building->float_data[n*4+1];
		cz += building->float_data[n*4+2];
	}
	cx /= p;
	cy /= p;
	cz /= p;
	for (int n=0;n<building->numPoints;n++) {
		building->float_data[n*4] -= cx;
		building->float_data[n*4+1] -= cy;
		building->float_data[n*4+2] -= cz;
	}
	for (int i=0;i<cloud.size();i++) {
		for (int n = 0; n < cloud[i]->numPoints; n++){
			cloud[i]->float_data[n*4] -= cx;
			cloud[i]->float_data[n*4+1] -= cy;
			cloud[i]->float_data[n*4+2] -= cz;
		}
	}
	for (int i=0;i<vertices.size();i++) {
		vertices[i].x -= cx;
		vertices[i].y -= cy;
		vertices[i].z -= cz;
	}
}

void drawBox(Point* p) {
	glBegin(GL_LINES);

	glVertex3d(p[0].x,p[0].y,p[0].z);
	glVertex3d(p[1].x,p[1].y,p[1].z);
	glVertex3d(p[0].x,p[0].y,p[0].z);
	glVertex3d(p[2].x,p[2].y,p[2].z);
	glVertex3d(p[1].x,p[1].y,p[1].z);
	glVertex3d(p[3].x,p[3].y,p[3].z);
	glVertex3d(p[2].x,p[2].y,p[2].z);
	glVertex3d(p[3].x,p[3].y,p[3].z);

	glVertex3d(p[0].x,p[0].y,p[0].z);
	glVertex3d(p[4].x,p[4].y,p[4].z);
	glVertex3d(p[1].x,p[1].y,p[1].z);
	glVertex3d(p[5].x,p[5].y,p[5].z);
	glVertex3d(p[2].x,p[2].y,p[2].z);
	glVertex3d(p[6].x,p[6].y,p[6].z);
	glVertex3d(p[3].x,p[3].y,p[3].z);
	glVertex3d(p[7].x,p[7].y,p[7].z);

	glVertex3d(p[4].x,p[4].y,p[4].z);
	glVertex3d(p[5].x,p[5].y,p[5].z);
	glVertex3d(p[4].x,p[4].y,p[4].z);
	glVertex3d(p[6].x,p[6].y,p[6].z);
	glVertex3d(p[5].x,p[5].y,p[5].z);
	glVertex3d(p[7].x,p[7].y,p[7].z);
	glVertex3d(p[6].x,p[6].y,p[6].z);
	glVertex3d(p[7].x,p[7].y,p[7].z);

	glEnd();
}

void draw() {
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glPushMatrix();
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(cameraX,cameraY,cameraZ,centerX,centerY,centerZ,upX,upY,upZ);

	glPointSize(1.0);
	glBegin(GL_POINTS);
	for (int i=0;i<cloud.size();i++) {
//		glColor3ub(255,0,0);
		glColor3ub(50,50,50);
		for (int n = 0; n < cloud[i]->numPoints; n++){
			glVertex3d(cloud[i]->float_data[n*4],cloud[i]->float_data[n*4+1],cloud[i]->float_data[n*4+2]);
		}
	}
	glColor3ub(50,50,50);
	for (int n=0;n<building->numPoints;n++) {
		int c = building->float_data[n*4+3];
		unsigned char b = c & 0xFF;
		unsigned char g = (c >> 8) & 0xFF;
		unsigned char r = (c >> 16) & 0xFF;
		glColor3ub(r,g,b);
		glVertex3d(building->float_data[n*4],building->float_data[n*4+1],building->float_data[n*4+2]);
	}
	glEnd();

	glLineWidth(1.0);
	glColor3ub(255,255,0);
	for (int i=0;i<boxes.size();i++) {
		drawBox(boxes[i].data());
	}

	glBegin(GL_TRIANGLES);
	glColor3ub(0, 0, 255);
	for (size_t i = 0; i < faces.size(); i++) {
//	for (size_t i = 0; i < 224; i++) {
		Point p1 = vertices[faces[i].id1];
		Point p2 = vertices[faces[i].id2];
		Point p3 = vertices[faces[i].id3];
		glVertex3f(p1.x, p1.y, p1.z);
		glVertex3f(p2.x, p2.y, p2.z);
		glVertex3f(p3.x, p3.y, p3.z);
	}
	glEnd();

	glFlush();
	SDL_GL_SwapBuffers();
	glPopMatrix();
	glPopAttrib();
}

int main(int argc,char* argv[]) {
	if (argc < 3) {
		printf("%s building.pcd [clusterFolder | element.obj]\n",argv[0]);
		return 1;
	}

	building = NewPCD(argv[1]);
	if (!building)
		return 1;

	size_t i=0;
	if (strncmp(".obj",argv[2]+strlen(argv[2])-4,4)==0) {
		loadOBJ(argv[2],&vertices,&faces);
		printf("Loaded %lu vertices %lu faces\n",vertices.size(),faces.size());
	} else {
		char buf[256];
		while (true) {
			sprintf(buf,"%s/%lu-cloud.pcd",argv[2],i);
			PCD* c = NewPCD(buf);
			if (!c)
				break;
			cloud.push_back(c);
			i++;
		}
		printf("Loaded %lu point clouds\n",cloud.size());
	}	
	centerMap();
	for (i=0;i<cloud.size();i++) {
		std::vector<Point> box;
		getPCA(cloud[i],&box);
		boxes.push_back(box);
	}

	SDL_Init(SDL_INIT_VIDEO);
	SDL_WM_SetCaption("Point Cloud", NULL);
	SDL_SetVideoMode(1800,1000, 32, SDL_OPENGL);
    glEnable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	glDisable(GL_CULL_FACE);
	glPolygonMode(GL_FRONT, GL_FILL);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(70,(double)640/480,1,1000);

	int interval = 10000;
	SDL_Event event;
	while (SDL_PollEvent(&event)); //clear event buffer
	draw();
	while (true) {
		while (SDL_PollEvent(&event)) {
			switch(event.type){
				case SDL_KEYDOWN:
					switch( event.key.keysym.sym ){
						case SDLK_LEFT:
						cameraX -= 1;
						break;
						case SDLK_RIGHT:
						cameraX += 1;
						break;
						case SDLK_UP:
						cameraZ += 1;
						break;
						case SDLK_DOWN:
						cameraZ -= 1;
						break;
						default:
						break;
					}
					draw();
					break;
				case SDL_MOUSEBUTTONDOWN:
					if (event.button.button == SDL_BUTTON_WHEELUP) {
						cameraX /= scrollSpeed;
						cameraY /= scrollSpeed;
						cameraZ /= scrollSpeed;
						draw();
					} else if (event.button.button == SDL_BUTTON_WHEELDOWN) {
						cameraX *= scrollSpeed;
						cameraY *= scrollSpeed;
						cameraZ *= scrollSpeed;
						draw();
					} else {
						mouseIndex = event.button.button == SDL_BUTTON_LEFT ? 1 : 2;
						previousX = event.button.x;
						previousY = event.button.y;
					}
					break;
				case SDL_MOUSEBUTTONUP:
					mouseIndex = 0;
					break;
				case SDL_MOUSEMOTION:
					if (mouseIndex == 1) {
						double rho = sqrt(cameraX*cameraX+cameraY*cameraY);
						double xstep = cameraY / rho;
						double ystep = -cameraX / rho;
						cameraX += 0.05 * (event.motion.x-previousX) * xstep;
						cameraY += 0.05 * (event.motion.x-previousX) * ystep;
						cameraZ += 0.05 * (event.motion.y-previousY);
						previousX = event.motion.x;
						previousY = event.motion.y;
						draw();
					}
					break;
				case SDL_QUIT:
					exit(0);
					break;
			}
		}
		usleep(interval);
	}

//	atexit(SQL_Quit);

	return 0;
}
