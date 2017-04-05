#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <dirent.h>
#include <unistd.h>
#include <SDL/SDL.h>
#include <GL/gl.h>
#include <GL/glu.h>
#define NUM_CAMERAS 8
#define RESOLUTION 800
#define ZFAR 100000
#define VERBOSE 0
#define NUM_COLORS 5
#define NUM_TEXTURE 10

struct Point {
	double x,y,z;
};
struct Vector {
	double x,y,z;
};
struct Triangle {
	size_t id1,id2,id3;
};
struct Color {
	unsigned char r,g,b;
};
struct Image {
	int width,height;
	unsigned char* data;
};

double cameraX=20,cameraY=20,cameraZ=10;
double centerX=0,centerY=0,centerZ=0;
double upX=0,upY=0,upZ=1;
int mouseIndex = 0;
int previousX,previousY;
double scrollSpeed = 1.1;
std::vector<Point> vertices;
std::vector<Triangle> faces;
std::vector<Vector> normals;
std::vector<Color> shading;
bool save_images;
Color previousColor;

Color palette[] = {
	{10,10,10}, //black
	{255,255,255}, //white
	{255,10,10}, //red
	{10,255,10}, //green
	{10,10,255}, //blue
};

bool loadPPM(char* name,Image *image) {
	char buffer[128];
	FILE* pgm = fopen(name,"r");
	if (!pgm) {
		printf("%s not found\n",name);
		return false;
	}
	fgets(buffer,128,pgm); //P5 or P6
	do {
		fgets(buffer,128,pgm);
	} while (buffer[0]=='#'); //remove comments
	char *c = buffer;
	image->width = strtol(c,&c,10);
	image->height = strtol(c,&c,10);
	fgets(buffer,128,pgm); //255
	if (image->data)
		delete[] image->data;
	image->data = new unsigned char[image->width*image->height*3];
	fread(image->data,1,image->width*image->height*3,pgm);
	//		for (int i=image->width*image->height;i>=0;i--) {
	//			unsigned char tmp = image->data[i*3+2];
	//			image->data[i*3+2] = image->data[i*3];
	//			image->data[i*3] = tmp;
	//		}
	return true;
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

Vector getTriangleNormal(Point p1, Point p2, Point p3) { //away from origin
	Point centroid = {
		(p1.x + p2.x + p3.x) / 3,
		(p1.y + p2.y + p3.y) / 3,
		(p1.z + p2.z + p3.z) / 3
	};
	Vector v1 = {p2.x-p1.x, p2.y-p1.y, p2.z-p1.z};
	Vector v2 = {p3.x-p1.x, p3.y-p1.y, p3.z-p1.z};
	Vector n = {
		v1.y*v2.z - v1.z*v2.y,
		v1.z*v2.x - v1.x*v2.z,
		v1.x*v2.y - v1.y*v2.x
	};
	double magn = sqrt(n.x*n.x + n.y*n.y + n.z*n.z);
	n.x /= magn;
	n.y /= magn;
	n.z /= magn;
	double dotProduct = centroid.x * n.x + centroid.y * n.y + centroid.z * n.z;
	if (dotProduct < 0) {
		n.x = -n.x;
		n.y = -n.y;
		n.z = -n.z;
	}
	return n;
}

void getShading(Color objectColor, Vector lightDir) {
	if (shading.size() == 0) {
		Color lightColor = {255,255,255};
		double ambientStrength = 0.5;

		for (size_t i=0;i<faces.size();i++) {
			Point p1 = vertices[faces[i].id1];
			Point p2 = vertices[faces[i].id2];
			Point p3 = vertices[faces[i].id3];
			Vector n = getTriangleNormal(p1,p2,p3);
			double diffuseStrength = lightDir.x*n.x + lightDir.y*n.y + lightDir.z*n.z;
			if (diffuseStrength < 0)
				diffuseStrength = 0;
			double strength = ambientStrength + diffuseStrength;
			if (strength > 1)
				strength = 1;
			Color c = {
				strength * lightColor.r * objectColor.r / 255,
				strength * lightColor.g * objectColor.g / 255,
				strength * lightColor.b * objectColor.b / 255,
			};
			shading.push_back(c);
			normals.push_back(n);
		}
	} else {
		for (size_t i=0;i<faces.size();i++) {
			shading[i].r = 1.0 * shading[i].r / previousColor.r * objectColor.r;
			shading[i].g = 1.0 * shading[i].g / previousColor.g * objectColor.g;
			shading[i].b = 1.0 * shading[i].b / previousColor.b * objectColor.b;
		}
	}
	previousColor = objectColor;
}

void writeImage(char* filename) {
	FILE* f = fopen(filename,"w");
	fprintf(f,"P6\n%d %d\n255\n",RESOLUTION,RESOLUTION);
	unsigned char* pixels = new unsigned char[RESOLUTION*RESOLUTION*3]();
	glReadPixels(0,0,RESOLUTION,RESOLUTION,GL_RGB,GL_UNSIGNED_BYTE,pixels);
	unsigned char* src = pixels + (RESOLUTION-1)*RESOLUTION*3;
	for (int i=0;i<RESOLUTION;i++) {
		fwrite(src,1,RESOLUTION*3,f);
		src -= RESOLUTION*3;
	}
	delete[] pixels;
#if VERBOSE
	printf("Wrote %d pixels to %s\n",RESOLUTION*RESOLUTION,filename);
#endif
	fclose(f);
}

void drawNoTexture() {
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glPushMatrix();
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(cameraX,cameraY,cameraZ,centerX,centerY,centerZ,upX,upY,upZ);

//	glLineWidth(1.0);
//	glColor3ub(200, 200, 200);
//	glBegin(GL_LINES);
//	for (size_t i = 0; i < faces.size(); i++) {
//		Point p1 = vertices[faces[i].id1];
//		Point p2 = vertices[faces[i].id2];
//		Point p3 = vertices[faces[i].id3];
//		glVertex3f(p1.x, p1.y, p1.z);
//		glVertex3f(p2.x, p2.y, p2.z);
//		glVertex3f(p2.x, p2.y, p2.z);
//		glVertex3f(p3.x, p3.y, p3.z);
//		glVertex3f(p3.x, p3.y, p3.z);
//		glVertex3f(p1.x, p1.y, p1.z);
//	}
//	glEnd();

	glBegin(GL_TRIANGLES);
	for (size_t i = 0; i < faces.size(); i++) {
		glColor3ub(shading[i].r, shading[i].g, shading[i].b);
		Point p1 = vertices[faces[i].id1];
		Point p2 = vertices[faces[i].id2];
		Point p3 = vertices[faces[i].id3];
		glVertex3f(p1.x, p1.y, p1.z);
		glVertex3f(p2.x, p2.y, p2.z);
		glVertex3f(p3.x, p3.y, p3.z);
	}
	glEnd();

	glFlush();
	if (!save_images)
		SDL_GL_SwapBuffers();
	glPopMatrix();
	glPopAttrib();
}

void drawTexture() {
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glPushMatrix();
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(cameraX,cameraY,cameraZ,centerX,centerY,centerZ,upX,upY,upZ);

	glBegin(GL_TRIANGLES);
	for (size_t i = 0; i < faces.size(); i++) {
		Point p1 = vertices[faces[i].id1];
		Point p2 = vertices[faces[i].id2];
		Point p3 = vertices[faces[i].id3];
		glTexCoord2d(0.0,0.0);
		glNormal3f(normals[i].x,normals[i].y,normals[i].z);
		glVertex3f(p1.x, p1.y, p1.z);
		glTexCoord2d(1.0,0.0);
		glNormal3f(normals[i].x,normals[i].y,normals[i].z);
		glVertex3f(p2.x, p2.y, p2.z);
		glTexCoord2d(0.5,1.0);
		glNormal3f(normals[i].x,normals[i].y,normals[i].z);
		glVertex3f(p3.x, p3.y, p3.z);
	}
	glEnd();

	glFlush();
	if (!save_images)
		SDL_GL_SwapBuffers();
	glPopMatrix();
	glPopAttrib();
}

int main(int argc, char* argv[]) {
	if (argc < 2) {
		printf("Usage: ./%s input.ply [outputFolder]\n",argv[0]);
		return 1;
	}

	char buffer[128];
	readPLY(argv[1],&vertices,&faces);
	save_images = argc == 3;
	srand(0);

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
#if VERBOSE
	printf("Bounding box: x:(%.2f %.2f) y:(%.2f %.2f) z:(%.2f %.2f)\n",minX,maxX,minY,maxY,minZ,maxZ);
#endif
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

	SDL_Init(SDL_INIT_VIDEO);
	SDL_WM_SetCaption("Render PLY", NULL);
	SDL_SetVideoMode(RESOLUTION,RESOLUTION, 32, SDL_OPENGL);
    glEnable(GL_DEPTH_TEST);
#if NUM_TEXTURE
	glEnable( GL_TEXTURE_2D );
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
#else
	glDisable(GL_LIGHTING);
#endif
	glDisable(GL_CULL_FACE);
	glPolygonMode(GL_FRONT, GL_FILL);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(70,1,1,ZFAR);
	cameraX = maxX - minX;
	cameraY = maxY - minY;
	cameraZ = maxZ - minZ;

	void (*draw)(void);
#if NUM_TEXTURE
	GLuint texture[NUM_TEXTURE];
	glGenTextures( NUM_TEXTURE, texture );
	Image img = {0,0,NULL};
	for (int i=0;i<NUM_TEXTURE;i++) {
		sprintf(buffer,"textures/%d.ppm",i);
		loadPPM(buffer,&img);
		glBindTexture( GL_TEXTURE_2D, texture[i] );
		gluBuild2DMipmaps( GL_TEXTURE_2D, 3, img.width, img.height, GL_RGB, GL_UNSIGNED_BYTE, img.data );
	}
	draw = drawTexture;
#else
	draw = drawNoTexture;
#endif

	if (save_images) {
		sprintf(buffer,"%s/labels.txt",argv[2]);
		FILE* labelFile = fopen(buffer,"r");
		for (int i=0;i<atoi(argv[1])+1;i++)
			fgets(buffer,128,labelFile);
		int classLabel = atoi(buffer);
		fclose(labelFile);
		sprintf(buffer,"%s/factors.txt",argv[2]);
		FILE* factorFile = fopen(buffer,"a");

		float rho = sqrt(cameraX*cameraX + cameraY*cameraY);
		float theta = atan2(cameraY, cameraX);
		for (int k = 0; k < NUM_CAMERAS; k++) {
			int c1 = rand() % (NUM_COLORS-1)+1, c2;
			do {
				c2 = rand() % NUM_COLORS;
			} while (c2 == c1);
			Color objectColor = palette[c1];
			Color bgColor = palette[c2];
			Vector lightDir = {0,0,1};
			getShading(objectColor,lightDir);
			glClearColor(1.0/255*bgColor.r,1.0/255*bgColor.g,1.0/255*bgColor.b,1);
#if NUM_TEXTURE
			int t = rand() % NUM_TEXTURE;
			glBindTexture(GL_TEXTURE_2D,texture[t]);
#endif

			cameraX = rho * cos(theta);
			cameraY = rho * sin(theta);
			theta += 2 * 3.14159265 / NUM_CAMERAS;
			draw();
			int n=0;
			while (true) {
				sprintf(buffer,"%s/%d.ppm",argv[2],n);
				FILE* f = fopen(buffer,"r");
				if (!f) {
					writeImage(buffer);
					break;
				}
				fclose(f);
				n++;
			}
#if NUM_TEXTURE
			fprintf(factorFile,"%d %d %d %d\n",classLabel,k,t,c2);
#else
			fprintf(factorFile,"%d %d %d %d\n",classLabel,k,c1,c2);
#endif
		}
		fclose(factorFile);
	} else {
		Color objectColor = palette[1];
		Color bgColor = palette[0];
		Vector lightDir = {0,0,1};
		getShading(objectColor,lightDir);
		glClearColor(1.0/255*bgColor.r,1.0/255*bgColor.g,1.0/255*bgColor.b,1);

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
	}

	return 0;
}
