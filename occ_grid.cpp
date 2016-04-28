#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <math.h>
#include <unistd.h>
#include <SDL/SDL.h>
#include <GL/gl.h>
#include <GL/glu.h>

#define GRAYSCALE 1

struct Dimensions {
	int length,width,height;
};

double cameraX=20,cameraY=0,cameraZ=40;
double centerX=0,centerY=0,centerZ=0;
double upX=0,upY=0,upZ=1;
int mouseIndex = 0;
int previousX,previousY;
double scrollSpeed = 1.01;
std::vector<float> gridValue;
Dimensions gridDimensions;

void colormap (float f,unsigned char *r,unsigned char *g,unsigned char *b) {
	*r=0;
	*g=0;
	*b=0;
	if (f<=0) {
		*b = 128;
	} else if (f <= 0.25) {
		*g = (unsigned char) f / 0.25 * 255;
		*b = (unsigned char) 128 * (1 - f / 0.25);
	} else if (f <= 0.5) {
		*g = 255;
		*r = (unsigned char) (f - 0.25) / 0.25 * 255;
	} else if (f <= 0.75) {
		*r = 255;
		*g = (unsigned char) 255 + (0.5 - f) / 0.25 * 127;
	} else if (f <= 1) {
		*r = 255;
		*g = (unsigned char) 128 * (1 - f) / 0.25;
	} else {
		*r = 255;
	}
}

void graymap(float f,unsigned char *r,unsigned char *g,unsigned char *b) {
	*r = *g = *b = (unsigned char) (f * 255);
}

void loadOccupancyGrid(char* fileName) {
	FILE* f = fopen(fileName,"r");
	if (!f) {
		printf("File not found: %s\n", fileName);
		return;
	}
	printf("Loaded %s\n",fileName);
	char buffer[128];
	char *c = buffer;
	fgets(buffer,128,f);
	gridDimensions.length = strtol(c,&c,10);
	gridDimensions.width = strtol(c,&c,10);
	gridDimensions.height = strtol(c,&c,10);
	int numPoints = gridDimensions.length * gridDimensions.width * gridDimensions.height;
	gridValue.resize(numPoints);
	float maxValue = 0;
	for (int i=0;i<numPoints;i++) {
		fscanf(f,"%f",gridValue.data()+i);
		if (gridValue[i] > maxValue)
			maxValue = gridValue[i];
	}
	for (int i=0;i<numPoints;i++) {
		gridValue[i] = gridValue[i] / maxValue;
	}
	fclose(f);
}

void drawCube(float x,float y,float z) {
	glBegin(GL_QUADS);
	glVertex3d(x,y,z);
	glVertex3d(x+1,y,z);
	glVertex3d(x+1,y+1,z);
	glVertex3d(x,y+1,z);
	glVertex3d(x,y,z);
	glVertex3d(x+1,y,z);
	glVertex3d(x+1,y,z+1);
	glVertex3d(x,y,z+1);
	glVertex3d(x,y,z);
	glVertex3d(x,y+1,z);
	glVertex3d(x,y+1,z+1);
	glVertex3d(x,y,z+1);

	glVertex3d(x,y,z+1);
	glVertex3d(x+1,y,z+1);
	glVertex3d(x+1,y+1,z+1);
	glVertex3d(x,y+1,z+1);
	glVertex3d(x,y+1,z);
	glVertex3d(x+1,y+1,z);
	glVertex3d(x+1,y+1,z+1);
	glVertex3d(x,y+1,z+1);
	glVertex3d(x+1,y,z);
	glVertex3d(x+1,y+1,z);
	glVertex3d(x+1,y+1,z+1);
	glVertex3d(x+1,y,z+1);
	glEnd();
}

void drawDot(float x,float y,float z) {
	glBegin(GL_POINTS);
	glVertex3d(x,y,z);
	glEnd();
}

void draw() {
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glPushMatrix();
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(cameraX,cameraY,cameraZ,centerX,centerY,centerZ,upX,upY,upZ);

	//draw ground
	int maxDimension = gridDimensions.length;
	if (gridDimensions.width > maxDimension)
		maxDimension = gridDimensions.width;
	if (gridDimensions.height > maxDimension)
		maxDimension = gridDimensions.height;
	maxDimension /= 2;
	glBegin(GL_QUADS);
	glColor3ub(100,100,100);
	glVertex3d(-maxDimension,-maxDimension,0);
	glVertex3d(maxDimension,-maxDimension,0);
	glVertex3d(maxDimension,maxDimension,0);
	glVertex3d(-maxDimension,maxDimension,0);
	glEnd();
	glLineWidth(5.0);
	glBegin(GL_LINES);
	glColor3ub(150,150,150);
	glVertex3d(-maxDimension,-maxDimension,0);
	glVertex3d(maxDimension,-maxDimension,0);
	glColor3ub(150,150,150);
	glVertex3d(-maxDimension,-maxDimension,0);
	glVertex3d(-maxDimension,maxDimension,0);
	glColor3ub(150,150,150);
	glVertex3d(-maxDimension,-maxDimension,0);
	glVertex3d(-maxDimension,-maxDimension,2*maxDimension);
	glEnd();

	//draw object
	unsigned char r,g,b;
	int xl,xr,xd,yl,yr,yd,zl,zr,zd;
	if (cameraX > 0) {
		xl=(1-gridDimensions.length)/2; xr=gridDimensions.length/2+1; xd=1;
	} else {
		xr=(1-gridDimensions.length)/2-1; xl=gridDimensions.length/2; xd=-1;
	}
	if (cameraY > 0) {
		yl=(1-gridDimensions.width)/2; yr=gridDimensions.width/2+1; yd=1;
	} else {
		yr=(1-gridDimensions.width)/2-1; yl=gridDimensions.width/2; yd=-1;
	}
	if (cameraZ > 0) {
		zl=0; zr=gridDimensions.height; zd=1;
	} else {
		zr=-1; zl=gridDimensions.height-1; zd=-1;
	}

	void (*draw3d)(float,float,float);
	if (maxDimension >= 15) {
		draw3d = drawDot;
		glPointSize(2.0);
	} else {
		draw3d = drawCube;
	}
	for (int z=zl;z!=zr;z+=zd) {
		for (int y=yl;y!=yr;y+=yd) {
			for (int x=xl;x!=xr;x+=xd) {
				int xi = x + (gridDimensions.length-1)/2;
				int yi = y + (gridDimensions.width-1)/2;
				int zi = z;
				int gridIndex = xi + yi*gridDimensions.length + zi*gridDimensions.length*gridDimensions.width;
				if (gridValue[gridIndex] > 0) {
#if GRAYSCALE
					graymap(gridValue[gridIndex],&r,&g,&b);
#else
					colormap(gridValue[gridIndex],&r,&g,&b);
#endif
					glColor3ub(r,g,b);
					draw3d(x,y,z);
				}
			}
		}
	}

	glFlush();
	SDL_GL_SwapBuffers();
	glPopMatrix();
	glPopAttrib();
}

int main(int argc,char* argv[]) {

	if (argc < 2) {
		printf("./occ_grid [0.og ... ]\n");
		return 1;
	}

	std::vector<char*> data;
	for (int i=1;i<argc;i++)
		data.push_back(argv[i]);
	int data_index=0;
	loadOccupancyGrid(data[data_index]);

	SDL_Init(SDL_INIT_VIDEO);
	SDL_WM_SetCaption("Point Cloud", NULL);
//	SDL_SetVideoMode(640,480, 32, SDL_OPENGL);
	SDL_SetVideoMode(1600,1200, 32, SDL_OPENGL);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(70,(double)640/480,1,1000);

	int interval = 5000;
	SDL_Event event= {};
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
						case 'n':
						if (data_index < data.size() - 1)
							data_index++;
						loadOccupancyGrid(data[data_index]);
						break;
						case 'b':
						if (data_index > 0)
							data_index--;
						loadOccupancyGrid(data[data_index]);
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
