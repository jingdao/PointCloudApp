#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <math.h>
#include <unistd.h>
#include <SDL/SDL.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include "lapack.h"
#include <ft2build.h>
#include FT_FREETYPE_H
#define IGNORE_ZEROS 1
#define INCLUDE_OUTLIERS 1
#define SHOW_PERCENTAGE 0

//double cameraX=0,cameraY=-5,cameraZ=3;
//double cameraX=245,cameraY=-223,cameraZ=201;
double cameraX=20,cameraY=20,cameraZ=20;
double centerX=0,centerY=0,centerZ=0;
double upX=0,upY=0,upZ=1;
#if SHOW_PERCENTAGE
	const int labelWidth=180;
#else 
	const int labelWidth=152;
#endif
const int labelHeight=30,fontpixels=20,grayLevel=50;
int mouseIndex = 0;
int previousX,previousY;
double scrollSpeed = 1.01;
unsigned char raster[labelWidth * labelHeight * 3];

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

int rChoice[] = {grayLevel,255,0,0,255,255,0};
int gChoice[] = {grayLevel,0,255,0,255,0,255};
int bChoice[] = {grayLevel,0,0,255,0,255,255};
std::vector<int> labels;
std::vector<PCD*> cloud;
PCD* outlier = NULL;
std::vector< std::vector<float> > box;
std::vector<char*> categories;
std::vector<char*> descriptions;
std::vector<float> score;
FT_Library ft;
FT_Face face;
FT_GlyphSlot glyph;
GLuint textures;

void draw();

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

void colormap(float f, unsigned char *r,unsigned char *g,unsigned char *b) {
	*r = 0;
	*g = 0;
	*b = 0;
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

void render_text(const char *text, unsigned char *data) {
	memset(data,0,labelWidth*labelHeight*3);
	const char *p;
	unsigned char color[] = {255,255,255};
	int x = 0;
	for(p = text; *p; p++) {
		if(FT_Load_Char(face, *p, FT_LOAD_RENDER))
			continue;
		unsigned char *src = glyph->bitmap.buffer;
		int k=0;
		for (int i=0;i<glyph->bitmap.rows;i++) {
			unsigned char *dest = data + ((glyph->bitmap_top - i + fontpixels/2) * labelWidth + x + glyph->bitmap_left)* 3;
			for (int j=0;j<glyph->bitmap.width;j++) {
				memset(dest,*src,3); // draw in grayscale
				//*dest = *src; //draw in red
				src++;
				dest+=3;
			}
		}
		x += glyph->bitmap_left + glyph->bitmap.width;
//		x += glyph->bitmap.width;
		if (x >= labelWidth)
			break;
	}
}

void drawLine(int i,int j,int k) {
	glBegin(GL_LINES);
	glVertex3d(box[i][j*3],box[i][j*3+1],box[i][j*3+2]);
	glVertex3d(box[i][k*3],box[i][k*3+1],box[i][k*3+2]);
	glEnd();
}

void drawText(int i) {
	glRasterPos3f(box[i][12],box[i][13],box[i][14]);
	render_text(descriptions[i],raster);
	glDrawPixels(labelWidth,labelHeight,GL_RGB,GL_UNSIGNED_BYTE,raster);
}

void draw() {
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glPushMatrix();

	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

//	glMatrixMode(GL_PROJECTION);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	gluLookAt(cameraX,cameraY,cameraZ,centerX,centerY,centerZ,upX,upY,upZ);

	glPointSize(1.0);
	glBegin(GL_POINTS);
	if (outlier) {
		glColor3ub(grayLevel,grayLevel,grayLevel);
		for (int n=0;n<outlier->numPoints;n++)
			glVertex3d(outlier->float_data[n*4],outlier->float_data[n*4+1],outlier->float_data[n*4+2]);
	}
	for (int i=0;i<cloud.size();i++) {
		if (!labels[i])
			glColor3ub(150,150,150);
		else
			glColor3ub(rChoice[labels[i]],gChoice[labels[i]],bChoice[labels[i]]);
		for (int n = 0; n < cloud[i]->numPoints; n++){
			glVertex3d(cloud[i]->float_data[n*4],cloud[i]->float_data[n*4+1],cloud[i]->float_data[n*4+2]);
		}
	}
	glEnd();

	glLineWidth(1.0);
	glColor3ub(255,255,255);
	for (int i=0;i<box.size();i++) {
		drawLine(i,0,1);
		drawLine(i,0,2);
		drawLine(i,1,3);
		drawLine(i,2,3);
		drawLine(i,0,4);
		drawLine(i,1,5);
		drawLine(i,2,6);
		drawLine(i,3,7);
		drawLine(i,4,5);
		drawLine(i,4,6);
		drawLine(i,5,7);
		drawLine(i,6,7);
#if IGNORE_ZEROS
		drawText(i);
#endif
	}

	glFlush();
	SDL_GL_SwapBuffers();

	glPopMatrix();
	glPopAttrib();
}

void centerMap() {
	int p=0;
	float cx=0,cy=0,cz=0;
	if (outlier) {
		for (int n=0;n<outlier->numPoints;n++) {
			p++;
			cx += outlier->float_data[n*4];
			cy += outlier->float_data[n*4+1];
			cz += outlier->float_data[n*4+2];
		}
	}
	for (int i=0;i<cloud.size();i++) {
		for (int n = 0; n < cloud[i]->numPoints; n++){
			p++;
			cx += cloud[i]->float_data[n*4];
			cy += cloud[i]->float_data[n*4+1];
			cz += cloud[i]->float_data[n*4+2];
		}
	}
	cx /= p;
	cy /= p;
	cz /= p;
	if (outlier) {
		for (int n=0;n<outlier->numPoints;n++) {
			outlier->float_data[n*4] -= cx;
			outlier->float_data[n*4+1] -= cy;
			outlier->float_data[n*4+2] -= cz;
		}
	}
	for (int i=0;i<cloud.size();i++) {
		for (int n = 0; n < cloud[i]->numPoints; n++){
			cloud[i]->float_data[n*4] -= cx;
			cloud[i]->float_data[n*4+1] -= cy;
			cloud[i]->float_data[n*4+2] -= cz;
		}
	}
	for (int i=0;i<box.size();i++) {
		for (int j=0;j<8;j++) {
			box[i][j*3] -= cx;
			box[i][j*3+1] -= cy;
			box[i][j*3+2] -= cz;
		}
	}
}

int main(int argc,char* argv[]) {
	if (argc < 2) {
		printf("./semantics clusterFolder\n");
		return 1;
	}

	char buf[128];
	char categories_buf[256];
	char desc_buf[1024];
	sprintf(buf,"%s/../labelCategory.txt",argv[1]);
	FILE* categoryFile = fopen(buf,"r");
	if (!categoryFile)
		return 1;
	char* c = categories_buf;
	while (fgets(buf,128,categoryFile)) {
		int l = strlen(buf);
		strncpy(c,buf+3,l-4);
		categories.push_back(c);
		c += l-4;
		*c++ = 0;
	}
	fclose(categoryFile);
	sprintf(buf,"%s/prediction.txt",argv[1]);
	FILE* labelFile = fopen(buf,"r");
	if (!labelFile)
		return 1;
	c = buf;
	char* d = desc_buf;
	while (fgets(buf,128,labelFile)) {
		int l = strtol(buf,&c,10);
		labels.push_back(l);
#if IGNORE_ZEROS
		if (l != 0) {
#else
		if (l >= 0) {
#endif
			double max_score = 0;
			while (*c != '\n') {
				double score = strtod(c,&c);
				max_score = score > max_score ? score : max_score;
			}
#if SHOW_PERCENTAGE
			int k = sprintf(d,"%s:%.0f%%",categories[l],max_score*100);
#else
			int k = sprintf(d,"%s",categories[l]);
#endif
			score.push_back(max_score);
			descriptions.push_back(d);
			d += k + 1;
		}
	}
	fclose(labelFile);

	for (int i=0;i<labels.size();i++) {
		sprintf(buf,"%s/%d-cloud.pcd",argv[1],i);
		PCD* c = NewPCD(buf);
		cloud.push_back(c);
#if IGNORE_ZEROS
		if (labels[i] != 0) {
#else
		if (labels[i] >= 0) {
#endif
			std::vector<float> currentBox;
			double lambda[3], v[9];
			double principalLengths[3], principalAxes[9], bbCenter[3];
			getPCA_XY(c,lambda,v);
			setPCA(c,lambda,v,principalLengths,principalAxes,bbCenter);
			for (int i=0;i<8;i++) {
				float coords[3];
				for (int j=0;j<3;j++) {
					coords[j] = bbCenter[j];
					for (int axis=0;axis<3;axis++) {
						float sign = (i & 1<<axis) ? 1 : -1;
						coords[j] += sign * principalLengths[axis] / 2 * principalAxes[axis*3 + j];
					}
					currentBox.push_back(coords[j]);
				}
			}
			box.push_back(currentBox);
		}
	}
	printf("Loaded %lu point clouds\n",cloud.size());
	printf("Loaded %lu bounding boxes\n",box.size());
#if INCLUDE_OUTLIERS
	sprintf(buf,"%s/outlier.pcd",argv[1]);
	outlier = NewPCD(buf);
#endif
	centerMap();

	FT_Init_FreeType(&ft);
	FT_New_Face(ft,"/usr/share/fonts/truetype/freefont/FreeSans.ttf",0,&face);
	FT_Set_Pixel_Sizes(face,fontpixels,fontpixels);
	glyph = face->glyph;

	SDL_Init(SDL_INIT_VIDEO);
	SDL_WM_SetCaption("Point Cloud", NULL);
	SDL_SetVideoMode(1800,1000, 32, SDL_OPENGL);
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
					/*	double rho = sqrt(cameraX*cameraX+cameraY*cameraY+cameraZ*cameraZ);
						double r = sqrt(cameraX*cameraX+cameraY*cameraY);
						double phi = asin(cameraZ/rho);
						double theta = asin(cameraY/r);
						theta += 0.01 * (event.motion.x-previousX);
						phi += 0.01 * (event.motion.y-previousY);
						cameraX = rho * cos(theta) * cos(phi);
						cameraY = rho * sin(theta) * cos(phi);
						cameraZ = rho * sin(phi);
						previousX = event.motion.x;
						previousY = event.motion.y;
						draw();
					} else if (mouseIndex == 2) {*/
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
