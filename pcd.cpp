#include "pcd.h"

PCD::PCD(const char* fileName) {
	numFields = 4;
	data = NULL;
	float_data = NULL;
	fields = NULL;
	numPoints = 0;
	capacity = 0;
	data_storage = NONE;
	kdtreeLeaves = NULL;
	kdtreeDepth = 0;
	FILE* f = fopen(fileName, "r");
	if (!f) {
		printf("File not found: %s\n", fileName);
		return;
	}
	char buf[256];
	int pointsParsed = 0;
	while (fgets(buf, 256, f)) {
		if (sscanf(buf, "POINTS %d", &numPoints) == 1) {
			capacity = numPoints;
			float_data = (float*)malloc(4 * numPoints * sizeof(float));
		}
		else if (sscanf(buf, "VIEWPOINT %f %f %f %f %f %f %f", &viewpoint[0], &viewpoint[1], &viewpoint[2],
			&viewpoint[3], &viewpoint[4], &viewpoint[5], &viewpoint[6]) == 7) {
			camera.position.r = 0;
			camera.focal_length = DEFAULT_FOCAL_LENGTH;
			camera.position.i = (float)viewpoint[0];
			camera.position.j = (float)viewpoint[1];
			camera.position.k = (float)viewpoint[2];
			camera.rotation.r = (float)viewpoint[3];
			camera.rotation.i = (float)viewpoint[4];
			camera.rotation.j = (float)viewpoint[5];
			camera.rotation.k = (float)viewpoint[6];
		} else if (strncmp(buf,"DATA ascii",10)==0) {
			data_storage = ASCII;
		} else if (strncmp(buf,"DATA binary_compressed",23)==0) {
			data_storage = BINARY;
			fread(float_data,sizeof(float),numPoints*4,f);
			break;
		}
		else if (data_storage == ASCII) {
			if (sscanf(buf, "%f %f %f %f", &float_data[pointsParsed * 4], &float_data[pointsParsed * 4 + 1],
				&float_data[pointsParsed * 4 + 2], &float_data[pointsParsed * 4 + 3]) == 4) {
				pointsParsed++;
			}
		}
	}
	fclose(f);
}

PCD::PCD(int size) {
	numFields = 4;
	data = NULL;
	float_data = (float*) malloc(size*4*sizeof(float));
	fields = NULL;
	numPoints = size;
	capacity = size;
	data_storage = ASCII;
	kdtreeLeaves = NULL;
	kdtreeDepth = 0;
}

PCD::~PCD() {
	if (fields) free(fields);
	if (data) free(data);
	if (float_data) free(float_data);
	if (kdtreeLeaves) free(kdtreeLeaves);
	for (size_t i=0;i<kdtreeBranches.size();i++)
		free(kdtreeBranches[i]);
}

bool PCD::expand(int size) {
	if (size < 0)
		return false;
	numPoints = size;
	if (size > capacity) {
		if (capacity == 0)
			capacity = size;
		else while (capacity < size)
			capacity *= 2;
		void* memblk  = realloc(float_data, capacity * 4 * sizeof(float));
		if (memblk) {
			float_data = (float*) memblk;
		} else {
			printf("Error: Cannot realloc %lu bytes\n",capacity * 4 * sizeof(float));
			return false;
		}
	}
	return true;
}

void PCD::drawLine(float x1,float y1,float z1,float x2,float y2,float z2,float rgb,float interval) {
	float distance = 0;
	int offset = numPoints * 4,numSteps;
	distance += (x2-x1) * (x2-x1);
	distance += (y2-y1) * (y2-y1);
	distance += (z2-z1) * (z2-z1);
	numSteps = (int) (distance/interval/interval);
	expand(numPoints + numSteps - 1);
	float* data = float_data + offset;

	float dx = (x2-x1) / numSteps;
	float dy = (y2-y1) / numSteps;
	float dz = (z2-z1) / numSteps;
	float x = x1+dx, y = y1+dy, z = z1+dz;
	for (int i=1;i<numSteps;i++) {
		*data++ = x;
		*data++ = y;
		*data++ = z;
		*data++ = rgb;
		x += dx;
		y += dy;
		z += dz;
	}
}

void PCD::drawAxis(float x,float y,float z,float* R,float length,float interval) {
	float I[9] = {1,0,0,0,1,0,0,0,1}; //identity matrix
	float *rotation = R ? R : I;
	float rgb[3] = {255<<16, 255<<8, 255};

	for (int i=0;i<3;i++) { // three axes
		float x2 = x + rotation[i] * length;
		float y2 = y + rotation[i + 3] * length;
		float z2 = z + rotation[i + 6] * length;
		drawLine(x,y,z,x2,y2,z2,rgb[i],interval);
	}
}

void PCD::loadDescriptor(const char* filename) {
	FILE* f = fopen(filename,"r");
	if (!f) {
		printf("File not found: %s\n", filename);
		return;
	}
	char buf[256];
	int numDescriptors=0,count=0;
	while (fgets(buf, 256, f)) {
		if (sscanf(buf, "POINTS %d", &numDescriptors) == 1) {
			if (numDescriptors!=numPoints) {
				printf("Number of descriptors does not match number of points!\n");
				return;
			}
		} else if (sscanf(buf,"COUNT %d",&count)==1) {
		} else if (strncmp(buf,"DATA ascii",10)==0) 
			break;
	}
	float val,sum,maxSum=0;
	char* c;
	float* values = new float[numDescriptors];
	for (int i=0;i<numDescriptors;i++) {
		fgets(buf,256,f);
		sum=0;
		c = buf;
		for (int j=0;j<count;j++) {
			val = strtod(c,&c);
			sum += val * val;
		}
		values[i] = sum;
		if (sum > maxSum) maxSum = sum;
	}
	for (int i=0;i<numDescriptors;i++) {
		values[i] /= maxSum;
		float_data[i * 4 + 3] = colormap(values[i]);
	}
}

void PCD::writeToPCD(const char* filename) {
	if (!filename)
		return;
	FILE* f = fopen(filename, "w");
	if (!f) {
		printf("Cannot write to file: %s\n", filename);
		return;
	}
	fprintf(f,"# .PCD v0.7 - Point Cloud Data file format\n"
	"VERSION 0.7\n"
	"FIELDS x y z rgb\n"
	"SIZE 4 4 4 4\n"
	"TYPE F F F I\n"
	"COUNT 1 1 1 1\n"
	"WIDTH %d\n"
	"HEIGHT 1\n"
	"VIEWPOINT 0 0 0 1 0 0 0\n"
	"POINTS %d\n"
	"DATA ascii\n",numPoints,numPoints);
	int k=0;
	float x,y,z,rgb;
	for (int i=0;i<numPoints;i++) {
		x=float_data[k++];
		y=float_data[k++];
		z=float_data[k++];
		rgb=float_data[k++];
//		fprintf(f,"%f %f %f %f\n",x,y,z,rgb);
		fprintf(f,"%f %f %f %d\n",x,y,z,(int)rgb);
	}
	fclose(f);
	printf("Wrote %d points to %s\n",numPoints,filename);
}

void PCD::writeToPLY(const char* filename) {
	if (!filename)
		return;
	FILE* f = fopen(filename, "w");
	if (!f) {
		printf("Cannot write to file: %s\n", filename);
		return;
	}
	fprintf(f,"ply\n"
	"format ascii 1.0\n"
	"element face 0\n"
	"property list uchar int vertex_indices\n"
	"element vertex %d\n"
	"property float x\n"
	"property float y\n"
	"property float z\n"
	"property uchar diffuse_red\n"
	"property uchar diffuse_green\n"
	"property uchar diffuse_blue\n"
	"end_header\n",numPoints);
	int k=0;
	float x,y,z;
	int rgb;
	for (int i=0;i<numPoints;i++) {
		x=float_data[k++];
		y=float_data[k++];
		z=float_data[k++];
		rgb=(int)float_data[k++];
		fprintf(f,"%f %f %f %d %d %d\n",x,y,z,(rgb>>16)&0xff,(rgb>>8)&0xff,rgb&0xff);
	}
	fclose(f);
	printf("Wrote %d points to %s\n",numPoints,filename);
}

void PCD::writeClustersToPCD(std::vector<std::vector<int>> *indices,const char* dir) {
	char buffer[256];
	for (size_t i=0;i<indices->size();i++) {
		snprintf(buffer,256,"%s/%lu.pcd",dir,i);
		std::vector<int> currentIndices = (*indices)[i];
		PCD* subcloud = extractIndices(&currentIndices);
		subcloud->writeToPCD(buffer);
		delete subcloud;
	}
}

PCD* PCD::LoadFromMatrix(const char* fileName) {
	FILE* f = fopen(fileName, "r");
	if (!f) {
		printf("File not found: %s\n", fileName);
		return NULL;
	}

	int dimensions = -1;
	std::vector<float> arr;
	char buf[256];
	while (fgets(buf,255,f)) {
		char* c = buf;
		while (*c) {
			if (*c<='9'&&*c>='0') {
				float f = (float) strtod(c,&c);
				arr.push_back(f);
			} else if (dimensions<0 && *c=='}') {
				dimensions = arr.size();
			} else if (*c=='}'&&*(c+1)=='}') {
				break;
			} else c++;
		}
		if (*c=='}'&&*(c+1)=='}')
			break;
	}

	fclose(f);
	int numPoints = arr.size()/dimensions;
	PCD* p = new PCD(numPoints);
	int i,j;
	for (i=0;i<numPoints;i++) {
		for (j=0;j<dimensions;j++) {
			p->float_data[i*4+j] = arr[i*dimensions+j];
		}
		p->float_data[i*4+3] = 255.0f;
	}

	return p;
}

PCD* PCD::LoadFromBundler(const char* fileName) {
	FILE* f = fopen(fileName,"r");
	if (!f) {
		printf("Input file not found\n");
		return 0;
	}

	const int buffer_size = 1024;
	int num_cameras,num_points;
	char buf[buffer_size], *ptr;
	float focal_length,distortion1,distortion2;
	float rotation[9];
	float translation[3];
	float x,y,z;
	float rgb_f;
	int r,g,b,rgb;
	float* camera_positions, *camera_rotations;

	fgets(buf,buffer_size-1,f); //header
	fgets(buf,buffer_size-1,f); //num_cameras, num_points
	ptr=buf;
	num_cameras = strtol(ptr,&ptr,10);
	num_points = strtol(ptr,&ptr,10);
	camera_positions = new float[num_cameras * 3];
	camera_rotations = new float[num_cameras * 9];
	for (int i=0;i<num_cameras;i++) {
		fgets(buf,buffer_size-1,f); //read camera information
		if (!sscanf(buf, "%f %f %f",&focal_length,&distortion1,&distortion2) == 3) {
			printf("Error parsing %s\n",fileName);
			return NULL;
		}
		fgets(buf,buffer_size-1,f);
		if (!sscanf(buf, "%f %f %f",rotation,rotation+1,rotation+2) == 3) {
			printf("Error parsing %s\n",fileName);
			return NULL;
		}
		fgets(buf,buffer_size-1,f);
		if (!sscanf(buf, "%f %f %f",rotation+3,rotation+4,rotation+5) == 3) {
			printf("Error parsing %s\n",fileName);
			return NULL;
		}
		fgets(buf,buffer_size-1,f);
		if (!sscanf(buf, "%f %f %f",rotation+6,rotation+7,rotation+8) == 3) {
			printf("Error parsing %s\n",fileName);
			return NULL;
		}
		fgets(buf,buffer_size-1,f);
		if (!sscanf(buf, "%f %f %f",translation,translation+1,translation+2) == 3) {
			printf("Error parsing %s\n",fileName);
			return NULL;
		}
		memcpy(camera_positions + i * 3, translation, 3 * sizeof(float));
		memcpy(camera_rotations + i * 9, rotation, 9 * sizeof(float));
	}

	PCD* pointcloud = new PCD(num_points);
	float* data = pointcloud->float_data;
	std::vector<float> *lineInfo = new std::vector<float>[num_cameras];
	for (int i = 0; i < num_points; ++i) {
		fgets(buf,buffer_size-1,f); //read X,Y,Z position
		if (!sscanf(buf, "%f %f %f",&x,&y,&z) == 3) {
			printf("Error parsing %s\n",fileName);
			return NULL;
		}
		fgets(buf,buffer_size-1,f); //read R,G,B color
		if (!sscanf(buf, "%d %d %d",&r,&g,&b) == 3) {
			printf("Error parsing %s\n",fileName);
			return NULL;
		}
		rgb = ((r<<16)|(g<<8)|b);
		fgets(buf,buffer_size-1,f); //read view list
		ptr = buf;
		int numViews = strtol(ptr,&ptr,10);
		int cameraIndex;
//		int keyIndex;
//		float viewX,viewY;
		for (int i=0;i<numViews;i++) {
			cameraIndex = strtol(ptr,&ptr,10);
//			keyIndex = strtol(ptr,&ptr,10);
//			viewX = strtod(ptr,&ptr);
//			viewY = strtod(ptr,&ptr);
			strtol(ptr,&ptr,10);
			strtol(ptr,&ptr,10);
			strtol(ptr,&ptr,10);
			lineInfo[cameraIndex].push_back(x);
			lineInfo[cameraIndex].push_back(y);
			lineInfo[cameraIndex].push_back(z);
			lineInfo[cameraIndex].push_back(camera_positions[cameraIndex * 3]);
			lineInfo[cameraIndex].push_back(camera_positions[cameraIndex * 3 + 1]);
			lineInfo[cameraIndex].push_back(camera_positions[cameraIndex * 3 + 2]);
			lineInfo[cameraIndex].push_back(colormap(1.0*cameraIndex/num_cameras));
		}
		rgb_f = (float) rgb;
		*data++ = x;
		*data++ = y;
		*data++ = z;
		*data++ = rgb_f;
	}

	//draw lines from camera to keypoints on point cloud
	PCD* lines = new PCD(0);
	const float interval = 2;
	const int maxLines = 10;
	for (int j=0;j<num_cameras;j++) {
		int numLines = lineInfo[j].size() / 7;
		for (int i=0;i<numLines&&i<maxLines;i++) {
			lines->drawLine(
				lineInfo[j][i * 7],
				lineInfo[j][i * 7 + 1],
				lineInfo[j][i * 7 + 2],
				lineInfo[j][i * 7 + 3],
				lineInfo[j][i * 7 + 4],
				lineInfo[j][i * 7 + 5],
				lineInfo[j][i * 7 + 6],
				interval
			);
		}
	}
	lines->drawAxis(0,0,0,NULL,10,1);
	printf("Wrote %d points\n",lines->numPoints);
	lines->writeToPCD("p2.pcd");
	delete lines;

	fclose(f);
	delete[] camera_positions;
	delete[] camera_rotations;
	delete[] lineInfo;
	return pointcloud;
}

PCD* PCD::LoadFromKITTI(const char* fileName) {
	FILE* f = fopen(fileName,"r");
	if (!f) {
		printf("File not found: %s\n", fileName);
		return NULL;
	}

	fseek (f, 0, SEEK_END);   // non-portable
	int size_bytes=ftell(f);
	int size = size_bytes / 4 / sizeof(float);

	PCD* p = new PCD(size);
	fseek(f,0,SEEK_SET);
	fread(p->float_data,sizeof(float),size*4,f);
	
	//map reflectance to color
	for (int i=0;i<size;i++) {
		p->float_data[i * 4 + 3] = colormap(p->float_data[i * 4 + 3]);
//		p->float_data[i * 4 + 3] = 255 << 16;
//		p->float_data[i * 4 + 3] = (int) (255 * p->float_data[i * 4 + 3]) << 16;
	}

	fclose (f);
	return p;
}

PCD* PCD::LoadFromCameraParameter(const char* dir) {

	DIR* d = opendir(dir);
	if (!d) {
		printf("Invalid directory %s\n",dir);
		return NULL;
	}
	struct dirent *dt;
	int numRead=0; 
	PCD* cloud = new PCD(0);
	char buffer[1024];

	int ndir = strlen(dir);
	strncpy(buffer,dir,ndir);
	char* buffer_c = buffer + ndir;
	*buffer_c++ = '/';

	double focal_length = 600, w = 640, h = 480;
//	double focal_length = 1801, w = 1824, h =2736;
	w = (1 - w) / 2;
	h = (1 - h) / 2;

	while ((dt = readdir(d))) {
		int l = strlen(dt->d_name);
		if (strncmp(dt->d_name + l - 4,".txt",4) == 0) {
			strcpy(buffer_c,dt->d_name);
			FILE* f = fopen(buffer,"r");
			if (f) {
				char line[256];
				double projection[12], *p = projection, translation[3];
				float x,y,z,tmp,rotation[9];

				fgets(line,255,f);
				if (strncmp(line,"CONTOUR",7)!=0) 
					break;
				while (fgets(line,255,f)) {
					if (sscanf(line,"%lf %lf %lf %lf",p,p+1,p+2,p+3) == 4) {
						p += 4;
					}
				}
				translation[2] = -projection[11];
				translation[0] = (projection[3] - translation[2] * w) / focal_length;
				translation[1] = (translation[2] * h - projection[7]) / focal_length;
				rotation[6] = -projection[8];
				rotation[7] = -projection[9];
				rotation[8] = -projection[10];
				rotation[0] = (projection[0] - rotation[6] * w) / focal_length;
				rotation[1] = (projection[1] - rotation[7] * w) / focal_length;
				rotation[2] = (projection[2] - rotation[8] * w) / focal_length;
				rotation[3] = (rotation[6] * h - projection[4]) / focal_length;
				rotation[4] = (rotation[7] * h - projection[5]) / focal_length;
				rotation[5] = (rotation[8] * h - projection[6]) / focal_length;
				//transpose R
				tmp = rotation[1]; rotation[1] = rotation[3]; rotation[3] = tmp;
				tmp = rotation[2]; rotation[2] = rotation[6]; rotation[6] = tmp;
				tmp = rotation[5]; rotation[5] = rotation[7]; rotation[7] = tmp;
				// C = - R' * T
				x = - (rotation[0]*translation[0]+rotation[1]*translation[1]+rotation[2]*translation[2]);
				y = - (rotation[3]*translation[0]+rotation[4]*translation[1]+rotation[5]*translation[2]);
				z = - (rotation[6]*translation[0]+rotation[7]*translation[1]+rotation[8]*translation[2]);
				cloud->drawAxis(x,y,z,rotation,3,1);
				numRead++;
				fclose(f);
			}
		}
	}
	printf("Loaded %d camera parameters\n",numRead);
	closedir(d);

	cloud->drawAxis(0,0,0,NULL,10,1);
	return cloud;
}

PCD* PCD::LoadFromPTS(const char* fileName) {
	FILE* f = fopen(fileName, "r");
	if (!f) {
		printf("File not found: %s\n", fileName);
		return NULL;
	}
	char buf[256];
	int pointsParsed = 0;
	PCD* pcd = new PCD(0);
	int r,g,b,d,rgb;
	while (fgets(buf, 256, f)) {
		int n = atoi(buf);
		printf("n = %d\n",n);
		if (!pcd->expand(n + pointsParsed))
			break;
		for (int i=0;i<n;i++) {
			fgets(buf,256,f);
			if (sscanf(buf, "%f %f %f %d %d %d %d", pcd->float_data + pointsParsed * 4, pcd->float_data + pointsParsed * 4 + 1,
				pcd->float_data + pointsParsed * 4 + 2, &d, &r, &g, &b) == 7) {
				rgb = (r<<16) | (g<<8) | b;
				pcd->float_data[pointsParsed * 4 + 3] = (float) rgb;
				pointsParsed++;
			} else {
				printf("Error parsing %s\n",fileName);
				printf("Line %d: %s\n",pointsParsed,buf);
				break;
			}
		}
	}
	fclose(f);
	return pcd;
}

float PCD::colormap(float f) {
	unsigned char r=0,g=0,b=0;
	if (f<=0) {
		b = 128;
	} else if (f <= 0.25) {
		g = (unsigned char) f / 0.25 * 255;
		b = (unsigned char) 128 * (1 - f / 0.25);
	} else if (f <= 0.5) {
		g = 255;
		r = (unsigned char) (f - 0.25) / 0.25 * 255;
	} else if (f <= 0.75) {
		r = 255;
		g = (unsigned char) 255 + (0.5 - f) / 0.25 * 127;
	} else if (f <= 1) {
		r = 255;
		g = (unsigned char) 128 * (1 - f) / 0.25;
	} else {
		r = 255;
	}
	int res = (r<<16)|(g<<8)|b;
	return (float) res;
}

PCD::Plane PCD::segmentPlane(int iter,float threshold,float inlierRatio) {
	Plane bestPlane = {0,0,0,0};
	int maxInliers = 0;
	int optimumInliers = inlierRatio * numPoints;
	for (int i=0;i<iter;i++) {
		//Pick 3 points
		int p0 = rand() % numPoints;
		int p1 = rand() % numPoints;
		int p2 = rand() % numPoints;
		float x[3] = {float_data[p0*4], float_data[p1*4], float_data[p2*4]};
		float y[3] = {float_data[p0*4+1], float_data[p1*4+1], float_data[p2*4+1]};
		float z[3] = {float_data[p0*4+2], float_data[p1*4+2], float_data[p2*4+2]};
		Plane currentPlane = {
			(y[1]-y[0])*(z[2]-z[0]) - (y[2]-y[0])*(z[1]-z[0]),
			(z[1]-z[0])*(x[2]-x[0]) - (z[2]-z[0])*(x[1]-x[0]),
			(x[1]-x[0])*(y[2]-y[0]) - (x[2]-x[0])*(y[1]-y[0]),
			0
		};
		currentPlane.d = -(currentPlane.a * x[0] + currentPlane.b * y[0] + currentPlane.c * z[0]);
		if (currentPlane.a == 0 && currentPlane.b == 0 && currentPlane.c ==0 )
			continue; //picked collinear points
		float distanceThreshold = threshold * sqrt(
			currentPlane.a * currentPlane.a + 
			currentPlane.b * currentPlane.b + 
			currentPlane.c * currentPlane.c
			);
		int numInliers = 0;
		for (int j=0;j<numPoints;j++) {
			if ( fabs( currentPlane.a * float_data[j * 4] +
				 currentPlane.b * float_data[j * 4 + 1] +
				 currentPlane.c * float_data[j * 4 + 2] +
				 currentPlane.d )
				 < distanceThreshold)	
				numInliers++;
		}
		if (numInliers > maxInliers) {
			maxInliers = numInliers;
			bestPlane = currentPlane;
			if (maxInliers > optimumInliers)
				break;
		}
//		printf("Current plane: %f %f %f %f %d\n",currentPlane.a,currentPlane.b,currentPlane.c,currentPlane.d,numInliers);
	}
	return bestPlane;
}

void PCD::filterPlane(std::vector<int> *ind,Plane p, float threshold,bool positiveOrNegative) {
	float distanceThreshold = threshold * sqrt( p.a * p.a + p.b * p.b + p.c * p.c);
	if (positiveOrNegative) {
		for (int j=0;j<numPoints;j++) {
			if ( fabs( p.a * float_data[j * 4] +
				 p.b * float_data[j * 4 + 1] +
				 p.c * float_data[j * 4 + 2] +
				 p.d )
				 < distanceThreshold)	
				ind->push_back(j);
		}
	} else {
		for (int j=0;j<numPoints;j++) {
			if ( fabs( p.a * float_data[j * 4] +
				 p.b * float_data[j * 4 + 1] +
				 p.c * float_data[j * 4 + 2] +
				 p.d )
				 > distanceThreshold)	
				ind->push_back(j);
		}
	}
}

PCD* PCD::extractIndices(std::vector<int> *ind) {
	PCD* p = new PCD(ind->size());
	for (size_t i=0;i<ind->size();i++) {
		p->float_data[i * 4] = float_data[(*ind)[i] * 4];
		p->float_data[i * 4 + 1] = float_data[(*ind)[i] * 4 + 1];
		p->float_data[i * 4 + 2] = float_data[(*ind)[i] * 4 + 2];
		p->float_data[i * 4 + 3] = float_data[(*ind)[i] * 4 + 3];
	}
	return p;
}

PCD::KdTreeNode* PCD::buildKdTreeNode(int* pointList, int n, int depth) {
	if (depth > kdtreeDepth)
		kdtreeDepth = depth;
	if (n == 1) //only one point
		return kdtreeLeaves + *pointList; 
	int dimension = depth % 3;

	//find median element
	int i1 = rand() % n, i2 = rand() % n, i3 = rand() % n;
	int medianIndex;
	float f1 = float_data[pointList[i1] * 4 + dimension];
	float f2 = float_data[pointList[i2] * 4 + dimension];
	float f3 = float_data[pointList[i3] * 4 + dimension];
	if (f1 >= f2 && f1 <= f3) medianIndex = i1;
	else if (f2 >= f1 && f2 <= f3) medianIndex = i2;
	else medianIndex = i3;
//	int medianIndex = rand() % n;
	int median = pointList[medianIndex];
	float medianValue = float_data[median * 4 + dimension];
	KdTreeNode* currentNode = (KdTreeNode*) calloc(1,sizeof(KdTreeNode));
	currentNode->value = medianValue;
	kdtreeBranches.push_back(currentNode);

	//partition around median
	int leftSize, *leftPoints = pointList;
	int rightSize, *rightPoints;
	int tmp;
	pointList[medianIndex] = pointList[n-1];
	int i=-1,j=-1;
	do {
		j++;
		if (float_data[pointList[j] * 4 + dimension] <= medianValue) {
			i++;
			tmp = pointList[i];
			pointList[i] = pointList[j];
			pointList[j] = tmp;
		}
	} while (j < n - 2);
	pointList[n-1] = pointList[i+1];
	pointList[i+1] = median;
	if (i + 1 == n - 1) //median is last element
		leftSize = i + 1;
	else 
		leftSize = i + 2;
	rightSize = n - leftSize;
	rightPoints = pointList + leftSize;
//	printf("n %d i,j %d %d left %d right %d median %d\n",n,i,j,leftSize,rightSize,medianIndex);

	//get left and right children
	if (leftSize > 0)
		currentNode->left = buildKdTreeNode(leftPoints,leftSize,depth+1);
	if (rightSize > 0)
		currentNode->right = buildKdTreeNode(rightPoints,rightSize,depth+1);

//	printf("depth %d left %d right %d median %d %f\n",depth,leftSize,rightSize,medianIndex,medianValue);
	return currentNode;
}

void PCD::buildKdTree() {
	if (kdtreeLeaves || numPoints == 0)
		return;
	kdtreeLeaves = (KdTreeNode*) calloc(numPoints, sizeof(KdTreeNode));
	int* pointList = (int*) malloc(numPoints * sizeof(int));
	for (int i=0;i<numPoints;i++) {
		pointList[i] = i;
		kdtreeLeaves[i].index = i;
		kdtreeLeaves[i].isLeaf = true;
	}
	kdtreeDepth = 0;
	kdtreeRoot = buildKdTreeNode(pointList,numPoints,0);
	free(pointList);
}

void PCD::searchKdTree(std::vector<int> *indices,float x,float y,float z,float distance) {
	buildKdTree();
	Cube range = {x-distance,y-distance,z-distance,x+distance,y+distance,z+distance};
	boundingBox = getBoundingBox();
	searchKdTreeNode(kdtreeRoot,indices,&range,boundingBox,0);
}

void PCD::searchKdTreeNode(KdTreeNode* node,std::vector<int> *indices,Cube* range,Cube cell,int depth){
	if (!node) return; //nothing here
	if (node->isLeaf) {
		float x = float_data[node->index * 4];
		float y = float_data[node->index * 4 + 1];
		float z = float_data[node->index * 4 + 2];
		if (x >= range->x1 && x <= range->x2 && 
			y >= range->y1 && y <= range->y2 && 
			z >= range->z1 && z <= range->z2)
			indices->push_back(node->index);
	} else {
		if (cell.x1 >= range->x1 && cell.x2 <= range->x2 &&
			cell.y1 >= range->y1 && cell.y2 <= range->y2 &&
			cell.z1 >= range->z1 && cell.z2 <= range->z2) {
			//region fully contained in range
			getIndexFromSubTree(node->left,indices);
			getIndexFromSubTree(node->right,indices);
		} else if ((cell.x2 >= range->x1 && cell.x1 <= range->x2) &&
			(cell.y2 >= range->y1 && cell.y1 <= range->y2) &&
			(cell.z2 >= range->z1 && cell.z1 <= range->z2)) {
			//region intersects range
			int dimension = depth % 3;
			Cube leftCell = cell;
			Cube rightCell = cell;
			//trim cells according to splitting plane
			switch (dimension) {
				case 0:
					leftCell.x2 = node->value;
					rightCell.x1 = node->value;
					break;
				case 1:
					leftCell.y2 = node->value;
					rightCell.y1 = node->value;
					break;
				case 2:
					leftCell.z2 = node->value;
					rightCell.z1 = node->value;
					break;
			}
			searchKdTreeNode(node->left,indices,range,leftCell,depth+1);
			searchKdTreeNode(node->right,indices,range,rightCell,depth+1);
		}
//		else {
//			printf("rejected\n");
//			printf("cell: %f %f %f %f %f %f\n",cell.x1,cell.y1,cell.z1,cell.x2,cell.y2,cell.z2);
//			printf("range: %f %f %f %f %f %f\n",range->x1,range->y1,range->z1,range->x2,range->y2,range->z2);
//		}
	}
}

void PCD::getIndexFromSubTree(KdTreeNode* node, std::vector<int> *indices) {
	if (!node) return;
	if (node->isLeaf)
		indices->push_back(node->index);
	else {
		getIndexFromSubTree(node->left,indices);
		getIndexFromSubTree(node->right,indices);
	}
}

PCD::Cube PCD::getBoundingBox() {
	float minX=float_data[0],maxX=float_data[0];
	float minY=float_data[1],maxY=float_data[1];
	float minZ=float_data[2],maxZ=float_data[2];
	for (int i=1;i<numPoints;i++) {
		if (float_data[i * 4] < minX) minX = float_data[i * 4];
		else if (float_data[i * 4] > maxX) maxX = float_data[i * 4];
		if (float_data[i * 4 + 1] < minY) minY = float_data[i * 4 + 1];
		else if (float_data[i * 4 + 1] > maxY) maxY = float_data[i * 4 + 1];
		if (float_data[i * 4 + 2] < minZ) minZ = float_data[i * 4 + 2];
		else if (float_data[i * 4 + 2] > maxZ) maxZ = float_data[i * 4 + 2];
	}
	Cube res = {minX, minY, minZ, maxX, maxY, maxZ};
	return res;
}

void PCD::euclideanCluster(std::vector<std::vector<int>> *indices,float distance,size_t minSize,size_t maxSize,size_t maxClusters) {
	bool* visited = new bool[numPoints]();
	for (int i=0;i<numPoints;i++) {
		if (visited[i]) continue;
		printf("%d %lu\n",i,indices->size());
		std::vector<int> P,Q,neighbors;
		Q.push_back(i);
		visited[i] = true;
		while (Q.size() > 0) {
			int p = Q[Q.size()-1];
			Q.pop_back();
			P.push_back(p);
			searchKdTree(&neighbors,float_data[p*4],float_data[p*4+1],float_data[p*4+2],distance);
			for (size_t j=0;j<neighbors.size();j++) {
				if (!visited[neighbors[j]]) {
					Q.push_back(neighbors[j]);
					visited[neighbors[j]] = true;
				}
			}
			neighbors.clear();
		}
		if (P.size() >= minSize && P.size() <= maxSize)
			indices->push_back(P);
		if (indices->size() >= maxClusters)
			break;
	}
	delete[] visited;
}
