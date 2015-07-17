#include "pcd.h"

//http://www.rcn.montana.edu/Resources/Converter.aspx
void wgs2utm(double lat, double lon, double *easting, double *northing) {
	double a = 6378137;
	double f = 1 / 298.257223563;
	double k0 = 0.9996;
	double b = a * (1 - f);
	double e = sqrt(1 - b * b / a / a);
	double phi = lat / 180 * M_PI;

	double utmz = 1 + floor((lon + 180) / 6);            // longitude to utm zone
	double zcm = 3 + 6 * (utmz - 1) - 180;                     // central meridian of a zone
	double esq = (1 - (b / a) * (b / a));
	double e0sq = e * e / (1 - pow(e, 2));
	double M = 0;

	double N = a / sqrt(1 - pow(e * sin(phi), 2));
	double T = pow(tan(phi), 2);
	double C = e0sq * pow(cos(phi), 2);
	double A = (lon - zcm) * M_PI / 180 * cos(phi);

	// calculate M (USGS style)
	M = phi * (1 - esq * (1.0 / 4 + esq * (3.0 / 64 + 5 * esq / 256)));
	M = M - sin(2 * phi) * (esq * (3.0 / 8 + esq * (3.0 / 32 + 45 * esq / 1024)));
	M = M + sin(4 * phi) * (esq * esq * (15.0 / 256 + esq * 45.0 / 1024));
	M = M - sin(6 * phi) * (esq * esq * esq * (35.0 / 3072));
	M = M * a;                                      //Arc length along standard meridian
	double M0 = 0;                                         // if another point of origin is used than the equator

	// now we are ready to calculate the UTM values...
	// first the easting
	double x = k0 * N * A * (1 + A * A * ((1 - T + C) / 6 + A * A * (5 - 18 * T + T * T + 72 * C - 58 * e0sq) / 120)); //Easting relative to CM
	x = x + 500000; // standard easting
	// now the northing
	double y = k0 * (M - M0 + N * tan(phi) * (A * A * (1.0 / 2 + A * A * ((5 - T + 9 * C + 4 * C * C) / 24 + A * A * (61 - 58 * T + T * T + 600 * C - 330 * e0sq) / 720))));    // first from the equator
	if (y < 0) {
		y = 10000000 + y;   // add in false northing if south of the equator
	}

	*easting = x;
	*northing = y;
}

PCD::PCD(const char* fileName) {
	numFields = 4;
	data = NULL;
	float_data = NULL;
	fields = NULL;
	numPoints = 0;
	capacity = 0;
	data_storage = NONE;
	kdtree = NULL;
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
	kdtree = NULL;
	camera = { {0,0,0,0}, {1,0,0,0}, 1};
}

PCD::~PCD() {
	if (fields) free(fields);
	if (data) free(data);
	if (float_data) free(float_data);
	if (kdtree) delete kdtree;
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
		snprintf(buffer,256,"%s/%lu-cloud.pcd",dir,i);
		std::vector<int> currentIndices = (*indices)[i];
		PCD* subcloud = extractIndices(&currentIndices);
		subcloud->writeToPCD(buffer);
		delete subcloud;
	}
}

void PCD::writeToOFF(const char* filename) {
	if (!filename)
		return;
	FILE* f = fopen(filename, "w");
	if (!f) {
		printf("Cannot write to file: %s\n", filename);
		return;
	}

	fprintf(f,"OFF\n8 6 12\n");
#define USE_OBB 1
#if USE_OBB
	Descriptor d(this);
	double lambda[3], v[9];
	d.getPCA_XY(lambda,v);
	d.setPCA(lambda,v);
	for (int i=0;i<8;i++) {
		float coords[3];
		for (int j=0;j<3;j++) {
			coords[j] = d.bbCenter[j];
			for (int axis=0;axis<3;axis++) {
				float sign = (i & 1<<axis) ? 1 : -1;
				coords[j] += sign * d.principalLengths[axis] / 2 * d.principalAxes[axis][j];
			}
			fprintf(f,"%f ",coords[j]);
		}
		fprintf(f,"\n");
	}
#else
	if (!kdtree) kdtree = new KdTree(this);
	KdTree::Cube bb = kdtree->getBoundingBox();
	fprintf(f,"%f %f %f\n",bb.x1,bb.y1,bb.z1);
	fprintf(f,"%f %f %f\n",bb.x2,bb.y1,bb.z1);
	fprintf(f,"%f %f %f\n",bb.x1,bb.y2,bb.z1);
	fprintf(f,"%f %f %f\n",bb.x2,bb.y2,bb.z1);
	fprintf(f,"%f %f %f\n",bb.x1,bb.y1,bb.z2);
	fprintf(f,"%f %f %f\n",bb.x2,bb.y1,bb.z2);
	fprintf(f,"%f %f %f\n",bb.x1,bb.y2,bb.z2);
	fprintf(f,"%f %f %f\n",bb.x2,bb.y2,bb.z2);
#endif
	fprintf(f,"4 0 1 3 2 1 0 0 0\n"
	"4 0 4 5 1 0 1 0 0\n"
	"4 0 2 6 4 1 0 0 0\n"
	"4 1 3 7 5 0 1 0 0\n"
	"4 2 6 7 3 1 0 0 0\n"
	"4 4 5 7 6 0 1 0 0\n");
	fclose(f);
}

void PCD::writeToOBJ(const char* filename) {
	if (!filename)
		return;
	FILE* f = fopen(filename, "w");
	if (!f) {
		printf("Cannot write to file: %s\n", filename);
		return;
	}

#if USE_OBB
	Descriptor d(this);
	double lambda[3], v[9];
	d.getPCA_XY(lambda,v);
	d.setPCA(lambda,v);
	for (int i=0;i<8;i++) {
		float coords[3];
		fprintf(f,"v ");
		for (int j=0;j<3;j++) {
			coords[j] = d.bbCenter[j];
			for (int axis=0;axis<3;axis++) {
				float sign = (i & 1<<axis) ? 1 : -1;
				coords[j] += sign * d.principalLengths[axis] / 2 * d.principalAxes[axis][j];
			}
			fprintf(f,"%f ",coords[j]);
		}
		fprintf(f,"\n");
	}
#else
	if (!kdtree) kdtree = new KdTree(this);
	KdTree::Cube bb = kdtree->getBoundingBox();
	fprintf(f,"%f %f %f\n",bb.x1,bb.y1,bb.z1);
	fprintf(f,"%f %f %f\n",bb.x2,bb.y1,bb.z1);
	fprintf(f,"%f %f %f\n",bb.x1,bb.y2,bb.z1);
	fprintf(f,"%f %f %f\n",bb.x2,bb.y2,bb.z1);
	fprintf(f,"%f %f %f\n",bb.x1,bb.y1,bb.z2);
	fprintf(f,"%f %f %f\n",bb.x2,bb.y1,bb.z2);
	fprintf(f,"%f %f %f\n",bb.x1,bb.y2,bb.z2);
	fprintf(f,"%f %f %f\n",bb.x2,bb.y2,bb.z2);
#endif
	fprintf(f,
	"l 1 2\n"
	"l 1 3\n"
	"l 2 4\n"
	"l 3 4\n"
	"l 1 5\n"
	"l 2 6\n"
	"l 3 7\n"
	"l 4 8\n"
	"l 5 6\n"
	"l 5 7\n"
	"l 6 8\n"
	"l 7 8\n");
	fclose(f);
}

PCD* PCD::LoadFromPLY(const char* fileName) {
	FILE* f = fopen(fileName, "r");
	if (!f) {
		printf("File not found: %s\n", fileName);
		return NULL;
	}
	char buf[256];
	int n;
	PCD* pcd;
	int r,g,b,rgb;
	while (fgets(buf, 256, f)) {
		if (sscanf(buf,"element vertex %d",&n)==1) {
			pcd = new PCD(n);
		} else if (strncmp(buf,"end_header",10)==0) {
			for (int i=0;i<n;i++) {
				fgets(buf,256,f);
				if (sscanf(buf, "%f %f %f %d %d %d", pcd->float_data + i * 4, pcd->float_data + i * 4 + 1,
					pcd->float_data + i * 4 + 2, &r, &g, &b) == 6) {
					rgb = (r<<16) | (g<<8) | b;
					pcd->float_data[i * 4 + 3] = (float) rgb;
				} else if (sscanf(buf, "%f %f %f", pcd->float_data + i * 4, pcd->float_data + i * 4 + 1,
					pcd->float_data + i * 4 + 2) == 3) {
					rgb = 255;
					pcd->float_data[i * 4 + 3] = (float) rgb;
				} else {
					printf("Error parsing %s\n",fileName);
					printf("Line %d: %s\n",i,buf);
					break;
				}
			}
			break;
		}
	}
	fclose(f);
	return pcd;

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
	PCD lines(0);
	const float interval = 2;
	const int maxLines = 10;
	for (int j=0;j<num_cameras;j++) {
		int numLines = lineInfo[j].size() / 7;
		for (int i=0;i<numLines&&i<maxLines;i++) {
			lines.drawLine(
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
	lines.drawAxis(0,0,0,NULL,10,1);
	printf("Wrote %d points\n",lines.numPoints);
	lines.writeToPCD("p2.pcd");

	fclose(f);
	delete[] camera_positions;
	delete[] camera_rotations;
	delete[] lineInfo;
	return pointcloud;
}

PCD* PCD::LoadFromKITTI(const char* fileName, const char* oxts) {
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
//		p->float_data[i * 4 + 3] = colormap(p->float_data[i * 4 + 3]);
		p->float_data[i * 4 + 3] = 255 << 16;
//		p->float_data[i * 4 + 3] = (int) (255 * p->float_data[i * 4 + 3]) << 16;
	}
	fclose (f);
	if (!oxts) return p;

	f = fopen(oxts,"r");
	if (!f) {
		printf("File not found: %s\n", fileName);
		return p;
	}
	char buf[512];
	fgets(buf,512,f);
	char* c = buf;

	double latitude = strtod(c,&c);
	double longitude = strtod(c,&c);
	double altitude = strtod(c,&c);
	double roll = strtod(c,&c);
	double pitch = strtod(c,&c);
	double yaw = strtod(c,&c);
	double E=0,N=0;
	wgs2utm(latitude,longitude,&E,&N);
	p->translate(E - 457840.938302, N - 5428721.722951, altitude - 113.966202);

	fclose(f);
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
//		printf("n = %d\n",n);
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

PCD* PCD::LoadFromCluster(const char* dir) {
	char buffer[1024];
	int ndir = strlen(dir);
	strncpy(buffer,dir,ndir);
	char* buffer_c = buffer + ndir;
	*buffer_c++ = '/';

	strcpy(buffer_c,"prediction.txt");
	FILE* labelFile = fopen(buffer,"r");
	std::vector<int> labels;
	char line[64];
	while (fgets(line,64,labelFile)) {
		labels.push_back(atoi(line));
	}
	fclose(labelFile);

	PCD* combined = new PCD(0);
	float colorChoice[] = {16777215, 255<<16, 255<<8, 255, (255<<16)|(255<<8), (255<<16)|255, (255<<8)|255};
	for (size_t i=0;i<labels.size();i++) {
		snprintf(buffer_c,64,"%lu-cloud.pcd",i);
		PCD cloud(buffer);
		int offset = combined->numPoints;
		combined->expand(offset + cloud.numPoints);
		for (int j = 0;j<cloud.numPoints;j++) {
			combined->float_data[(offset+j) * 4] = cloud.float_data[j * 4];
			combined->float_data[(offset+j) * 4 + 1] = cloud.float_data[j * 4 + 1];
			combined->float_data[(offset+j) * 4 + 2] = cloud.float_data[j * 4 + 2];
			combined->float_data[(offset+j) * 4 + 3] = colorChoice[labels[i]];
		}
	}

	return combined;
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

PCD::Quaternion PCD::quaternionFromAngle(float rx,float ry, float rz) {
	float r,i,j,k;
	float r1 = rz * M_PI / 180;
	float r2 = ry * M_PI / 180;
	float r3 = rx * M_PI / 180;
	float M[9];
	M[0] = cos(r2)*cos(r3);
	M[1] = -cos(r2) * sin(r3);
	M[2] = sin(r2);
	M[3] = cos(r1) * sin(r3) + cos(r3) * sin(r1) * sin(r2);
	M[4] = cos(r1) * cos(r3) - sin(r1) * sin(r2) * sin(r3);
	M[5] = -cos(r2) * sin(r1);
	M[6] = sin(r1) * sin(r3) - cos(r1) * cos(r3) * sin(r2);
	M[7] = cos(r3) * sin(r1) + cos(r1) * sin(r2) * sin(r3);
	M[8] = cos(r1) * cos(r2);

	double tr = M[0] + M[4] + M[8];

	if (tr > 0) { 
		float S = sqrt(tr+1.0) * 2;
		r = (float)(0.25 * S);
		i = (float)((M[7] - M[5]) / S);
		j = (float)((M[2] - M[6]) / S);
		k = (float)((M[3] - M[1]) / S);
	} else if ((M[0] > M[4])&&(M[0]> M[8])) {
		float S = sqrt(1.0 + M[0]- M[4] - M[8]) * 2;
		i = (float)(0.25 * S);
		r = (float)((M[7] - M[5]) / S);
		k = (float)((M[2] + M[6]) / S);
		j = (float)((M[3] + M[1]) / S);
	} else if (M[4]>M[8]) {
		float S = sqrt(1.0 - M[0] + M[4] - M[8]) * 2; 
		j = (float)(0.25 * S);
		k = (float)((M[7] + M[5]) / S);
		r = (float)((M[2] - M[6]) / S);
		i = (float)((M[3] + M[1]) / S); 
	} else {
		float S = sqrt(1.0 - M[0] - M[4] + M[8]) * 2; 
		k = (float)(0.25 * S);
		j = (float)((M[7] + M[5]) / S);
		i = (float)((M[2] + M[6]) / S);
		r = (float)((M[3] - M[1]) / S); 
	}
	Quaternion res = {r,i,j,k};
	return res;
}

PCD::Quaternion PCD::quaternionMult(Quaternion qa, Quaternion qb) {
	Quaternion qc;
	qc.r = qa.r*qb.r - qa.i*qb.i - qa.j*qb.j - qa.k*qb.k;
	qc.i = qa.r*qb.i + qa.i*qb.r + qa.j*qb.k - qa.k*qb.j;
	qc.j = qa.r*qb.j + qa.j*qb.r + qa.k*qb.i - qa.i*qb.k;
	qc.k = qa.r*qb.k + qa.k*qb.r + qa.i*qb.j - qa.j*qb.i;
	return qc;
}

PCD::Quaternion PCD::quaternionInv(Quaternion q) {
	Quaternion qc;
	qc.r = q.r;
	qc.i = -q.i;
	qc.j = -q.j;
	qc.k = -q.k;
	return qc;
}

void PCD::rotate(Quaternion q) {
	for (int i=0;i<numPoints;i++) {
		Quaternion vec = {0, float_data[i*4], float_data[i*4+1], float_data[i*4+2]};
		vec = quaternionMult(quaternionMult(q, vec), quaternionInv(q));
		float_data[i * 4] = vec.i;
		float_data[i * 4 + 1] = vec.j;
		float_data[i * 4 + 2] = vec.k;
	}
}

void PCD::translate(float x, float y, float z) {
	for (int i=0;i<numPoints;i++) {
		float_data[i * 4] += x;
		float_data[i * 4 + 1] += y;
		float_data[i * 4 + 2] += z;
	}
}

//indices: subcloud under consideration, complete cloud if null
PCD::Quaternion PCD::getCentroid(std::vector<int> *indices) {
	float sumX=0, sumY=0, sumZ=0;
	if (indices) {
		for (size_t i=0;i<indices->size();i++) {
			sumX += float_data[(*indices)[i] * 4];
			sumY += float_data[(*indices)[i] * 4 + 1];
			sumZ += float_data[(*indices)[i] * 4 + 2];
		}
	} else {
		for (int i=0;i<numPoints;i++) {
			sumX += float_data[i * 4];
			sumY += float_data[i * 4 + 1];
			sumZ += float_data[i * 4 + 2];
		}
	}
	Quaternion q = {0, sumX / numPoints, sumY / numPoints, sumZ / numPoints};
	return q;
}

int PCD::get3DProjection(unsigned int* pixels, int screenWidth, int screenHeight, int pointSize, bool zoomToFit) {
	int i, j, k;
	Quaternion q, newPoint;
	std::vector<float> x_data;
	std::vector<float> y_data;
	std::vector<unsigned int> rgb_data;
	q.r = 0;
	printf("Using camera: (%f %f %f %f %f %f %f)\n", camera.position.i, camera.position.j, camera.position.k,
		camera.rotation.r, camera.rotation.i, camera.rotation.j, camera.rotation.k);
	for (i = 0; i < numPoints; i++) {
		//first transform point cloud to camera coordinates
		q.i = float_data[i * 4] - camera.position.i;
		q.j = float_data[i * 4 + 1] - camera.position.j;
		q.k = float_data[i * 4 + 2] - camera.position.k;
		newPoint = quaternionMult(quaternionMult(camera.rotation, q), quaternionInv(camera.rotation));
		//next map points to image plane based on camera focal length
		if ( - newPoint.k > camera.focal_length) { //if not behind camera
			x_data.push_back(newPoint.i * camera.focal_length / - newPoint.k);
			y_data.push_back(newPoint.j * camera.focal_length / - newPoint.k);
			rgb_data.push_back((unsigned int) float_data[i * 4 + 3]);
		}
	}
	//lastly plot these points on xy plane of bitmap
	int pixelsDrawn = 0;
	float maxX = -1;
	float maxY = -1;
	float x, y, scaleX, scaleY;
	float margin = 0.05f;
	if (zoomToFit){
		for (size_t n = 0; n < x_data.size(); n++) {
			if (fabs(x_data[n]) > maxX) maxX = fabs(x_data[n]);
			if (fabs(y_data[n]) > maxY) maxY = fabs(y_data[n]);
		}
	} else {
		maxX = 1;
		maxY = 1;
	}
//	printf("display rectangle: %f x %f\n", maxX, maxY);
	if (maxX > maxY) {
		scaleX = (screenWidth - 1) * (1 - margin * 2) / (maxX * 2);
		scaleY = (screenHeight - 1) * (1 - margin * 2) / (maxX * 2);
	}
	else {
		scaleX = (screenWidth - 1) * (1 - margin * 2) / (maxY * 2);
		scaleY = (screenHeight - 1) * (1 - margin * 2) / (maxY * 2);
	}
	for (size_t n = 0; n < x_data.size(); n++) {
		x = x_data[n] * scaleX + screenWidth / 2;
		y = y_data[n] * scaleY + screenHeight / 2;
		x -= (pointSize - 1) / 2;
		y -= (pointSize - 1) / 2;
		for (j = 0; j < pointSize; j++) {
			for (k = 0; k < pointSize; k++) {
				if (x+j>=0 && x+j<screenWidth && y+k>=0 && y+k<screenHeight) {
					pixelsDrawn++;
					pixels[((int)y+k)*screenWidth + ((int)x+j)] = rgb_data[n];
				}
			}
		}
	}
	return pixelsDrawn;
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

void PCD::euclideanClusteringKDTree(std::vector<std::vector<int>> *indices,float distance,size_t minSize,size_t maxSize,size_t maxClusters) {
	if (!kdtree) kdtree = new KdTree(this);
	bool* visited = new bool[numPoints]();
	for (int i=0;i<numPoints;i++) {
		if (visited[i]) continue;
//		printf("%d %lu\n",i,indices->size());
		std::vector<int> P,Q,neighbors;
		Q.push_back(i);
		visited[i] = true;
		while (Q.size() > 0) {
			int p = Q[Q.size()-1];
			Q.pop_back();
			P.push_back(p);
			kdtree->search(&neighbors,float_data[p*4],float_data[p*4+1],float_data[p*4+2],distance);
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

struct Ind {
	float val;
	int index;
};

int compare(const void* v1,const void* v2) {
	Ind* i1 = (Ind*) v1;
	Ind* i2 = (Ind*) v2;
	if (i1->val < i2->val)
		return -1;
	else if (i1->val > i2->val)
		return 1;
	else
		return 0;
}

void* closestBSearch(void* key,void* base,int num,int size,int (*compar)(const void*,const void*),bool lowerOrUpper) {
	char* A = (char*) base;
	int l = 0, h = num-1, mid=(l+h)/2, diff;
	if (lowerOrUpper) { //lower
		while (l < h) {
			diff = compar(key,A + mid * size);
			if (diff > 0)
				l = mid + 1;
			else if (diff <= 0)
				h = mid;
			mid = (l + h) / 2;
		}
	} else { //upper
		while (l < h - 1) {
			diff = compar(key,A + mid * size);
			if (diff >= 0)
				l = mid;
			else if (diff < 0)
				h = mid - 1;
			mid = (l + h) / 2;
		}
		if (l == h - 1) //edge case
			if (compar(key, A + h * size) >= 0)
				mid++;
	}
	return A + mid * size;
}

void PCD::euclideanClustering(std::vector<std::vector<int>> *indices,float distance,size_t minSize,size_t maxSize,size_t maxClusters) {
	Ind* byX = new Ind[numPoints];
	Ind* byY = new Ind[numPoints];
	Ind* byZ = new Ind[numPoints];
	for (int i=0;i<numPoints;i++) {
		byX[i] = {float_data[i*4], i};
		byY[i] = {float_data[i*4 + 1], i};
		byZ[i] = {float_data[i*4 + 2], i};
	}
	qsort(byX,numPoints,sizeof(Ind),compare);
	qsort(byY,numPoints,sizeof(Ind),compare);
	qsort(byZ,numPoints,sizeof(Ind),compare);
	bool* visited = new bool[numPoints]();
	for (int i=0;i<numPoints;i++) {
		if (visited[i]) continue;
		std::vector<int> P,Q,neighbors;
		Q.push_back(i);
		visited[i] = true;
		while (Q.size() > 0) {
			int p = Q[Q.size()-1];
			Q.pop_back();
			P.push_back(p);
			Ind key,*xl,*xh,*yl,*yh,*zl,*zh;
			key.val = float_data[p*4] - distance;
			xl = (Ind*) closestBSearch(&key,byX,numPoints,sizeof(Ind),compare,true);
			key.val = float_data[p*4] + distance;
			xh = (Ind*) closestBSearch(&key,byX,numPoints,sizeof(Ind),compare,false);
			key.val = float_data[p*4+1] - distance;
			yl = (Ind*) closestBSearch(&key,byY,numPoints,sizeof(Ind),compare,true);
			key.val = float_data[p*4+1] + distance;
			yh = (Ind*) closestBSearch(&key,byY,numPoints,sizeof(Ind),compare,false);
			key.val = float_data[p*4+2] - distance;
			zl = (Ind*) closestBSearch(&key,byZ,numPoints,sizeof(Ind),compare,true);
			key.val = float_data[p*4+2] + distance;
			zh = (Ind*) closestBSearch(&key,byZ,numPoints,sizeof(Ind),compare,false);
			HashTable ht;
			int closeX,closeY,closeZ;
			for (Ind* k = xl; k <= xh; k++) ht.insert(k->index,&closeX);
			for (Ind* k = yl; k <= yh; k++) if (ht.remove(k->index)) ht.insert(k->index,&closeY);
			for (Ind* k = zl; k <= zh; k++) if (ht.remove(k->index) == &closeY) ht.insert(k->index,&closeZ);
			for (int k=0; k<ht.size; k++) {
				if (ht.entries[k] == &closeZ) {
					neighbors.push_back(ht.keys[k]);
				}
			}
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
	delete[] byX;
	delete[] byY;
	delete[] byZ;
	delete[] visited;
}
