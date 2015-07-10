#include "descriptor.h"

Descriptor::Descriptor() {
	this->pointcloud = NULL;
	this->data = NULL;
}

Descriptor::Descriptor(const char* filename) {
	FILE* f = fopen(filename,"r");
	if (!f) {
		printf("File not found: %s\n", filename);
		return;
	}
	char buf[256];
	while (fgets(buf, 256, f)) {
		if (sscanf(buf, "POINTS %d", &numPoints) == 1) {
		} else if (sscanf(buf,"COUNT %d",&dimension)==1) {
		} else if (strncmp(buf,"DATA ascii",10)==0) 
			break;
	}
	float val;
	char* c;
	data = new float[numPoints * dimension];
	for (int i=0;i<numPoints;i++) {
		fgets(buf,256,f);
		c = buf;
		for (int j=0;j<dimension;j++) {
			val = strtod(c,&c);
			data[i * dimension + j] = val;
		}
	}
	fclose(f);
}

Descriptor::Descriptor(PCD* p) {
	numPoints = 1;
	this->pointcloud = p;
	this->data = NULL;
}

Descriptor::~Descriptor() {
	if (data) delete[] data;
}

Descriptor* Descriptor::LoadFromDir(const char* dir) {
	Descriptor* res = new Descriptor();

	char buffer[1024];
	int ndir = strlen(dir);
	strncpy(buffer,dir,ndir);
	char* buffer_c = buffer + ndir;
	*buffer_c++ = '/';

	std::vector<double> measurements;
	int j;
	for (j=0;;j++) {
		snprintf(buffer_c,64,"%d-cloud.pcd",j);
		PCD cloud(buffer);
		if (cloud.numPoints==0)
			break;
		Descriptor desc(&cloud);
		double lambda[3],v[9];
		desc.getPCA_XY(lambda,v);
		desc.setPCA(lambda,v);
		measurements.push_back(desc.principalLengths[0]); //length;
		measurements.push_back(desc.principalLengths[1]); //width;
		measurements.push_back(desc.principalLengths[2]); //height;
	}

	//set data
	res->numPoints = j;
	res->dimension = 3;
	res->data = new float[res->numPoints * res->dimension];
	for (size_t i=0;i<measurements.size();i++) {
		res->data[i] = measurements[i];
	}

	return res;
}

void Descriptor::writeToPCD(const char* filename) {
	if (!filename)
		return;
	FILE* f = fopen(filename, "w");
	if (!f) {
		printf("Cannot write to file: %s\n", filename);
		return;
	}
	fprintf(f,"# .PCD v0.7 - Point Cloud Data file format\n"
	"VERSION 0.7\n"
	"FIELDS descriptor\n"
	"SIZE 4\n"
	"TYPE F\n"
	"COUNT %d\n"
	"WIDTH %d\n"
	"HEIGHT 1\n"
	"POINTS %d\n"
	"DATA ascii\n",dimension,numPoints,numPoints);
	for (int i=0;i<numPoints;i++) {
		for (int j=0;j<dimension;j++)
			fprintf(f,"%f ",data[i*dimension+j]);
		fprintf(f,"\n");
	}
	fclose(f);
	printf("Wrote %d points to %s\n",numPoints,filename);
}

void Descriptor::getPCA(double* lambda_real,double* v) {
	double cov[9] = {}; //column major
	PCD::Quaternion center = pointcloud->getCentroid(NULL);
	for (int j=0;j<pointcloud->numPoints;j++) {
		float deltaP[3] = {
			pointcloud->float_data[j * 4] - center.i,
			pointcloud->float_data[j * 4 + 1] - center.j,
			pointcloud->float_data[j * 4 + 2] - center.k,
		};
		cov[0] += deltaP[0] * deltaP[0];
		cov[1] += deltaP[1] * deltaP[0];
		cov[2] += deltaP[2] * deltaP[0];
		cov[3] += deltaP[0] * deltaP[1];
		cov[4] += deltaP[1] * deltaP[1];
		cov[5] += deltaP[2] * deltaP[1];
		cov[6] += deltaP[0] * deltaP[2];
		cov[7] += deltaP[1] * deltaP[2];
		cov[8] += deltaP[2] * deltaP[2];
	}
	//compute PCA
	double lambda_imag[3];
	Normal::eigenvalue(3,cov,lambda_real,lambda_imag,v);
}

void Descriptor::setPCA(double* lambda_real, double* v) {
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
	this->bbCenter[0] = cx;
	this->bbCenter[1] = cy;
	this->bbCenter[2] = cz;
	for (int i=0;i<3;i++) {
		this->principalLengths[i] = maxScale[i] - minScale[i];
		for (int j=0;j<3;j++)
			this->principalAxes[i][j] = v[i * 3 + j];
	}
}

void Descriptor::getPCA_XY(double* lambda_real,double* v) {
	double cov[4] = {}; //column major
	PCD::Quaternion center = pointcloud->getCentroid(NULL);
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
	Normal::eigenvalue(2,cov,lambda_real,lambda_imag,v);
	//complete 3D structure
	lambda_real[2] = -1;
	v[4] = v[3];
	v[3] = v[2];
	v[2] = v[5] = v[6] = v[7] = 0;
	v[8] = 1;
}

void Descriptor::kMeansClustering(std::vector<std::vector<int>> *indices) {
	int k = 10;
//	int k = (int) sqrt(numPoints / 2);
	//initialize
	for (int i=0;i<k;i++) {
		std::vector<int> v;
		indices->push_back(v);
	}
	int* associations = new int[numPoints]();
	float* means = new float[k * dimension];
	int iterations = 0;

	//choose initial means
	for (int i=0;i<k;i++) {
		int index = rand() % numPoints;
		memcpy(means + i * dimension, data + index * dimension, dimension * sizeof(float));
	}

	bool update = true;
	while (update) {
		update = false;
		for (int i=0;i<k;i++)
			(*indices)[i].clear();

		//assignment step
		for (int i=0;i<numPoints;i++) {
			//find nearest mean
			int minIndex;
			float minDistance;
			for (int j=0;j<k;j++) {
				float distance = 0;
				for (int l=0;l<dimension;l++)
					distance += (data[i*dimension+l] - means[j*dimension+l]) * (data[i*dimension+l] - means[j*dimension+l]);
				if (j==0 || distance < minDistance) {
					minIndex = j;
					minDistance = distance;
				}
			}
			if (minIndex != associations[i])
				update = true;
			associations[i] = minIndex;
			(*indices)[minIndex].push_back(i);
		}

		//update step
		for (int i=0;i<k;i++) {
			for (int j=0;j<dimension;j++) {
				means[i * dimension + j] = 0;
				size_t clusterSize = (*indices)[i].size();
				for (size_t l=0;l<clusterSize;l++) {
					means[i * dimension + j] += data[(*indices)[i][l] * dimension + j];
				}
				means[i * dimension + j] /= clusterSize;
			}
		}
		iterations++;
	}

//	printf("K-Means: Found %lu clusters (%d iterations)\n",indices->size(),iterations);

	delete[] associations;
	delete[] means;
}
