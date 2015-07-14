#include "pcd.h"
#include "descriptor.h"
#include "normal.h"
#include <string>

int numCommandLineArgs;
char** commandLineArgs;

void testEigenvalue() {
	int n=4;
	double A[4*4] = {0.35,0.09,-0.44,0.25,
					 0.45,0.07,-0.33,-0.32,
					 -0.14,-0.54,-0.03,-0.13,
					 -0.17,0.35,0.17,0.11};
	double lambda_real[4] = {};
	double lambda_imag[4] = {};
	double v[4*4] = {};
	Normal::eigenvalue(n,A,lambda_real,lambda_imag,v);
	printf("lambda: ");
	for (int i=0;i<n;i++)
		printf("%f+%fi ",lambda_real[i],lambda_imag[i]);
	printf("\n");
	for (int i=0;i<n;i++) {
		printf("v%d: ",i);
		for (int j=0;j<n;j++) {
			if (lambda_imag[i]==0)
				printf("%f ",v[i*n+j]);
			else if (lambda_imag[i]>0)
				printf("%f+%fi ",v[i*n+j],v[(i+1)*n+j]);
			else
				printf("%f+%fi ",v[(i-1)*n+j],-v[i*n+j]);

		}
		printf("\n");
	}
}

void testSVD() {
	int m=6,n=4,l=(m>n?n:m);
	double A[6*4] = {2.27,0.28,-0.48,1.07,-2.35,0.62,
					-1.54,-1.67,-3.09,1.22,2.93,-7.39,
					1.15,0.94,0.99,0.79,-1.45,1.03,
					-1.94,-0.78,-0.21,0.63,2.30,-2.57};
	double U[6*6] = {};
	double S[4] = {};
	double VT[4*4] = {};
	Normal::svd(m,n,A,U,S,VT);
	printf("S: ");
	for (int i=0;i<l;i++)
		printf("%f ",S[i]);
	printf("\n");
	for (int i=0;i<m;i++) {
		printf("u%d: ",i);
		for (int j=0;j<m;j++) {
			printf("%f ",U[i*m+j]);
		}
		printf("\n");
	}
	for (int i=0;i<n;i++) {
		printf("v%d: ",i);
		for (int j=0;j<n;j++) {
			printf("%f ",VT[i*n+j]);
		}
		printf("\n");
	}
}

void fileConversion(char* in, char* out) {
	PCD* p,*q;
	int l = strlen(in);
	if (strncmp(in + l - 4,".pcd",4) == 0) {
		p = new PCD(in);
	} else if (strncmp(in + l - 4,".ply",4) == 0) {
		p = PCD::LoadFromPLY(in);
	} else if (strncmp(in + l - 4,".bin",4) == 0) {
		p = PCD::LoadFromKITTI(in,NULL);
		printf("Loaded %s (%s, %d points)\n",in,
			p->data_storage==PCD::ASCII ? "ascii" : "binary",
			p->numPoints);

		float segThreshold = 0.1, filterThreshold = 0.5, clusterThreshold = 0.4;
		
		PCD::Plane  coefficients = p->segmentPlane(10000,segThreshold,0.5);
		printf("Plane Coefficients: %f %f %f %f\n",coefficients.a,coefficients.b,coefficients.c,coefficients.d);
		std::vector<int> ind;
		p->filterPlane(&ind,coefficients,filterThreshold,false);
		q = p->extractIndices(&ind);
		char originalPath[256];
		snprintf(originalPath,255,"%s/original.pcd",out);
		q->writeToPCD(originalPath);

//		q = p;
		q->kdtree = new KdTree(q);
		printf("kdtree depth: %d\n",q->kdtree->kdtreeDepth);
		std::vector<std::vector<int>> clusters;
		q->euclideanClusteringKDTree(&clusters,clusterThreshold,100,100000,200);
		q->writeClustersToPCD(&clusters,out);

		delete p;
		delete q;
		return;
	} else if (strncmp(in + l - 4,".pts",4) == 0) {
		p = PCD::LoadFromPTS(in);
	} else if (DIR* d = opendir(in)){
		closedir(d);
		l = strlen(out);
		if (strncmp(out + l - 14,"descriptor.pcd",14) == 0) {
			Descriptor* d = Descriptor::LoadFromDir(in);
			d->writeToPCD(out);
			delete d;
			return;
		} else {
			p = PCD::LoadFromCluster(in);
		}
	}


	printf("Loaded %s (%s, %d points)\n",in,
		p->data_storage==PCD::ASCII ? "ascii" : "binary",
		p->numPoints);

	l = strlen(out);
	if (strncmp(out + l - 4,".pcd",4) == 0)
		p->writeToPCD(out);
	else if (strncmp(out + l - 4,".ply",4) == 0)
		p->writeToPLY(out);
	else if (strncmp(out + l - 4,".off",4) == 0)
		p->writeToOFF(out);
	else if (strncmp(out + l - 4,".obj",4) == 0)
		p->writeToOBJ(out);
	delete p;
}

int main(int argc,char* argv[]) {
	srand(time(NULL));
//	srand(0);
	numCommandLineArgs = argc;
	commandLineArgs = argv;

	if (argc >= 3)
		fileConversion(argv[1], argv[2]);

//	printf("Loaded %s (%d points)\n",f, d->numPoints);
//	std::vector<std::vector<int>> clusters;
//	d->kMeansClustering(&clusters);
//	for (size_t i=0;i<clusters.size();i++) {
//		for (size_t j=0;j<clusters[i].size();j++) {
//			printf("%d.pcd ",clusters[i][j]);
//		}
//		printf("\n");
//	}

//	p->loadDescriptor(argv[3]);

//	delete d;
}
