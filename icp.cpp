#include "pcd.h"
#include "lapack.h"
#include <string>

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

void svd(int M,int N,double* A,double *U, double* S, double* VT) {

	int info, lwork=5*(M>N?N:M);
	double* work = new double[lwork];
	char jobu = 'A', jobvt = 'A';
	dgesvd_(&jobu, &jobvt, &M, &N, A, &M, 
	     S, U, &M, VT, &N, work, &lwork, &info);
//	printf("info: %d\n",info);
//	printf("optimal: %f\n",work[0]);
	if (info!=0) {
		printf("Error in subroutine dgesvd_ (info=%d)\n",info);
	}
	delete[] work;
}

void testEigenvalue() {
	int n=4;
	double A[n*n] = {0.35,0.09,-0.44,0.25,
					 0.45,0.07,-0.33,-0.32,
					 -0.14,-0.54,-0.03,-0.13,
					 -0.17,0.35,0.17,0.11};
	double lambda_real[n] = {};
	double lambda_imag[n] = {};
	double v[n*n] = {};
	eigenvalue(n,A,lambda_real,lambda_imag,v);
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
	double A[m*n] = {2.27,0.28,-0.48,1.07,-2.35,0.62,
					-1.54,-1.67,-3.09,1.22,2.93,-7.39,
					1.15,0.94,0.99,0.79,-1.45,1.03,
					-1.94,-0.78,-0.21,0.63,2.30,-2.57};
	double U[m*m] = {};
	double S[l] = {};
	double VT[n*n] = {};
	svd(m,n,A,U,S,VT);
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

int main(int argc,char* argv[]) {
	srand(time(NULL));
//	srand(0);

	std::string file("test_pcd.pcd");
	const char* f = file.c_str();
	if (argc>=2)
		f = argv[1];
	PCD* p = new PCD(f),*q;
	printf("Loaded %s (%s, %d points)\n",f,
		p->data_storage==PCD::ASCII ? "ascii" : "binary",
		p->numPoints);
	
//	PCD::Plane  coefficients = p->segmentPlane(10000,0.1,0.4);
//	printf("Plane Coefficients: %f %f %f %f\n",coefficients.a,coefficients.b,coefficients.c,coefficients.d);
//	std::vector<int> ind;
//	p->filterPlane(&ind,coefficients,0.1,false);
//	q = p->extractIndices(&ind);

	q = p;
	q->kdtree = new KdTree(q);
	printf("kdtree depth: %d\n",q->kdtree->kdtreeDepth);
//	std::vector<std::vector<int>> clusters;
//	q->euclideanCluster(&clusters,0.5,1000,50000,50);
//	q->writeClustersToPCD(&clusters,"clusters");

//	p->loadDescriptor(argv[3]);

	if (argc>=3) {
		int l = strlen(argv[2]);
		if (strncmp(argv[2] + l - 4,".pcd",4) == 0)
			p->writeToPCD(argv[2]);
		else if (strncmp(argv[2] + l - 4,".ply",4) == 0)
			p->writeToPLY(argv[2]);
	}

	delete p;
//	delete q;
}
