#include "normal.h"

Normal::Normal(PCD* pointcloud, float radius) {
	if (!pointcloud || pointcloud->numPoints == 0)
		return;

	//set up data structures
	cloud = pointcloud;
	norm = new float[pointcloud->numPoints * 3];
	curvature = new float[pointcloud->numPoints]();
	if (!pointcloud->kdtree)
		pointcloud->kdtree = new KdTree(pointcloud);
	
	for (int i=0;i<pointcloud->numPoints;i++) {
		float P[3] = {
			pointcloud->float_data[i * 4],
			pointcloud->float_data[i * 4 + 1],
			pointcloud->float_data[i * 4 + 2]
		};
		//form covariance matrix
		std::vector<int> neighbors;
		pointcloud->kdtree->search(&neighbors,P[0],P[1],P[2],radius);
		PCD::Quaternion center = pointcloud->getCentroid(&neighbors);
		double cov[9] = {}; //column major
		for (size_t j=0;j<neighbors.size();j++) {
			float deltaP[3] = {
				pointcloud->float_data[neighbors[j] * 4] - center.i,
				pointcloud->float_data[neighbors[j] * 4 + 1] - center.j,
				pointcloud->float_data[neighbors[j] * 4 + 2] - center.k,
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
		double lambda_real[3], lambda_imag[3], v[9];
		eigenvalue(3,cov,lambda_real,lambda_imag,v);
		//get normal and curvature
		double minLambda;
		int minIndex;
		for (int j=0;j<3;j++) {
			curvature[i] += lambda_real[j];
			if (j==0 || lambda_real[j] < minLambda) {
				minIndex = j;
				minLambda = lambda_real[j];
			}
		}
		curvature[i] = minLambda / curvature[i];
		double N[3] = {
			v[minIndex * 3],
			v[minIndex * 3 + 1],
			v[minIndex * 3 + 2]
		};
		//normal points to origin?
		if (P[0] * N[0] + P[1] * N[1] + P[2] * N[2] < 0) {
			norm[i * 3] = N[0];
			norm[i * 3 + 1] = N[1];
			norm[i * 3 + 2] = N[2];
		} else {
			norm[i * 3] = -N[0];
			norm[i * 3 + 1] = -N[1];
			norm[i * 3 + 2] = -N[2];
		}
//		printf("%4f,%4f,%4f: %lu neighbors %4f,%4f,%4f %4f\n",
//			P[0],P[1],P[2],neighbors.size(),N[0],N[1],N[2],curvature[i]);
	}
}

Normal::~Normal() {
	if (norm) delete[] norm;
	if (curvature) delete[] curvature;
}

//matrices are column major
void Normal::eigenvalue(int N,double* A,double* lambda_real,double* lambda_imag,double* v) {

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

//matrices are column major
void Normal::svd(int M,int N,double* A,double *U, double* S, double* VT) {

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

