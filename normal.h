#pragma once
#include <math.h>
#include <vector>
#include "pcd.h"
#include "lapack.h"

class PCD;

class Normal {
	public:
		PCD* cloud;
		float* norm;
		float* curvature;

		Normal(PCD* pointcloud, float radius);
		~Normal();
		
		static void eigenvalue(int N,double* A,double* lambda_real,double* lambda_imag,double* v);
		static void svd(int M,int N,double* A,double *U, double* S, double* VT);
};
