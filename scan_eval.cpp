#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>

struct Point {
	double x,y,z;
};

double getCompleteness(Point p,std::vector<Point> *scans) {
	double d = 0;
	for (size_t j=0;j<scans->size();j++) {
		Point q = (*scans)[j];
		double e = (p.x-q.x)*(p.x-q.x) + (p.y-q.y)*(p.y-q.y);
		d += 1 / sqrt(e);
	}
	return d;
}

void histogram(std::vector<double> x,std::vector<double> y,double* xh,double* yh,int numBins) {
	double xmin=x[0],xmax=x[0];
	for (size_t i=1;i<x.size();i++) {
		if (x[i] < xmin)
			xmin = x[i];
		if (x[i] > xmax)
			xmax = x[i];
	}
	double xstep = (xmax-xmin)/numBins;
	for (int i=0;i<numBins;i++)
		xh[i] = xmin + i * xstep;
	int* tp = new int[numBins]();
	int* fp = new int[numBins]();
	for (size_t i=0;i<x.size();i++) {
		int idx = (int)((x[i]-xmin)/xstep);
		if (y[i] > 0)
			tp[idx]++;
		else
			fp[idx]++;
	}
	for (int i=0;i<numBins;i++) {
//		int total = tp[i]+fp[i];
//		if (total > 0)
//			yh[i] = 100.0 * tp[i] / total;
		yh[i] = tp[i];
	}
	delete[] tp;
	delete[] fp;	
}

int main(int argc,char* argv[]) {
	if (argc < 2) {
		printf("./scan_eval clusterFolder\n");
		return 1;
	}

	char buffer[128];
	sprintf(buffer,"%s/scan_origin.txt",argv[1]);
	FILE* scan_origin = fopen(buffer,"r");
	if (!scan_origin) {
		printf("Cannot open %s\n",buffer);
		return 1;
	}
	sprintf(buffer,"%s/completeness.txt",argv[1]);
	FILE* completeness = fopen(buffer,"w");
	if (!completeness) {
		printf("Cannot open %s\n",buffer);
		return 1;
	}
	sprintf(buffer,"%s/labels.txt",argv[1]);
	FILE* labelFile = fopen(buffer,"r");
	if (!labelFile)
		return 1;
	sprintf(buffer,"%s/prediction.txt",argv[1]);
	FILE* predictionFile = fopen(buffer,"r");
	if (!predictionFile)
		return 1;

	std::vector<Point> scans;
	std::vector<Point> clusters;
	std::vector<int> labels;
	std::vector<int> predictions;
	std::vector<double> score;

	while (fgets(buffer,128,labelFile)) {
		int l = strtol(buffer,NULL,10);
		labels.push_back(l);
	}
	while (fgets(buffer,128,predictionFile)) {
		char* c;
		int l = strtol(buffer,&c,10);
		predictions.push_back(l);
		double max_score = 0;
		while (*c != '\n') {
			double score = strtod(c,&c);
			max_score = score > max_score ? score : max_score;
		}
		score.push_back(max_score);
	}
	while (fgets(buffer,128,scan_origin)) {
		Point p;
		char* c = buffer;
		p.x = strtod(c,&c);
		p.y = strtod(c,&c);
		p.z = strtod(c,&c);
		scans.push_back(p);
	}

	for (size_t i=0;i<labels.size();i++) {
		sprintf(buffer,"%s/%lu-cloud.pcd",argv[1],i);
		FILE* f = fopen(buffer,"r");
		if (!f) {
			printf("Cannot open %s\n",buffer);
			return 1;
		}
		while (fgets(buffer,128,f)) {
			if (strncmp(buffer,"DATA",4)==0)
				break;
		}
		double sumX=0,sumY=0,sumZ=0;
		double x,y,z;
		int count=0;
		while (fscanf(f,"%lf %lf %lf",&x,&y,&z)==3) {
			sumX += x;
			sumY += y;
			sumZ += z;
			count++;
		}
		Point p = {sumX/count, sumY/count, sumZ/count};
		clusters.push_back(p);
		fclose(f);
	}

	std::vector<double> x,y;
	for (size_t i=0;i<labels.size();i++) {
		x.push_back(getCompleteness(clusters[i],&scans));
		y.push_back(labels[i]==predictions[i]);
	}
	int numBins = 8;
	double* xh = new double[numBins]();
	double* yh = new double[numBins]();
	histogram(x,y,xh,yh,numBins);
	for (int i=0;i<numBins;i++) {
		fprintf(completeness,"%f %f\n",xh[i],yh[i]);
	}
	delete[] xh;
	delete[] yh;

	fclose(scan_origin);
	fclose(completeness);
	fclose(labelFile);
	fclose(predictionFile);
}
