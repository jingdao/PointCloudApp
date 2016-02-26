#include <pcl/io/pcd_io.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/point_types.h>
#include <pcl/features/normal_3d.h>
#include <vector>
#define NUM_BINS 256

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

std::vector<float> histeq(std::vector<float> val) {
	int count[NUM_BINS];
	memset(count,0,NUM_BINS*sizeof(int));
	int totalCount=0;
	for (size_t i=0;i<val.size();i++) {
		if (val[i] > 0) {
			count[(int)(val[i] * (NUM_BINS-1))]++;
			totalCount++;
		}
	}
	for (int i=1;i<NUM_BINS;i++)
		count[i] += count[i-1];
	std::vector<float> result;
	for (size_t i=0;i<val.size();i++) {
		float t = 1.0 * count[(int)(val[i] * (NUM_BINS-1))] / totalCount;
		result.push_back(t);
	}
	for (size_t i=0;i<NUM_BINS;i++) {
		printf("%d\n",count[i]);
	}
	return result;
}

int main(int argc, char** argv) {
	pcl::PointCloud<pcl::PointXYZ>::Ptr object(new pcl::PointCloud<pcl::PointXYZ>);

	// Read a PCD file from disk.
	if (pcl::io::loadPCDFile<pcl::PointXYZ>(argv[1], *object) != 0)
	{
		return -1;
	}
	char outputFile[128];
	sprintf(outputFile,"%s-curvature.pcd",argv[1]);

	pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
	kdtree.setInputCloud (object);
	int K=50;
	std::vector<int> pointIdxNKNSearch(K);
	std::vector<float> pointNKNSquaredDistance(K);
	std::vector<float> curvatures;
	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> ne;
	float maxc,minc;

	for (size_t i=0;i<object->points.size();i++) {
		float nx,ny,nz,curv;
		pcl::PointXYZ searchPoint = object->points[i];
		kdtree.nearestKSearch (searchPoint, K, pointIdxNKNSearch, pointNKNSquaredDistance);
		ne.computePointNormal(*object, pointIdxNKNSearch, nx,ny,nz, curv);
		curvatures.push_back(curv);
		if (i==0 || curv < minc)
			minc = curv;
		if (i==0 || curv > maxc)
			maxc = curv;
	}

	for (size_t i=0;i<curvatures.size();i++) {
		curvatures[i] = (curvatures[i] - minc) / (maxc - minc);
	} 
	curvatures = histeq(curvatures);

	pcl::PointCloud<pcl::PointXYZRGB>::Ptr output (new pcl::PointCloud<pcl::PointXYZRGB>);
	output->width = object->points.size();
	output->height = 1;
	for (size_t i=0;i<object->points.size();i++) {
		pcl::PointXYZRGB p;
		p.x = object->points[i].x;
		p.y = object->points[i].y;
		p.z = object->points[i].z;
		colormap(curvatures[i],&p.r,&p.g,&p.b);
		output->points.push_back(p);
	}
	pcl::PCDWriter writer;
	writer.write<pcl::PointXYZRGB> (outputFile, *output, false); //*
}
