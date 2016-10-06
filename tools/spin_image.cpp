#include <stdio.h>
#include <iostream>
#include <vector>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/features/spin_image.h>
#include <pcl/features/normal_3d.h>
#define NUM_BIN 50

void voxel_grid(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud) {
	float minX = cloud->points[0].x;
	float maxX = cloud->points[0].x;
	float minY = cloud->points[0].y;
	float maxY = cloud->points[0].y;
	float minZ = cloud->points[0].z;
	float maxZ = cloud->points[0].z;
	std::vector<pcl::PointXYZ> v;
	for (size_t i=1;i<cloud->points.size();i++) {
		pcl::PointXYZ p = cloud->points[i];
		if (p.x < minX) minX = p.x;
		if (p.x > maxX) maxX = p.x;
		if (p.y < minY) minY = p.y;
		if (p.y > maxY) maxY = p.y;
		if (p.z < minZ) minZ = p.z;
		if (p.z > maxZ) maxZ = p.z;
		v.push_back(p);
	}
	float midX = (minX + maxX) / 2;
	float midY = (minY + maxY) / 2;
	float midZ = (minZ + maxZ) / 2;
	float extent = maxX - minX;
	if (maxY - minY > extent) extent = maxY - minY;
	if (maxZ - minZ > extent) extent = maxZ - minZ;
	bool *** occupied = new bool**[NUM_BIN];
	for (int i=0;i<NUM_BIN;i++) {
		occupied[i] = new bool*[NUM_BIN];
		for (int j=0;j<NUM_BIN;j++) {
			occupied[i][j] = new bool[NUM_BIN]();
		}
	}
	cloud->points.clear();
	for (size_t i=0;i<v.size();i++) {
		pcl::PointXYZ p = v[i];
		float xf = (p.x - midX) / extent;
		float yf = (p.y - midY) / extent;
		float zf = (p.z - midZ) / extent;
		int xi = xf>=0.5 ? NUM_BIN-1 : int((xf + 0.5) * NUM_BIN);
		int yi = yf>=0.5 ? NUM_BIN-1 : int((yf + 0.5) * NUM_BIN);
		int zi = zf>=0.5 ? NUM_BIN-1 : int((zf + 0.5) * NUM_BIN);
		if (!occupied[xi][yi][zi]) {
			occupied[xi][yi][zi] = true;
			p.x = xf;
			p.y = yf;
			p.z = yf;
			cloud->points.push_back(p);
		}
	}
	cloud->width = cloud->points.size();
}

int main (int argc, char** argv)
{
	std::string fileName = argv[1];
	std::cout << "Reading " << fileName << std::endl;
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PCDWriter writer;

	if (pcl::io::loadPCDFile <pcl::PointXYZ> (fileName.c_str(), *cloud) == -1)
	// load the file
	{
	PCL_ERROR ("Couldn't read file");
	return (-1);
	}
	voxel_grid(cloud);
//	writer.write<pcl::PointXYZ> ("grid.pcd", *cloud, false);
//	std::cout << "Loaded " << cloud->points.size() << " points." << std::endl;

	// Compute the normals
	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> normalEstimation;
	normalEstimation.setInputCloud (cloud);

	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ>);
	normalEstimation.setSearchMethod (tree);

	pcl::PointCloud<pcl::Normal>::Ptr normals (new pcl::PointCloud< pcl::Normal>);
	normalEstimation.setRadiusSearch (0.2);
	normalEstimation.compute (*normals);

	// Setup the spin images computation
	typedef pcl::Histogram<153> SpinImage;
	pcl::SpinImageEstimation<pcl::PointXYZ, pcl::Normal, SpinImage> spinImageEstimation(8, 0.5, 16);
	spinImageEstimation.setInputCloud (cloud);
	spinImageEstimation.setInputNormals (normals);

	// Use the same KdTree from the normal estimation
	spinImageEstimation.setSearchMethod (tree);
	pcl::PointCloud<SpinImage>::Ptr spinImages(new pcl::PointCloud<SpinImage>);
	spinImageEstimation.setRadiusSearch (0.5);
//	spinImageEstimation.setKSearch(10);

	// Actually compute the spin images
	spinImageEstimation.compute (*spinImages);
//	std::cout << "SI output points.size (): " << spinImages->points.size () << std::endl;

	// Display and retrieve the spin image descriptor vector for the first point.
//	SpinImage descriptor = spinImages->points[0];
//	std::cout << descriptor << std::endl;
	char outputFile[128];
	sprintf(outputFile,"%s-spin.pcd",argv[1]);
	FILE* f = fopen(outputFile,"w");
	fprintf(f,"# .PCD v0.7 - Point Cloud Data file format\n"
	"VERSION 0.7\n"
	"FIELDS SpinImage\n"
	"SIZE 4\n"
	"TYPE F\n"
	"COUNT 153\n"
	"WIDTH %lu\n"
	"HEIGHT 1\n"
	"VIEWPOINT 0 0 0 1 0 0 0\n"
	"POINTS %lu\n"
	"DATA ascii\n",spinImages->points.size(),spinImages->points.size());
	for (size_t i=0;i<spinImages->points.size();i++) {
		for (int j=0;j<153;j++)
			fprintf(f,"%f ",spinImages->points[i].histogram[j]);
		fprintf(f,"\n");
	}
	fclose(f);
	printf("%s: saved %lu descriptors\n",outputFile,spinImages->points.size());

	return 0;
}
