#include <pcl/point_types.h>
#include <pcl/features/fpfh.h>
#include <pcl/io/pcd_io.h>
#include <pcl/features/normal_3d.h>
#define NUM_BIN 20

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

int main(int argc, char** argv) {
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::Normal>::Ptr normals (new pcl::PointCloud<pcl::Normal> ());
	if (pcl::io::loadPCDFile<pcl::PointXYZ>(argv[1], *cloud) != 0)
		return -1;
	voxel_grid(cloud);
	char outputFile[128];
	sprintf(outputFile,"%s-fpfh.pcd",argv[1]);

	// Estimate the normals.
	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> normalEstimation;
	normalEstimation.setInputCloud(cloud);
	normalEstimation.setRadiusSearch(0.2);
	pcl::search::KdTree<pcl::PointXYZ>::Ptr kdtree(new pcl::search::KdTree<pcl::PointXYZ>);
	normalEstimation.setSearchMethod(kdtree);
	normalEstimation.compute(*normals);

	// Create the FPFH estimation class, and pass the input dataset+normals to it
	pcl::FPFHEstimation<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh;
	fpfh.setInputCloud (cloud);
	fpfh.setInputNormals (normals);
	// alternatively, if cloud is of tpe PointNormal, do fpfh.setInputNormals (cloud);
	// Create an empty kdtree representation, and pass it to the FPFH estimation object.
	// Its content will be filled inside the object, based on the given input dataset (as no other search surface is given).
	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ>);
	fpfh.setSearchMethod (tree);

	// Output datasets
	pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhs (new pcl::PointCloud<pcl::FPFHSignature33> ());
	// Use all neighbors in a sphere of radius 5cm
	// IMPORTANT: the radius used here has to be larger than the radius used to estimate the surface normals!!!
	fpfh.setRadiusSearch (0.5);
	// Compute the features
	fpfh.compute (*fpfhs);

	pcl::PCDWriter writer;
	writer.write<pcl::FPFHSignature33> (outputFile, *fpfhs, false);
	printf("%s: saved %lu descriptors\n",outputFile,fpfhs->points.size());

	// fpfhs->points.size () should have the same size as the input cloud->points.size ()*
}
