#include <pcl/point_types.h>
#include <pcl/features/usc.h>
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
	sprintf(outputFile,"%s-usc.pcd",argv[1]);

	// Estimate the normals.
	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> normalEstimation;
	normalEstimation.setInputCloud(cloud);
	normalEstimation.setRadiusSearch(0.2);
	pcl::search::KdTree<pcl::PointXYZ>::Ptr kdtree(new pcl::search::KdTree<pcl::PointXYZ>);
	normalEstimation.setSearchMethod(kdtree);
	normalEstimation.compute(*normals);

	pcl::UniqueShapeContext<pcl::PointXYZ, pcl::UniqueShapeContext1960> uscd;
	pcl::PointCloud<pcl::UniqueShapeContext1960>::Ptr desc (new pcl::PointCloud<pcl::UniqueShapeContext1960> ());
	uscd.setInputCloud (cloud);
	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ>);
	uscd.setSearchMethod (tree);
	uscd.setRadiusSearch (0.5);
	uscd.setMinimalRadius (0.05);
	uscd.setPointDensityRadius (0.1);
	uscd.setLocalRadius (0.5);
	uscd.compute (*desc);

	pcl::PCDWriter writer;
	writer.write<pcl::UniqueShapeContext1960> (outputFile, *desc, false);
	printf("%s: saved %lu descriptors\n",outputFile,desc->points.size());

	// fpfhs->points.size () should have the same size as the input cloud->points.size ()*
}
