#include <pcl/io/pcd_io.h>
#include <pcl/features/normal_3d.h>
#include <pcl/features/our_cvfh.h>
 
int main(int argc, char** argv) {
	// Cloud for storing the object.
	pcl::PointCloud<pcl::PointXYZ>::Ptr object(new pcl::PointCloud<pcl::PointXYZ>);
	// Object for storing the normals.
	pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);
	// Object for storing the OUR-CVFH descriptors.
	pcl::PointCloud<pcl::VFHSignature308>::Ptr descriptors(new pcl::PointCloud<pcl::VFHSignature308>);

	// Note: you should have performed preprocessing to cluster out the object
	// from the cloud, and save it to this individual file.

	// Read a PCD file from disk.
	if (pcl::io::loadPCDFile<pcl::PointXYZ>(argv[1], *object) != 0)
	{
	return -1;
	}
	char outputFile[128];
	sprintf(outputFile,"%s-esf.pcd",argv[1]);

	// Estimate the normals.
	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> normalEstimation;
	normalEstimation.setInputCloud(object);
	normalEstimation.setRadiusSearch(0.03);
	pcl::search::KdTree<pcl::PointXYZ>::Ptr kdtree(new pcl::search::KdTree<pcl::PointXYZ>);
	normalEstimation.setSearchMethod(kdtree);
	normalEstimation.compute(*normals);

	// OUR-CVFH estimation object.
	pcl::OURCVFHEstimation<pcl::PointXYZ, pcl::Normal, pcl::VFHSignature308> ourcvfh;
	ourcvfh.setInputCloud(object);
	ourcvfh.setInputNormals(normals);
	ourcvfh.setSearchMethod(kdtree);
	ourcvfh.setEPSAngleThreshold(5.0 / 180.0 * M_PI); // 5 degrees.
	ourcvfh.setCurvatureThreshold(1.0);
	ourcvfh.setNormalizeBins(false);
	// Set the minimum axis ratio between the SGURF axes. At the disambiguation phase,
	// this will decide if additional Reference Frames need to be created, if ambiguous.
	ourcvfh.setAxisRatio(0.8);

	ourcvfh.compute(*descriptors);
	if (descriptors->points.size() > 0) {
		pcl::PCDWriter writer;
		writer.write<pcl::VFHSignature308> (outputFile, *descriptors, false);
	}
}
