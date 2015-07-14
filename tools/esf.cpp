#include <pcl/io/pcd_io.h>
#include <pcl/features/esf.h>

typedef pcl::ESFSignature640 DescriptorType;

int
main(int argc, char** argv)
{
	pcl::PointCloud<pcl::PointXYZ>::Ptr object(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<DescriptorType>::Ptr descriptors(new pcl::PointCloud<DescriptorType>);

	// Read a PCD file from disk.
	if (pcl::io::loadPCDFile<pcl::PointXYZ>(argv[1], *object) != 0)
	{
		return -1;
	}
	char* outputFile = argv[2];

	pcl::ESFEstimation<pcl::PointXYZ, pcl::ESFSignature640> esf;
	esf.setInputCloud(object);
	esf.compute(*descriptors);

	if (descriptors->points.size() > 0) {
		pcl::PCDWriter writer;
		writer.write<DescriptorType> (outputFile, *descriptors, false); //*
	}
}
