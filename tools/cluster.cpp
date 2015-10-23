#include <pcl/ModelCoefficients.h>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/features/normal_3d.h>
#include <pcl/kdtree/kdtree.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/segmentation/extract_clusters.h>


int main (int argc, char** argv)
{
	// Read in the cloud data
	pcl::PCDReader reader;
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZRGB>), cloud_f (new pcl::PointCloud<pcl::PointXYZRGB>);
	reader.read (argv[1], *cloud);
	std::cout << "PointCloud before filtering has: " << cloud->points.size () << " data points." << std::endl; //*

	// Create the filtering object: downsample the dataset using a leaf size of 1cm
//	pcl::VoxelGrid<pcl::PointXYZRGB> vg;
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_filtered =cloud;
//	vg.setInputCloud (cloud);
//	float lsize = 0.05;
//	vg.setLeafSize (lsize,lsize,lsize);
//	vg.filter (*cloud_filtered);
//	std::cout << "PointCloud after filtering has: " << cloud_filtered->points.size ()  << " data points." << std::endl; //*

	//parse command line options
	int segIterations=100000,minCluster=500,maxCluster=100000,maxClusters=200;
	float segThreshold=1,segRatio=0.5,clusterTolerance=0.3;
	if (argc >= 3) segIterations = atoi(argv[2]);
	if (argc >= 4) segThreshold = strtod(argv[3],NULL);
	if (argc >= 5) segRatio = strtod(argv[4],NULL);
	if (argc >= 6) minCluster = atoi(argv[5]);
	if (argc >= 7) maxCluster = atoi(argv[6]);
	if (argc >= 8) maxClusters = atoi(argv[7]);
	if (argc >= 9) clusterTolerance = strtod(argv[8],NULL);
	printf("segmentation: %d iter %f thres %f ratio\n",segIterations,segThreshold,segRatio);
	printf("clustering: %d min %d max %d total %f thres\n",minCluster,maxCluster,maxClusters,clusterTolerance);

	// Create the segmentation object for the planar model and set all the parameters
	pcl::SACSegmentation<pcl::PointXYZRGB> seg;
	pcl::PointIndices::Ptr inliers (new pcl::PointIndices);
	pcl::ModelCoefficients::Ptr coefficients (new pcl::ModelCoefficients);
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_plane (new pcl::PointCloud<pcl::PointXYZRGB> ());
	pcl::PCDWriter writer;
	seg.setOptimizeCoefficients (true);
	seg.setModelType (pcl::SACMODEL_PLANE);
	seg.setMethodType (pcl::SAC_RANSAC);
	seg.setMaxIterations (segIterations);
	seg.setDistanceThreshold (segThreshold);

	int i=0, nr_points = (int) cloud_filtered->points.size ();
//	while (cloud_filtered->points.size () > segRatio * nr_points)
	for (i=0;i<1;i++)
	{
		// Segment the largest planar component from the remaining cloud
		seg.setInputCloud (cloud_filtered);
		seg.segment (*inliers, *coefficients);
		if (inliers->indices.size () == 0)
		{
		std::cout << "Could not estimate a planar model for the given dataset." << std::endl;
		break;
		}

		// Extract the planar inliers from the input cloud
		pcl::ExtractIndices<pcl::PointXYZRGB> extract;
		extract.setInputCloud (cloud_filtered);
		extract.setIndices (inliers);
		extract.setNegative (false);

		// Get the points associated with the planar surface
		extract.filter (*cloud_plane);
		std::cout << "PointCloud representing the planar component: " << cloud_plane->points.size () << " data points." << std::endl;
		std::cout << coefficients->values[0] << " " << coefficients->values[1] << " " << coefficients->values[2] << " " << coefficients->values[3] << std::endl;

		// Remove the planar inliers, extract the rest
		extract.setNegative (true);
		extract.filter (*cloud_f);
		*cloud_filtered = *cloud_f;
	}

	std::cout << "PointCloud after filtering has: " << cloud_filtered->points.size ()  << " data points." << std::endl; //*
	writer.write<pcl::PointXYZRGB> ("original.pcd", *cloud_filtered, false); //*
	// Creating the KdTree object for the search method of the extraction
	pcl::search::KdTree<pcl::PointXYZRGB>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZRGB>);
	tree->setInputCloud (cloud_filtered);

	std::vector<pcl::PointIndices> cluster_indices;
	pcl::EuclideanClusterExtraction<pcl::PointXYZRGB> ec;
	ec.setClusterTolerance (clusterTolerance);
	ec.setMinClusterSize (minCluster);
	ec.setMaxClusterSize (maxCluster);
	ec.setSearchMethod (tree);
	ec.setInputCloud (cloud_filtered);
	ec.extract (cluster_indices);

	int j = 0;
	for (std::vector<pcl::PointIndices>::const_iterator it = cluster_indices.begin ();
	 	j<maxClusters && it != cluster_indices.end (); ++it)
	{
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_cluster (new pcl::PointCloud<pcl::PointXYZRGB>);
		for (std::vector<int>::const_iterator pit = it->indices.begin (); pit != it->indices.end (); ++pit)
		cloud_cluster->points.push_back (cloud_filtered->points[*pit]); //*
		cloud_cluster->width = cloud_cluster->points.size ();
		cloud_cluster->height = 1;
		cloud_cluster->is_dense = true;

		std::cout << "PointCloud representing the Cluster: " << cloud_cluster->points.size () << " data points." << std::endl;
		std::stringstream ss;
		ss << j << "-cloud.pcd";
		writer.write<pcl::PointXYZRGB> (ss.str (), *cloud_cluster, false); //*
		j++;
	}

	return (0);
}
