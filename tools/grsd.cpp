#include <pcl/io/pcd_io.h>
#include <pcl/features/normal_3d.h>
#include <pcl/filters/voxel_grid.h>
#include "grsd.h"
#define SAVE_RSD 0

char* fn;
 
template <typename PointInT, typename PointNT, typename PointOutT> int
pcl::GRSDEstimation<PointInT, PointNT, PointOutT>::getSimpleType (float min_radius, float max_radius,
                                                                  double min_radius_plane,
                                                                  double max_radius_noise,
                                                                  double min_radius_cylinder,
                                                                  double max_min_radius_diff)
{
  if (min_radius > min_radius_plane)
    return (1); // plane
  else if (max_radius > min_radius_cylinder)
    return (2); // cylinder (rim)
  else if (min_radius < max_radius_noise)
    return (0); // noise/corner
  else if (max_radius - min_radius < max_min_radius_diff)
    return (3); // sphere/corner
  else
    return (4); // edge
}

//////////////////////////////////////////////////////////////////////////////////////////////
template <typename PointInT, typename PointNT, typename PointOutT> void
pcl::GRSDEstimation<PointInT, PointNT, PointOutT>::computeFeature (PointCloudOut &output)
{
  // Check if search_radius_ was set
  if (width_ < 0)
  {
    PCL_ERROR ("[pcl::%s::computeFeature] A voxel cell width needs to be set!\n", getClassName ().c_str ());
    output.width = output.height = 0;
    output.points.clear ();
    return;
  }

  // Create the voxel grid
  PointCloudInPtr cloud_downsampled (new PointCloudIn());
  pcl::VoxelGrid<PointInT> grid;
  grid.setLeafSize (width_, width_, width_);
  grid.setInputCloud (input_);
  grid.setSaveLeafLayout (true); // TODO maybe avoid this using nearest neighbor search
  grid.filter (*cloud_downsampled);

  // Compute RSD
  pcl::PointCloud<pcl::PrincipalRadiiRSD>::Ptr radii (new pcl::PointCloud<pcl::PrincipalRadiiRSD>());
  pcl::RSDEstimation<PointInT, PointNT, pcl::PrincipalRadiiRSD> rsd;
  rsd.setInputCloud (cloud_downsampled);
  rsd.setSearchSurface (input_);
  rsd.setInputNormals (normals_);
  rsd.setRadiusSearch (std::max (search_radius_, std::sqrt (3.0) * width_ / 2));
//  pcl::KdTree<PointInT>::Ptr tree = boost::make_shared<pcl::KdTreeFLANN<PointInT> >();
//  tree->setInputCloud(input_);
//  rsd.setSearchMethod(tree);
  rsd.compute (*radii);
  
  // Save the type of each point
  int NR_CLASS = 5; // TODO make this nicer
  std::vector<int> types (radii->points.size ());
  for (size_t idx = 0; idx < radii->points.size (); ++idx)
    types[idx] = getSimpleType (radii->points[idx].r_min, radii->points[idx].r_max);

  // Get the transitions between surface types between neighbors of occupied cells
  Eigen::MatrixXi transition_matrix = Eigen::MatrixXi::Zero (NR_CLASS + 1, NR_CLASS + 1);
  for (size_t idx = 0; idx < cloud_downsampled->points.size (); ++idx)
  {
    int source_type = types[idx];
    std::vector<int> neighbors = grid.getNeighborCentroidIndices (cloud_downsampled->points[idx], relative_coordinates_all_);
    for (unsigned id_n = 0; id_n < neighbors.size (); id_n++)
    {
      int neighbor_type;
      if (neighbors[id_n] == -1) // empty
        neighbor_type = NR_CLASS;
      else
        neighbor_type = types[neighbors[id_n]];
      transition_matrix (source_type, neighbor_type)++;
    }
  }

  // Save feature values
  output.points.resize (1);
  output.height = output.width = 1;
  int nrf = 0;
  for (int i = 0; i < NR_CLASS + 1; i++)
    for (int j = i; j < NR_CLASS + 1; j++)
      output.points[0].histogram[nrf++] = transition_matrix (i, j) + transition_matrix (j, i);

	//save as colored cloud
#if SAVE_RSD
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_colored(new pcl::PointCloud<pcl::PointXYZRGB>);
	for (size_t i=0;i<cloud_downsampled->points.size();i++) {
		pcl::PointXYZ p = cloud_downsampled->points[i];
		pcl::PointXYZRGB q;
		q.x = p.x;
		q.y = p.y;
		q.z = p.z;
		if (types[i] == 0) {
			q.r=255; q.g=0; q.b=0;
		} else if (types[i] == 3) {
			q.r=0; q.g=0; q.b=255;
		} else if (types[i] == 2) {
			q.r=0; q.g=255; q.b=0;
		} else if (types[i] == 1) {
			q.r=255; q.g=255; q.b=0;
		} else {
			q.r=255; q.g=0; q.b=0;
		}
		cloud_colored->points.push_back(q);
	}
	cloud_colored->width = cloud_downsampled->points.size();
	cloud_colored->height = 1;
	pcl::PCDWriter writer;
	char outputFile[128];
	sprintf(outputFile,"%s-colored.pcd",fn);
	writer.write<pcl::PointXYZRGB> (outputFile, *cloud_colored, false);
	printf("Saved to %s\n",outputFile);
#endif
}

int main(int argc, char** argv) {
	// Object for storing the point cloud.
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
	// Object for storing the normals.
	pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);
	// Object for storing the GRSD descriptors for each point.
	pcl::PointCloud<pcl::GRSDSignature21>::Ptr descriptors(new pcl::PointCloud<pcl::GRSDSignature21>());

	// Read a PCD file from disk.
	if (pcl::io::loadPCDFile<pcl::PointXYZ>(argv[1], *cloud) != 0) {
		return -1;
	}
	char outputFile[128];
	fn = argv[1];
	sprintf(outputFile,"%s-esf.pcd",fn);

	// Note: you would usually perform downsampling now. It has been omitted here
	// for simplicity, but be aware that computation can take a long time.
	double lsize = 0.5;
//	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filtered(new pcl::PointCloud<pcl::PointXYZ>);
//	pcl::VoxelGrid<pcl::PointXYZ> vg;
//	vg.setInputCloud(cloud);
//	vg.setLeafSize(lsize,lsize,lsize);
//	vg.filter(*cloud_filtered);
//	cloud = cloud_filtered;

	// Estimate the normals.
	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> normalEstimation;
	normalEstimation.setInputCloud(cloud);
	normalEstimation.setRadiusSearch(lsize);
	pcl::search::KdTree<pcl::PointXYZ>::Ptr kdtree(new pcl::search::KdTree<pcl::PointXYZ>);
	normalEstimation.setSearchMethod(kdtree);
	normalEstimation.compute(*normals);

	// GRSD estimation object.
	pcl::GRSDEstimation<pcl::PointXYZ, pcl::Normal, pcl::GRSDSignature21> grsd;
	grsd.setInputCloud(cloud);
	grsd.setInputNormals(normals);
	grsd.setSearchMethod(kdtree);
	// Search radius, to look for neighbors. Note: the value given here has to be
	// larger than the radius used to estimate the normals.
	grsd.setRadiusSearch(lsize * 1.7);

	grsd.compute(*descriptors);
	if (descriptors->points.size() > 0) {
		pcl::PCDWriter writer;
		writer.write<pcl::GRSDSignature21> (outputFile, *descriptors, false);
//		printf("Saved to %s\n",outputFile);
	}
}
