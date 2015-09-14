#include <iostream>

#include <boost/thread/thread.hpp>
#include <pcl/common/common_headers.h>
#include <pcl/features/normal_3d.h>
#include <pcl/io/pcd_io.h>
#include <pcl/visualization/pcl_visualizer.h>

int idx = 0;
int v1(0);
int v2(0);
pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud1(new pcl::PointCloud<pcl::PointXYZRGB>);
pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud2(new pcl::PointCloud<pcl::PointXYZRGB>);
char buf[256];
char* buffer_c;

void keyboardEventOccurred (const pcl::visualization::KeyboardEvent &event, void* viewer_void) {
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer = 
		*static_cast<boost::shared_ptr<pcl::visualization::PCLVisualizer> *> (viewer_void);
	if (event.keyDown ()) {
		if (event.getKeySym() == "v") {
			std::cout << '0' <<std::endl;
		} else if (event.getKeySym() == "b") {
			std::cout << '1' <<std::endl;
		} else if (event.getKeySym() == "n") {
			std::cout << '2' <<std::endl;
		} else if (event.getKeySym() == "m") {
			std::cout << '3' <<std::endl;
		} else if (event.getKeySym() == "comma") {
			std::cout << '4' <<std::endl;
		} else if (event.getKeySym() == "period") {
			std::cout << '5' <<std::endl;
		} else if (event.getKeySym() == "slash") {
			std::cout << '6' <<std::endl;
		} else if (event.getKeySym() == "semicolon") {
			std::cout << '7' <<std::endl;
		} else if (event.getKeySym() == "apostrophe") {
			std::cout << '8' <<std::endl;
		}
//		std::cout << event.getKeySym() << std::endl;
		idx++;
		viewer->removeShape("v1 text");
		viewer->removeShape("sample cloud1");
		viewer->removeShape("sample cloud3");
		char text[256];
		snprintf(text,255,"idx: %d",idx);
		viewer->addText(text, 10, 10, "v1 text", v1);
		pcl::PCDReader reader;
		snprintf(buffer_c,64,"%d-cloud.pcd",idx);
		if (reader.read(buf, *cloud1)!=0)
			viewer->close();
		else {
			pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZRGB> red_color(cloud1, 255, 0, 0);
			viewer->addPointCloud<pcl::PointXYZRGB> (cloud1, red_color, "sample cloud1", v1);
			viewer->addPointCloud<pcl::PointXYZRGB> (cloud1, red_color, "sample cloud3", v2);
		}
	}
}

int main(int argc, char* argv[]) {
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer (new pcl::visualization::PCLVisualizer ("3D Viewer"));
//	viewer->initCameraParameters ();

	viewer->createViewPort(0.0, 0.0, 0.5, 1.0, v1);
	viewer->setBackgroundColor (0, 0, 0, v1);
	snprintf(buf,255,"idx: %d",idx);
	viewer->addText(buf, 10, 10, "v1 text", v1);
//	pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb(cloud);
//	viewer->addPointCloud<pcl::PointXYZRGB> (cloud, rgb, "sample cloud1", v1);

	viewer->createViewPort(0.5, 0.0, 1.0, 1.0, v2);
	viewer->setBackgroundColor (0.3, 0.3, 0.3, v2);
//	viewer->addText("Radius: 0.1", 10, 10, "v2 text", v2);
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZRGB> single_color(cloud, 0, 255, 0);
//	viewer->addPointCloud<pcl::PointXYZRGB> (cloud, single_color, "sample cloud2", v2);
	
	char* dir = argv[1];
	int ndir = strlen(dir);
	strncpy(buf,dir,ndir);
	buffer_c = buf + ndir;
	*buffer_c++ = '/';

	pcl::PCDReader reader;
	snprintf(buffer_c,64,"%d-cloud.pcd",idx);
	reader.read(buf, *cloud1);
	snprintf(buffer_c,64,"original.pcd");
	reader.read(buf, *cloud2);
	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZRGB> red_color(cloud1, 255, 0, 0);
	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZRGB> green_color(cloud2, 0, 255, 0);
	viewer->addPointCloud<pcl::PointXYZRGB> (cloud2, green_color, "sample cloud2", v1);
	viewer->addPointCloud<pcl::PointXYZRGB> (cloud1, red_color, "sample cloud1", v1);
	viewer->addPointCloud<pcl::PointXYZRGB> (cloud1, red_color, "sample cloud3", v2);

	viewer->registerKeyboardCallback (keyboardEventOccurred, (void*)&viewer);

	while (!viewer->wasStopped ())
	{
		viewer->spinOnce (100);
		boost::this_thread::sleep (boost::posix_time::microseconds (100000));
	}
}
