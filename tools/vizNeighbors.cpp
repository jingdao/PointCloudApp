#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>
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
char cadPath[256];
char* cad_c;
char samplePath[256];
char* sample_c;
std::vector<int> label;
std::vector<int> match;
std::vector<double> dist;

void keyboardEventOccurred (const pcl::visualization::KeyboardEvent &event, void* viewer_void) {
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer = 
		*static_cast<boost::shared_ptr<pcl::visualization::PCLVisualizer> *> (viewer_void);
	if (event.keyDown ()) {
		idx++;
		if (idx >= label.size())
			viewer->close();
		else {
			viewer->removeShape("v1 text");
			viewer->removeShape("sample cloud1");
			viewer->removeShape("sample cloud2");
			char text[256];
			snprintf(text,255,"dist: %f",dist[idx]);
			viewer->addText(text, 50, 50,20,1,1,1, "v1 text", v1);
			pcl::PCDReader reader;
			snprintf(sample_c,64,"%d-cloud.pcd",idx);
			reader.read(samplePath, *cloud1);
			snprintf(cad_c,64,"%d-cloud.pcd",match[idx]);
			reader.read(cadPath, *cloud2);
			snprintf(buffer_c,64,"%d-cloud.pcd",idx);
			pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZRGB> red_color(cloud1, 255, 0, 0);
			pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZRGB> green_color(cloud2, 0, 255, 0);
			viewer->addPointCloud<pcl::PointXYZRGB> (cloud1, red_color, "sample cloud1", v1);
			viewer->addPointCloud<pcl::PointXYZRGB> (cloud2, green_color, "sample cloud2", v2);
		}
	}
}

int main(int argc, char* argv[]) {
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer (new pcl::visualization::PCLVisualizer ("3D Viewer"));
//	viewer->initCameraParameters ();

	viewer->createViewPort(0.0, 0.0, 0.5, 1.0, v1);
	viewer->setBackgroundColor (0, 0, 0, v1);
//	pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb(cloud);
//	viewer->addPointCloud<pcl::PointXYZRGB> (cloud, rgb, "sample cloud1", v1);

	viewer->createViewPort(0.5, 0.0, 1.0, 1.0, v2);
	viewer->setBackgroundColor (0.3, 0.3, 0.3, v2);
//	viewer->addText("Radius: 0.1", 10, 10, "v2 text", v2);
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZRGB> single_color(cloud, 0, 255, 0);
//	viewer->addPointCloud<pcl::PointXYZRGB> (cloud, single_color, "sample cloud2", v2);

	FILE* neighborFile = fopen("neighbors.txt","r");
	while (fgets(buf,256,neighborFile)) {
		buffer_c = buf;
		label.push_back(strtol(buffer_c,&buffer_c,10));
		buffer_c++;
		match.push_back(strtol(buffer_c,&buffer_c,10));
		buffer_c++;
		dist.push_back(strtod(buffer_c,&buffer_c));
	}
	fclose(neighborFile);
	printf("Processed %lu samples\n",label.size());
	
	char* cadDir = "/home/jd/Documents/PointCloudApp/cloud/psb/";
	int ndir = strlen(cadDir);
	strncpy(cadPath,cadDir,ndir);
	cad_c = cadPath + ndir;
	char* sampleDir = "clusters07/";
	ndir = strlen(sampleDir);
	strncpy(samplePath,sampleDir,ndir);
	sample_c = samplePath + ndir;

	pcl::PCDReader reader;
	snprintf(sample_c,64,"%d-cloud.pcd",idx);
	reader.read(samplePath, *cloud1);
	snprintf(cad_c,64,"%d-cloud.pcd",match[idx]);
	reader.read(cadPath, *cloud2);
	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZRGB> red_color(cloud1, 255, 0, 0);
	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZRGB> green_color(cloud2, 0, 255, 0);
	viewer->addPointCloud<pcl::PointXYZRGB> (cloud1, red_color, "sample cloud1", v1);
	viewer->addPointCloud<pcl::PointXYZRGB> (cloud2, green_color, "sample cloud2", v2);
	snprintf(buf,255,"dist: %f",dist[idx]);
	viewer->addText(buf, 50, 50,20,255,255,255, "v1 text", v1);

	viewer->registerKeyboardCallback (keyboardEventOccurred, (void*)&viewer);

	while (!viewer->wasStopped ())
	{
		viewer->spinOnce (100);
		boost::this_thread::sleep (boost::posix_time::microseconds (100000));
	}
}
