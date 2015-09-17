#include <stdio.h>
#include <vector>
#include <math.h>
#include "pcd.h"

//http://stackoverflow.com/questions/4279478/largest-circle-inside-a-non-convex-polygon

#define GRID_LENGTH 2048

struct Point {
	float x,y;
};

struct Coordinates {
	int x,y;
};

struct Vertex {
	int left,right;
};

typedef std::vector<Point> Polygon;

void writeToPBM(char* filename, bool* arr, int width, int height) {
	FILE* f = fopen(filename,"w");
	if (!f) return;
	fprintf(f,"P1\n");
	fprintf(f,"%d %d\n",width,height);
	for (int i=0;i<height;i++) {
		for (int j=0;j<width;j++) {
			if (*arr)
				fprintf(f,"1 ");
			else
				fprintf(f,"0 ");
			arr++;
		}
		fprintf(f,"\n");
	}
	fclose(f);

}

void writeComponentsToPBM(std::vector<Polygon> polygons,int width,int height) {
	char buffer[64];
	for (size_t i=0;i<polygons.size();i++) {
		snprintf(buffer,64,"component%lu.pbm",i);
		bool* componentImage = new bool[width*height]();
		for (size_t j=0;j<polygons[i].size();j++) {
			Point p = polygons[i][j];
			int k = (int)(p.y * width + p.x);
			componentImage[k] = true;
		}
		writeToPBM(buffer,componentImage,width,height);
		delete componentImage;
	}
}

void erosion(bool* input, bool* output, int width, int height) {
	for (int i=0;i<height;i++) {
		for (int j=0;j<width;j++) {
			if (i==0 || i==height-1 || j==0 || j==width-1) {
				*output = false;
			} else if (*input && *(input-1) && *(input+1) && *(input-width) && *(input+width)) {
				*output = true;
			} else {
				*output = false;
			}
			input++;
			output++;
		}
	}
}

void dilation(bool* input, bool* output, int width, int height) {
	for (int i=0;i<height;i++) {
		for (int j=0;j<width;j++) {
			if (i==0 || i==height-1 || j==0 || j==width-1) {
				*output = *input;
			} else if (*input || *(input-1) || *(input+1) || *(input-width) || *(input+width)) {
				*output = true;
			} else {
				*output = false;
			}
			input++;
			output++;
		}
	}
}

void opening(int k,bool* input, bool* output, int width, int height) {
	bool* tmp = new bool[width*height];
	memcpy(tmp,input,width*height*sizeof(bool));
	for (int i=0;i<k;i++) {
		erosion(tmp,output,width,height);
		memcpy(tmp,output,width*height*sizeof(bool));
	}
	for (int i=0;i<k;i++) {
		dilation(tmp,output,width,height);
		memcpy(tmp,output,width*height*sizeof(bool));
	}
	delete[] tmp;
}

void closing(int k,bool* input, bool* output, int width, int height) {
	bool* tmp = new bool[width*height];
	memcpy(tmp,input,width*height*sizeof(bool));
	for (int i=0;i<k;i++) {
		dilation(tmp,output,width,height);
		memcpy(tmp,output,width*height*sizeof(bool));
	}
	for (int i=0;i<k;i++) {
		erosion(tmp,output,width,height);
		memcpy(tmp,output,width*height*sizeof(bool));
	}
	delete[] tmp;
}

void sobel(bool* input,bool* output, int width, int height) {
	for (int i=0;i<height;i++) {
		for (int j=0;j<width;j++) {
			if (i==0 || i==height-1 || j==0 || j==width-1) {
				*output = false;
			} else {
				int gx = 0, gy = 0;
				if (*(input - width - 1)) { gx-=1; gy-=1;}
				if (*(input - width)) gx-=2;
				if (*(input - width + 1)) { gx-=1; gy+=1;}
				if (*(input - 1)) gy-=2;
				if (*(input + 1)) gy+=2;
				if (*(input + width - 1)) { gx+=1; gy-=1;}
				if (*(input + width)) gx+=2;
				if (*(input + width + 1)) { gx+=1; gy+=1;}
				if (abs(gx) + abs(gy) > 4)
					*output = true;
				else
					*output = false;
			}
			input++;
			output++;
		}
	}
}

int connected_components(bool* input,int* labels, int width,int height) {
	int idx = 1;
	for (int i=0;i<height;i++) {
		for (int j=0;j<width;j++) {
			if (input[i*width+j])
				continue;
			if (labels[i*width+j])
				continue;
			Coordinates p = {j,i};
			std::vector<Coordinates> stack;
			stack.push_back(p);
			while (!stack.empty()) {
				Coordinates q = stack.back();
				stack.pop_back();
				if (q.x<0 || q.x>=width || q.y<0 || q.y>=height)
					continue;
				int k = (q.y*width+q.x);
				if (input[k])
					continue;
				if (labels[k])
					continue;
				labels[k]=idx;
				Coordinates r1 = {q.x-1,q.y-1};
				Coordinates r2 = {q.x-1,q.y};
				Coordinates r3 = {q.x-1,q.y+1};
				Coordinates r4 = {q.x,q.y-1};
				Coordinates r5 = {q.x,q.y+1};
				Coordinates r6 = {q.x+1,q.y-1};
				Coordinates r7 = {q.x+1,q.y};
				Coordinates r8 = {q.x+1,q.y+1};
				stack.push_back(r1);
				stack.push_back(r2);
				stack.push_back(r3);
				stack.push_back(r4);
				stack.push_back(r5);
				stack.push_back(r6);
				stack.push_back(r7);
				stack.push_back(r8);
			}
			idx++;
		}
	}
	return idx-1;
}

std::vector<Polygon> form_polygons(int numComponents,bool* input,int* labels,int width,int height) {
	std::vector<Polygon> polygons;
	for (int i=2;i<=numComponents;i++) { //ignore background
		Polygon component;
		bool* arr_ptr = input;
		int* label_ptr = labels;
		for (int j=0;j<height;j++) {
			for (int k=0;k<width;k++) {
				if (*arr_ptr && *label_ptr==i) {
					Point p = {(float)k,(float)j};
					component.push_back(p);
				}
				arr_ptr++;
				label_ptr++;
			}
		}
		polygons.push_back(component);
	}
	return polygons;
}

Polygon convex_hull(Polygon polygon) { //reassign pts to sorted order on hull
	Polygon hull;
	Point leftmost = polygon[0];
	size_t leftIndex = 0;
	for (size_t i=1;i<polygon.size();i++) {
		if (polygon[i].x < leftmost.x) {
			leftmost = polygon[i];
			leftIndex = i;
		}
	}
	size_t pointOnHull = leftIndex;
	while (true) {
		hull.push_back(polygon[pointOnHull]);
		size_t endIndex = 0;
		Point endpoint = polygon[0];
		for (size_t j=1;j<polygon.size();j++) {
			if (j == pointOnHull)
				continue;
			else if (endIndex == pointOnHull) {
				endpoint = polygon[j];
				endIndex = j;
			} else {
				//check if polygon[j] is left of line from pointOnHull to endpoint
				double theta1 = atan2(endpoint.y-polygon[pointOnHull].y,endpoint.x-polygon[pointOnHull].x);
				double theta2 = atan2(polygon[j].y-polygon[pointOnHull].y,polygon[j].x-polygon[pointOnHull].x);
				while (theta2 < theta1) theta2 += 2*M_PI;
				if (theta2 - theta1 < M_PI) {
					endpoint = polygon[j];
					endIndex = j;
				}
			}
		}
		pointOnHull = endIndex;
		if (endIndex == leftIndex)
			break;
	}
	return hull;
}

std::vector<Polygon> marching_squares(bool* input, int width, int height) {
	bool* input_ptr = input;
	for (int i=1;i<height;i++) {
		for (int j=1;j<width;j++) {
			int lookup_key = 0;
			lookup_key |= *(input_ptr - width - 1);
			lookup_key |= *(input_ptr - width) << 1;
			lookup_key |= *(input_ptr - 1) << 2;
			lookup_key |= *(input_ptr) << 3;
		}
	}
	std::vector<Polygon> contour;
	return contour;
}

float getPolygonDiameter(Polygon *polygon) {
	float maxDistance = 0;
	for (size_t i=0;i<polygon->size();i++) {
		for (size_t j=i+1;j<polygon->size();j++) {
			float d = 0;
			d += ((*polygon)[i].x - (*polygon)[j].x) * ((*polygon)[i].x - (*polygon)[j].x);
			d += ((*polygon)[i].y - (*polygon)[j].y) * ((*polygon)[i].y - (*polygon)[j].y);
			if (d > maxDistance)
				maxDistance = d;
		}
	}
	return sqrt(maxDistance);
}

int main(int argc, char* argv[]) {

	if (argc < 3) {
		printf("./hole_detector input.pcd output.pbm\n");
		return 1;
	}

	PCD cloud(argv[1]);
	printf("Loaded %s (%d points)\n",argv[1],cloud.numPoints);

	float minX=cloud.float_data[0],maxX=cloud.float_data[0];
	float minY=cloud.float_data[1],maxY=cloud.float_data[1];
	for (int i=1;i<cloud.numPoints;i++) {
		if (cloud.float_data[i * 4] < minX) minX = cloud.float_data[i * 4];
		else if (cloud.float_data[i * 4] > maxX) maxX = cloud.float_data[i * 4];
		if (cloud.float_data[i * 4 + 2] < minY) minY = cloud.float_data[i * 4 + 2];
		else if (cloud.float_data[i * 4 + 2] > maxY) maxY = cloud.float_data[i * 4 + 2];
	}
	float xstep = (maxX - minX) / GRID_LENGTH;
	float ystep = (maxY - minY) / GRID_LENGTH;
	float step;
	int xgrid, ygrid;
	printf("x: %f %f %f\n",minX,maxX,xstep);
	printf("y: %f %f %f\n",minY,maxY,ystep);

	if (xstep > ystep) {
		step = ystep;
		ygrid = GRID_LENGTH;
		xgrid = GRID_LENGTH / ystep * xstep + 1;
	} else {
		step = xstep;
		xgrid = GRID_LENGTH;
		ygrid = GRID_LENGTH / xstep * ystep + 1;
	}

	printf("Grid size %dx%d (%f,%f : %f)\n",xgrid,ygrid,xstep,ystep,step);

	bool* arr = new bool[xgrid * ygrid]();
	bool* arr2 = new bool[xgrid * ygrid]();
	int* labels = new int[xgrid * ygrid]();
	for (int i=0;i<cloud.numPoints;i++) {
		int xcoord = (int) ((cloud.float_data[i * 4] - minX) / step);
		int ycoord = (int) ((cloud.float_data[i * 4 + 2] - minY) / step);
		if (xcoord >= xgrid || ycoord >= ygrid) {
			printf("Error at coords: %f %f %d %d\n",cloud.float_data[i * 4], cloud.float_data[i * 4 + 2],xcoord,ycoord);
			break;
		}
		arr[xcoord + ycoord * xgrid] = true;
	}

	writeToPBM("pointcloud.pbm",arr,xgrid,ygrid);
	//smoothing
	closing(3,arr,arr2,xgrid,ygrid);
	writeToPBM("smoothed.pbm",arr2,xgrid,ygrid);

	//edge detection, labelling, and polygon formation
	sobel(arr2,arr,xgrid,ygrid);
	writeToPBM("edges.pbm",arr,xgrid,ygrid);
	int numComponents = connected_components(arr2,labels,xgrid,ygrid);
	std::vector<Polygon> polygons = form_polygons(numComponents,arr,labels,xgrid,ygrid);
	writeComponentsToPBM(polygons,xgrid,ygrid);
	for (size_t i=0;i<polygons.size();i++) {
		polygons[i] = convex_hull(polygons[i]);
		float diameter = getPolygonDiameter(&(polygons[i]));
		diameter *= step / 1000;
		printf("Diameter: %fm (%lu points)\n",diameter,polygons[i].size());
	}

	writeToPBM(argv[2],arr,xgrid,ygrid);

	delete[] arr;
	delete[] arr2; 
	delete[] labels;
}
