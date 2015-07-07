CXX       = g++
CXXFLAGS  = -O2 -std=c++11
INCLUDE_DIRS = /usr/local/include/pcl-1.7 /usr/local/include/eigen3
INC = $(INCLUDE_DIRS:%=-I%)
LIB_DIRS = -L/usr/local/lib
LDFLAGS = -lpcl_segmentation -lpcl_io -lpcl_common -lboost_system
LDFLAGS2 = -lpcl_io -lpcl_common -lpcl_features -lpcl_search -lpcl_keypoints -lpcl_recognition -lpcl_visualization -lpcl_kdtree -lboost_system -lvtkCommonDataModel-6.2 -lvtkCommonCore-6.2 -lvtkCommonMath-6.2 -lvtkRenderingCore-6.2
debug: CXXFLAGS = -ggdb3 -O0 -Wall -std=c++11
SRC  = planar_segmentation.cpp
EXE = icp
OBJ = $(patsubst %.cpp,%.o,$(SRC))

#all: $(EXE) bundle2pcd
all: $(EXE)

#$(EXE): $(OBJ)
#	$(CXX) $(LIB_DIRS) -o $(EXE) $(OBJ) $(LDFLAGS)

bundle2pcd: bundle2pcd.cpp
	$(CXX) $(INC) $(LIB_DIRS) -o $@ $< $(LDFLAGS)

hall_segment: hall_segment.cpp
	$(CXX) $(INC) $(LIB_DIRS) -o $@ $< $(LDFLAGS)

correspondence_grouping: correspondence_grouping.cpp
	$(CXX) $(INC) -I/usr/local/include/vtk-6.2/ $(LIB_DIRS) -o $@ $< $(LDFLAGS2)
	
icp: icp.o pcd.o pcd.h kdtree.o kdtree.h descriptor.h descriptor.o normal.h normal.o hashtable.h hashtable.o
	$(CXX) -o $@ icp.o pcd.o kdtree.o descriptor.o normal.o hashtable.o -llapack

fpf2jpg: fpf2jpg.cpp
	$(CXX) -o $@ $< -ljpeg

Viewer: Viewer.cpp Viewer.h
	$(CXX) -I/usr/local/include/eigen3 -I/usr/include/SDL -o $@ $< -lSDL -lGL -lGLU 

sdlViewer: sdlViewer.o sdlViewer.h pcd.o pcd.h kdtree.o kdtree.h
	$(CXX) -ggdb3 -o $@ sdlViewer.o pcd.o kdtree.o -lSDL

debug: icp

single: all

clean:
	rm -f $(EXE) $(OBJ) *.o .depend

.PHONY: clean

#.cpp.o:
#	$(CXX) $(CXXFLAGS) $(INC) -c $<

.depend:
	$(CXX) $(CXXFLAGS) $(INC) -MM *.cpp > .depend

-include .depend
