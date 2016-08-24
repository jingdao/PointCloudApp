CXX       = g++
CXXFLAGS  = -O2 -std=c++11
INCLUDE_DIRS = /usr/local/include/pcl-1.7 /usr/local/include/eigen3
INC = $(INCLUDE_DIRS:%=-I%)
LIB_DIRS = -L/usr/local/lib
LDFLAGS = -lpcl_segmentation -lpcl_io -lpcl_common -lboost_system
LDFLAGS2 = -lpcl_io -lpcl_common -lpcl_features -lpcl_search -lpcl_keypoints -lpcl_recognition -lpcl_visualization -lpcl_kdtree -lboost_system -lvtkCommonDataModel-6.2 -lvtkCommonCore-6.2 -lvtkCommonMath-6.2 -lvtkRenderingCore-6.2
debug: CXXFLAGS = -ggdb3 -O0 -Wall -std=c++11
SRC  = planar_segmentation.cpp
EXE = main
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
	
main: main.o pcd.o pcd.h kdtree.o kdtree.h descriptor.h descriptor.o normal.h normal.o hashtable.h hashtable.o
	$(CXX) -o $@ main.o pcd.o kdtree.o descriptor.o normal.o hashtable.o -llapack

fpf2jpg: fpf2jpg.cpp
	$(CXX) -o $@ $< -ljpeg

Viewer: Viewer.cpp Viewer.h
	$(CXX) -I/usr/local/include/eigen3 -I/usr/include/SDL -o $@ $< -lSDL -lGL -lGLU 

sdlViewer: sdlViewer.o sdlViewer.h pcd.o pcd.h kdtree.o kdtree.h
	$(CXX) -ggdb3 -o $@ sdlViewer.o pcd.o kdtree.o -lSDL
	
hashtable: hashtable.h hashtable.o
	$(CXX) -ggdb3 -o $@ hashtable.o

hole_detector: hole_detector.cpp pcd.cpp pcd.h kdtree.cpp kdtree.h descriptor.h descriptor.cpp normal.h normal.cpp hashtable.h hashtable.cpp
	$(CXX) -std=c++11 -ggdb3 -o $@ hole_detector.cpp pcd.cpp kdtree.cpp descriptor.cpp normal.cpp hashtable.cpp -llapack

ray_tracing: ray_tracing.cpp
	$(CXX) -O2 -o $@ ray_tracing.cpp

hole_tracing: hole_tracing.cpp
	$(CXX) -O2 -o $@ hole_tracing.cpp

semantics: semantics.cpp
	$(CXX) -ggdb3 -o $@ semantics.cpp -lSDL -lGL -lGLU -llapack -lfreetype -I/usr/include/freetype2

rd_desc: rd_desc.cpp
	$(CXX) -ggdb3 -o $@ $<

hpcd: hpcd.cpp
	$(CXX) -O3 -Wall -o $@ $<

scan_eval: scan_eval.cpp
	$(CXX) -ggdb3 -Wall -o $@ $<
	
pad3d: pad3d.cpp
	$(CXX) -ggdb3 -Wall -o $@ $< -llapack

occ_grid: occ_grid.cpp
	$(CXX) -ggdb3 -o $@ $< -lSDL -lGL -lGLU 

size_filter: size_filter.cpp
	$(CXX) -ggdb3 -o $@ $< -llapack

zbuffer: zbuffer.cpp
	$(CXX) -ggdb3 -o $@ $< -lGLU -lOSMesa

bow: bow.cpp
	$(CXX) -ggdb3 -o $@ $<

debug: main

single: all

clean:
	rm -f $(EXE) $(OBJ) *.o .depend

.PHONY: clean

#.cpp.o:
#	$(CXX) $(CXXFLAGS) $(INC) -c $<

.depend:
	$(CXX) $(CXXFLAGS) $(INC) -MM *.cpp > .depend

-include .depend
