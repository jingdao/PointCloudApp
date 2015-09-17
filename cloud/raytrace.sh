#!/bin/bash

inputdir=$1
view_result=false

for f in $inputdir/*.ply
do
	dir=${f%.*}
	mkdir $dir
	rm $dir/*-cloud.pcd
	../ray_tracing $f $dir
done

if $view_result
then
	for f in $inputdir/*.ply
	do
		dir=${f%.*}
		numCloud=`ls $dir | wc -l`
		for i in `seq 0 $((numCloud-1))`
		do
			pcl_viewer $dir/$i-cloud.pcd
		done
	done
fi
