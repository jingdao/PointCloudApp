#!/bin/bash

start=0000000000
end=0000000800
increment=50
velodyne_dir=/home/jd/Desktop/kitti2/velodyne_points/data
svm_dir=/home/jd/Downloads/libsvm-3.20
cluster_dir=/home/jd/Desktop/kitti2
compute_cluster=false
write_labels=false

if $compute_cluster
then
	for i in `seq -w $start $increment $end`
	do
		mkdir $cluster_dir/clusters$i
		rm $cluster_dir/clusters$i/*-cloud.pcd
		#read velodyne data and cluster
		../main $velodyne_dir/$i.bin $cluster_dir/clusters$i/ &
	done
	wait
fi

if $write_labels
then
	for i in `seq -w $start $increment $end`
	do
		echo "Processing clusters$i ..."

		#user input for labels
		../tools/viz $cluster_dir/clusters$i > $cluster_dir/clusters$i/labels.txt

		#compute descriptors
		for j in $cluster_dir/clusters$i/*-cloud.pcd
		do
			../tools/esf $j $j-esf.pcd
		done

		#write labels
		./writeLabels.py $cluster_dir/clusters$i/
	done
fi

#visualize data
for i in `seq -w $start $increment $end`
do
	cp $cluster_dir/clusters$i/labels.txt $cluster_dir/clusters$i/prediction.txt
	../main $cluster_dir/clusters$i/ $cluster_dir/clusters$i/combined.pcd
	pcl_viewer $cluster_dir/clusters$i/combined.pcd
done
