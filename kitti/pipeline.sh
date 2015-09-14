#!/bin/bash

start=0000000060
end=0000000200
increment=20
velodyne_dir=/home/jd/Documents/PointCloudApp/cloud/2011_09_28_drive_0002/velodyne_points/data
svm_dir=/home/jd/Downloads/libsvm-3.20
cluster_dir=/home/jd/Documents/PointCloudApp/cloud/2011_09_28_drive_0002/
compute_cluster=true
write_labels=false
viz_data=false

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

if $viz_data
then
	#visualize data
	for i in `seq -w $start $increment $end`
	do
		cp $cluster_dir/clusters$i/labels.txt $cluster_dir/clusters$i/prediction.txt
		../main $cluster_dir/clusters$i/ $cluster_dir/clusters$i/combined.pcd
		pcl_viewer $cluster_dir/clusters$i/combined.pcd
	done
fi
