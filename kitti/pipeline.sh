#!/bin/bash

start=0000000000
end=0000000010
velodyne_dir=/home/jd/Downloads/kitti_vision/velodyne_points/data
svm_dir=/home/jd/Downloads/libsvm-3.20

for i in `seq -w $start $end`
do
	mkdir clusters$i
	rm clusters$i/*-cloud.pcd
	#read velodyne data and cluster
	../main $velodyne_dir/$i.bin clusters$i/ &
done
wait

for i in `seq -w $start $end`
do
	#user input for labels
	../tools/viz clusters$i > clusters$i/labels.txt

	#compute descriptors
	for j in clusters$i/*-cloud.pcd
	do
		../tools/esf $j $j-esf.pcd
	done

	#write labels
	./writeLabels.py clusters$i/
	#../main clusters$i/ clusters$i/descriptor.pcd
	#./writeDescriptors.py clusters$i/

done

#scale data
$svm_dir/svm-scale -l 0 -u 1 -s range.txt clusters$start/svmdata.txt > clusters$start/svmdata-scale.txt
#train classifier
$svm_dir/svm-train -t 0 clusters$start/svmdata-scale.txt model.txt

#calculate output
for i in `seq -w $start $end`
do
	$svm_dir/svm-scale -r range.txt clusters$i/svmdata.txt > clusters$i/svmdata-scale.txt
	$svm_dir/svm-predict clusters$i/svmdata-scale.txt model.txt clusters$i/prediction.txt
	../main clusters$i/ clusters$i/combined.pcd
done

#view output
for i in `seq -w $start $end`
do
	pcl_viewer clusters$i/combined.pcd
done
