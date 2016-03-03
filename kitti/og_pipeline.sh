#!/bin/bash

compact_equipment_dir=/home/jd/Documents/PointCloudApp/cloud/compact_equipment
seq_out="8 1 8"
output_dir=$compact_equipment_dir

rm $output_dir/svm_prediction.txt
for i in `seq -w $seq_out`
do
	rm $output_dir/clusters$i/prediction.txt
	numLines=`cat $output_dir/clusters$i/labels.txt | wc -l`
	for n in `seq 0 $((numLines-1))`
	do
		j=$output_dir/clusters$i/$n-cloud.pcd
		~/Documents/PointCloudApp/esfc $j
		~/Documents/PointCloudApp/kitti/part_feature.py $j.og $output_dir/clusters$i/prediction.txt
	done
done

echo "labels 1 2 3 4 5" | cat $output_dir/$clusterDir/prediction.txt > $output_dir/svm_prediction.txt

~/Documents/PointCloudApp/kitti/countScore.py $output_dir $seq_out -i -c -m


