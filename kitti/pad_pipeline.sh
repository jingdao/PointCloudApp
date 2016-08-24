#!/bin/bash

input_dir=/home/jd/Documents/PointCloudApp/cloud/assets/crane
#output_dir=/home/jd/Documents/PointCloudApp/cloud/compact_equipment/clusters8
output_dir=/home/jd/Documents/PointCloudApp/cloud/review/test
#cad_dir=/home/jd/Documents/PointCloudApp/cloud/caterpillar/mixed_pose
cad_dir=/home/jd/Documents/PointCloudApp/cloud/review/train
script_dir=/home/jd/Documents/PointCloudApp/kitti
svm_dir=/home/jd/Downloads/libsvm-3.20
#DESC_EXE=/home/jd/Documents/PointCloudApp/pad3d
DESC_EXE=/home/jd/Documents/PointCloudApp/tools/esf
DESC=`basename $DESC_EXE`
SIZE=/home/jd/Documents/PointCloudApp/size_filter
computeDescriptors=true
computeFilter=false

if $computeFilter
then
	rm $output_dir/*.pcd
	$SIZE $cad_dir $input_dir $output_dir
fi

if $computeDescriptors
then
	rm $cad_dir/*-$DESC.pcd
	for j in $cad_dir/*-cloud.pcd
	do
		$DESC_EXE $j
	done
	$script_dir/writeLabels.py $cad_dir $DESC

	rm $output_dir/*-$DESC.pcd
	for j in $output_dir/*-cloud.pcd
	do
		$DESC_EXE $j
	done
	$script_dir/writeLabels.py $output_dir $DESC
fi

rm $output_dir/svm_train_data.txt $output_dir/svm_test_data.txt
cp $cad_dir/svmdata.txt $output_dir/svm_train_data.txt
cp $output_dir/svmdata.txt $output_dir/svm_test_data.txt
$svm_dir/svm-scale -l 0 -u 1 -s $output_dir/range.txt $output_dir/svm_train_data.txt > $output_dir/svm_train_scaled.txt 2>/dev/null
$svm_dir/svm-scale -r $output_dir/range.txt $output_dir/svm_test_data.txt > $output_dir/svm_test_scaled.txt 2>/dev/null
$script_dir/svc.py $output_dir
cat $output_dir/svm_prediction.txt | tail -n +2 > $output_dir/prediction.txt

$script_dir/countScore.py $output_dir
