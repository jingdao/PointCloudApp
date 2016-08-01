#!/bin/bash

input_dir=/home/jd/Documents/PointCloudApp/cloud/assets/crane
output_dir=/home/jd/Documents/PointCloudApp/cloud/assets/crane2
cad_dir=/home/jd/Documents/PointCloudApp/cloud/caterpillar/truck/mixed
script_dir=/home/jd/Documents/PointCloudApp/kitti
svm_dir=/home/jd/Downloads/libsvm-3.20
PAD=/home/jd/Documents/PointCloudApp/pad3d
SIZE=/home/jd/Documents/PointCloudApp/size_filter
computeDescriptors=false
computeFilter=false

if $computeFilter
then
	rm $output_dir/*.pcd
	$SIZE $cad_dir $input_dir $output_dir
fi

if $computeDescriptors
then
	rm $cad_dir/*-pad.pcd
	for j in $cad_dir/*-cloud.pcd
	do
		$PAD $j
	done
	$script_dir/writeLabels.py $cad_dir

	rm $output_dir/*-pad.pcd
	for j in $output_dir/*-cloud.pcd
	do
		$PAD $j
	done
	$script_dir/writeLabels.py $output_dir
fi

rm $output_dir/svm_train_data.txt $output_dir/svm_test_data.txt
cp $cad_dir/svmdata.txt $output_dir/svm_train_data.txt
cp $output_dir/svmdata.txt $output_dir/svm_test_data.txt
$svm_dir/svm-scale -l 0 -u 1 -s $output_dir/range.txt $output_dir/svm_train_data.txt > $output_dir/svm_train_scaled.txt 2>/dev/null
$svm_dir/svm-scale -r $output_dir/range.txt $output_dir/svm_test_data.txt > $output_dir/svm_test_scaled.txt 2>/dev/null
$script_dir/svc.py $output_dir
cat $output_dir/svm_prediction.txt | tail -n +2 > $output_dir/prediction.txt

$script_dir/countScore.py $output_dir
