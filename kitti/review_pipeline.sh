#!/bin/bash

train_dir=/home/jd/Documents/PointCloudApp/cloud/review/train
test_home=/home/jd/Documents/PointCloudApp/cloud/review
#test_mod="test_res1 test_res05 test_res01 test_res005 test_res001 test_noise1 test_noise05 test_noise01 test_noise005 test_noise001 test_heavy test_partial test_complete"
test_mod="test"

script_dir=/home/jd/Documents/PointCloudApp/kitti
svm_dir=/home/jd/Downloads/libsvm-3.20
base_dir=/home/jd/Documents/PointCloudApp
#DESC_LIST="spin fpfh shot usc esf vfh ../pad3d"
DESC_LIST="spin fpfh shot usc"
LOCAL_DESC="spin fpfh shot usc"
declare -A K_PARAM
K_PARAM[spin]=8
K_PARAM[fpfh]=8
K_PARAM[shot]=11
K_PARAM[usc]=4
GLOBAL_DESC="esf vfh pad3d"
computeTrainDescriptors=false
computeTestDescriptors=false

if $computeTrainDescriptors
then
	for DESC_EXE in $DESC_LIST
	do
		DESC=`basename $DESC_EXE`
		rm $train_dir/*-$DESC.pcd
		for j in $train_dir/*-cloud.pcd
		do
			$base_dir/tools/$DESC_EXE $j >/dev/null
		done
	done
fi

for t in $test_mod
do
	test_dir=$test_home/$t
	echo $test_dir
	if $computeTestDescriptors
	then
		for DESC_EXE in $DESC_LIST
		do
			DESC=`basename $DESC_EXE`
			rm $test_dir/*-$DESC.pcd
			for j in $test_dir/*-cloud.pcd
			do
				$base_dir/tools/$DESC_EXE $j >/dev/null
			done
		done
	fi

	for DESC in $LOCAL_DESC
	do
		for c in `seq 2 5`
		do
			$base_dir/bow $train_dir $test_dir $DESC -n ${K_PARAM[$DESC]} -c $c
		done
	done

#	for DESC in $GLOBAL_DESC
#	do
#		$script_dir/writeLabels.py $train_dir $DESC 200 >/dev/null
#		$script_dir/writeLabels.py $test_dir $DESC 54 >/dev/null
#		cp $train_dir/svmdata.txt $test_dir/svm_train_data.txt
#		cp $test_dir/svmdata.txt $test_dir/svm_test_data.txt
#		$svm_dir/svm-scale -l 0 -u 1 -s $test_dir/range.txt $test_dir/svm_train_data.txt > $test_dir/svm_train_scaled.txt 2>/dev/null
#		$svm_dir/svm-scale -r $test_dir/range.txt $test_dir/svm_test_data.txt > $test_dir/svm_test_scaled.txt 2>/dev/null
#		echo $DESC `$script_dir/svc.py $test_dir`
#		cat $test_dir/svm_prediction.txt | tail -n +2 > $test_dir/prediction.txt
#		$script_dir/countScore.py $test_dir
#	done
done

