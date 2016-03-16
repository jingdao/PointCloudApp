#!/bin/bash

compact_equipment_dir=/home/jd/Documents/PointCloudApp/cloud/compact_equipment
cad_dir=/home/jd/Documents/PointCloudApp/cloud/caterpillar/
seq_in="8 1 8"
seq_out="8 1 8"
input_dir=$cad_dir
output_dir=$compact_equipment_dir
clusterDir="clusters$(echo $seq_out | awk '{print $1;}')"
svm_dir=/home/jd/Downloads/libsvm-3.20
direct_clf=false
scaleOption=true
computeDescriptors=true

if $computeDescriptors
then
	if ! $direct_clf
	then
		for i in `seq -w $seq_in`
		do
			numLines=`cat $input_dir/clusters$i/labels.txt | wc -l`
			for n in `seq 0 $((numLines-1))`
			do
				j=$input_dir/clusters$i/$n-cloud.pcd
				~/Documents/PointCloudApp/esfc $j
				~/Documents/PointCloudApp/kitti/part_feature.py $j.og $j-esf.pcd 
			done
			~/Documents/PointCloudApp/kitti/writeLabels.py $input_dir/clusters$i
		done
	fi

	rm $output_dir/svm_prediction.txt
	for i in `seq -w $seq_out`
	do
		rm $output_dir/clusters$i/prediction.txt
		numLines=`cat $output_dir/clusters$i/labels.txt | wc -l`
		for n in `seq 0 $((numLines-1))`
		do
			j=$output_dir/clusters$i/$n-cloud.pcd
			~/Documents/PointCloudApp/esfc $j
			if $direct_clf
			then
				~/Documents/PointCloudApp/kitti/part_feature.py $j.og $output_dir/clusters$i/prediction.txt
			else
				~/Documents/PointCloudApp/kitti/part_feature.py $j.og $j-esf.pcd
			fi
		done
		if ! $direct_clf
		then
			~/Documents/PointCloudApp/kitti/writeLabels.py $output_dir/clusters$i
		fi
	done
fi

if $direct_clf
then
	echo "labels 1 2 3 4 5" | cat $output_dir/$clusterDir/prediction.txt > $output_dir/svm_prediction.txt
else
	rm $output_dir/svm_train_data.txt $output_dir/svm_test_data.txt
	for f in `seq -w $seq_in`
	do
		cat $input_dir/clusters$f/svmdata.txt >> $output_dir/svm_train_data.txt
	done
	for f in `seq -w $seq_out` 
	do
		cat $output_dir/clusters$f/svmdata.txt >> $output_dir/svm_test_data.txt
	done
	if $scaleOption
	then
		$svm_dir/svm-scale -l 0 -u 1 -s $output_dir/range.txt $output_dir/svm_train_data.txt > $output_dir/svm_train_scaled.txt 2>/dev/null
		$svm_dir/svm-scale -r $output_dir/range.txt $output_dir/svm_test_data.txt > $output_dir/svm_test_scaled.txt 2>/dev/null
	else
		cp $output_dir/svm_train_data.txt $output_dir/svm_test_scaled.txt
		cp $output_dir/svm_test_data.txt $output_dir/svm_test_scaled.txt
	fi
#	~/Documents/PointCloudApp/kitti/knn.py 10 $output_dir/svm_train_scaled.txt $output_dir/svm_test_scaled.txt $output_dir/svm_prediction.txt
	~/Documents/PointCloudApp/kitti/svc.py $output_dir
#	~/Documents/PointCloudApp/kitti/lda.py $output_dir
#	~/Documents/PointCloudApp/kitti/logistic.py $output_dir
	cat $output_dir/svm_prediction.txt | tail -n +2 > $output_dir/$clusterDir/prediction.txt
fi

~/Documents/PointCloudApp/kitti/countScore.py $output_dir $seq_out -i -c -m


