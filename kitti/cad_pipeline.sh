#!/bin/bash

start=0000000000
end=0000000000
svm_dir=/home/jd/Downloads/libsvm-3.20
script_dir=/home/jd/Documents/PointCloudApp/kitti
cad_dir=/home/jd/Documents/PointCloudApp/cloud/cad

computeDescriptors=true
scaleOption=true
parameterOption=true
ignoreZeroOption=false
knnOption=true
kernel=0
svm_type=0
k_parameter=3

#compute descriptors
if $computeDescriptors
then
	for j in $cad_dir/*-cloud.pcd
	do
		$script_dir/../tools/esf $j $j-esf.pcd
	done
	$script_dir/writeLabels.py $cad_dir
fi

#divide into training and testing data
rm svm_train_data.txt svm_test_data.txt
for f in `seq -w $start $end` 
do
	cat clusters$f/svmdata.txt >> svm_test_data.txt
done
cp $cad_dir/svmdata.txt svm_train_data.txt

#filter data
if $ignoreZeroOption
then
	cat svm_train_data.txt | grep -v ^0 > svm_train_filtered.txt
	mv svm_train_filtered.txt svm_train_data.txt
	cat svm_test_data.txt | grep -v ^0 > svm_test_filtered.txt
	mv svm_test_filtered.txt svm_test_data.txt
fi

#check data
if ! $svm_dir/tools/checkdata.py svm_train_data.txt || ! $svm_dir/tools/checkdata.py svm_test_data.txt
then
	exit 1
fi
echo "Initialize SVM: train=$(wc -l < svm_train_data.txt) test=$(wc -l < svm_test_data.txt) ..."

#scale data
if $scaleOption
then
	$svm_dir/svm-scale -l 0 -u 1 -s range.txt svm_train_data.txt > svm_train_scaled.txt
	$svm_dir/svm-scale -r range.txt svm_test_data.txt > svm_test_scaled.txt
else
	cp svm_train_data.txt svm_test_scaled.txt
	cp svm_test_data.txt svm_test_scaled.txt
fi

if $knnOption
then
	#use K Nearest Neighbors algorithm
	$script_dir/knn.py $k_parameter svm_train_scaled.txt svm_test_scaled.txt svm_prediction.txt
else
	#train classifier
	if $parameterOption
	then
		parameters=`$svm_dir/tools/grid.py -gnuplot null -out null -t $kernel svm_train_scaled.txt | tail -1`
		c_param=`echo $parameters | cut -d' ' -f1`
		g_param=`echo $parameters | cut -d' ' -f2`
		echo "SVM Train: c=$c_param g=$g_param ..."
		$svm_dir/svm-train -c $c_param -g $g_param -s $svm_type -t $kernel -b 1 svm_train_scaled.txt model.txt > /dev/null
	else
		$svm_dir/svm-train -s $svm_type -t $kernel -b 1 svm_train_scaled.txt model.txt > /dev/null
	fi

	#calculate output
	echo "SVM Test: (Saving to svm_prediction.txt) ..."
	$svm_dir/svm-predict -b 1 svm_test_scaled.txt model.txt svm_prediction.txt
fi

#save predicted values to individual directories
i=0
j=1 #first line is label
while read line
do
	predictions[i]="$line"
	((i++))
done < svm_prediction.txt
for f in `seq -w $start $end`
do
	rm clusters$f/prediction.txt
	numLines=`wc -l < clusters$f/labels.txt`
	read -a labels <<< `cat clusters$f/labels.txt`
	echo "Updating clusters$f/prediction.txt ..."
	for l in ${labels[@]}
	do
		if $ignoreZeroOption && [ "$l" -eq 0 ]
		then
			echo 0 >> clusters$f/prediction.txt
		else
			echo ${predictions[j]} >> clusters$f/prediction.txt
			((j++))
		fi
	done
	$script_dir/../main clusters$f/ clusters$f/combined.pcd
done

#count score
if $ignoreZeroOption
then
	$script_dir/countScore.py . -i -c
else
	$script_dir/countScore.py . -c
fi
