#!/bin/bash

#seq_in="0000000000 50 0000000800"
#seq_in="0000000000 5 0000000120"
#seq_out="0000000000 50 0000000800"
#seq_out="0000000000 5 0000000120"
#seq_in="7 1 7"
#seq_in="5 1 5"
#seq_out="6 1 6"
#seq_out="8 1 8"
#seq_out="1 1 1"
seq_out="11 1 11"
kitti1_dir=/home/jd/Documents/PointCloudApp/cloud/2011_09_29_drive_0071
kitti2_dir=/home/jd/Documents/PointCloudApp/cloud/2011_09_26_drive_0106
kitti3_dir=/home/jd/Documents/PointCloudApp/cloud/2011_09_28_drive_0002
china_dir=/home/jd/Documents/PointCloudApp/cloud/china
big_equipment_dir=/home/jd/Documents/PointCloudApp/cloud/big_equipment
compact_equipment_dir=/home/jd/Documents/PointCloudApp/cloud/compact_equipment
wheel_track_dir=/home/jd/Documents/PointCloudApp/cloud/wheel_vs_track
toy_excavator_dir=/home/jd/Documents/PointCloudApp/cloud/toy_excavator
atlantic_dr_dir=/home/jd/Documents/PointCloudApp/cloud/atlantic_dr
input_dir=$kitti2_dir
output_dir=$atlantic_dr_dir
svm_dir=/home/jd/Downloads/libsvm-3.20
script_dir=/home/jd/Documents/PointCloudApp/kitti
cad_dir=/home/jd/Documents/PointCloudApp/cloud/caterpillar/mixed_pose
ESFO=/home/jd/Documents/PointCloudApp/tools/esf
VFH=/home/jd/Documents/PointCloudApp/tools/vfh
OURCVFH=/home/jd/Documents/PointCloudApp/tools/our_cvfh
GRSD=/home/jd/Documents/PointCloudApp/tools/grsd
ESFC=/home/jd/Documents/PointCloudApp/esfc
ESF=$ESFC

useCAD=true
computeDescriptors=true
scaleOption=true
parameterOption=true
ignoreZeroOption=false
combineAssembly=false
assemblies="{3,4}"
#assemblies="3"

knnOption=false
ldaOption=false
lrOption=false
bayesOption=false
svcOption=true
adaboostOption=false
decisionTreeOption=false
savePointCloud=false
kernel=0
svm_type=0
k_parameter=10

if $useCAD
then
	input_dir=$cad_dir
	if $computeDescriptors
	then
		rm $cad_dir/*-part.pcd
		rm $cad_dir/*-esf.pcd
		if $combineAssembly
		then
			for j in $cad_dir/*-cloud.pcd
			do
				$ESFC $j
			done
			for j in $cad_dir/*-part.pcd
			do
				$ESF $j
			done
			numLines=`cat $cad_dir/labels.txt | wc -l`
			for j in `seq 0 $((numLines-1))`
			do
				eval ~/Documents/PointCloudApp/cloud/concat_esf.py $cad_dir/$j-cloud.pcd-$assemblies-part.pcd-esf.pcd $cad_dir/$j-cloud.pcd-esf.pcd
			done
		else
			for j in $cad_dir/*-cloud.pcd
			do
				$ESF $j
			done
		fi
		$script_dir/writeLabels.py $cad_dir
	fi
else
	if $computeDescriptors
	then
		for i in `seq -w $seq_in`
		do
			rm $input_dir/clusters$i/*-part.pcd
			rm $input_dir/clusters$i/*-esf.pcd
			if $combineAssembly
			then
				for j in $input_dir/clusters$i/*-cloud.pcd
				do
					$ESFC $j
				done
				for j in $input_dir/clusters$i/*-part.pcd
				do
					$ESF $j
				done
				numLines=`cat $input_dir/clusters$i/labels.txt | wc -l`
				for j in `seq 0 $((numLines-1))`
				do
					eval ~/Documents/PointCloudApp/cloud/concat_esf.py $input_dir/clusters$i/$j-cloud.pcd-$assemblies-part.pcd-esf.pcd $input_dir/clusters$i/$j-cloud.pcd-esf.pcd
				done
			else
				for j in $input_dir/clusters$i/*-cloud.pcd
				do
					$ESF $j
				done
			fi
			$script_dir/writeLabels.py $input_dir/clusters$i
		done
	fi
fi

#compute descriptors
if $computeDescriptors
then
	for i in `seq -w $seq_out`
	do
		rm $output_dir/clusters$i/*-part.pcd
		rm $output_dir/clusters$i/*-esf.pcd
		if $combineAssembly
		then
			for j in $output_dir/clusters$i/*-cloud.pcd
			do
				$ESFC $j
			done
			for j in $output_dir/clusters$i/*-part.pcd
			do
				$ESF $j
			done
			numLines=`cat $output_dir/clusters$i/labels.txt | wc -l`
			for j in `seq 0 $((numLines-1))`
			do
				eval ~/Documents/PointCloudApp/cloud/concat_esf.py $output_dir/clusters$i/$j-cloud.pcd-$assemblies-part.pcd-esf.pcd $output_dir/clusters$i/$j-cloud.pcd-esf.pcd
			done
		else
			for j in $output_dir/clusters$i/*-cloud.pcd
			do
				$ESF $j
			done
		fi
		$script_dir/writeLabels.py $output_dir/clusters$i
	done
fi

#divide into training and testing data
rm $output_dir/svm_train_data.txt $output_dir/svm_test_data.txt
for f in `seq -w $seq_out` 
do
	cat $output_dir/clusters$f/svmdata.txt >> $output_dir/svm_test_data.txt
done
if $useCAD
then
	cp $cad_dir/svmdata.txt $output_dir/svm_train_data.txt
else
	for f in `seq -w $seq_in`
	do
		cat $input_dir/clusters$f/svmdata.txt >> $output_dir/svm_train_data.txt
	done
fi

#filter data
if $ignoreZeroOption
then
	cat $output_dir/svm_train_data.txt | grep -v ^0 > $output_dir/svm_train_filtered.txt
	mv $output_dir/svm_train_filtered.txt $output_dir/svm_train_data.txt
	cat $output_dir/svm_test_data.txt | grep -v ^0 > $output_dir/svm_test_filtered.txt
	mv $output_dir/svm_test_filtered.txt $output_dir/svm_test_data.txt
fi

#check data
if ! $svm_dir/tools/checkdata.py $output_dir/svm_train_data.txt || ! $svm_dir/tools/checkdata.py $output_dir/svm_test_data.txt
then
	exit 1
fi
echo "Initialize SVM: train=$(wc -l < $output_dir/svm_train_data.txt) test=$(wc -l < $output_dir/svm_test_data.txt) ..."

#scale data
if $scaleOption
then
	if [ $ESF == $VFH ] || [ $ESF == $GRSD ]
	then
		$svm_dir/svm-scale -l 0 -u 1 -s $output_dir/range.txt $output_dir/svm_train_data.txt > $output_dir/svm_train_scaled.txt 2>/dev/null
		$svm_dir/svm-scale -l 0 -u 1 -s $output_dir/range.txt $output_dir/svm_test_data.txt > $output_dir/svm_test_scaled.txt 2>/dev/null
	else
		$svm_dir/svm-scale -l 0 -u 1 -s $output_dir/range.txt $output_dir/svm_train_data.txt > $output_dir/svm_train_scaled.txt 2>/dev/null
		$svm_dir/svm-scale -r $output_dir/range.txt $output_dir/svm_test_data.txt > $output_dir/svm_test_scaled.txt 2>/dev/null
	fi
else
	cp $output_dir/svm_train_data.txt $output_dir/svm_test_scaled.txt
	cp $output_dir/svm_test_data.txt $output_dir/svm_test_scaled.txt
fi

if $knnOption
then
	#use K Nearest Neighbors algorithm
	$script_dir/knn.py $k_parameter $output_dir/svm_train_scaled.txt $output_dir/svm_test_scaled.txt $output_dir/svm_prediction.txt
elif $ldaOption
then
	$script_dir/lda.py $output_dir
elif $lrOption
then
	$script_dir/logistic.py $output_dir
elif $bayesOption
then
	$script_dir/bayes.py $output_dir
elif $svcOption
then
	$script_dir/svc.py $output_dir
elif $adaboostOption
then
	$script_dir/adaboost.py $output_dir
elif $decisiontreeOption
then
	$script_dir/decision_tree.py $output_dir
else
	#train classifier
	if $parameterOption
	then
		parameters=`$svm_dir/tools/grid.py -gnuplot null -out null -t $kernel $output_dir/svm_train_scaled.txt | tail -1`
		c_param=`echo $parameters | cut -d' ' -f1`
		g_param=`echo $parameters | cut -d' ' -f2`
		echo "SVM Train: c=$c_param g=$g_param ..."
		$svm_dir/svm-train -c $c_param -g $g_param -s $svm_type -t $kernel -b 1 $output_dir/svm_train_scaled.txt $output_dir/model.txt > /dev/null
	else
		$svm_dir/svm-train -s $svm_type -t $kernel -b 1 $output_dir/svm_train_scaled.txt $output_dir/model.txt > /dev/null
	fi

	#calculate output
	echo "SVM Test: (Saving to $output_dir/svm_prediction.txt) ..."
	$svm_dir/svm-predict -b 1 $output_dir/svm_test_scaled.txt $output_dir/model.txt $output_dir/svm_prediction.txt
fi

#save predicted values to individual directories
i=0
j=1 #first line is label
while read line
do
	predictions[i]="$line"
	((i++))
done < $output_dir/svm_prediction.txt
for f in `seq -w $seq_out`
do
	rm $output_dir/clusters$f/prediction.txt
	numLines=`wc -l < $output_dir/clusters$f/labels.txt`
	read -a labels <<< `cat $output_dir/clusters$f/labels.txt`
	#echo "Updating $output_dir/clusters$f/prediction.txt ..."
	for l in ${labels[@]}
	do
		if $ignoreZeroOption && [ "$l" -eq 0 ]
		then
			echo 0 >> $output_dir/clusters$f/prediction.txt
		else
			echo ${predictions[j]} >> $output_dir/clusters$f/prediction.txt
			((j++))
		fi
	done
	if $savePointCloud
	then
		$script_dir/../main $output_dir/clusters$f/ $output_dir/clusters$f/combined.pcd
	fi
done

#count score
if $ignoreZeroOption
then
	$script_dir/countScore.py $output_dir $seq_out -i -c
else
	$script_dir/countScore.py $output_dir $seq_out -c
fi
