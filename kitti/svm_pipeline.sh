#!/bin/bash

start=0000000000
end=0000000010
svm_dir=/home/jd/Downloads/libsvm-3.20

scaleOption=true
parameterOption=true
ignoreZeroOption=true
kernel=2

#combine data
rm svm_aggregate.txt
for i in `seq -w $start $end`
do
	cat clusters$i/svmdata.txt >> svm_aggregate.txt
done

#check data
if ! $svm_dir/tools/checkdata.py svm_aggregate.txt
then
	exit 1
fi

#filter data
if $ignoreZeroOption
then
	cat svm_aggregate.txt | grep -v ^0 > svm_aggregate_filtered.txt
	mv svm_aggregate_filtered.txt svm_aggregate.txt
fi

#divide into training and testing data (20/80)
totalLines=`wc -l svm_aggregate.txt | cut -d' ' -f1`
subset_size=$((totalLines / 5))
$svm_dir/tools/subset.py svm_aggregate.txt $subset_size svm_train_data.txt svm_test_data.txt
echo "Initialize SVM: train=$subset_size test=$((totalLines-subset_size)) ..."

 
#scale data
if $scaleOption
then
	$svm_dir/svm-scale -l 0 -u 1 -s range.txt svm_train_data.txt > svm_train_scaled.txt
	$svm_dir/svm-scale -r range.txt svm_test_data.txt > svm_test_scaled.txt
else
	cp svm_train_data.txt svm_test_scaled.txt
	cp svm_test_data.txt svm_test_scaled.txt
fi

#train classifier
if $parameterOption
then
	parameters=`$svm_dir/tools/grid.py -gnuplot null -out null -t $kernel svm_train_scaled.txt | tail -1`
	c_param=`echo $parameters | cut -d' ' -f1`
	g_param=`echo $parameters | cut -d' ' -f2`
	echo "SVM Train: c=$c_param g=$g_param ..."
	$svm_dir/svm-train -c $c_param -g $g_param -t $kernel svm_train_scaled.txt model.txt
else
	$svm_dir/svm-train -t $kernel svm_train_scaled.txt model.txt
fi

#calculate output
echo "SVM Test: (Saving to svm_prediction.txt) ..."
$svm_dir/svm-predict svm_test_scaled.txt model.txt svm_prediction.txt

#save predicted values to individual directories
i=0
j=0
read -a arr <<< `seq $start $end`
rm clusters${arr[i]}/prediction.txt
numLines=`wc -l < clusters${arr[i]}/labels.txt`
read -a labels <<< `cat clusters${arr[i]}/labels.txt`
while read line
do
	if $ignoreZeroOption && [ "${labels[j]}" -eq 0 ]
	then
		echo 0 >> clusters${arr[i]}/prediction.txt
	else
		echo $line >> clusters${arr[i]}/prediction.txt
	fi
	j=$((j+1))
	if [ "$j" -eq "$numLines" ]
	then
		i=$((i+1))
		j=0
		rm clusters${arr[i]}/prediction.txt
		numLines=`wc -l < clusters${arr[i]}/labels.txt`
		read -a labels <<< `cat clusters${arr[i]}/labels.txt`
	fi
done < svm_prediction.txt

