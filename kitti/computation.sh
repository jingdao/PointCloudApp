#!/bin/bash

test_dir=/home/jd/Documents/PointCloudApp/cloud/review/test
base_dir=/home/jd/Documents/PointCloudApp
DESC_LIST="spin fpfh shot usc esf vfh ../pad3d"
TIMEFORMAT='%U'

for DESC_EXE in $DESC_LIST
do
	DESC=`basename $DESC_EXE`
	rm $test_dir/tmp-$DESC.txt
done

for i in $test_dir/*-cloud.pcd
do
	l=`cat $i | wc -l`
	l=$((l-11))
	for DESC_EXE in $DESC_LIST
	do
		DESC=`basename $DESC_EXE`
		t=$((time $base_dir/tools/$DESC_EXE $i) 2>&1 1>/dev/null)
		printf "%7d %8.3f\n" "$l" "$t" >> $test_dir/tmp-$DESC.txt
	done
done

for DESC_EXE in $DESC_LIST
do
	DESC=`basename $DESC_EXE`
	cat $test_dir/tmp-$DESC.txt | sort -n > $test_dir/timing-$DESC.txt
done
