#!/bin/bash

for l in truth/*
do
	for c in 0.6 0.7 0.8
	do
		echo $l $c
		rm detection/*-cloud.pcd
		~/Documents/PointCloudApp/lbp3d combined.pcd $l/0-cloud.pcd detection -c $c > /dev/null
		python get_accuracy.py detection/ $l
	done
done

for l in truth/*
do
	for c in 0.6 0.7 0.8
	do
		echo $l $c
		rm detection/*-cloud.pcd
		~/Documents/PointCloudApp/lbp3d original.pcd $l/0-cloud.pcd detection -c $c > /dev/null
		python get_accuracy.py detection/ $l
	done
done
