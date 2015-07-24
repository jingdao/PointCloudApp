#!/bin/bash

inputFile=$1
outputDir=$2
baseDir=`pwd`

rm $outputDir/*

#Set parameters
leafX=0.05
leafY=0.05
leafZ=0.05
segThreshold=0.1
clusterThreshold=0.5
minCluster=500
maxCluster=200000

echo "leaf: $leafX $leafY $leafZ" > $outputDir/param.txt
echo "segmentation: $segThreshold" >> $outputDir/param.txt
echo "clustering: $minCluster-$maxCluster $clusterThreshold" >> $outputDir/param.txt

#Downsample
pcl_voxel_grid $inputFile $outputDir/filtered.pcd -leaf $leafX,$leafY,$leafZ

#Segmentation and Clustering
cd $outputDir
$baseDir/../tools/cluster $baseDir/$outputDir/filtered.pcd 100000 $segThreshold 0.5 $minCluster $maxCluster 200 $clusterThreshold 
#../main $outputDir/$inputFile-filtered.pcd $outputDir
cd $baseDir

#write to PLY
pcl_convert_pcd_ascii_binary $outputDir/filtered.pcd $outputDir/filtered.pcd 0
../main $outputDir/filtered.pcd $outputDir/filtered.ply

#write OBJ
for i in $outputDir/*-cloud.pcd
do
	../main $i $i.obj
done

exit

#User input for labels
../tools/viz $outputDir > $outputDir/labels.txt

#Compute descriptors
for j in $outputDir/*-cloud.pcd
do
	../tools/esf $j $j-esf.pcd
done

#Write labels
../kitti/writeLabels.py $outputDir/

#Combine point cloud
cp $outputDir/labels.txt $outputDir/prediction.txt
../main $outputDir $outputDir/combined.pcd

#View result
pcl_viewer $outputDir/combined.pcd
