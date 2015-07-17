#!/bin/bash

inputFile=$1
outputDir=$2
baseDir=`pwd`

#Set parameters
leafX=200
leafY=200
leafZ=200
segThreshold=200
clusterThreshold=500
minCluster=100
maxCluster=100000

echo "leaf: $leafX $leafY $leafZ" > $outputDir/param.txt
echo "segmentation: $segThreshold" >> $outputDir/param.txt
echo "clustering: $minCluster-$maxCluster $clusterThreshold" >> $outputDir/param.txt

#Downsample
pcl_voxel_grid $inputFile $outputDir/filtered.pcd -leaf $leafX,$leafY,$leafZ

#Segmentation and Clustering
rm $outputDir/*-cloud.pcd
cd $outputDir
$baseDir/../tools/cluster $baseDir/$outputDir/filtered.pcd 100000 $segThreshold 0.5 $minCluster $maxCluster 200 $clusterThreshold 
#../main $outputDir/$inputFile-filtered.pcd $outputDir
cd $baseDir

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
