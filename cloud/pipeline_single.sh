#!/bin/bash

inputFile=$1
outputDir=$2
baseDir=`pwd`

#Set parameters
leafX=0.05
leafY=0.05
leafZ=0.05
segThreshold=0.4
clusterThreshold=0.4
minCluster=100
maxCluster=200000
downsample=false
subsample=false

echo "leaf: $leafX $leafY $leafZ" > $outputDir/param.txt
echo "segmentation: $segThreshold" >> $outputDir/param.txt
echo "clustering: $minCluster-$maxCluster $clusterThreshold" >> $outputDir/param.txt
if $subsample
then
	minX=-20
	maxX=20
	minY=-30
	maxY=0
	minZ=238
	maxZ=248
	echo "x range: $minX $maxX" >> $outputDir/param.txt
	echo "y range: $minY $maxY" >> $outputDir/param.txt
	echo "z range: $minZ $maxZ" >> $outputDir/param.txt
fi

if $downsample
then
	#Downsample
	if $subsample
	then
		pcl_passthrough_filter $inputFile $outputDir/filteredX.pcd -field x -min $minX -max $maxX
		pcl_passthrough_filter $outputDir/filteredX.pcd $outputDir/filteredY.pcd -field y -min $minY -max $maxY
		pcl_passthrough_filter $outputDir/filteredY.pcd $outputDir/filteredZ.pcd -field z -min $minZ -max $maxZ
		pcl_voxel_grid $outputDir/filteredZ.pcd $outputDir/filtered.pcd -leaf $leafX,$leafY,$leafZ
	else
		pcl_voxel_grid $inputFile $outputDir/filtered.pcd -leaf $leafX,$leafY,$leafZ
	fi
fi

#Segmentation and Clustering
rm $outputDir/*-cloud.pcd
cd $outputDir
$baseDir/../tools/cluster $baseDir/$outputDir/filtered.pcd 100000 $segThreshold 0.5 $minCluster $maxCluster 200 $clusterThreshold 
#../main $outputDir/$inputFile-filtered.pcd $outputDir
cd $baseDir

#write to PLY
pcl_convert_pcd_ascii_binary $outputDir/filtered.pcd $outputDir/filtered.pcd 0
../main $outputDir/filtered.pcd $outputDir/filtered.ply

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
