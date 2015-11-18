#!/bin/bash
rm *.pcd
../tools/cluster ../cloud/kitti.pcd 10000 0.186 0.5 10 100000 200 0.186 -80 80 -80 80 -10 3 0.186
