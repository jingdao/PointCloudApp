#!/bin/bash
rm *.pcd
../tools/cluster ../cloud/kitti.pcd 10000 0.14 0.5 50 100000 200 0.14 -80 80 -80 30 -10 3 0.14
