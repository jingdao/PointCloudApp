#!/bin/bash
rm *.pcd
../hpcd ../cloud/kitti.pcd ./ -b -80 -80 -10 80 80 3 -n 70
