#!/bin/bash
set -o nounset
set -o errexit
set -o posix
set -o pipefail
export LC_ALL=C

if [ $# != 4 ]; then
  echo "Initialize an incremental distance tree"
  echo "#1: Output directory"
  echo "#2: grid_min (> 0)"
  echo "#3: hybridness_min (> 1)"
  echo "#4: dissim_boundary (> 0)"
  exit 1
fi

inc=$1
grid_min=$2
hybridness_min=$3
dissim_boundary=$4

if [ $grid_min -le 0 ]; then
  exit 1
fi


mkdir $inc

cp /dev/null $inc/tree
cp /dev/null $inc/dissim
cp /dev/null $inc/leaf

mkdir $inc/new
mkdir $inc/search
mkdir $inc/outlier

echo "1" > $inc/version

mkdir $inc/hist

echo $hybridness_min  > $inc/hybridness_min
echo $grid_min        > $inc/grid_min
echo $dissim_boundary > $inc/dissim_boundary

cp /dev/null $inc/runlog

echo "exit 1" > $inc/request2dissim.sh
chmod a+x $inc/request2dissim.sh

echo "exit 1" > $inc/objects_in_tree.sh
chmod a+x $inc/objects_in_tree.sh

echo "exit 1" > $inc/request_closest.sh
chmod a+x $inc/request_closest.sh

echo "exit 1" > $inc/hybrid2db.sh
chmod a+x $inc/hybrid2db.sh

echo "exit 1" > $inc/db2hybrid.sh
chmod a+x $inc/db2hybrid.sh

echo "exit 1" > $inc/db2unhybrid.sh
chmod a+x $inc/db2unhybrid.sh
