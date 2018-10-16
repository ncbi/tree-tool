#!/bin/bash
source bash_common.sh
if [ $# -ne 4 ]; then
  echo "Initialize an incremental distance tree"
  echo "#1: Output directory"
  echo "#2: grid_min (> 0)"
  echo "#3: hybridness_min (> 1); 0 - no hybrids"
  echo "#4: dissim_boundary (> 0)"
  exit 1
fi

INC=$1
GRID_MIN=$2
HYBRIDNESS_MIN=$3
DISSIM_BOUNDARY=$4


if [ $GRID_MIN -le 0 ]; then
  exit 1
fi


mkdir $INC

cp /dev/null $INC/dissim
cp /dev/null $INC/leaf
cp /dev/null $INC/tree

mkdir $INC/new
mkdir $INC/outlier
mkdir $INC/search

echo "1" > $INC/version

mkdir $INC/hist

echo $GRID_MIN > $INC/GRID_MIN
cp /dev/null $INC/runlog


if [ $HYBRIDNESS_MIN != 0 ]; then
	echo $DISSIM_BOUNDARY > $INC/DISSIM_BOUNDARY
	echo $HYBRIDNESS_MIN  > $INC/HYBRIDNESS_MIN
fi


function create_script
{
	name=$1.sh
	echo "exit 1" > $INC/$name
	chmod a+x $INC/$name
	echo "Set $INC/$name !"
}
create_script objects_in_tree
create_script request2dissim
create_script request_closest
if [ $HYBRIDNESS_MIN != 0 ]; then
 #create_script db2hybrid
	create_script db2unhybrid
	create_script hybrid2db
fi


