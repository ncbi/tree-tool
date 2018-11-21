#!/bin/bash
source bash_common.sh
if [ $# -ne 5 ]; then
  echo "Initialize an incremental distance tree"
  echo "#1: Output directory"
  echo "#2: grid_min (> 0)"
  echo "#3: hybridness_min (> 1); 0 - no hybrids"
  echo "#4: dissim_boundary (> 0 or NAN)"
  echo "#5: genogroup_barrier (> 0 or NAN)"
  exit 1
fi
INC=$1
GRID_MIN=$2
HYBRIDNESS_MIN=$3
DISSIM_BOUNDARY=$4
GENOGROUP_BARRIER=$5


if [ $GRID_MIN -le 0 ]; then
  exit 1
fi


mkdir $INC

cp /dev/null $INC/dissim
cp /dev/null $INC/leaf
cp /dev/null $INC/tree

mkdir $INC/new
mkdir $INC/hybrid
mkdir $INC/search

echo "1" > $INC/version

mkdir $INC/hist

echo $GRID_MIN > $INC/grid_min
cp /dev/null $INC/runlog

echo $DISSIM_BOUNDARY > $INC/dissim_boundary
echo $GENOGROUP_BARRIER > $INC/genogroup_barrier

if [ $HYBRIDNESS_MIN != 0 ]; then
	echo $HYBRIDNESS_MIN > $INC/hybridness_min
fi


function create_script
{
	name=$1.sh
	echo "exit 1" > $INC/$name
	chmod a+x $INC/$name
	echo "Set $INC/$name !"
}
create_script objects_in_tree
create_script outlier2db
create_script request2dissim
create_script request_closest
if [ $HYBRIDNESS_MIN != 0 ]; then
	create_script db2unhybrid
	create_script hybrid2db
fi
if [ $DISSIM_BOUNDARY != "NAN" ]; then
  create_script genogroup2db
fi


