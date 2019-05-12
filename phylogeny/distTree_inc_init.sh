#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 7 ]; then
  echo "Initialize the main directory of an incremental distance tree"
  echo "#1: Output directory"
  echo "#2: grid_min (> 0)"
  echo "#3: variance: lin|linExp"
  echo "#4: hybridness_min (> 1); 0 - no hybrids"
  echo "#5: dissim_boundary (> 0 or NAN)"
  echo "#6: genogroup_barrier (> 0 or NAN)"
  echo "#7: complete path to phen/ or ''"
  exit 1
fi
INC=$1
GRID_MIN=$2
VARIANCE="$3"
HYBRIDNESS_MIN=$4
DISSIM_BOUNDARY=$5
GENOGROUP_BARRIER=$6
PHEN=$7


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

echo "$VARIANCE"        > $INC/variance
echo $DISSIM_BOUNDARY   > $INC/dissim_boundary
echo $GENOGROUP_BARRIER > $INC/genogroup_barrier
echo $HYBRIDNESS_MIN    > $INC/hybridness_min


function create_script
{
	NAME=$1.sh
	#
	cp $THIS/distTree_inc/$NAME $INC/
	echo "Implement $INC/$NAME !"
}
create_script objects_in_tree
create_script outlier2db
create_script request2dissim
create_script request_closest
create_script qc
if [ $HYBRIDNESS_MIN != 0 ]; then
	create_script db2unhybrid
	create_script hybrid2db
fi
if [ "$DISSIM_BOUNDARY" != "NAN" ]; then
  create_script genogroup2db
fi


if [ $PHEN ]; then
  ln -s $PHEN $INC/phen
fi