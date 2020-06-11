#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 8 ]; then
  echo "Initialize an incremental distance tree directory"
  echo "#1: output directory"
  echo "#2: grid_min (> 0): min. number of dissimilarity requests to use GRID"
  echo "#3: variance: lin|linExp|..."
  echo "#4: hybridness_min (> 1); 0 - no hybrids"
  echo "#5: dissim_boundary (> 0 or NAN)"
  echo "#6: genogroup_barrier (> 0 or NAN)"
  echo "#7: complete path to 'phenotype' directory or ''"
  echo "#8: phen_large (0/1): 1 - files in #7 directory are grouped into subdirectories named file2hash(<file name>)"
  exit 1
fi
INC=$1
GRID_MIN=$2
VARIANCE="$3"
HYBRIDNESS_MIN=$4
DISSIM_BOUNDARY=$5
GENOGROUP_BARRIER=$6
PHEN=$7
PHEN_LARGE=$8


if [ $GRID_MIN -le 0 ]; then
  echo "Bad GRID_MIN"
  exit 1
fi


mkdir $INC

touch $INC/dissim
touch $INC/leaf
touch $INC/tree

mkdir $INC/new
mkdir $INC/search

echo "1" > $INC/version

mkdir $INC/hist

echo $GRID_MIN > $INC/grid_min
touch $INC/runlog

echo "$VARIANCE"        > $INC/variance
echo $DISSIM_BOUNDARY   > $INC/dissim_boundary
echo $GENOGROUP_BARRIER > $INC/genogroup_barrier
echo $HYBRIDNESS_MIN    > $INC/hybridness_min


function create_script
{
	NAME=$1.sh
	#
	cp $THIS/distTree_inc/$NAME $INC/
	echo -e "${YELLOW}Implement $INC/$NAME !${NOCOLOR}"
}
create_script objects_in_tree
create_script outlier2db
create_script request2dissim
create_script request_closest
create_script qc
if false; then  # deprecated
  if [ $HYBRIDNESS_MIN != 0 ]; then
  	create_script db2unhybrid
  	create_script hybrid2db
  fi
fi
if [ "$DISSIM_BOUNDARY" != "NAN" ]; then
  create_script genogroup2db
fi


if [ $PHEN ]; then
  ln -s $PHEN $INC/phen
  echo $PHEN_LARGE > $INC/phen_large
fi
