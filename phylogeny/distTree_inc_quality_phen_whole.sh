#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 4 ]; then
  echo "Phenotypic quality of an incremental tree: initial tree"
  echo "Requires: no added outliers"
  echo "#1: incremental tree data structure"
  echo "#2: target tree version"
  echo "#3: version increment (>= 1)"
  echo "#4: Directory for output trees: #2.tree, #2.makeDistTree, #2.makeFeatureTree, #2.qual"
  exit 1
fi
INC=$1
VER_MAX=$2
VER_STEP=$3
OUT_DIR=$4


if [ $VER_MAX < 1 ]; then
  exit 1
fi
if [ $VER_STEP < 1 ]; then  
  exit 1
fi


TMP=`mktemp`
#echo $TMP  


cp $INC/hist/tree.1 $TMP.tree

$THIS/tree2obj.sh $INC/hist/tree.1 > $TMP.init


i=1
while [ $i -le $VER_MAX ]
do
	$THIS/distTree_inc_quality_phen_step.sh $INC $TMP.init $i $TMP $OUT_DIR
	i = $(( $i + $VER_STEP ))
done


rm -fr $TMP*
