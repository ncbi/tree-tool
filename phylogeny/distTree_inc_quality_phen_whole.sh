#!/bin/bash
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


if [ $2 < 1 ]; then
  exit 1
fi
if [ $3 < 1 ]; then  
  exit 1
fi


TMP=`mktemp`
#echo $TMP  


cp $1/hist/tree.1 $TMP.tree

$THIS/tree2obj.sh $1/hist/tree.1 > $TMP.init


i=1
while [ $i <= $2 ]
do
	$THIS/distTree_inc_quality_phen_step.sh $1 $TMP.init $i $TMP $4
	i = $(( $i + $3 ))
done


rm -fr $TMP*
