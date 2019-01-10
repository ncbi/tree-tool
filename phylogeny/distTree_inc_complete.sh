#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Build a distance tree with complete pair-wise dissimilarity matrix using the incremental tree data structure"
  echo "#1: Incremental distance tree directory"
  echo "#2: List of objects"
  exit 1
fi


if [ -s $1/tree ]; then
  echo "$1/tree must be empty"
  exit 1
fi

N=`ls $1/new/ | wc -l`
if [ $N -gt 0 ]; then
  echo "$1/new/ must be empty"
  exit 1
fi

N=`ls $1/hybrid/ | wc -l`
if [ $N -gt 0 ]; then
  echo "$1/hybrid/ must be empty"
  exit 1
fi


sort.sh $2

if [ ! -s $1/dissim ]; then
	$THIS/../list2pairs $2 > $1/dissim_request
	$THIS/distTree_inc_request2dissim.sh $1 $1/dissim_request $1/dissim
	rm $1/dissim_request
fi


echo ""
echo "data.dm ..."
$THIS/../dm/pairs2attr2 $1/dissim 1 cons 6 -distance > $1/data.dm

echo ""
echo "Tree ..."
HYBRID=""
HYBRIDNESS_MIN=`cat $1/hybridness_min`
if [ "$HYBRIDNESS_MIN" != 0 ]; then
  DISSIM_BOUNDARY=`cat $1/dissim_boundary`
	HYBRID="-hybrid_parent_pairs $1/hybrid_parent_pairs  -delete_hybrids $1/hybrid.new  -delete_all_hybrids  -hybridness_min $HYBRIDNESS_MIN  -dissim_boundary $DISSIM_BOUNDARY"
fi
$THIS/makeDistTree  -threads 5  -data $1/data  -dissim cons  -optimize  $HYBRID \
  -output_tree $1/tree  -output_feature_tree $1/_feature_tree \
  > $1/hist/makeDistTree-complete.1
rm $1/data.dm

echo ""
echo "Database ..."
$1/objects_in_tree.sh $2 1

$THIS/distTree_inc_hybrid.sh $1 
$THIS/distTree_inc_unhybrid.sh $1


if [ -e $1/phen ]; then
  echo ""
  echo "Quality ..."
	$THIS/makeFeatureTree  -input_tree $1/_feature_tree  -features $1/phen  -nominal_singleton_is_optional  -prefer_gain  -output_core $1/_core  -qual $1/_qual > $1/hist/makeFeatureTree.1
	rm $1/_core
	rm $1/_qual
	grep ' !' $1/hist/makeFeatureTree.1
fi
rm $1/_feature_tree


