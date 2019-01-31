#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Compute a complete pair-wise dissimilarity matrix and build a distance tree using the incremental tree data structure"
  echo "#1: Incremental distance tree directory"
  echo "#2: List of objects"
  echo "Output: #1/, data.dm"
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


$THIS/../sort.sh $2

if [ ! -s $1/dissim ]; then
	$THIS/../list2pairs $2 > $1/dissim_request
	$THIS/distTree_inc_request2dissim.sh $1 $1/dissim_request $1/dissim.raw
	rm $1/dissim_request
	cat $1/dissim.raw | grep -vwi nan | grep -vwi inf > $1/dissim
	rm $1/dissim.raw
fi


echo ""
echo "data.dm ..."
$THIS/../dm/pairs2attr2 $1/dissim 1 cons 6 -distance > data.dm

echo ""
echo "Tree ..."
HYBRID=""
HYBRIDNESS_MIN=`cat $1/hybridness_min`
if [ "$HYBRIDNESS_MIN" != 0 ]; then
  DISSIM_BOUNDARY=`cat $1/dissim_boundary`
	HYBRID="-hybrid_parent_pairs $1/hybrid_parent_pairs  -delete_hybrids $1/hybrid.new  -delete_all_hybrids  -hybridness_min $HYBRIDNESS_MIN  -dissim_boundary $DISSIM_BOUNDARY"
fi
$THIS/makeDistTree  -threads 5  -data data  -dissim cons  -optimize  $HYBRID  -output_tree $1/tree  > $1/hist/makeDistTree-complete.1
#rm $1/data.dm

echo ""
echo "Database ..."
$1/objects_in_tree.sh $2 1

if [ "$HYBRIDNESS_MIN" != 0 ]; then
  echo ""
  echo "Hybrid ..."
	$THIS/distTree_inc_hybrid.sh $1
  echo "Unhybrid ..."
  $THIS/distTree_inc_unhybrid.sh $1
fi


if [ -e $1/phen ]; then
  echo ""
  echo "Quality ..."
  $THIS/tree_quality_phen.sh $1/tree "" $1/phen > $1/hist/tree_quality_phen.1
	grep ' !' $1/hist/tree_quality_phen.1
fi


echo ""
date


