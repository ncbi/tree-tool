#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Compute a complete pair-wise dissimilarity matrix and build a distance tree using the incremental tree data structure"
  echo "#1: incremental distance tree directory"
  echo "#2: sorted list of objects"
  echo "Output: #1/, #1/../data.dm"
  exit 1
fi
INC=$1
OBJS=$2


#set -x


#if false; then  
section "QC $INC/"
if [ -s $INC/tree ]; then
  error "$INC/tree must be empty"
fi

N=`distTree_inc_new_list.sh $INC | wc -l`
if [ $N -gt 0 ]; then
  error "$INC/new/ must be empty"
fi

sort -c $OBJS


section "Computing dissimilarities"
$THIS/../list2pairs $OBJS > $INC/dissim_request
$THIS/distTree_inc_request2dissim.sh $INC $INC/dissim_request $INC/dissim.raw
rm $INC/dissim_request
cat $INC/dissim.raw | grep -vwi "nan$" | grep -vwi "inf$" > $INC/dissim
rm $INC/dissim.raw

section "data.dm"
$THIS/../dm/pairs2dm $INC/dissim 1 "dissim" 6 -distance > $INC/../data.dm
echo "nan:"
set +o errexit
grep -wic 'nan$' $INC/../data.dm
set -o errexit

$THIS/../dm/dm2objs $INC/../data -noprogress | sort > $INC/tree.list

SERVER=`cat $INC/server`

if [ $SERVER ]; then
  section "Outliers"
  $THIS/../setMinus $OBJS $INC/tree.list > $INC/outlier-alien
  wc -l $INC/outlier-alien
  $THIS/../trav $INC/outlier-alien "$INC/outlier2db.sh %f alien"
  rm $INC/outlier-alien
fi


HYBRIDNESS_MIN=`cat $INC/hybridness_min`
if [ $HYBRIDNESS_MIN != 0 ]; then
  section "distTriangle"
 #cat data.dm | sed 's/nan/inf/g' > $INC/data1.dm
  mkdir $INC/clust
  $THIS/../dm/distTriangle $INC/../data "dissim"  -clustering_dir $INC/clust  -hybridness_min $HYBRIDNESS_MIN  -hybrid $INC/hybrid.new
  N=`ls $INC/clust/ | wc -l`
  if [ $N -gt 1 ]; then
    error "# Clusters: $N"
  fi
  mv $INC/clust/1/data.dm $INC/../data.dm
  rm -r $INC/clust/
  if [ $SERVER ]; then
    section "Hybrid"
  	$THIS/distTree_inc_hybrid.sh $INC 
  fi
fi


section "Tree"
HYBRID=""
if [ $HYBRIDNESS_MIN != 0 ]; then
  DISSIM_BOUNDARY=`cat $INC/dissim_boundary`
	HYBRID="-hybrid_parent_pairs $INC/hybrid_parent_pairs  -delete_hybrids $INC/hybrid.new  -hybridness_min $HYBRIDNESS_MIN  -dissim_boundary $DISSIM_BOUNDARY"
fi
VARIANCE=`cat $INC/variance`
$THIS/makeDistTree  -threads 5  -data $INC/../data  -dissim_attr "dissim"  -variance $VARIANCE  -optimize  -subgraph_iter_max 10  $HYBRID  -output_tree $INC/tree  > $INC/hist/makeDistTree-complete.1

if [ $SERVER ]; then
  section "Database"
  $INC/objects_in_tree.sh $INC/tree.list 1
  if [ $HYBRIDNESS_MIN != 0 ]; then
    section "Hybrid"
  	$THIS/distTree_inc_hybrid.sh $INC
  fi
fi

rm $INC/tree.list

cp $INC/tree $INC/hist/tree.1


if [ -e $INC/phen ]; then
  section "Quality"
  LARGE=0
  if [ -e $INC/large ]; then
    LARGE=1
  fi
  $THIS/tree_quality_phen.sh $INC/tree "" $INC/phen $LARGE 0 "" > $INC/hist/tree_quality_phen.1
	grep 'V !$' $INC/hist/tree_quality_phen.1
fi


echo ""
date


