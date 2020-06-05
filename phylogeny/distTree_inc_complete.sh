#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Compute a complete pair-wise dissimilarity matrix and build a distance tree using the incremental tree data structure"
  echo "#1: incremental distance tree directory"
  echo "#2: list of objects"
  echo "#3: update database (0/1)"
  echo "Output: #1/, #1/../data.dm"
  exit 1
fi
INC=$1
OBJS=$2
DB=$3


#if false; then  
if [ -s $INC/tree ]; then
  echo "$INC/tree must be empty"
  exit 1
fi

N=`ls $INC/new/ | wc -l`
if [ $N -gt 0 ]; then
  echo "$INC/new/ must be empty"
  exit 1
fi


$THIS/../sort.sh $OBJS

$THIS/../list2pairs $OBJS > $INC/dissim_request
$THIS/distTree_inc_request2dissim.sh $INC $INC/dissim_request $INC/dissim.raw
rm $INC/dissim_request
cat $INC/dissim.raw | grep -vwi "nan" | grep -vwi "inf" > $INC/dissim
rm $INC/dissim.raw

echo ""
echo "data.dm ..."
$THIS/../dm/pairs2attr2 $INC/dissim 1 "cons" 6 -distance > $INC/../data.dm
echo "nan:"
set +o errexit
grep -wic 'nan' $INC/../data.dm
set -o errexit

echo ""
$THIS/../dm/dm2objs $INC/../data -noprogress | sort > $INC/tree.list
if [ $DB == 1 ]; then
  $THIS/../setMinus $OBJS $INC/tree.list > $INC/outlier-alien
  wc -l $INC/outlier-alien
  $THIS/../trav $INC/outlier-alien "$INC/outlier2db.sh %f alien"
  rm $INC/outlier-alien
fi


HYBRIDNESS_MIN=`cat $INC/hybridness_min`
if [ $HYBRIDNESS_MIN != 0 ]; then
  echo ""
  echo "distTriangle ..."
 #cat data.dm | sed 's/nan/inf/g' > $INC/data1.dm
  mkdir $INC/clust
  $THIS/../dm/distTriangle $INC/../data "cons"  -clustering_dir $INC/clust  -hybridness_min $HYBRIDNESS_MIN  -hybrid $INC/hybrid.new
  N=`ls $INC/clust/ | wc -l`
  if [ $N -gt 1 ]; then
    echo "# Clusters: $N"
    exit 1
  fi
  mv $INC/clust/1/data.dm $INC/../data.dm
  rm -r $INC/clust/
  $THIS/../dm/dm2objs $INC/../data | sort > $INC/tree.list  
  if [ $DB == 1 ]; then
    echo ""
    echo "Hybrid ..."
  	$THIS/distTree_inc_hybrid.sh $INC 
   #echo "Unhybrid ..."
   #$THIS/distTree_inc_unhybrid.sh $INC 
  fi
fi


echo ""
echo "Tree ..."
HYBRID=""
if [ $HYBRIDNESS_MIN != 0 ]; then
  DISSIM_BOUNDARY=`cat $INC/dissim_boundary`
	HYBRID="-hybrid_parent_pairs $INC/hybrid_parent_pairs  -delete_hybrids $INC/hybrid.new  -hybridness_min $HYBRIDNESS_MIN  -dissim_boundary $DISSIM_BOUNDARY"
fi
VARIANCE=`cat $INC/variance`
$THIS/makeDistTree  -threads 5  -data $INC/../data  -dissim_attr "cons"  -variance $VARIANCE  -optimize  -subgraph_iter_max 10  $HYBRID  -output_tree $INC/tree  > $INC/hist/makeDistTree-complete.1

if [ $DB == 1 ]; then
  echo ""
  echo "Database ..."
  $INC/objects_in_tree.sh $INC/tree.list 1
  if [ $HYBRIDNESS_MIN != 0 ]; then
    echo ""
    echo "Hybrid ..."
  	$THIS/distTree_inc_hybrid.sh $INC
   #echo "Unhybrid ..."
   #$THIS/distTree_inc_unhybrid.sh $INC
  fi
fi

rm $INC/tree.list

cp $INC/tree $INC/hist/tree.1


if [ -e $INC/phen ]; then
  echo ""
  echo "Quality ..."
	PHEN_LARGE=`cat $INC/phen_large`
  $THIS/tree_quality_phen.sh $INC/tree "" $INC/phen $PHEN_LARGE 0 "" > $INC/hist/tree_quality_phen.1
	grep ' !' $INC/hist/tree_quality_phen.1
fi


echo ""
date


