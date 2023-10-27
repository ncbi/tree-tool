#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Compute a complete pair-wise dissimilarity matrix and build a distance tree using the incremental tree data structure"
  echo "#1: incremental distance tree directory"
  echo "#2: sorted and distinct list of objects"
  echo "output: #1/, #1/../data.dm"
  exit 1
fi
INC=$1
OBJS=$2


#set -x


VER=`cat $INC/version`
if [ $VER -ne 1 ]; then
  error "version must be 1"
fi

SERVER=`cat $INC/server`


#if false; then  
section "QC $INC/"
if [ -s $INC/tree ]; then
  error "$INC/tree must be empty"
fi

if [ -s $INC/dissim ]; then
  error "$INC/dissim must be empty"
fi

if [ -s $INC/dissim.bad ]; then
  error "$INC/dissim.bad must be empty"
fi

if [ -s $INC/indiscern ]; then
  error "$INC/indiscern must be empty"
fi

TMP=`mktemp`
$THIS/distTree_inc_new_list.sh $INC > $TMP || (cat $TMP && exit 1)
N=`cat $TMP | wc -l`
if [ $N -gt 0 ]; then
  error "$INC/new/ must be empty"
fi
rm $TMP

if [ -s $INC/runlog ]; then
  error "$INC/runlog must be empty"
fi

sort -cu $OBJS


section "Computing dissimilarities"
$THIS/../list2pairs $OBJS > $INC/dissim_request
$THIS/distTree_inc_request2dissim.sh $INC $INC/dissim_request $INC/dissim
rm $INC/dissim_request
$THIS/distTree_inc_dissim2indiscern.sh $INC $INC/dissim

sort -k3,3g $INC/dissim > $TMP.dissim.sort
LONGEST=`tail -1 $TMP.dissim.sort | cut -f 3`
echo "Longest dissimilarity: $LONGEST"
wc -l $INC/dissim
warning "# Pairs with longest dissimilarity:"
cut -f 3 $INC/dissim | grep -cx $LONGEST

section "data.dm"
$THIS/../dm/conversion/pairs2dm $INC/dissim 1 "dissim" 6 -distance > $INC/../data.dm
warning "nan:"
grep -wic 'nan' $INC/../data.dm || true

$THIS/../dm/dm2objs $INC/../data -noprogress | sort > $INC/tree.list

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
  	$THIS/../dm/dm2subset $INC/../data $INC/hist/hybrid-indiscern.$VER -exclude > $INC/data.dm
  	mv $INC/data.dm $INC/../
  fi
fi

super_section "Tree"
HYBRID=""
if [ $HYBRIDNESS_MIN != 0 ]; then
  DISSIM_BOUNDARY=`cat $INC/dissim_boundary`
	HYBRID="-hybrid_parent_pairs $INC/hybrid_parent_pairs  -delete_hybrids $INC/hybrid.new  -hybridness_min $HYBRIDNESS_MIN  -dissim_boundary $DISSIM_BOUNDARY"
fi
VARIANCE=`cat $INC/variance`
$THIS/makeDistTree  -threads 5  -data $INC/../data  -dissim_attr "dissim"  -variance $VARIANCE  -optimize  -subgraph_iter_max 10  $HYBRID  -output_tree $INC/tree  -output_tree_tmp $INC/tree.tmp > $INC/hist/makeDistTree-complete.1
rm $INC/tree.tmp

if [ $SERVER ]; then
  section "Database"
  $THIS/tree2obj.sh $INC/tree > $INC/tree.list
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


