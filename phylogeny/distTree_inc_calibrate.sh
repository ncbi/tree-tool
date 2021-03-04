#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "#1: incremental distance tree directory"
  echo "#2: dissimilarity power (> 0)"
  echo "#3: dissim_coeff (> 0)"
  exit 1
fi
INC=$1
DISSIM_POWER=$2
DISSIM_COEFF=$3


$THIS/../check_tmp.sh
TMP=`mktemp`
echo $TMP 


section "Building tree ..."
VARIANCE=`cat $INC/variance`
$THIS/makeDistTree  -threads 15  -data $INC/  -variance $VARIANCE  -dissim_coeff $DISSIM_COEFF  -dissim_power $DISSIM_POWER \
  -optimize  -skip_len  -reinsert  -subgraph_iter_max 5  -noqual  -output_feature_tree $TMP.feature_tree 

section "Evaluating tree ..."
LARGE=""
if [ -e $INC/large ]; then
  LARGE="-large"
fi
$THIS/makeFeatureTree  -threads 15  -input_tree $TMP.feature_tree  -features $INC/phen  $LARGE \
  -prefer_gain  -nominal_singleton_is_optional  -qual $TMP.qual   -output_core $TMP.core  -save_mem


rm $TMP*  
