#!/bin/bash
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


TMP=`mktemp`
echo $TMP 


VARIANCE=`cat $INC/variance`

$THIS/makeDistTree  -threads 15  -data $INC/  -variance $VARIANCE  -dissim_coeff $DISSIM_COEFF  -dissim_power $DISSIM_POWER \
  -optimize  -skip_len  -reinsert  -subgraph_iter_max 5  -noqual  -output_feature_tree $TMP.feature_tree 

echo ""
echo ""
$THIS/makeFeatureTree  -threads 15  -input_tree $TMP.feature_tree  -features $INC/phen  \
  -prefer_gain  -nominal_singleton_is_optional  -qual $TMP.qual   -output_core $TMP.core  -save_mem


rm $TMP*  