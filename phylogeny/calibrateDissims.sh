#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 6 ]; then
  echo "#1: input .dm-file without '.dm'"
  echo "#2: dissimilarity attribute in #1"
  echo "#3: parameters of makeDistTree after -variance"
  echo "#4: max. optimization iterations (> 0, normally 20)"
  echo "#5: phen/"
  echo "#6: phen is large (0/1)"
  exit 1
fi
INPUT=$1
DISSIM="$2"
VARIANCE_PAR="$3"
ITER_MAX=$4
PHEN=$5
PHEN_LARGE=$6


TMP=`mktemp`
echo $TMP 


echo ""
echo ""
$THIS/makeDistTree  -threads 5  -data $INPUT  -dissim_attr $DISSIM  -variance $VARIANCE_PAR \
     -optimize  -subgraph_iter_max $ITER_MAX  -noqual  -output_feature_tree $TMP.feature_tree

echo ""
echo ""
LARGE=""
if [ $PHEN_LARGE -eq 1 ]; then
  LARGE="-large"
fi
$THIS/makeFeatureTree  -input_tree $TMP.feature_tree  -features $PHEN  $LARGE  -prefer_gain  -nominal_singleton_is_optional  -qual $TMP.qual


rm $TMP*  
