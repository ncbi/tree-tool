#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Input: #1.dm with dissimilarity 'cons'"
  echo "#1: dissimilarity data"
  echo "#2: phen/"
  echo "#3: dissim_coeff"
  exit 1
fi
DATA=$1
PHEN=$2
DISSIM_COEFF=$3


echo ""
echo "dissim_coeff = $DISSIM_COEFF"


TMP=`mktemp`


$THIS/makeDistTree  -threads 5  -data $DATA  -dissim cons  -dissim_coeff $DISSIM_COEFF  -optimize  -noqual  -output_tree $TMP.tree  -output_feature_tree $TMP.feature_tree

echo ""
echo ""
$THIS/makeFeatureTree  -input_tree $TMP.feature_tree  -features $PHEN  -nominal_singleton_is_optional  -output_core $TMP.core  -qual $TMP.qual


rm -f $TMP*


