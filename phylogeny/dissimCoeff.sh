#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 5 ]; then
  echo "Input: #1.dm with dissimilarity 'cons'"
  echo "#1: dissimilarity data"
  echo "#2: phen/"
  echo "#3: dissim_power"
  echo "#4: dissim_coeff"
  echo "#5: variance"
  exit 1
fi
DATA=$1
PHEN=$2
DISSIM_POWER=$3
DISSIM_COEFF=$4
VARIANCE=$5


echo ""
echo "dissim_power = $DISSIM_POWER"
echo "dissim_coeff = $DISSIM_COEFF"
echo "varaince     = $VARIANCE"


TMP=`mktemp`


$THIS/makeDistTree  -threads 5  -data $DATA  -dissim cons  -dissim_power $DISSIM_POWER  -dissim_coeff $DISSIM_COEFF  -variance $VARIANCE  \
  -optimize  -subgraph_iter_max 10  -noqual  -output_tree $TMP.tree  -output_feature_tree $TMP.feature_tree

echo ""
echo ""
$THIS/makeFeatureTree  -input_tree $TMP.feature_tree  -features $PHEN  -nominal_singleton_is_optional  -qual $TMP.qual


rm -f $TMP*


