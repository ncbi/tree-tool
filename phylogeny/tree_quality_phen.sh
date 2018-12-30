#!/bin/bash
source bash_common.sh
if [ $# -ne 2 ]; then
  echo "Phenotypic quality of a distance tree"
  echo "#1: distance tree"
  echo "#2: phen/"
  exit 1
fi
INC=$1
PHEN=$2


TMP=`mktemp`
echo $TMP

makeDistTree  -threads 15  -input_tree $INC  -output_feature_tree $TMP.feature_tree  > $TMP.distTree
makeFeatureTree  -threads 15  -input_tree $TMP.feature_tree  \
  -features $PHEN  -nominal_singleton_is_optional  -output_core $TMP.core  -qual $TMP.qual 

rm -f $TMP*
