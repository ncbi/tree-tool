#!/bin/bash
source bash_common.sh
if [ $# -ne 2 ]; then
  echo "Quality of the target list of objects in the current tree"
  echo "Print the result of makeFeatureTree"
  echo "Input: #1/phen/"
  echo "#1: incremental distance tree directory"
  echo "#2: target list of objects"
  exit 1
fi
INC=$1
TARGET=$2


TMP=`mktemp`

tree2obj.sh $INC/tree > $TMP.cur
setMinus $TMP.cur $TARGET > $TMP.del
makeDistTree  -threads 15  -input_tree $INC/tree  -delete $TMP.del  -check_delete  -output_feature_tree $TMP.feature_tree >& /dev/null
makeFeatureTree  -threads 15  -input_tree $TMP.feature_tree  -features $INC/phen  -prefer_gain  -output_core $TMP.core  -qual $TMP.qual

rm $TMP*

