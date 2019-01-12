#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Phenotypic quality of a distance tree, find root"
  echo "#1: distance tree"
  echo "#2: target list of objects | '' - all"
  echo "#3: phen/"
  exit 1
fi
TREE=$1
TARGET="$2"
PHEN=$3


TMP=`mktemp`
echo $TMP


DELETE=""
if [ $TARGET ]; then
  tree2obj.sh $TREE > $TMP.cur
  $THIS/../setMinus $TMP.cur $TARGET > $TMP.del
  DELETE="-delete $TMP.del  -check_delete"
fi

$THIS/makeDistTree  -threads 15  -input_tree $TREE  $DELETE  -output_feature_tree $TMP.feature_tree > $TMP.distTree

$THIS/makeFeatureTree  -threads 15  -input_tree $TMP.feature_tree  -features $PHEN  \
  -prefer_gain  -nominal_singleton_is_optional  -output_core $TMP.core  -qual $TMP.qual


rm -f $TMP*
