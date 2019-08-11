#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 5 ]; then
  echo "#1: input .dm-file without '.dm'"
  echo "#2: dissimilarity attribute in #1"
  echo "#3: dissimilarity parameters of makeDistTree"
  echo "#4: variance parameters of makeDistTree (without '-variance ')"
  echo "#5: phen/"
  exit 1
fi
INPUT=$1
DISSIM=$2
DISSIM_PAR="$3"
VARIANCE_PAR="$4"
PHEN=$5


TMP=`mktemp`
#echo $TMP 


$THIS/../dm/dm2objs $INPUT | sort > $TMP.list
ls $PHEN > $TMP.phen
$THIS/../setMinus $TMP.list $TMP.phen > $TMP.delete
if [ -s $TMP.delete ]; then
  wc -l $TMP.delete
  $THIS/../dm/dm2subset $INPUT $TMP.delete -exclude > $TMP.dm
  INPUT=$TMP.dm
fi

echo ""
echo ""
$THIS/makeDistTree  -threads 5  -data $INPUT  -dissim_attr $DISSIM \
  $DISSIM_PAR  -variance $VARIANCE_PAR \
  -optimize  -subgraph_iter_max 20  -noqual  -output_feature_tree $TMP.feature_tree

echo ""
echo ""
$THIS/makeFeatureTree  -input_tree $TMP.feature_tree  -features $PHEN  -prefer_gain  -nominal_singleton_is_optional  -qual $TMP.qual


rm $TMP*  
