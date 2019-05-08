#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 6 ]; then
  echo "#1: input .dm-file without '.dm'"
  echo "#2: dissimilarity attribute in #1"
  echo "#3: dissimilarity power (> 0)"
  echo "#4: dissim_coeff: >0 <=> variance linExp"
  echo "#5: variance power (NAN <=> variance = linExp)"
  echo "#6: phen/"
  exit 1
fi
INPUT=$1
DISSIM=$2
DISSIM_POWER=$3
DISSIM_COEFF=$4
VAR_POWER=$5
PHEN=$6


TMP=`mktemp`
#echo $TMP 


VARIANCE=linExp
DISSIM_COEFF_OPTION="-dissim_coeff $DISSIM_COEFF"
if [ $DISSIM_COEFF == 0 ]; then
  VARIANCE="pow  -variance_power $VAR_POWER"
  DISSIM_COEFF_OPTION=""
fi

$THIS/../dm/dm2objs $INPUT | sort > $TMP.list
ls $PHEN > $TMP.phen
$THIS/../setMinus $TMP.list $TMP.phen > $TMP.delete

echo ""
echo ""
$THIS/makeDistTree  -threads 5  -data $INPUT  -dissim_attr $DISSIM  -delete $TMP.delete  \
  -variance $VARIANCE  $DISSIM_COEFF_OPTION  -dissim_power $DISSIM_POWER  \
  -optimize  -subgraph_iter_max 20  -noqual  -output_feature_tree $TMP.feature_tree

echo ""
echo ""
$THIS/makeFeatureTree  -input_tree $TMP.feature_tree  -features $PHEN  -prefer_gain  -nominal_singleton_is_optional  -qual $TMP.qual


rm $TMP*  
