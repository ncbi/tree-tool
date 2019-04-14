#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 6 ]; then
  echo "Input: #3-calibrate.dm"
  echo "Output: hmm-univ.stat"
  echo "#1: power"
  echo "#2: outlierSEs"
  echo "#3: phen/"
  echo "#4: Input .dm-file without '.dm'"
  echo "#5: delete hybrids (0/1)"
  echo "#6: dissim_coeff: 0 - variance = lin, else linVar"
  exit 1
fi
POWER=$1
OUTLIER_SES=$2
PHEN=$3
INPUT=$4
DELETE_HYBRIDS=$5
DISSIM_COEFF=$6


echo ""
echo "power = $POWER"
echo "outlier_SEs = $OUTLIER_SES"


TMP=`mktemp`
#echo $TMP 


$THIS/../dm/positiveAverage $INPUT $POWER $OUTLIER_SES  -output_dissim $TMP.pairs > hmm-univ.stat
tail -n +5 $TMP.pairs.dm | sed 's/-/ /1' > $TMP.pairs
$THIS/../dm/pairs2attr2 $TMP.pairs 1 cons 6  -distance > $TMP.dm


echo ""
echo ""

HYBRID=""
if [ $DELETE_HYBRIDS -eq 1 ]; then
  HYBRID="-hybrid_parent_pairs hybrid_parent_pairs  -delete_hybrids hybrid"
fi

VARIANCE=linExp
DISSIM_COEFF_OPTION="-dissim_coeff $DISSIM_COEFF"
if [ $DISSIM_COEFF == 0 ]; then
  VARIANCE=lin
  DISSIM_COEFF_OPTION=""
fi


$THIS/../dm/dm2objs $TMP > $TMP.list
ls $PHEN > $TMP.phen
$THIS/../setMinus $TMP.list $TMP.phen > $TMP.delete


echo ""
echo ""
$THIS/makeDistTree  -threads 5  -data $TMP  -dissim_attr cons  -variance $VARIANCE  $DISSIM_COEFF_OPTION  -delete $TMP.delete  -optimize  -subgraph_iter_max 20  $HYBRID  -noqual  -output_feature_tree $TMP.feature_tree


echo ""
echo ""
$THIS/makeFeatureTree  -input_tree $TMP.feature_tree  -features $PHEN  -prefer_gain  -nominal_singleton_is_optional  -qual $TMP.qual


rm $TMP*  
