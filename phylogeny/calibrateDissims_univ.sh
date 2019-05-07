#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 6 ]; then
  echo "Output: #1-out.dm, hmm-univ.stat"
  echo "#1: Input .dm-file without '.dm'"
  echo "#2: delete hybrids (0/1)"
  echo "#3: power"
  echo "#4: outlierSEs"
  echo "#5: dissim_coeff: 0 - variance = lin, else linVar"
  echo "#6: phen/"
  exit 1
fi
INPUT=$1
DELETE_HYBRIDS=$2
POWER=$3
OUTLIER_SES=$4
DISSIM_COEFF=$5
PHEN=$6


echo ""
echo "power = $POWER"
echo "outlier_SEs = $OUTLIER_SES"


TMP=`mktemp`
#echo $TMP 


$THIS/../dm/positiveAverage $INPUT $POWER $OUTLIER_SES  -output_dissim $TMP.pairs > hmm-univ.stat
tail -n +5 $TMP.pairs.dm | sed 's/-/ /1' > $TMP.pairs
$THIS/../dm/pairs2attr2 $TMP.pairs 1 cons 6  -distance > $INPUT-out.dm


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


$THIS/../dm/dm2objs $INPUT-out > $TMP.list
ls $PHEN > $TMP.phen
$THIS/../setMinus $TMP.list $TMP.phen > $TMP.delete


echo ""
echo ""
$THIS/makeDistTree  -threads 5  -data $INPUT-out  -dissim_attr cons  -variance $VARIANCE  $DISSIM_COEFF_OPTION  -delete $TMP.delete  -optimize  -subgraph_iter_max 20  $HYBRID  -noqual  -output_feature_tree $TMP.feature_tree


echo ""
echo ""
$THIS/makeFeatureTree  -input_tree $TMP.feature_tree  -features $PHEN  -prefer_gain  -nominal_singleton_is_optional  -qual $TMP.qual


rm $TMP*  
