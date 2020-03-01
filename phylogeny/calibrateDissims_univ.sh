#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 7 ]; then
  echo "Output: hmm-univ.stat, positiveAverage.out"
  echo "#1: input .dm-file without '.dm'"
 #echo "#2: file with universal attrbutes or '' if all attributes are universal"
 #echo "#3: weight factor of universal attributes, >= 1.0"
  echo "#2: delete hybrids (0/1)"
  echo "#3: dissimilarity power (> 0)"
  echo "#4: outlierSEs"
  echo "#5: dissim_coeff: >0 <=> variance linExp"
  echo "#6: variance power (NAN <=> variance = linExp)"
  echo "#7: phen/ (large storage)"
  exit 1
fi
INPUT=$1
#UNIV="$2"
#UNIV_WEIGHT=$3
DELETE_HYBRIDS=$2
DISSIM_POWER=$3
OUTLIER_SES=$4
DISSIM_COEFF=$5
VAR_POWER=$6
PHEN=$7


echo ""
#echo "univ_weight = $UNIV_WEIGHT"
echo "power = $DISSIM_POWER"
echo "outlier_SEs = $OUTLIER_SES"


TMP=`mktemp`
echo $TMP 


#UNIV_PAR=""
#if [ "$UNIV" ]; then
#  UNIV_PAR="-universal $UNIV"
#fi
# PAR
$THIS/../dm/positiveAverage $INPUT $DISSIM_POWER $OUTLIER_SES hmm-univ.stat  -output_dissim $TMP.pairs > positiveAverage.out
tail -n +5 $TMP.pairs.dm | sed 's/-/ /1' > $TMP.pairs
$THIS/../dm/pairs2attr2 $TMP.pairs 1 cons 6  -distance > $TMP.dm
#cp $TMP.dm data.dm 


echo ""
echo ""

HYBRID=""
if [ $DELETE_HYBRIDS -eq 1 ]; then
  HYBRID="-hybrid_parent_pairs hybrid_parent_pairs  -delete_hybrids hybrid"
fi

VARIANCE=linExp
DISSIM_COEFF_OPTION="-dissim_coeff $DISSIM_COEFF"
if [ $DISSIM_COEFF == 0 ]; then
  VARIANCE="pow  -variance_power $VAR_POWER"
  DISSIM_COEFF_OPTION=""
fi

echo ""
echo ""
$THIS/makeDistTree  -threads 5  -data $TMP  -dissim_attr cons  -variance $VARIANCE  $DISSIM_COEFF_OPTION  -optimize  -subgraph_iter_max 10  $HYBRID  -noqual  -output_feature_tree $TMP.feature_tree  

echo ""
echo ""
$THIS/makeFeatureTree  -input_tree $TMP.feature_tree  -features $PHEN  -large  -prefer_gain  -nominal_singleton_is_optional  -qual $TMP.qual


rm $TMP*  
