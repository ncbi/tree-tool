#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 9 ]; then
  echo "Output: hmm-univ.stat, positiveAverage.out, data.dm, tree"
  echo "#1: input .dm-file without '.dm' created by univ_separate.sh"
  echo "#2: delete hybrids (0/1)"
  echo "#3: dissimilarity power (> 0)"
  echo "#4: outlierSEs"
  echo "#5: dissim_coeff: >=0 (>0 <=> variance = linExp)"
  echo "#6: variance_power (NAN <=> variance = linExp)"
  echo "#7: variance_dissim (0/1)"
  echo "#8: phen/"
  echo "#9: phen/ is large (0/1)"
  exit 1
fi
INPUT=$1
DELETE_HYBRIDS=$2
DISSIM_POWER=$3
OUTLIER_SES=$4
DISSIM_COEFF=$5
VAR_POWER=$6
VARIANCE_DISSIM=$7
PHEN=$8
LARGE=$9


TMP=`mktemp`
echo $TMP 


section "Estimating hmm-univ.stat"
$THIS/../dm/positiveAverage $INPUT $DISSIM_POWER $OUTLIER_SES hmm-univ.stat  -ignoreZero  -iter_max 30  -output_dissim $TMP.pairs > positiveAverage.out
tail -n +5 $TMP.pairs.dm | sed 's/-/ /1' > $TMP.pairs
$THIS/../dm/pairs2dm $TMP.pairs 1 "cons" 6  -distance > data.dm


section "Building tree"
VARIANCE="linExp"
DISSIM_COEFF_OPTION="-dissim_coeff $DISSIM_COEFF"
if [ $DISSIM_COEFF == 0 ]; then
  VARIANCE="pow  -variance_power $VAR_POWER"
  DISSIM_COEFF_OPTION=""
fi
if [ $VARIANCE_DISSIM == 1 ]; then
  VARIANCE="$VARIANCE -variance_dissim"
fi

HYBRID=""
if [ $DELETE_HYBRIDS -eq 1 ]; then
  HYBRID="-hybrid_parent_pairs hybrid_parent_pairs  -delete_hybrids hybrid"
fi

$THIS/makeDistTree  -threads 5  -data data  -dissim_attr "cons"  -variance $VARIANCE  $DISSIM_COEFF_OPTION  -optimize  -subgraph_iter_max 10  $HYBRID  -noqual  -output_tree tree  -output_feature_tree $TMP.feature_tree  


section "Evaluating tree"
LARGE_PAR=""
if [ $LARGE == 1 ]; then
  LARGE_PAR="-large"
fi
$THIS/makeFeatureTree  -input_tree $TMP.feature_tree  -features $PHEN  $LARGE_PAR  -nominal_singleton_is_optional  -qual $TMP.qual


rm $TMP*  
