#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 8 ]; then
  echo "Output: hmm-univ.stat, positiveAverage.out, data.dm, tree"
  echo "#1: input .dm-file without '.dm' created by univ_separate.sh"
  echo "#2: phen/"
  echo "#3: phen/ is large (0/1)"
 #echo "#2: delete hybrids (0/1)"
  echo "#4: dissimilarity power (> 0)"
  echo "#5: outlierSEs"
 #echo "#5: ignoreZero (0/1)"
  echo "#6: dissim_coeff: >=0 (>0 <=> variance = linExp)"
  echo "#7: variance_power (NAN <=> variance = linExp)"
  echo "#8: variance_dissim (0/1)"
  exit 1
fi
INPUT=$1
PHEN=$2
LARGE=$3
DISSIM_POWER=$4
OUTLIER_SES=$5
DISSIM_COEFF=$6
VAR_POWER=$7
VARIANCE_DISSIM=$8


DELETE_HYBRIDS=0
IGNORE_ZERO=1


section "Estimating hmm-univ.stat"
TMP=`mktemp`
comment $TMP 

# $TMP.pairs
IGNORE_ZERO_PAR=""
if [ $IGNORE_ZERO == 1 ]; then
  IGNORE_ZERO_PAR="-ignoreZero"
fi
$THIS/../dm/positiveAverage $INPUT $DISSIM_POWER $OUTLIER_SES hmm-univ.stat  $IGNORE_ZERO_PAR  -iter_max 30  -output_dissim $TMP.pairs > positiveAverage.out
if [ -s positiveAverage.out ]; then
  cat positiveAverage.out
  exit 1
fi
rm positiveAverage.out

tail -n +5 $TMP.pairs.dm | sed 's/-/ /1' > $TMP.pairs
$THIS/../dm/conversion/pairs2dm $TMP.pairs 1 "cons" 6  -distance > data.dm

rm $TMP*  


section "Parameters after -variance"
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

PARAM="$VARIANCE  $DISSIM_COEFF_OPTION  $HYBRID"
warning "$PARAM"
echo ""
$THIS/calibrateDissims.sh data "cons" "$PARAM" 10 $PHEN $LARGE

