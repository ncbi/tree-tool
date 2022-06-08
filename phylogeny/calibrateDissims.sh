#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 6 ]; then
  echo "output: tree"
  echo "#1: input .dm-file without '.dm'"
  echo "#2: dissimilarity attribute in #1"
  echo "#3: parameters of makeDistTree after -variance"
  echo "#4: max. optimization iterations (> 0, normally 10)"
  echo "#5: phen/"
  echo "#6: large directores (0/1)"
  exit 1
fi
INPUT=$1
DISSIM="$2"
VARIANCE_PAR="$3"
ITER_MAX=$4
PHEN=$5
LARGE=$6


TMP=`mktemp`
echo $TMP 


section "Bulding tree"
$THIS/makeDistTree  -threads 5  -data $INPUT  -dissim_attr $DISSIM  -variance $VARIANCE_PAR \
     -optimize  -subgraph_iter_max $ITER_MAX  -noqual  -output_feature_tree $TMP.feature_tree  -output_tree tree

section "Evaluating tree"
LARGE_PAR=""
if [ $LARGE -eq 1 ]; then
  LARGE_PAR="-large"
fi
$THIS/makeFeatureTree  -input_tree $TMP.feature_tree  -features $PHEN  $LARGE_PAR  -nominal_singleton_is_optional  -qual $TMP.qual
  # -prefer_gain  


rm $TMP*  
