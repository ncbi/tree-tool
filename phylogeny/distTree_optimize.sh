#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 7 ]; then
  echo "Optimize a distance tree and evaluate"
  echo "#1: input .dm-file without '.dm'"
  echo "#2: input tree"
  echo "#3: target list of objects | '' - all"
  echo "#4: parametsrs after -variance (non-empty string)"
  echo "#5: phen/"
  echo "#6: phen/ is a large directory (0/1)"
  echo "#7: output tree"
  exit 1
fi
DATA=$1
IN_TREE=$2
TARGET="$3"
VARIANCE="$4"
PHEN=$5
LARGE=$6
OUT_TREE=$7


section "Tree"
$THIS/makeDistTree  -data $DATA  -input_tree $IN_TREE  -variance $VARIANCE  -optimize  -subgraph_iter_max 5  -output_tree $OUT_TREE

super_section "Quality"
$THIS/tree_quality_phen.sh $OUT_TREE "$TARGET" $PHEN $LARGE 1 ""
