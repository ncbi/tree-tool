#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 4 ]; then
  echo "Optimize a distance tree and evaluate"
  echo "#1: incremental distance tree directory"
  echo "#2: number of optimization iterations"
  echo "#3: parameters after -variance"
  echo "#4: output tree"
  echo "Time: O(n log^4(n))"
  exit 1
fi
INC=$1
ITER_MAX=$2
PAR="$3"
OUT_TREE=$4


if [ ! -e $INC/phen ]; then
  echo "No $INC/phen"
  exit 1
fi

$THIS/makeDistTree  -threads 15  -data $INC/  -variance $PAR  -optimize  -skip_len  -subgraph_iter_max $ITER_MAX  -output_tree $OUT_TREE
  # -reinsert  

echo ""
echo "Quality ..."
PHEN_LARGE=`cat $INC/phen_large`
$THIS/tree_quality_phen.sh $OUT_TREE "" $INC/phen $PHEN_LARGE 1 ""

