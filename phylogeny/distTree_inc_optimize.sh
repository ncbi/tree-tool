#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 4 ]; then
  echo "Optimize a distance tree and evaluate"
  echo "#1: incremental distance tree directory (no updating)"
  echo "#2: target list of objects | '' - all"
  echo "#3: new contents of #1/variance (non-empty string)"
  echo "#4: output tree"
  exit 1
fi
INC=$1
TARGET="$2"
VARIANCE="$3"
OUT_TREE=$4


section "Tree"
THREADS=$( file2var $INC/threads 15 )
$THIS/makeDistTree  -threads $THREADS  -data $INC/  -variance $VARIANCE  -optimize  -skip_len  -subgraph_iter_max 5  -output_tree $OUT_TREE

if [ -e $INC/phen ]; then
  super_section "Quality"
  LARGE=0
  if [ -e $INC/large ]; then
    LARGE=1
  fi
  $THIS/tree_quality_phen.sh $OUT_TREE "$TARGET" $INC/phen $LARGE 1 ""
fi
