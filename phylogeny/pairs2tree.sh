#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 6 ]; then
  echo "Given a dissimilarity build a tree and evaluate the quality"
  echo "#1: file with dissimilarities, line format: <obj1> <obj2> <dissim>"
  echo "#2: directory with reference classifications or ''"
    # #2 is small (one-level)
  echo "#3: parameters after '-variance' in makeDistTree "
  echo "#4: object list to measure quality for or ''"
  echo "#5: output tree file"
  echo "#6: output tree format: ASNT, newick, dm"
  exit 1
fi
DISSIM=$1
PHEN="$2"
VARIANCE=$3
EVAL_LIST="$4"
TREE=$5
FORMAT=$6


TMP=`mktemp`


section "Preparing data"
$THIS/../dm/pairs2dm $DISSIM 1 symbet 8 -distance -qc > $TMP.dm

section "Building tree"
$THIS/makeDistTree  -data $TMP  -dissim_attr symbet  -variance $VARIANCE  -optimize  -subgraph_iter_max 5  -output_tree $TMP.tree  -threads 10

section "Creating $TREE"
$THIS/printDistTree $TMP.tree -qc  -order  -format $FORMAT > $TREE

if [ $PHEN ]; then
  section "Calculating miscongruence of $TREE with $PHEN"
  $THIS/tree_quality_phen.sh $TMP.tree "$EVAL_LIST" $PHEN 0 1 ''
fi
 

rm -f $TMP*


