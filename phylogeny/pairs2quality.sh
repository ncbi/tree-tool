#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 5 ]; then
  echo "Evaluate the quality of a dissimilarity"
  echo "#1: file with dissimilarities, line format: <obj1> <obj2> <dissim>"
  echo "#2: directory with reference classifications"
    # #2 is small
  echo "#3: parameters after '-variance' in makeDistTree "
  echo "#4: output tree file"
  echo "#5: output tree format: ASNT, newick, dm"
  exit 1
fi
DISSIM=$1
PHEN=$2
VARIANCE=$3
TREE=$4
FORMAT=$5


TMP=`mktemp`


section "Preparing data ..."
$THIS/../dm/pairs2dm $DISSIM 1 symbet 8 -distance -qc > $TMP.dm

section "Building tree ..."
$THIS/makeDistTree  -data $TMP  -dissim_attr symbet  -variance $VARIANCE  -optimize  -subgraph_iter_max 5  -output_tree $TMP.tree  -threads 10

section "Calculating miscongruence with $PHEN ..."
$THIS/tree_quality_phen.sh $TMP.tree "" $PHEN 0 1 ''

$THIS/printDistTree $TMP.tree -qc  -order  -format $FORMAT > $TREE
 

rm -f $TMP*


