#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 4 ]; then
  echo "Set names to a distance tree by phenotypes"
  echo "Output: core, qual, gain_nodes, disagreement_nodes, disagreement_nodes.txt"
  echo "#1: distance tree"
  echo "#2: target list of objects | '' - all"
  echo "#3: phen/"
  echo "#4: large directories (0/1)"
  exit 1
fi
TREE=$1
TARGET="$2"
PHEN=$3
LARGE=$4


TMP=$( mktemp )
#echo $TMP


DELETE=""
if [ "$TARGET" ]; then
  sort -cu $TARGET
  $THIS/tree2obj.sh $TREE > $TMP.cur
  $THIS/../setMinus $TMP.cur $TARGET > $TMP.del
  DELETE="-delete $TMP.del  -check_delete"
fi
$THIS/makeDistTree  -threads 15  -input_tree $TREE  $DELETE  -noqual  -output_feature_tree $TMP.feature_tree  > $TMP.distTree

LARGE_PAR=""
if [ $LARGE == 1 ]; then
  LARGE_PAR="-large"
fi
$THIS/makeFeatureTree  -threads 15  -input_tree $TMP.feature_tree  -features $PHEN  $LARGE_PAR    -save_mem  -prefer_gain \
  -qual qual  -gain_nodes gain_nodes  -disagreement_nodes disagreement_nodes
  # removing -prefer_gain can put a rare low rank too high in the tree
  
cut -f 1 disagreement_nodes | sort | uniq -c | sort -n -k 1 -r > disagreement_nodes.txt
set +o errexit
cat disagreement_nodes.txt | grep -v ":" | grep -v ' 1 ' > disagreement_objects
set -o errexit
echo ""
wc -l disagreement_nodes.txt


rm -f $TMP*
