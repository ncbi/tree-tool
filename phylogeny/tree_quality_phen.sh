#!/bin/bash
source bash_common.sh
if [ $# -ne 2 ]; then
  echo "Phenotypic quality of a distance tree"
  echo "Output: core, qual, gain_nodes, disagreement_nodes, disagreement_nodes.txt"
  echo "#1: distance tree"
  echo "#2: phen/"
  exit 1
fi


TMP=`mktemp`
echo $TMP


makeDistTree  -threads 15  -input_tree $1  -output_feature_tree $TMP.feature_tree  > $TMP.distTree

makeFeatureTree  -threads 15  -input_tree $TMP.feature_tree  -features $2  -prefer_gain  -output_core core \
  -qual qual \
  -gain_nodes gain_nodes \
  -disagreement_nodes disagreement_nodes
cut -f 1 disagreement_nodes | sort | uniq -c | sort -n -k 1 -r > disagreement_nodes.txt
echo ""
wc -l disagreement_nodes.txt


rm -f $TMP*
