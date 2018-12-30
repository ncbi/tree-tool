#!/bin/bash
source bash_common.sh
if [ $# -ne 2 ]; then
  echo "Input: sm1K.combo, dissim_scale"
  echo "#1: coeff (> 0)"
  echo "#2: phen/"
  exit 1
fi


combine_dissims sm1K.combo dissim_scale  -coeff $1 > sm1K.dissim

pairs2attr2 sm1K.dissim 1 cons 6  -distance > sm1K.dm
rm sm1K.dissim

echo ""
distTree.sh sm1K cons
rm sm1K.dm

echo ""
echo "Tree quality ..."
makeDistTree  -input_tree sm1K.tree  -output_feature_tree $TMP.feature_tree  > $TMP.distTree
makeFeatureTree  -input_tree $TMP.feature_tree  \
  -features $2  -nominal_singleton_is_optional  -prefer_gain  -output_core core \
  -qual qual \
  -gain_nodes gain_nodes \
  -disagreement_nodes disagreement_nodes

