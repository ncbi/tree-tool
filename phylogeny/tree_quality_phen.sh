#!/bin/csh -f

if ($# != 2) then
  echo "Phenotypic quality of a distance tree"
  echo "Output: core, qual, gain_nodes, disagreement_nodes, disagreement_nodes.txt"
  echo "#1: distance tree"
  echo "#2: phen/"
  exit 1
endif


set tmp = `mktemp`
if ($?) exit 1
echo $tmp


makeDistTree  -input_tree $1  -output_feature_tree $tmp.feature_tree  > $tmp.distTree
if ($?) exit 1

makeFeatureTree  -input_tree $tmp.feature_tree  -features $2  -output_core core \
  -qual qual \
  -gain_nodes gain_nodes \
  -disagreement_nodes disagreement_nodes
if ($?) exit 1
cut -f 1 disagreement_nodes | sort | uniq -c | sort -n -k 1 -r > disagreement_nodes.txt
if ($?) exit 1
echo ""
wc -l disagreement_nodes.txt


rm -f $tmp*
