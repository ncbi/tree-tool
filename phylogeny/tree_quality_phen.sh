#!/bin/csh -f

if ($# != 2) then
  echo "Phenotypic quality of a distance tree"
  echo "Output: core, qual, gain_nodes"
  echo "#1: distance tree"
  echo "#2: phen/"
  exit 1
endif


set tmp = `mktemp`
if ($?) exit 1


makeDistTree  -input_tree $1  -output_feature_tree $tmp.feature_tree > /dev/null
if ($?) exit 1

makeFeatureTree  -input_tree $tmp.feature_tree  -features $2  -output_core core  -qual qual  -gain_nodes gain_nodes
if ($?) exit 1


rm -f $tmp*
