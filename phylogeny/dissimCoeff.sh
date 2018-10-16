#!/bin/bash
source bash_common.sh
if [ $# -ne 3 ]; then
  echo "Input: #3.dm"
  echo "#1: dissim_coeff"
  echo "#2: phen/"
  echo "#3: dissimilarity data"
  exit 1
fi


echo ""
echo "dissim_coeff = $1"

makeDistTree  -threads 5  -data $3  -dissim cons  -dissim_coeff $1  -optimize  -output_tree tree  -output_feature_tree _feature_tree

echo ""
echo ""
makeFeatureTree  -input_tree _feature_tree  -features $2  -output_core _core  -qual _qual

