#!/bin/bash
source bash_common.sh
if [ $# -ne 5 ]; then
  echo "Input: #3-calibrate.dm"
  echo "#1: power"
  echo "#2: phen/"
  echo "#3: Input file prefix"
  echo "#4: delete hybrids (0/1)"
  echo "#5: dissim_coeff (default: 1)"
  exit 1
fi


echo ""
echo "power = $1"

calibrateDissims $3-calibrate $1  -output_dissim $3.pairs  > hmm-univ.stat

tail -n +5 $3.pairs.dm | sed 's/-/ /1' > $3.pairs
rm $3.pairs.dm

pairs2attr2 $3.pairs 1 cons 6  -distance > $3.dm
rm $3.pairs

hybrid=""
if [ $4 -eq 1 ]; then
  hybrid="-hybrid_parent_pairs hybrid_parent_pairs  -delete_hybrids hybrid  -delete_all_hybrids"
fi

echo ""
echo ""
makeDistTree  -threads 5  -data $3  -dissim cons  -dissim_coeff $5  -optimize  $hybrid  -noqual  -output_tree tree  -output_feature_tree _feature_tree

echo ""
echo ""
makeFeatureTree  -input_tree _feature_tree  -features $2  -nominal_singleton_is_optional  -output_core _core  -qual _qual

