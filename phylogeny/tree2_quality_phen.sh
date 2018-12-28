#!/bin/bash
source bash_common.sh
if [ $# -ne 3 ]; then
  echo "Phenotypic quality comparison of 2 distance trees"
  echo "#1: distance tree 1"
  echo "#2: distance tree 2"
  echo "#3: phen/"
  exit 1
fi
T1=$1
T2=$2
PHEN=$3


TMP=`mktemp`
#echo $TMP  


tree2obj.sh $T1 > $TMP.list1
tree2obj.sh $T2 > $TMP.list2

wc -l $TMP.list1
wc -l $TMP.list2

setMinus $TMP.list1 $TMP.list2 > $TMP.list1-del
setMinus $TMP.list2 $TMP.list1 > $TMP.list2-del

makeDistTree -input_tree $T1  -delete $TMP.list1-del  -check_delete  -output_tree $TMP.tree1  -output_feature_tree $TMP.feature_tree1 > $TMP.makeDistTree1
makeDistTree -input_tree $T2  -delete $TMP.list2-del  -check_delete  -output_tree $TMP.tree2  -output_feature_tree $TMP.feature_tree2 > $TMP.makeDistTree2

echo ""
makeFeatureTree  -input_tree $TMP.feature_tree1  -features $PHEN  -nominal_singleton_is_optional  -prefer_gain  -output_core $TMP.core1  -qual $TMP.qual1  -gain_nodes $TMP.gain_nodes1  -disagreement_nodes $TMP.disagreement_nodes1
echo ""
makeFeatureTree  -input_tree $TMP.feature_tree2  -features $PHEN  -nominal_singleton_is_optional  -prefer_gain  -output_core $TMP.core2  -qual $TMP.qual2  -gain_nodes $TMP.gain_nodes2  -disagreement_nodes $TMP.disagreement_nodes2


rm -f $TMP*
