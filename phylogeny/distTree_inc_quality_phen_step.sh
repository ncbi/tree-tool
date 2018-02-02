#!/bin/csh -f

if ($# != 4) then
  echo "Phenotypic quality of an incremental tree vs. #4.tree"
  echo "Requires: no adding outliers"
  echo "#1: incremental tree"
  echo "#2: test objects"
  echo "#3: target tree version"
  echo "#4: temporary file prefix"
  exit 1
endif


if (! -e $1/phen) then
  echo "Directory $1/phen does not exist"
  exit 1
endif


set tmp = $4


echo ""
echo "Target: $3"

tree2obj.sh $1/hist/tree.$3 > $tmp.target
if ($?) exit 1

# QC
setMinus $2 $tmp.target > $tmp.zero
if ($?) exit 1
if (! -z $tmp.zero)  exit 1

setMinus $tmp.target $2 > $tmp.remove
if ($?) exit 1

makeDistTree  -input_tree $1/hist/tree.$3  -remove $tmp.remove  -output_tree $tmp.trees/$3.tree  -output_feature_tree $tmp.feature_tree > /dev/null
if ($?) exit 1

makeFeatureTree  -input_tree $tmp.feature_tree  -features $1/phen  -output_core $tmp.core  -qual $tmp.qual | grep "Feature agreement:"
if ($?) exit 1

echo "match+:"
compareTrees $tmp.tree $tmp.trees/$3.tree  -frequency none | grep -c '^match+'
if ($?) exit 1
