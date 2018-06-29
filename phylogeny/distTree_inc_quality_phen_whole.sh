#!/bin/csh -f

if ($# != 4) then
  echo "Phenotypic quality of an incremental tree: initial tree"
  echo "Requires: no added outliers"
  echo "#1: incremental tree data structure"
  echo "#2: target tree version"
  echo "#3: version increment (>= 1)"
  echo "#4: Directory for output trees: #2.tree, #2.makeDistTree, #2.makeFeatureTree, #2.qual"
  exit 1
endif


if ($2 < 1)  exit 1
if ($3 < 1)  exit 1


set tmp = `mktemp`
if ($?) exit 1
echo $tmp  # ??


cp $1/hist/tree.1 $tmp.tree
if ($?) exit 1

tree2obj.sh $1/hist/tree.1 > $tmp.init
if ($?) exit 1


if ($?) exit 1
set i = 1
while ($i <= $2)
	distTree_inc_quality_phen_step.sh $1 $tmp.init $i $tmp $4
	if ($?) exit 1
	@ i = $i + $3
end


rm -fr $tmp*
