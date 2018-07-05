#!/bin/csh -f

if ($# != 5) then
  echo "Phenotypic quality of an incremental tree vs. #4.tree"
  echo "Requires: no adding outliers"
  echo "#1: incremental tree data structure"
  echo "#2: test objects"
  echo "#3: target tree version|'opt'"
  echo "#4: temporary file prefix"
  echo "#5: Directory for output trees"
  exit 1
endif


if (! -e $1/phen) then
  echo "Directory $1/phen does not exist"
  exit 1
endif



set tmp = $4


if (! -e $5/) then
  echo "Directory $5/ does not exist"
  exit 1
endif

echo ""
echo "Target: $3"

set INTREE = $1/hist/tree.$3
if ("$3" == opt)  set INTREE = $1/tree.opt

tree2obj.sh $INTREE > $tmp.target
if ($?) exit 1

# QC
setMinus $2 $tmp.target > $tmp.zero
if ($?) exit 1
#if (! -z $tmp.zero)  exit 1
set N = `wc -l $tmp.zero`
echo "# Outliers deleted: $N[1]"

setMinus $tmp.target $2 > $tmp.delete
if ($?) exit 1

makeDistTree  -input_tree $INTREE  -delete $tmp.delete  -output_tree $5/$3.tree  -output_feature_tree $tmp.feature_tree > $5/$3.makeDistTree
if ($?) exit 1

makeFeatureTree  -input_tree $tmp.feature_tree  -features $1/phen  -output_core $tmp.core  -qual $5/$3.qual > $5/$3.makeFeatureTree
if ($?) exit 1
grep ' !' $5/$3.makeFeatureTree
if ($?) exit 1

echo "match+:"
compareTrees $tmp.tree $5/$3.tree  -frequency none | grep -c '^match+'
if ($?) exit 1
