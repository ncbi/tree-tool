#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 5 ]; then
  echo "Phenotypic quality of an incremental tree vs. #4.tree"
  echo "Requires: no adding outliers"
  echo "#1: incremental tree data structure"
  echo "#2: test objects"
  echo "#3: target tree version|'opt'"
  echo "#4: temporary file prefix"
  echo "#5: Directory for output trees"
  exit 1
fi


if [ ! -e $1/phen ]; then
  echo "Directory $1/phen does not exist"
  exit 1
fi


TMP=$4


if [ ! -e $5/ ]; then
  echo "Directory $5/ does not exist"
  exit 1
fi

echo ""
echo "Target: $3"

INTREE=$1/hist/tree.$3
if [ "$3" == opt ]; then
  INTREE=$1/tree.opt
fi

$THIS/tree2obj.sh $INTREE > $TMP.target

# QC
$THIS/../setMinus $2 $TMP.target > $TMP.zero
#if (! -z $TMP.zero)  exit 1
N=`cat $TMP.zero | wc -l`
echo "# Outliers deleted: $N"

$THIS/../setMinus $TMP.target $2 > $TMP.delete
$THIS/makeDistTree  -input_tree $INTREE  -delete $TMP.delete  -output_tree $5/$3.tree  -output_feature_tree $TMP.feature_tree > $5/$3.makeDistTree
$THIS/makeFeatureTree  -input_tree $TMP.feature_tree  -features $1/phen  -nominal_singleton_is_optional  -qual $5/$3.qual > $5/$3.makeFeatureTree
grep ' !' $5/$3.makeFeatureTree

echo "match+:"
$THIS/compareTrees $TMP.tree $5/$3.tree  -frequency none | grep -c '^match+'


rm -f $TMP*
