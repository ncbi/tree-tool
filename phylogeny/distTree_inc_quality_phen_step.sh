#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 5 ]; then
  echo "Phenotypic quality of an incremental tree vs. #4.tree"
  echo "Requires: no adding outliers"
  echo "#1: incremental tree data structure"
  echo "#2: test objects"
  echo "#3: target tree version|'opt'"
  echo "#4: temporary file prefix"
  echo "#5: directory for output trees"
  exit 1
fi
INC=$1
OBJ_LIST=$2
TTV=$3
TMP=$4
OUTDIR=$5


if [ ! -e $INC/phen ]; then
  error "Directory $INC/phen does not exist"
fi


if [ ! -d $OUTDIR ]; then
  error "Directory $OUTDIR/ does not exist"
fi

section "Target: $TTV"

INTREE=$INC/hist/tree.$TTV
if [ "$TTV" == "opt" ]; then
  INTREE=$INC/tree.opt
fi

$THIS/tree2obj.sh $INTREE > $TMP.target

# QC
$THIS/../setMinus $OBJ_LIST $TMP.target > $TMP.zero
#if (! -z $TMP.zero)  exit 1
N=$( < $TMP.zero  wc -l )
echo "# Outliers deleted: $N"

$THIS/../setMinus $TMP.target $OBJ_LIST > $TMP.delete
$THIS/makeDistTree  -input_tree $INTREE  -delete $TMP.delete  -output_tree $OUTDIR/$TTV.tree  -output_feature_tree $TMP.feature_tree > $OUTDIR/$TTV.makeDistTree
$THIS/makeFeatureTree  -input_tree $TMP.feature_tree  -features $INC/phen  -nominal_singleton_is_optional  -qual $OUTDIR/$TTV.qual > $OUTDIR/$TTV.makeFeatureTree
grep 'V !$' $OUTDIR/$TTV.makeFeatureTree

echo "match+:"
$THIS/compareTrees $TMP.tree $OUTDIR/$TTV.tree  -frequency none | grep -c '^match+'


rm -f $TMP*
