#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Evaluate quality of an incremental tree by subsampling"
  echo "Requires: no added outliers"
  echo "#1: incremental tree"
  exit 1
fi


TMP=`mktemp`


grep '^# Discernible leaves:'  $1/hist/makeDistTree.1
echo ""


$THIS/tree2obj.sh $1/hist/tree.1 > $TMP.init

echo "First tree  ..."
$THIS/../list2pairs $TMP.init > $TMP.pairs
$THIS/distTree_inc_request2dissim.sh $1 $TMP.pairs $TMP.dissim
$THIS/../dm/pairs2attr2 $TMP.dissim 1 cons 6 -distance > $TMP.dm
cp $1/hist/tree.1 $TMP.tree
$THIS/subsample.sh $TMP cons distTree.sh 1

echo ""
echo "Last tree  ..."
$THIS/tree2obj.sh $1/tree.opt > $TMP.opt
$THIS/../setMinus $TMP.init $TMP.opt > $TMP.zero
if [ -s $TMP.zero ]; then
  ls -laF $TMP.zero
  exit 1
fi
$THIS/../setMinus $TMP.opt $TMP.init > $TMP.delete
$THIS/makeDistTree  -input_tree $1/tree.opt  -delete $TMP.delete  -output_tree $TMP-last.tree  > /dev/null
mv $TMP.dm $TMP-last.dm
$THIS/subsample.sh $TMP-last cons distTree.sh 1


rm -r $TMP*
