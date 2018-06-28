#!/bin/csh -f

if ($# != 1) then
  echo "Evaluate quality of an incremental tree by subsampling"
  echo "Requires: no added outliers"
  echo "#1: incremental tree"
  exit 1
endif


set tmp = $1/tmp
# `mktemp` cannot be used by the grid


grep '^# Discernible leaves:'  $1/hist/makeDistTree.1
if ($?) exit 1
echo ""


tree2obj.sh $1/hist/tree.1 > $tmp.init
if ($?) exit 1


echo "First tree  ..."
list2pairs $tmp.init > $tmp.pairs
if ($?) exit 1

distTree_inc_request2dissim.sh $1 $tmp.pairs $tmp.dissim
if ($?) exit 1

pairs2attr2 $tmp.dissim 1 cons 6 -distance > $tmp.dm
if ($?) exit 1

cp $1/hist/tree.1 $tmp.tree
if ($?) exit 1

subsample.sh $tmp cons distTree.sh 1
if ($?) exit 1


echo ""
echo "Last tree  ..."
tree2obj.sh $1/tree.opt > $tmp.opt
if ($?) exit 1

setMinus $tmp.init $tmp.opt > $tmp.zero
if ($?) exit 1
if (! -z $tmp.zero) then
  ls -laF $tmp.zero
  exit 1
endif
setMinus $tmp.opt $tmp.init > $tmp.remove
if ($?) exit 1

makeDistTree  -input_tree $1/tree.opt  -delete $tmp.remove  -output_tree $tmp-last.tree  > /dev/null
if ($?) exit 1

mv $tmp.dm $tmp-last.dm
if ($?) exit 1

subsample.sh $tmp-last cons distTree.sh 1
if ($?) exit 1


rm -fr $tmp*
