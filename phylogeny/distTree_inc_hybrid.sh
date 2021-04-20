#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Process #1/hybrid.new, #1/hybrid_parent_pairs (optional)"
  echo "#1: incremental distance tree directory"
  exit 1
fi
INC=$1


VER=`cat $INC/version`


if [ -e $INC/hybrid_parent_pairs ]; then
  mv $INC/hybrid_parent_pairs $INC/hist/hybrid_parent_pairs.$VER
fi


if [ ! -e $INC/hybrid.new ]; then
  exit 0
fi
wc -l $INC/hybrid.new

if [ ! -s $INC/hybrid.new ]; then
  rm $INC/hybrid.new
  exit 0
fi

# $INC/hybrid-indiscern
$THIS/hybrid2list.sh $INC/hybrid.new > $INC/hybrid-indiscern.raw
$THIS/distTree_inc_expand_indiscern.sh $INC $INC/hybrid-indiscern.raw > $INC/hybrid-indiscern
rm $INC/hybrid-indiscern.raw

$INC/objects_in_tree.sh $INC/hybrid-indiscern "null"
$THIS/../trav $INC/hybrid-indiscern "$INC/outlier2db.sh %f auto_hybrid"  

mv $INC/hybrid-indiscern $INC/hist/hybrid-indiscern.$VER
mv $INC/hybrid.new $INC/hist/hybrid.$VER
