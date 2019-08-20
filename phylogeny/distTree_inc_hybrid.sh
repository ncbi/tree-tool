#!/bin/bash
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

#$INC/hybrid2db.sh $INC/hybrid.new
$THIS/hybrid2list.sh $INC/hybrid.new > $INC/hybrid.add
$INC/objects_in_tree.sh $INC/hybrid.add null
  # was: 0
$THIS/../trav $INC/hybrid.add "$INC/outlier2db.sh %f auto_hybrid"  
rm $INC/hybrid.add

mv $INC/hybrid.new $INC/hist/hybrid.$VER
