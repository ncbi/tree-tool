#!/bin/bash
source bash_common.sh
if [ $# -ne 1 ]; then
  echo "Unhybrid objects"
  echo "#1: incremental distance tree directory"
  exit 1
fi
INC=$1


if [ -e $INC/unhybrid ]; then
  exit 1
fi
$INC/db2unhybrid.sh $INC/unhybrid

VER=`cat $INC/version`

if [ -s $INC/unhybrid ]; then
  echo ""
	uniq.sh $INC/unhybrid
  wc -l $INC/unhybrid
  trav $INC/unhybrid "mv $INC/outlier/%f $INC/new/"
  $INC/objects_in_tree.sh $INC/unhybrid null
	mv $INC/unhybrid $INC/hist/unhybrid.$VER
else
  rm $INC/unhybrid
fi


