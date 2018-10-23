#!/bin/bash
source bash_common.sh
if [ $# -ne 2 ]; then
  echo "Process #1/hybrid"
  echo "#1: incremental distance tree directory"
  echo "#2: tree version"
  exit 1
fi


if [ -e $1/unhybrid ]; then
  exit 1
fi


$1/db2unhybrid.sh $1/unhybrid

if [ -s $1/unhybrid ]; then
  echo ""
	uniq.sh $1/unhybrid
  wc -l $1/unhybrid
  trav $1/unhybrid "rm -f $1/outlier/%f"
  trav $1/unhybrid "cp /dev/null $1/new/%f"
  $1/objects_in_tree.sh $1/unhybrid null
	mv $1/unhybrid $1/hist/unhybrid.$2
else
  rm $1/unhybrid
fi

