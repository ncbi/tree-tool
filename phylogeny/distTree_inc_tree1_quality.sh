#!/bin/bash
source bash_common.sh
if [ $# -ne 1 ]; then
  echo "Quality of the initial tree distance tree data structure"
  echo "#1: incremental distance tree directory"
  exit 1
fi


if [ -e $1/phen ]; then
  echo ""
  echo "Quality of the initial tree ..."
  VER=`cat $1/version`
	tree2obj.sh $1/hist/tree.1 > $1/_init.list
  distTree_inc_target_quality.sh $1 $1/_init.list > $1/hist/makeFeatureTree-tree1.$VER
  rm $1/_init.list
	grep ' !' $1/hist/makeFeatureTree-tree1.$VER
fi


