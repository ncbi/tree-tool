#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Quality of the initial tree distance tree data structure"
  echo "#1: incremental distance tree directory"
  exit 1
fi


if [ -e $1/phen ]; then
  echo ""
  echo "Quality of the initial tree ..."
  VER=`cat $1/version`
	$THIS/tree2obj.sh $1/hist/tree.1 > $1/_init.list
  $THIS/tree_quality_phen.sh $1/tree $1/_init.list $1/phen > $1/hist/makeFeatureTree-tree1.$VER
  rm $1/_init.list
	grep ' !' $1/hist/makeFeatureTree-tree1.$VER
fi


