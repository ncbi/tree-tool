#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Quality of the initial tree distance tree data structure"
  echo "#1: Incremental distance tree directory"
  exit 1
fi
INC=$1


if [ -e $INC/phen ]; then
  echo ""
  echo "Quality of the initial tree ..."
  VER=`cat $INC/version`
	$THIS/tree2obj.sh $INC/hist/tree.1 > $INC/_init.list
	PHEN_LARGE=`cat $INC/phen_large`
  $THIS/tree_quality_phen.sh $INC/tree $INC/_init.list $INC/phen $PHEN_LARGE 0 "" > $INC/hist/makeFeatureTree-tree1.$VER
  rm $INC/_init.list
	grep ' !' $INC/hist/makeFeatureTree-tree1.$VER
fi


