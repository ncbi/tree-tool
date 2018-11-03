#!/bin/bash
source bash_common.sh
if [ $# -ne 1 ]; then
  echo "Quality of the initial tree distance tree data structure"
  echo "#1: incremental distance tree directory"
  exit 1
fi


if [ -e $1/phen ]; then
  echo ""
  echo "Quality ..."
  VER=`cat $1/version`
	tree2obj.sh $1/hist/tree.1 > $1/_init.list
	tree2obj.sh $1/tree > $1/_cur.list
	setMinus $1/_cur.list $1/_init.list >  $1/_delete.list
	setMinus $1/_init.list $1/_cur.list >> $1/_delete.list
	rm $1/_init.list 
	rm $1/_cur.list
	makeDistTree  -input_tree $1/tree  -delete $1/_delete.list  -output_feature_tree $1/_feature_tree >& /dev/null
	rm $1/_delete.list
	makeFeatureTree  -input_tree $1/_feature_tree  -prefer_gain  -features $1/phen  -output_core $1/_core  -qual $1/_qual > $1/hist/makeFeatureTree-tree1.$VER
	rm $1/_feature_tree
	rm $1/_core
	rm $1/_qual
	grep ' !' $1/hist/makeFeatureTree-tree1.$VER
fi


