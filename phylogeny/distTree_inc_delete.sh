#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Delete a list of objects from an incremental distance tree"
  echo "Update: #1/"
  echo "Delete: #2"
  echo "#1: Incremental distance tree directory"
  echo "#2: List of objects to delete; to be moved into #1/hist/"
  exit 1
fi
INC=$1
DEL=$2


VER=`cat $INC/version`

cp $INC/tree $INC/hist/tree.$VER
if [ $VER -gt 1 ]; then
  gzip $INC/hist/tree.$VER
fi

VER=$(( $VER + 1 ))
echo $VER > $INC/version

VARIANCE=`cat $INC/variance`

section "Building tree ..."
# Cf. distTree_inc_new.sh
$THIS/makeDistTree  -threads 15  -data $INC/  -variance $VARIANCE \
  -delete $DEL  \
  -optimize  -skip_len  -subgraph_iter_max 1 \
  -noqual \
  -output_dissim $INC/dissim.new \
  -output_tree $INC/tree.new > $INC/hist/makeDistTree-delete.$VER
  
mv $INC/tree.new $INC/tree

wc -l $INC/dissim
wc -l $INC/dissim.new
mv $INC/dissim.new $INC/dissim

section "Database ..."
$INC/objects_in_tree.sh $DEL null

$THIS/distTree_inc_new_list.sh $INC > $INC/new.list
$THIS/../setIntersect.sh $DEL $INC/new.list 0 > $INC/new-del.list
rm $INC/new.list
$THIS/distTree_inc_new_cmd.sh $INC "rm" $INC/new-del.list
rm $INC/new-del.list

mv $DEL $INC/hist/delete.$VER


section "QC ..."
$INC/qc.sh go  


$THIS/distTree_inc_tree1_quality.sh $INC


echo ""
date

