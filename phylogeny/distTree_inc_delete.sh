#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Delete a list of objects from an incremental distance tree"
  echo "Update: #1/"
  echo "#1: Incremental distance tree directory"
  echo "#2: List of objects to delete"
  exit 1
fi
INC=$1
DEL=$2


VER=`cat $INC/version`

cp $INC/tree $INC/hist/tree.$VER
gzip $INC/hist/tree.$VER

VER=$(( $VER + 1 ))
echo $VER > $INC/version

VARIANCE=`cat $INC/variance`

# Cf. distTree_inc_new.sh
$THIS/makeDistTree  -threads 15  -data $INC/  -variance $VARIANCE \
  -delete $DEL  \
  -optimize  -skip_len  -subgraph_iter_max 1 \
  -noqual \
  -output_tree $INC/tree.new > $INC/hist/makeDistTree-delete.$VER
mv $INC/tree.new $INC/tree

echo ""
$INC/objects_in_tree.sh $DEL null
$THIS/../trav $DEL "rm -f $INC/new/%f"

cp $DEL $INC/hist/delete.$VER


echo ""
echo "QC ..."
$INC/qc.sh $INC


$THIS/distTree_inc_tree1_quality.sh $INC


echo ""
date

