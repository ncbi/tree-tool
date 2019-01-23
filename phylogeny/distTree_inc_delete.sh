#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Delete a list of objects from an incremental distance tree"
  echo "Update: #1/"
  echo "#1: incremental distance tree directory"
  echo "#2: List of objects to delete"
  exit 1
fi


VER=`cat $1/version`

cp $1/tree $1/hist/tree.$VER
gzip $1/hist/tree.$VER

VER=$(( $VER + 1 ))
echo $VER > $1/version

VARIANCE=`cat $1/variance`

# Cf. distTree_inc_new.sh
$THIS/makeDistTree  -threads 15  -data $1/  -variance $VARIANCE \
  -delete $2  \
  -optimize  -skip_len  -subgraph_iter_max 1 \
  -noqual \
  -output_tree $1/tree.new > $1/hist/makeDistTree-delete.$VER
mv $1/tree.new $1/tree

echo ""
$1/objects_in_tree.sh $2 null
$THIS/../trav $2 "rm -f $1/hybrid/%f"
$THIS/../trav $2 "rm -f $1/new/%f"

cp $2 $1/hist/delete.$VER


$THIS/distTree_inc_tree1_quality.sh $1
