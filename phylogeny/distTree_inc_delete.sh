#!/bin/bash
source bash_common.sh
if [ $# -ne 2 ]; then
  echo "Delete a list of objects from an incremental distance tree"
  echo "Update: #1/"
  echo "#1: incremental distance tree directory"
  echo "#2: List of objects to delete"
  exit 1
fi


VER=`cat $1/version`

cp $1/tree $1/hist/tree.$VER

VER=$(( $VER + 1 ))
echo $VER > $1/version

# Cf. distTree_inc_new.sh
makeDistTree  -threads 15  -data $1/  -variance lin \
  -delete $2  \
  -optimize  -skip_len  -subgraph_iter_max 1 \
  -noqual \
  -output_tree $1/tree.new > $1/hist/makeDistTree.$VER
mv $1/tree.new $1/tree

echo ""
$1/objects_in_tree.sh $2 null
trav $2 "rm -f $1/outlier/%f"
trav $2 "rm -f $1/new/%f"

cp $2 $1/hist/delete.$VER


distTree_inc_tree1_quality.sh $1
