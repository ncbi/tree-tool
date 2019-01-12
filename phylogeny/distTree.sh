#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Make a distance tree"
  echo "Output: #1.tree, #1.makeDistTree"
  echo "#1: Input .dm file without .dm"
  echo "#2: Dissimilarity attribute in #1"
  exit 1
fi


# PAR
MDS=0  # 1
VARIANCE=linExp # lin
SPARSE=""  # -sparse  
WHOLE=""  # -whole  


input_tree=""
if [ $MDS == 1 ]; then
  echo "mdsTree.sh ..."
  $THIS/../dm/mdsTree.sh $1 $2 2 >& /dev/null
  input_tree="-input_tree $1.dir/"
fi

$THIS/makeDistTree  $input_tree  -data $1  -dissim $2  -variance $VARIANCE  $SPARSE  -optimize  $WHOLE  -output_tree $1.tree  
if [ $MDS == 1 ]; then
  rm -r $1.dir/
fi

echo ""
$THIS/makeDistTree  -input_tree $1.tree  -data $1  -dissim $2  -variance $VARIANCE  $SPARSE > $1.makeDistTree

