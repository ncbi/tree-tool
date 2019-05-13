#!/bin/bash 
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Build a tree using protein hashes"
  echo "#1: Directory with protein hashes"
  echo "#2: Output Newick file"
  exit 1
fi
IN=$1
NEWICK=$2


TMP=`mktemp`
#echo $TMP


ls $IN/ > $TMP.list

$THIS/hash2dissim $TMP.list $IN $TMP -intersection_min 10

echo ""
$THIS/makeDistTree  -data $TMP  -dissim_attr cons  -variance linExp  -optimize  -output_tree $TMP.tree
$THIS/printDistTree $TMP.tree > $NEWICK
 

rm $TMP*
