#!/bin/csh -f

if ($# != 2) then
  echo "Build a distance tree with complete pair-wise dissimilarity matrix using the incteremntal tree data structure"
  echo "#1: incremental distance tree directory"
  echo "#2: List of objects"
  exit 1
endif


if (! -z $1/tree) then
  echo "$1/tree must be empty"
  exit 1
endif

set N = `ls $1/new/ | wc -l`
if ($N[1]) then
  echo "$1/new/ must be empty"
  exit 1
endif


list2pairs $2 > $1/dissim_request
if ($?) exit 1

distTree_inc_request2dissim.sh $1 $1/dissim_request $1/dissim
if ($?) exit 1
rm $1/dissim_request

pairs2attr2 $1/dissim 1 cons 6 -distance > $1/data.dm
if ($?) exit 1

echo ""
echo "Tree ..."
makeDistTree  -data $1/data  -dissim cons  -optimize  -reroot  -root_topological  -output_tree $1/tree > $1/hist/makeDistTree.1
if ($?) exit 1
  # -remove_outliers $1/outlier.add  
#trav -noprogress $1/outlier.add "cp /dev/null $1/outlier/%f"
#if ($?) exit 1
rm $1/data.dm


echo ""
echo "Database ..."
$1/objects_in_tree.sh $2 1
if ($?) exit 1
#$1/objects_in_tree.sh $1/outlier.add 0
#if ($?) exit 1
#rm $1/outlier.add


