#!/bin/csh -f

if ($# != 2) then
  echo "Make a distance tree"
  echo "Output: #1.tree, #1.makeDistTree"
  echo "#1: Input .dm file without .dm"
  echo "#2: Distance attribute in #1"
  exit 1
endif



#mdsTree.sh $1 $2 2 >& /dev/null
#if ($?) exit 1

date

set variance = linExp 
makeDistTree  ""  -data $1  -dissim $2  -topology  -variance $variance  -output_tree $1.tree  
# $1.dir/
if ($?) exit 1
#rm -r $1.dir/
if ($?) exit 1
makeDistTree  $1.tree  -data $1  -dissim $2  -variance $variance  > $1.makeDistTree
if ($?) exit 1

date
