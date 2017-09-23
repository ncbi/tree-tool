#!/bin/csh -f

if ($# != 2) then
  echo "Make a distance tree"
  echo "Output: #1.tree, #1.makeDistTree"
  echo "#1: Input .dm file without .dm"
  echo "#2: Distance attribute in #1"
  exit 1
endif


# PAR
set mds         = 0  # 1
set variance    = linExp 
set sparse_init = -sparse_init  # ""
set whole       = ""  # -whole  


set input_tree = ""
if ($mds) then
  mdsTree.sh $1 $2 2 >& /dev/null
  if ($?) exit 1
  set input_tree = "-input_tree $1.dir/"
endif

makeDistTree  $input_tree  -data $1  -dissim $2  -variance $variance  $sparse_init  -topology  $whole  -output_tree $1.tree  
if ($?) exit 1
if ($mds) then
  rm -r $1.dir/
  if ($?) exit 1
endif

echo ""
makeDistTree  -input_tree $1.tree  -data $1  -dissim $2  -variance $variance  $sparse_init > $1.makeDistTree
if ($?) exit 1

