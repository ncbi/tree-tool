#!/bin/csh -f

if ($# != 3) then
  echo "Subsampling bootstrap"
  echo "#1: Input: #1.tree, #1.trees/. Output: #1.bootstrap-#2-#3"
  echo "#2: # replicas"
  echo "#3: Frequent nodes filtering: none|directed|undirected"
  exit 1
endif


echo "$3 ..."

set N = `ls $1.trees/ | wc -l`
if ($N[1] != $2) then
  echo "# Files: $N[1]  Should be: $2"
  exit 1
endif
trav -step 1 $1.trees "compareTrees $1.tree %d/%f -frequency $3" | grep '^match[+-]' | sort | uniq -c > $1.bootstrap-$2-$3
if ($?) exit 1
bootstrapReport  $1.bootstrap-$2-$3  -replicas $2  
if ($?) exit 1


