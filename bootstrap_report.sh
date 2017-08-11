#!/bin/csh -f

if ($# != 2) then
  echo "Subsampling bootstrap"
  echo "#1: Input: #1.tree, #1.trees/. Output: #1.bootstrap"
  echo "#2: # replicas"
  exit 1
endif


trav $1.trees "compareTrees $1.tree %d/%f" | grep '^match[+-] ' | sort | uniq -c > $1.bootstrap-$2
if ($?) exit 1
bootstrapReport  $1.bootstrap-$2  -replicas $2  
if ($?) exit 1


