#!/bin/csh -f

if ($# != 1) then
  echo "Subsampling bootstrap"
  echo "#1: Input: #1.tree, #1.trees/. Output: #1.bootstrap"
  exit 1
endif


trav $1.trees "compareTrees $1.tree %d/%f" | grep '^match[+-] ' | sort | uniq -c > $1.bootstrap
if ($?) exit 1
#bootstrapReport $1.bootstrap -replicas 100 | sort -k 4 -n -r
bootstrapReport  $1.bootstrap  -replicas 100  
if ($?) exit 1


