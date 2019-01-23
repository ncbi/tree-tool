#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Sampling report"
  echo "#1: Input: #1.tree, #1.trees/. Output: #1.sampling-#2-#3"
  echo "#2: # replicas"
  echo "#3: Frequent nodes filtering: none|directed|undirected"
  exit 1
fi


echo "$3 ..."

N=`ls $1.trees/ | wc -l`
if [ $N -ne $2 ]; then
  echo "# Files: $N  Should be: $2"
  exit 1
fi
$THIS/../trav -step 1 $1.trees "$THIS/compareTrees $1.tree %d/%f -frequency $3" | grep '^match[+-]' | sort | uniq -c > $1.sampling-$2-$3
$THIS/sampleReport  $1.sampling-$2-$3  -replicas $2  


