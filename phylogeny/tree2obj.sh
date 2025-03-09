#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Print objects of a distance tree, sorted"
  echo "#1: Distance tree"
  exit 1
fi
TREE=$1


$THIS/tree2obj_raw.sh $TREE | sort

