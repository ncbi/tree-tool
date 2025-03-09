#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Print objects sorted by criterion descending"
  echo "#1: distance tree"
  exit 1
fi
TREE=$1


grep 'leaf_error=' $TREE | sed 's/:.* leaf_error=/#/1' | sed 's/^ *//1' | sed 's/ \+[^ ]*$//1' | tr '#' ' ' | sort -k 2 -r -g

