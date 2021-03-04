#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Print objects sorted by criterion descending"
  echo "#1: Distance tree"
  exit 1
fi

grep 'leaf_error=' $1 | sed 's/:.* leaf_error=/#/1' | sed 's/^ *//1' | sed 's/ \+[^ ]*$//1' | tr '#' ' ' | sort -k 2 -r -g

