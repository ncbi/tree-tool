#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Print leaf statistics for an incremental distance tree"
  echo "#1: distance tree data"
  exit 1
fi

cat $1/hist/leaf* | sort -k 6 -g

