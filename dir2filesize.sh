#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/bash_common.sh
if [ $# -ne 1 ]; then
  echo "Sort files by size"
  echo "#1: Input directory"
  exit 1
fi
DIR=$1


ls -laF $DIR | tail -n +4 | sort -k 5 -n
