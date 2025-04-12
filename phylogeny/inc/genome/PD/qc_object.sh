#!/bin/bash --noprofile
#source bash_common.sh
if [ $# -ne 1 ]; then
  echo "QC an object"
  echo "#1: file or directory with object data"
  exit 1
fi
FD=$1


warning "$0 is not implemented"
