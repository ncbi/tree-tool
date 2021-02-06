#!/bin/bash
THIS=`dirname $0`
source /home/brovervv/code/cpp/bash_common.sh
if [ $# -ne 1 ]; then
  echo "QC an object"
  echo "#1: file or directory with object data"
  exit 1
fi
FD=$1


exit 0
