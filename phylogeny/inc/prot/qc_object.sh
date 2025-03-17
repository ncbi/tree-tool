#!/bin/bash --noprofile
THIS=$( dirname $0 )
source CPP_DIR/bash_common.sh
if [ $# -ne 1 ]; then
  echo "QC an object"
  echo "#1: file or directory with object data"
  exit 1
fi
OBJ=$1


INC=$THIS
CPP_DIR/check_file.sh $INC/../seq/$OBJ 1
