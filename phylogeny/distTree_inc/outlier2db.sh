#!/bin/bash --noprofile
source bash_common.sh
if [ $# -ne 2 ]; then
  echo "Record an outlier in a database"
  echo "#1: object id"
  echo "#2: outlier type"
  exit 1
fi
OBJ=$1
OUTLIER=$2


error "$0 is not implemented"
