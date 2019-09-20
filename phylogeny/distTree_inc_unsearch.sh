#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Restore #1/search/#2/"
  echo "#1: incremental distance tree directory"
  echo "#2: new object"
  exit 1
fi
INC=$1
OBJ=$2


touch $INC/new/$OBJ
rm -r $INC/search/$OBJ/

